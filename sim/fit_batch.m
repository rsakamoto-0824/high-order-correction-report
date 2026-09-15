function theta = fit_batch(A, Y, method, isPenalized, lambdaL1, lambdaL2, cfg)
%FIT_BATCH 同じ設計行列を使う多数のShotを一括で回帰する。
%   A           : 標準化済み設計行列（マーク数×項数、切片列を含む）
%   Y           : 測定値 [nm]（マーク数×Shot数）
%   method      : 回帰手法の構造体（run_all の build_methods で作る）
%                   loss    : 'OLS' | 'LAD' | 'Huber'
%                   penalty : 'none' | 'L2' | 'L1'
%   isPenalized : 罰則をかける項なら true（項数×1）。線形成分にもかける場合はすべて true
%   theta       : 回帰係数（項数×Shot数）
%
%   目的関数（n: マーク数、r_i: 残差、s: 残差のMADから求めた尺度、c: Huberの定数）
%     損失  OLS   : (1/2n) Σ r_i^2
%           LAD   : (1/n)  Σ s |r_i|
%           Huber : (1/n)  Σ s^2 ρ_c(r_i / s)
%     罰則  L2    : (λ2/2) Σ_{罰則項} θ_j^2
%           L1    :  λ1    Σ_{罰則項} |θ_j|
%   LADとHuberは、重み付き最小二乗 (1/2n) Σ w_i r_i^2 を、重みを更新しながら解き直して（IRLS）最小化する。
%   LADの損失に s を掛けるのは、残差が s 程度の点の重みを1程度にして、罰則の効き方をOLS・Huberとそろえるため。
if ~ismember(method.loss, {'OLS', 'LAD', 'Huber'}) || ~ismember(method.penalty, {'none', 'L2', 'L1'})
    error('fit_batch:unknownMethod', '未対応の回帰手法です: %s', method.name);
end
[markCount, shotCount] = size(Y);
if size(A, 1) ~= markCount
    error('fit_batch:sizeMismatch', '設計行列と測定値のマーク数が一致しません。');
end
if numel(isPenalized) ~= size(A, 2)
    error('fit_batch:maskMismatch', '罰則の対象を表すベクトルの長さが項数と一致しません。');
end

weights = ones(markCount, shotCount);
theta = solve_weighted(A, Y, weights, method.penalty, isPenalized, lambdaL1, lambdaL2, [], cfg);
if strcmp(method.loss, 'OLS')
    return;
end

% LAD・Huber: 残差の大きいマークほど重みを下げて解き直す（IRLS）。
% 尺度 s は最初の irlsScaleUpdateIter 回だけ更新し、その後は固定する。更新を続けるとMADの値が
% 隣り合う残差の間で行き来して収束しないShotが出るため。s を固定した後は凸な問題になり収束に向かう。
% 係数の変化が irlsTolNm 未満になったShotから反復を止める。
scale = zeros(1, shotCount);
isActive = true(1, shotCount);
for iteration = 1:cfg.irlsMaxIter
    active = find(isActive);
    residual = Y(:, active) - A * theta(:, active);
    if iteration <= cfg.irlsScaleUpdateIter
        scale(active) = max(cfg.madToSigma * median(abs(residual - median(residual, 1)), 1), eps);
    end
    weights = robust_weights(residual, scale(active), method.loss, cfg);
    thetaNew = solve_weighted(A, Y(:, active), weights, method.penalty, isPenalized, lambdaL1, lambdaL2, ...
        theta(:, active), cfg);
    change = max(abs(thetaNew - theta(:, active)), [], 1);
    theta(:, active) = thetaNew;
    if iteration > cfg.irlsScaleUpdateIter
        isActive(active(change < cfg.irlsTolNm)) = false;
    end
    if ~any(isActive)
        return;
    end
end
warning('fit_batch:irlsNotConverged', 'IRLS（%s）が %d 回で収束しなかったShotが %d 個あります（最大変化 %.2e nm）。', ...
    method.name, cfg.irlsMaxIter, sum(isActive), max(change));
end

function weights = robust_weights(residual, scale, loss, cfg)
% IRLSの重み。scale は Shotごとの尺度 s（1×Shot数）。
absResidual = abs(residual);
switch loss
    case 'LAD'
        % w = s/|r|。残差がほぼ0の点で重みが発散しないように |r| の下限を設ける
        weights = scale ./ max(absResidual, cfg.ladMinResidualRatio * scale);
    case 'Huber'
        % w = min(1, c s/|r|)
        weights = min(1, cfg.huberTuning * scale ./ max(absResidual, eps));
end
end

function theta = solve_weighted(A, Y, weights, penaltyType, isPenalized, lambdaL1, lambdaL2, thetaStart, cfg)
% 重み付き最小二乗（＋罰則）を解く。thetaStart はL1の座標降下法の初期値（空なら最小二乗解から始める）。
termCount = size(A, 2);
switch penaltyType
    case 'L2'
        % L2の罰則は「罰則項の係数を0に近づける観測行」を設計行列の下に加えて表す
        penaltyRows = sqrt(size(A, 1) * lambdaL2) * diag(double(isPenalized(:)));
        theta = weighted_least_squares(A, Y, weights, penaltyRows);
    case 'L1'
        if isempty(thetaStart)
            thetaStart = weighted_least_squares(A, Y, weights, zeros(0, termCount));
        end
        theta = weighted_lasso(A, Y, weights, isPenalized, lambdaL1, thetaStart, cfg);
    otherwise
        theta = weighted_least_squares(A, Y, weights, zeros(0, termCount));
end
end

function theta = weighted_least_squares(A, Y, weights, penaltyRows)
% 数値的に安定なQR分解（バックスラッシュ）でShotごとに解く
shotCount = size(Y, 2);
theta = zeros(size(A, 2), shotCount);
sqrtWeights = sqrt(weights);
zeroRows = zeros(size(penaltyRows, 1), 1);
for s = 1:shotCount
    theta(:, s) = [sqrtWeights(:, s) .* A; penaltyRows] \ [sqrtWeights(:, s) .* Y(:, s); zeroRows];
end
end

function thetaOut = weighted_lasso(A, Y, weights, isPenalized, lambda, theta, cfg)
% 座標降下法（グラム行列形式）。複数Shotをまとめて更新し、収束したShotから計算を外す。
% 収束の遅い一部のShotのために、全Shotの反復が続くのを避けるため。
[markCount, termCount] = size(A);
shotCount = size(Y, 2);
gram = zeros(termCount, termCount, shotCount);
rhs = zeros(termCount, shotCount);
for s = 1:shotCount
    weightedA = weights(:, s) .* A;
    gram(:, :, s) = (A' * weightedA) / markCount;
    rhs(:, s) = (weightedA' * Y(:, s)) / markCount;
end
diagonal = zeros(termCount, shotCount);
for j = 1:termCount
    diagonal(j, :) = reshape(gram(j, j, :), 1, shotCount);
end
gradientTerm = zeros(termCount, shotCount);   % G * theta
for s = 1:shotCount
    gradientTerm(:, s) = gram(:, :, s) * theta(:, s);
end

thetaOut = theta;
shotIndex = 1:shotCount;   % 作業中の配列の各列が、元の何番目のShotか
for sweep = 1:cfg.cdMaxSweep
    shotChange = zeros(1, numel(shotIndex));
    for j = 1:termCount
        rho = rhs(j, :) - gradientTerm(j, :) + diagonal(j, :) .* theta(j, :);
        if isPenalized(j)
            updated = sign(rho) .* max(abs(rho) - lambda, 0) ./ diagonal(j, :);   % ソフトしきい値処理
        else
            updated = rho ./ diagonal(j, :);
        end
        delta = updated - theta(j, :);
        if any(delta)
            theta(j, :) = updated;
            gradientTerm = gradientTerm + reshape(gram(:, j, :), termCount, numel(shotIndex)) .* delta;
            shotChange = max(shotChange, abs(delta) .* sqrt(diagonal(j, :)));
        end
    end
    isDone = shotChange < cfg.cdTolNm;
    if any(isDone)
        % 収束したShotの係数を確定し、以後の計算から外す
        thetaOut(:, shotIndex(isDone)) = theta(:, isDone);
        shotIndex    = shotIndex(~isDone);
        theta        = theta(:, ~isDone);
        gram         = gram(:, :, ~isDone);
        rhs          = rhs(:, ~isDone);
        diagonal     = diagonal(:, ~isDone);
        gradientTerm = gradientTerm(:, ~isDone);
        if isempty(shotIndex)
            return;
        end
    end
end
thetaOut(:, shotIndex) = theta;
warning('fit_batch:notConverged', '座標降下法が %d 回で収束しなかったShotが %d 個あります（最大変化 %.2e）。', ...
    cfg.cdMaxSweep, numel(shotIndex), max(shotChange(~isDone)));
end

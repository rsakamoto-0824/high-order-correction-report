function theta = fit_batch(A, Y, method, isPenalized, lambdaL1, lambdaL2, cfg)
%FIT_BATCH 同じ設計行列を使う多数のShotを一括で回帰する。
%   A           : 標準化済み設計行列（マーク数×項数、切片列を含む）
%   Y           : 測定値（マーク数×Shot数）
%   method      : 'OLS' | 'Huber' | 'L2' | 'L1' | 'Huber+L2' | 'Huber+L1'
%   isPenalized : 罰則をかける項なら true（項数×1の論理ベクトル）
%   theta       : 回帰係数（項数×Shot数）
%
%   目的関数（n: マーク数、w_i: IRLSの重み、非ロバストでは w_i = 1）
%     OLS/Huber : (1/2n) Σ w_i r_i^2
%     L2        : (1/2n) Σ w_i r_i^2 + (λ2/2) Σ_{罰則項} θ_j^2
%     L1        : (1/2n) Σ w_i r_i^2 +  λ1    Σ_{罰則項} |θ_j|
[penaltyType, isRobust] = parse_method(method);
[markCount, shotCount] = size(Y);
if size(A, 1) ~= markCount
    error('fit_batch:sizeMismatch', '設計行列と測定値のマーク数が一致しません。');
end

weights = ones(markCount, shotCount);
theta = solve_weighted(A, Y, weights, penaltyType, isPenalized, lambdaL1, lambdaL2, cfg);
if ~isRobust
    return;
end

% Huber M推定: 残差の大きいマークほど重みを下げて解き直す（IRLS）
for iteration = 1:cfg.irlsMaxIter
    residual = Y - A * theta;
    scale = cfg.madToSigma * median(abs(residual - median(residual, 1)), 1);
    scale = max(scale, eps);
    standardizedResidual = abs(residual) ./ (cfg.huberTuning * scale);
    weights = min(1, 1 ./ max(standardizedResidual, eps));
    thetaNew = solve_weighted(A, Y, weights, penaltyType, isPenalized, lambdaL1, lambdaL2, cfg);
    change = max(abs(thetaNew - theta), [], 'all');
    theta = thetaNew;
    if change < cfg.irlsTolNm
        break;
    end
end
end

function [penaltyType, isRobust] = parse_method(method)
switch method
    case 'OLS',      penaltyType = 'none'; isRobust = false;
    case 'Huber',    penaltyType = 'none'; isRobust = true;
    case 'L2',       penaltyType = 'L2';   isRobust = false;
    case 'L1',       penaltyType = 'L1';   isRobust = false;
    case 'Huber+L2', penaltyType = 'L2';   isRobust = true;
    case 'Huber+L1', penaltyType = 'L1';   isRobust = true;
    otherwise
        error('fit_batch:unknownMethod', '未対応の回帰手法です: %s', method);
end
end

function theta = solve_weighted(A, Y, weights, penaltyType, isPenalized, lambdaL1, lambdaL2, cfg)
[markCount, termCount] = size(A);
shotCount = size(Y, 2);
theta = zeros(termCount, shotCount);
sqrtWeights = sqrt(weights);

% 最小二乗とL2は、数値的に安定なQR分解（バックスラッシュ）でShotごとに解く
penaltyRows = sqrt(markCount * lambdaL2) * diag(double(isPenalized(:)));
for s = 1:shotCount
    weightedA = sqrtWeights(:, s) .* A;
    weightedY = sqrtWeights(:, s) .* Y(:, s);
    if strcmp(penaltyType, 'L2')
        theta(:, s) = [weightedA; penaltyRows] \ [weightedY; zeros(termCount, 1)];
    else
        theta(:, s) = weightedA \ weightedY;
    end
end
if ~strcmp(penaltyType, 'L1')
    return;
end

% L1は最小二乗解を初期値にして座標降下法で解く
gram = zeros(termCount, termCount, shotCount);
rhs = zeros(termCount, shotCount);
for s = 1:shotCount
    weightedA = weights(:, s) .* A;
    gram(:, :, s) = (A' * weightedA) / markCount;
    rhs(:, s) = (weightedA' * Y(:, s)) / markCount;
end
theta = coordinate_descent_l1(gram, rhs, isPenalized, lambdaL1, theta, cfg);
end

function theta = coordinate_descent_l1(gram, rhs, isPenalized, lambda, theta, cfg)
% 全Shotを同時に更新する座標降下法（グラム行列形式）
[termCount, ~, shotCount] = size(gram);
diagonal = zeros(termCount, shotCount);
for j = 1:termCount
    diagonal(j, :) = reshape(gram(j, j, :), 1, shotCount);
end
gradientTerm = zeros(termCount, shotCount);   % G * theta
for s = 1:shotCount
    gradientTerm(:, s) = gram(:, :, s) * theta(:, s);
end

for sweep = 1:cfg.cdMaxSweep
    maxChange = 0;
    for j = 1:termCount
        rho = rhs(j, :) - gradientTerm(j, :) + diagonal(j, :) .* theta(j, :);
        if isPenalized(j)
            updated = sign(rho) .* max(abs(rho) - lambda, 0) ./ diagonal(j, :);
        else
            updated = rho ./ diagonal(j, :);
        end
        delta = updated - theta(j, :);
        if any(delta)
            theta(j, :) = updated;
            gradientTerm = gradientTerm + reshape(gram(:, j, :), termCount, shotCount) .* delta;
            maxChange = max(maxChange, max(abs(delta) .* sqrt(diagonal(j, :))));
        end
    end
    if maxChange < cfg.cdTolNm
        return;
    end
end
warning('fit_batch:notConverged', '座標降下法が %d 回で収束しませんでした（最大変化 %.2e）。', ...
    cfg.cdMaxSweep, maxChange);
end

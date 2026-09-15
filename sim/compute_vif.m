function vif = compute_vif(A, exponents)
%COMPUTE_VIF 各項の分散拡大係数 VIF_j = 1 / (1 - R_j^2) を計算する。
%   R_j^2 は、項 j を他の項（切片を含む）で回帰したときの決定係数。
%   切片（指数 [0 0]）はVIFの対象外なので NaN を返す。
%   VIFはマーク配置と項の組だけで決まり、測定値には依存しない。
isIntercept = all(exponents == 0, 2);
% mm座標のべき乗は列ごとに桁が大きく違うため、各列を最大値で割って数値計算を安定させる。
% VIFは列を定数倍しても変わらない。
A = A ./ max(abs(A), [], 1);
termCount = size(A, 2);
vif = nan(termCount, 1);
for j = find(~isIntercept)'
    others = setdiff(1:termCount, j);
    target = A(:, j);
    fitted = A(:, others) * (A(:, others) \ target);
    totalSumSquares = sum((target - mean(target)).^2);
    residualSumSquares = sum((target - fitted).^2);
    rSquared = 1 - residualSumSquares / totalSumSquares;
    vif(j) = 1 / max(1 - rSquared, eps);
end
end

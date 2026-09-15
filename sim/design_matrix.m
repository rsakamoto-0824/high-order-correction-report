function A = design_matrix(x, y, exponents)
%DESIGN_MATRIX 各点の座標から、項 x^m y^n を並べた設計行列を作る。
%   x, y      : Shot中心を原点とする座標 [mm] の列ベクトル
%   exponents : 項の指数 [m n]（項数×2）
%   A         : 点数×項数
if numel(x) ~= numel(y)
    error('design_matrix:sizeMismatch', 'x と y の点数が一致しません。');
end
A = zeros(numel(x), size(exponents, 1));
for j = 1:size(exponents, 1)
    A(:, j) = x(:) .^ exponents(j, 1) .* y(:) .^ exponents(j, 2);
end
end

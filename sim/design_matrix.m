function A = design_matrix(xNorm, yNorm, exponents)
%DESIGN_MATRIX 各点の座標から、項 x^m y^n を並べた設計行列を作る。
%   xNorm, yNorm : 正規化座標（Shot端で±1）の列ベクトル
%   exponents    : 項の指数 [m n]（項数×2）
%   A            : 点数×項数
if numel(xNorm) ~= numel(yNorm)
    error('design_matrix:sizeMismatch', 'x と y の点数が一致しません。');
end
A = zeros(numel(xNorm), size(exponents, 1));
for j = 1:size(exponents, 1)
    A(:, j) = xNorm(:) .^ exponents(j, 1) .* yNorm(:) .^ exponents(j, 2);
end
end

function labels = monomial_label(exponents)
%MONOMIAL_LABEL 指数 [m n] から x^m y^n の表示文字列（TeX形式）を作る。
%   例: [0 0] -> '1', [2 1] -> 'x^{2}y'
labels = cell(size(exponents, 1), 1);
for i = 1:size(exponents, 1)
    labels{i} = [power_text('x', exponents(i, 1)) power_text('y', exponents(i, 2))];
    if isempty(labels{i})
        labels{i} = '1';
    end
end
end

function text = power_text(symbol, power)
if power == 0
    text = '';
elseif power == 1
    text = symbol;
else
    text = sprintf('%s^{%d}', symbol, power);
end
end

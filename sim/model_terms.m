function model = model_terms(modelName)
%MODEL_TERMS Shot高次補正演算式の項を返す。
%   各軸（x: dx の式, y: dy の式）について、項 x^m y^n の指数 [m n] と名前を持つ。
%   'ASML'  : iHOPC k1〜k20（Huang et al., Proc. SPIE 7272 (2009) の式）
%   'Nikon' : Shot線形6項 ＋ Shot高次係数21項。
%             高次係数名 SHOTFACmn を x^m y^n の係数と解釈している（装置仕様で要確認）。
%   'Truth' : 真の歪みを作るための項（ASMLとNikonの項の和集合）

switch modelName
    case 'ASML'
        model.x.exp  = [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 3 0; 2 1; 1 2; 0 3];
        model.x.name = {'k1', 'k3', 'k5', 'k7', 'k9', 'k11', 'k13', 'k15', 'k17', 'k19'};
        model.y.exp  = [0 0; 0 1; 1 0; 0 2; 1 1; 2 0; 0 3; 1 2; 2 1; 3 0];
        model.y.name = {'k2', 'k4', 'k6', 'k8', 'k10', 'k12', 'k14', 'k16', 'k18', 'k20'};
    case 'Nikon'
        model.x.exp  = [0 0; 1 0; 0 1; 0 2; 0 3; 0 4; 0 5; 0 6; 2 0; 3 0; 1 1; 1 2];
        model.y.exp  = [0 0; 0 1; 1 0; 1 1; 0 2; 1 2; 0 3; 1 3; 0 4; 1 4; 0 5; 1 5; 0 6; 2 0; 2 1];
        model.x.name = nikon_names(model.x.exp, 'X');
        model.y.name = nikon_names(model.y.exp, 'Y');
    case 'Truth'
        asml  = model_terms('ASML');
        nikon = model_terms('Nikon');
        model.x.exp  = unique([asml.x.exp; nikon.x.exp], 'rows', 'stable');
        model.y.exp  = unique([asml.y.exp; nikon.y.exp], 'rows', 'stable');
        model.x.name = monomial_label(model.x.exp)';
        model.y.name = monomial_label(model.y.exp)';
    otherwise
        error('model_terms:unknownModel', '未対応のモデル名です: %s', modelName);
end
model.label = modelName;
end

function names = nikon_names(exponents, axisSuffix)
% 1次以下は Shot線形成分（OFS: 並進、SCL: 倍率、ROT: 回転）、2次以上は SHOTFACmn と呼ぶ
names = cell(1, size(exponents, 1));
for i = 1:size(exponents, 1)
    m = exponents(i, 1);
    n = exponents(i, 2);
    isScaling = (axisSuffix == 'X' && m == 1 && n == 0) || (axisSuffix == 'Y' && m == 0 && n == 1);
    if m + n == 0
        names{i} = ['OFS' axisSuffix];
    elseif m + n == 1 && isScaling
        names{i} = ['SCL' axisSuffix];
    elseif m + n == 1
        names{i} = ['ROT' axisSuffix];
    else
        names{i} = sprintf('SHOTFAC%d%d%s', m, n, axisSuffix);
    end
end
end

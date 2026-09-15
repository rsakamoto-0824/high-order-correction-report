function model = model_terms(modelName)
%MODEL_TERMS Shot高次補正演算式の項を返す。
%   各軸（x: dx の式, y: dy の式）について、項 x^m y^n の指数 [m n] と係数番号を持つ。
%   x, y はShot中心を原点とするマーク座標 [mm]。
%   'ASML', 'Nikon' : 露光機メーカーごとのShot高次補正演算式（前提条件で決められた式）。
%                     項・係数番号・変数・符号を変更しない。メーカー間で式を共通化したり、
%                     記載のない項を補ったりもしない。
%   'Truth'         : 真の歪みを作るための項（ASML式とニコン式の項の和集合）

switch modelName
    case 'ASML'
        % dx = k1 + k3*x + k5*y + k7*x^2 + k11*y^2 + k13*x^3 + k19*y^3
        model.x.exp  = [0 0; 1 0; 0 1; 2 0; 0 2; 3 0; 0 3];
        model.x.name = {'k1', 'k3', 'k5', 'k7', 'k11', 'k13', 'k19'};
        % dy = k2 + k4*y + k6*x + k8*y^2 + k10*x*y + k12*x^2 + k14*y^3 + k16*x*y^2
        model.y.exp  = [0 0; 0 1; 1 0; 0 2; 1 1; 2 0; 0 3; 1 2];
        model.y.name = {'k2', 'k4', 'k6', 'k8', 'k10', 'k12', 'k14', 'k16'};
    case 'Nikon'
        % dx = k1 + k3*x + k5*y + k7*x^2 + k9*x*y + k11*y^2 + k13*x^3 + k17*x*y^2 + k19*y^3
        %      + k29*y^4 + k41*y^5 + k55*y^6
        model.x.exp  = [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 3 0; 1 2; 0 3; 0 4; 0 5; 0 6];
        model.x.name = {'k1', 'k3', 'k5', 'k7', 'k9', 'k11', 'k13', 'k17', 'k19', 'k29', 'k41', 'k55'};
        % dy = k2 + k4*y + k6*x + k8*y^2 + k10*x*y + k12*x^2 + k14*y^3 + k16*x*y^2 + k18*x^2*y
        %      + k22*y^4 + k24*x*y^3 + k32*y^5 + k34*x*y^4 + k44*y^6 + k46*x*y^5
        model.y.exp  = [0 0; 0 1; 1 0; 0 2; 1 1; 2 0; 0 3; 1 2; 2 1; 0 4; 1 3; 0 5; 1 4; 0 6; 1 5];
        model.y.name = {'k2', 'k4', 'k6', 'k8', 'k10', 'k12', 'k14', 'k16', 'k18', ...
            'k22', 'k24', 'k32', 'k34', 'k44', 'k46'};
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
if numel(model.x.name) ~= size(model.x.exp, 1) || numel(model.y.name) ~= size(model.y.exp, 1)
    error('model_terms:nameMismatch', '%s の項の数と係数名の数が一致しません。', modelName);
end
model.label = modelName;
end

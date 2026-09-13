function data = generate_shots(cfg, truthModel, markX, markY, evalX, evalY)
%GENERATE_SHOTS 高次歪みを持つShotを乱数で生成し、測定値を作る。
%   座標は正規化済み（Shot端で±1）。各軸（x, y）について次を返す（点数×Shot数）。
%     trueAtMarks     : マーク位置での真の歪み [nm]
%     trueAtEval      : 評価グリッドでの真の歪み [nm]
%     measuredClean   : 真値 ＋ 偶然誤差
%     measuredOutlier : 真値 ＋ 偶然誤差 ＋ 異常値
%   2つの測定値は同じ偶然誤差を共有するので、異常値の有無だけを比べられる。
shotCount = cfg.shotCount;
markCount = numel(markX);

% 異常値はマーク単位で発生し、そのマークのX・Y両方の測定値に乗る
isOutlierMark = rand(markCount, shotCount) < cfg.outlierProbability;
data.isOutlierMark = isOutlierMark;

axisNames = {'x', 'y'};
for a = 1:numel(axisNames)
    axisName = axisNames{a};
    exponents = truthModel.(axisName).exp;
    order = sum(exponents, 2);
    sigma = cfg.truthSigmaNmByOrder(order + 1);
    activeProbability = cfg.truthActiveProbByOrder(order + 1);
    termCount = size(exponents, 1);

    isActive = rand(termCount, shotCount) < activeProbability(:);
    coefficient = randn(termCount, shotCount) .* sigma(:) .* isActive;

    trueAtMarks = design_matrix(markX, markY, exponents) * coefficient;
    trueAtEval  = design_matrix(evalX, evalY, exponents) * coefficient;

    noise = cfg.noiseSigmaNm * randn(markCount, shotCount);
    outlierSize = cfg.outlierMinNm + (cfg.outlierMaxNm - cfg.outlierMinNm) * rand(markCount, shotCount);
    outlierSign = 2 * (rand(markCount, shotCount) < 0.5) - 1;
    outlier = isOutlierMark .* outlierSign .* outlierSize;

    data.(axisName).coefficient     = coefficient;
    data.(axisName).trueAtMarks     = trueAtMarks;
    data.(axisName).trueAtEval      = trueAtEval;
    data.(axisName).measuredClean   = trueAtMarks + noise;
    data.(axisName).measuredOutlier = trueAtMarks + noise + outlier;
end
end

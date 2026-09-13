function run_all(shotCountOverride)
%RUN_ALL Shot高次補正の回帰シミュレーションを実行し、結果を results/ に保存する。
%   run_all      : simulation_config の条件（1000 Shot）で実行する
%   run_all(50)  : Shot数を減らして動作確認する
%
%   手順
%     1. ASML・Nikonの補正式と、VIFで項を削減した補正式を作る
%     2. 高次歪みを持つShotを乱数で生成し、測定誤差（偶然誤差・異常値）を付ける
%     3. 6種類の回帰（OLS / Huber / L2 / L1 / Huber+L2 / Huber+L1）で補正値を求める
%     4. Shot全面の密グリッドで「真の歪み − 補正値」を評価する

cfg = simulation_config();
if nargin >= 1
    if ~isscalar(shotCountOverride) || shotCountOverride < 1
        error('run_all:invalidShotCount', 'Shot数は1以上の整数で指定してください。');
    end
    cfg.shotCount = round(shotCountOverride);
end
make_folder(cfg.resultDir);
make_folder(cfg.figureDir);
rng(cfg.randomSeed, 'twister');
totalTimer = tic;

%% 1. 座標と補正式
geometry = build_geometry(cfg);
[variants, vifTable] = build_model_variants(cfg, geometry);
writetable(vifTable, fullfile(cfg.resultDir, 'vif_table.csv'));

%% 2. Shot生成
truthModel = model_terms('Truth');
data = generate_shots(cfg, truthModel, geometry.markX, geometry.markY, geometry.evalX, geometry.evalY);

%% 3-4. 回帰と評価
methodNames   = {'OLS', 'Huber', 'L2', 'L1', 'Huber+L2', 'Huber+L1'};
scenarioNames = {'outlier', 'clean'};
measuredField = {'measuredOutlier', 'measuredClean'};
axisNames     = {'x', 'y'};

rmsTrue     = cell(numel(scenarioNames), numel(variants), numel(methodNames), numel(axisNames));
rmsApparent = cell(size(rmsTrue));
summaryRows = {};
for sc = 1:numel(scenarioNames)
    for v = 1:numel(variants)
        for a = 1:numel(axisNames)
            axisName = axisNames{a};
            design = variants(v).(axisName);
            measured = data.(axisName).(measuredField{sc});
            for m = 1:numel(methodNames)
                theta = fit_batch(design.Astd, measured, methodNames{m}, design.isPenalized, ...
                    cfg.lambdaL1, cfg.lambdaL2, cfg);
                trueResidual = data.(axisName).trueAtEval - design.Aeval * theta;
                rmsTrue{sc, v, m, a} = sqrt(mean(trueResidual.^2, 1));
                rmsApparent{sc, v, m, a} = sqrt(mean((measured - design.Astd * theta).^2, 1));
                summaryRows(end + 1, :) = {scenarioNames{sc}, variants(v).name, methodNames{m}, ...
                    upper(axisName), size(design.exponents, 1), ...
                    median(rmsTrue{sc, v, m, a}), prctile(rmsTrue{sc, v, m, a}, 95), ...
                    abs(mean(trueResidual(:))) + 3 * std(trueResidual(:)), ...
                    median(rmsApparent{sc, v, m, a})}; %#ok<AGROW>
                fprintf('[%6.1fs] %-7s %-9s %-8s d%s  median RMS = %.3f nm\n', toc(totalTimer), ...
                    scenarioNames{sc}, variants(v).name, methodNames{m}, axisName, median(rmsTrue{sc, v, m, a}));
            end
        end
    end
end
summaryTable = cell2table(summaryRows, 'VariableNames', {'Scenario', 'Model', 'Method', 'Axis', ...
    'TermCount', 'MedianRmsNm', 'P95RmsNm', 'MeanPlus3SigmaNm', 'MedianApparentRmsNm'});
writetable(summaryTable, fullfile(cfg.resultDir, 'summary.csv'));

%% 参考値: 補正なし、および各補正式で表現できる限界（雑音なし・全面データで当てはめた残差）
[referenceTable, rawRms, floorRms] = compute_reference(data, variants, axisNames);
writetable(referenceTable, fullfile(cfg.resultDir, 'reference.csv'));

%% 対になった比較（同じShotどうしで比べる）
comparisonTable = compare_conditions(rmsTrue, variants, methodNames, scenarioNames);
writetable(comparisonTable, fullfile(cfg.resultDir, 'comparisons.csv'));

%% λの感度解析（異常値ありのシナリオ、ロバスト＋正則化のみ）
sweepTable = sweep_lambda(cfg, data, variants, axisNames);
writetable(sweepTable, fullfile(cfg.resultDir, 'lambda_sweep.csv'));

%% 保存
save(fullfile(cfg.resultDir, 'simulation_results.mat'), 'cfg', 'geometry', 'variants', 'data', ...
    'methodNames', 'scenarioNames', 'axisNames', 'rmsTrue', 'rmsApparent', 'rawRms', 'floorRms', ...
    'summaryTable', 'referenceTable', 'comparisonTable', 'sweepTable', 'vifTable');
fprintf('完了: %.1f 秒。結果は %s に保存しました。\n', toc(totalTimer), cfg.resultDir);
end

% ------------------------------------------------------------------------
function make_folder(folderPath)
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
end

function geometry = build_geometry(cfg)
% マーク位置と評価グリッドを正規化座標（Shot端で±1）で作る
[markXgrid, markYgrid] = meshgrid(cfg.markXmm, cfg.markYmm);
geometry.markX = markXgrid(:) / cfg.shotHalfWidthMm;
geometry.markY = markYgrid(:) / cfg.shotHalfHeightMm;
evalAxis = linspace(-1, 1, cfg.evalGridCount);
[evalXgrid, evalYgrid] = meshgrid(evalAxis, evalAxis);
geometry.evalX = evalXgrid(:);
geometry.evalY = evalYgrid(:);
end

function [variants, vifTable] = build_model_variants(cfg, geometry)
% ASML・Nikonの補正式と、それぞれのVIF削減版（計4種類）を作る
baseNames = {'ASML', 'Nikon'};
axisNames = {'x', 'y'};
variants = struct('name', {}, 'x', {}, 'y', {});
vifRows = {};
for b = 1:numel(baseNames)
    base = model_terms(baseNames{b});
    fullDesign = struct();
    reducedDesign = struct();
    for a = 1:numel(axisNames)
        axisName = axisNames{a};
        exponents = base.(axisName).exp;
        names = base.(axisName).name;
        Araw = design_matrix(geometry.markX, geometry.markY, exponents);
        vifFull = compute_vif(Araw, exponents);
        [keepMask, removalLog] = reduce_terms_by_vif(Araw, exponents, cfg.vifThreshold);
        vifAfter = nan(size(vifFull));
        vifAfter(keepMask) = compute_vif(Araw(:, keepMask), exponents(keepMask, :));
        removalOrder = zeros(size(exponents, 1), 1);
        removalOrder(removalLog(:, 1)) = 1:size(removalLog, 1);
        labels = monomial_label(exponents);
        for t = 1:size(exponents, 1)
            vifRows(end + 1, :) = {baseNames{b}, upper(axisName), names{t}, labels{t}, ...
                sum(exponents(t, :)), vifFull(t), keepMask(t), removalOrder(t), vifAfter(t)}; %#ok<AGROW>
        end
        fullDesign.(axisName) = make_design(exponents, names, geometry);
        reducedDesign.(axisName) = make_design(exponents(keepMask, :), names(keepMask), geometry);
    end
    variants(end + 1) = struct('name', baseNames{b}, 'x', fullDesign.x, 'y', fullDesign.y); %#ok<AGROW>
    variants(end + 1) = struct('name', [baseNames{b} '-VIF'], 'x', reducedDesign.x, 'y', reducedDesign.y); %#ok<AGROW>
end
vifTable = cell2table(vifRows, 'VariableNames', ...
    {'Model', 'Axis', 'Term', 'Monomial', 'Order', 'VIF', 'Kept', 'RemovalOrder', 'VIFAfter'});
end

function design = make_design(exponents, names, geometry)
% 回帰用に列を標準化する（切片以外を平均0・標準偏差1）。評価グリッドにも同じ変換をかける。
Araw = design_matrix(geometry.markX, geometry.markY, exponents);
isIntercept = all(exponents == 0, 2)';
center = mean(Araw, 1);
scale = std(Araw, 1, 1);
center(isIntercept) = 0;
scale(isIntercept) = 1;
design.exponents   = exponents;
design.names       = names;
design.Astd        = (Araw - center) ./ scale;
design.AevalRaw    = design_matrix(geometry.evalX, geometry.evalY, exponents);
design.Aeval       = (design.AevalRaw - center) ./ scale;
design.isPenalized = sum(exponents, 2) >= 2;   % 並進・倍率・回転には罰則をかけない
end

function [referenceTable, rawRms, floorRms] = compute_reference(data, variants, axisNames)
rows = {};
rawRms = cell(1, numel(axisNames));
floorRms = cell(numel(variants), numel(axisNames));
for a = 1:numel(axisNames)
    axisName = axisNames{a};
    trueAtEval = data.(axisName).trueAtEval;
    rawRms{a} = sqrt(mean(trueAtEval.^2, 1));
    rows(end + 1, :) = {'NoCorrection', upper(axisName), median(rawRms{a}), prctile(rawRms{a}, 95)}; %#ok<AGROW>
    for v = 1:numel(variants)
        AevalRaw = variants(v).(axisName).AevalRaw;
        bestResidual = trueAtEval - AevalRaw * (AevalRaw \ trueAtEval);
        floorRms{v, a} = sqrt(mean(bestResidual.^2, 1));
        rows(end + 1, :) = {['Floor_' variants(v).name], upper(axisName), ...
            median(floorRms{v, a}), prctile(floorRms{v, a}, 95)}; %#ok<AGROW>
    end
end
referenceTable = cell2table(rows, 'VariableNames', {'Condition', 'Axis', 'MedianRmsNm', 'P95RmsNm'});
end

function comparisonTable = compare_conditions(rmsTrue, variants, methodNames, scenarioNames)
% 条件Aと条件Bを同じShotで比べる。X・Yの残差RMSは二乗平均で1つにまとめる。
% 列: 比較の種類, 条件A, 条件B, Aの中央値, Bの中央値, 中央値の低減率, AがBより良いShotの割合
variantNames = {variants.name};
pairs = {
    'C1_VIFbeforeL1',  'ASML-VIF',  'L1',       'ASML',      'L1'
    'C1_VIFbeforeL1',  'ASML-VIF',  'Huber+L1', 'ASML',      'Huber+L1'
    'C1_VIFbeforeL1',  'Nikon-VIF', 'L1',       'Nikon',     'L1'
    'C1_VIFbeforeL1',  'Nikon-VIF', 'Huber+L1', 'Nikon',     'Huber+L1'
    'C2_L2vsL1_VIF',   'ASML-VIF',  'L2',       'ASML-VIF',  'L1'
    'C2_L2vsL1_VIF',   'ASML-VIF',  'Huber+L2', 'ASML-VIF',  'Huber+L1'
    'C2_L2vsL1_VIF',   'Nikon-VIF', 'L2',       'Nikon-VIF', 'L1'
    'C2_L2vsL1_VIF',   'Nikon-VIF', 'Huber+L2', 'Nikon-VIF', 'Huber+L1'
    'C2_L2vsL1_Full',  'ASML',      'Huber+L2', 'ASML',      'Huber+L1'
    'C2_L2vsL1_Full',  'Nikon',     'Huber+L2', 'Nikon',     'Huber+L1'
    'C3_Robust',       'ASML',      'Huber',    'ASML',      'OLS'
    'C3_Robust',       'ASML-VIF',  'Huber',    'ASML-VIF',  'OLS'
    'C3_Robust',       'Nikon',     'Huber',    'Nikon',     'OLS'
    'C3_Robust',       'Nikon-VIF', 'Huber',    'Nikon-VIF', 'OLS'
    'C3_Robust',       'ASML-VIF',  'Huber+L2', 'ASML-VIF',  'L2'
    'C3_Robust',       'Nikon-VIF', 'Huber+L2', 'Nikon-VIF', 'L2'
    'C4_VIFonly',      'ASML-VIF',  'OLS',      'ASML',      'OLS'
    'C4_VIFonly',      'Nikon-VIF', 'OLS',      'Nikon',     'OLS'
    'C4_VIFonly',      'ASML-VIF',  'Huber',    'ASML',      'Huber'
    'C4_VIFonly',      'Nikon-VIF', 'Huber',    'Nikon',     'Huber'
    };
rows = {};
for sc = 1:numel(scenarioNames)
    for p = 1:size(pairs, 1)
        rmsA = combined_rms(rmsTrue, sc, find(strcmp(variantNames, pairs{p, 2})), ...
            find(strcmp(methodNames, pairs{p, 3})));
        rmsB = combined_rms(rmsTrue, sc, find(strcmp(variantNames, pairs{p, 4})), ...
            find(strcmp(methodNames, pairs{p, 5})));
        rows(end + 1, :) = {scenarioNames{sc}, pairs{p, 1}, ...
            [pairs{p, 2} ' ' pairs{p, 3}], [pairs{p, 4} ' ' pairs{p, 5}], ...
            median(rmsA), median(rmsB), 100 * (1 - median(rmsA) / median(rmsB)), ...
            100 * mean(rmsA < rmsB)}; %#ok<AGROW>
    end
end
comparisonTable = cell2table(rows, 'VariableNames', {'Scenario', 'Group', 'ConditionA', 'ConditionB', ...
    'MedianA_nm', 'MedianB_nm', 'MedianReductionPercent', 'WinRateA_percent'});
end

function rmsCombined = combined_rms(rmsTrue, scenarioIndex, variantIndex, methodIndex)
rmsX = rmsTrue{scenarioIndex, variantIndex, methodIndex, 1};
rmsY = rmsTrue{scenarioIndex, variantIndex, methodIndex, 2};
rmsCombined = sqrt((rmsX.^2 + rmsY.^2) / 2);
end

function sweepTable = sweep_lambda(cfg, data, variants, axisNames)
methodNames = {'Huber+L1', 'Huber+L2'};
rows = {};
for scaleValue = cfg.lambdaScaleList
    for v = 1:numel(variants)
        for m = 1:numel(methodNames)
            rmsByAxis = cell(1, numel(axisNames));
            for a = 1:numel(axisNames)
                axisName = axisNames{a};
                design = variants(v).(axisName);
                theta = fit_batch(design.Astd, data.(axisName).measuredOutlier, methodNames{m}, ...
                    design.isPenalized, cfg.lambdaL1 * scaleValue, cfg.lambdaL2 * scaleValue, cfg);
                residual = data.(axisName).trueAtEval - design.Aeval * theta;
                rmsByAxis{a} = sqrt(mean(residual.^2, 1));
            end
            rmsCombined = sqrt((rmsByAxis{1}.^2 + rmsByAxis{2}.^2) / 2);
            rows(end + 1, :) = {scaleValue, variants(v).name, methodNames{m}, ...
                median(rmsCombined), prctile(rmsCombined, 95)}; %#ok<AGROW>
            fprintf('  sweep x%-5.2f %-9s %-8s median = %.3f nm\n', scaleValue, variants(v).name, ...
                methodNames{m}, median(rmsCombined));
        end
    end
end
sweepTable = cell2table(rows, 'VariableNames', {'LambdaScale', 'Model', 'Method', 'MedianRmsNm', 'P95RmsNm'});
end

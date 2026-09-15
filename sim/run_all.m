function run_all(shotCountOverride)
%RUN_ALL Shot高次補正の回帰シミュレーションを実行し、結果を results/ に保存する。
%   run_all      : simulation_config の条件（1000 Shot）で実行する
%   run_all(50)  : Shot数を減らして動作確認する
%
%   手順
%     1. ASML式・ニコン式と、VIFで項を削減した補正式を作る（削除する項がなければ削減版は作らない）
%     2. 高次歪みを持つShotを乱数で生成し、測定誤差（偶然誤差・異常値）を付ける
%     3. 15種類の回帰で補正値を求める
%        損失関数（OLS / LAD / Huber）× 罰則（なし / L2 / L1）× 罰則の対象（高次項のみ / 全項）
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
methods       = build_methods();
methodNames   = {methods.name};
scenarioNames = {'outlier', 'clean'};
measuredField = {'measuredOutlier', 'measuredClean'};
axisNames     = {'x', 'y'};

rmsTrue     = cell(numel(scenarioNames), numel(variants), numel(methods), numel(axisNames));
rmsApparent = cell(size(rmsTrue));
summaryRows = {};
for sc = 1:numel(scenarioNames)
    for v = 1:numel(variants)
        for a = 1:numel(axisNames)
            axisName = axisNames{a};
            design = variants(v).(axisName);
            measured = data.(axisName).(measuredField{sc});
            for m = 1:numel(methods)
                isPenalized = penalty_mask(design, methods(m).scope);
                theta = fit_batch(design.Astd, measured, methods(m), isPenalized, cfg.lambdaL1, cfg.lambdaL2, cfg);
                trueResidual = data.(axisName).trueAtEval - design.Aeval * theta;
                rmsTrue{sc, v, m, a} = sqrt(mean(trueResidual.^2, 1));
                rmsApparent{sc, v, m, a} = sqrt(mean((measured - design.Astd * theta).^2, 1));
                summaryRows(end + 1, :) = {scenarioNames{sc}, variants(v).name, methods(m).name, ...
                    methods(m).loss, methods(m).penalty, methods(m).scope, upper(axisName), ...
                    size(design.exponents, 1), median(rmsTrue{sc, v, m, a}), prctile(rmsTrue{sc, v, m, a}, 95), ...
                    abs(mean(trueResidual(:))) + 3 * std(trueResidual(:)), ...
                    median(rmsApparent{sc, v, m, a})}; %#ok<AGROW>
                fprintf('[%7.1fs] %-7s %-9s %-13s d%s  median RMS = %.3f nm\n', toc(totalTimer), ...
                    scenarioNames{sc}, variants(v).name, methods(m).name, axisName, median(rmsTrue{sc, v, m, a}));
            end
        end
    end
end
summaryTable = cell2table(summaryRows, 'VariableNames', {'Scenario', 'Model', 'Method', 'Loss', 'Penalty', ...
    'Scope', 'Axis', 'TermCount', 'MedianRmsNm', 'P95RmsNm', 'MeanPlus3SigmaNm', 'MedianApparentRmsNm'});
writetable(summaryTable, fullfile(cfg.resultDir, 'summary.csv'));

%% dx と dy をまとめた要約（レポートの表で使う）
combinedTable = summarize_combined(rmsTrue, rmsApparent, variants, methods, scenarioNames);
writetable(combinedTable, fullfile(cfg.resultDir, 'summary_combined.csv'));

%% 参考値: 補正なし、および各補正式で表現できる限界（雑音なし・全面データで当てはめた残差）
[referenceTable, rawRms, floorRms] = compute_reference(data, variants, axisNames);
writetable(referenceTable, fullfile(cfg.resultDir, 'reference.csv'));

%% 対になった比較（同じShotどうしで比べる）
comparisonTable = compare_conditions(rmsTrue, variants, methods, scenarioNames);
writetable(comparisonTable, fullfile(cfg.resultDir, 'comparisons.csv'));

%% λの感度解析（異常値ありのシナリオ、Huber回帰＋罰則）
sweepTable = sweep_lambda(cfg, data, variants, axisNames);
writetable(sweepTable, fullfile(cfg.resultDir, 'lambda_sweep.csv'));

%% 保存
save(fullfile(cfg.resultDir, 'simulation_results.mat'), 'cfg', 'geometry', 'variants', 'data', ...
    'methods', 'methodNames', 'scenarioNames', 'axisNames', 'rmsTrue', 'rmsApparent', 'rawRms', 'floorRms', ...
    'summaryTable', 'combinedTable', 'referenceTable', 'comparisonTable', 'sweepTable', 'vifTable');
fprintf('完了: %.1f 秒。結果は %s に保存しました。\n', toc(totalTimer), cfg.resultDir);
end

% ------------------------------------------------------------------------
function make_folder(folderPath)
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
end

function methods = build_methods()
% 回帰手法の一覧を作る。損失関数ごとに「罰則なし、L2・L1を高次項のみ、L2・L1を全項」の5種類。
lossNames = {'OLS', 'LAD', 'Huber'};
settings = {'none', 'none'; 'L2', 'HO'; 'L1', 'HO'; 'L2', 'All'; 'L1', 'All'};   % {罰則, 罰則の対象}
methods = struct('name', {}, 'loss', {}, 'penalty', {}, 'scope', {});
for l = 1:numel(lossNames)
    for s = 1:size(settings, 1)
        methods(end + 1) = struct('name', method_name(lossNames{l}, settings{s, 1}, settings{s, 2}), ...
            'loss', lossNames{l}, 'penalty', settings{s, 1}, 'scope', settings{s, 2}); %#ok<AGROW>
    end
end
end

function name = method_name(loss, penalty, scope)
% 手法名。罰則なしは損失関数名だけ（例: 'Huber'）、罰則ありは対象も付ける（例: 'Huber+L2(HO)'）
% HO: 高次項（2次以上）のみに罰則、All: 切片・線形成分を含む全項に罰則
if strcmp(penalty, 'none')
    name = loss;
else
    name = sprintf('%s+%s(%s)', loss, penalty, scope);
end
end

function isPenalized = penalty_mask(design, scope)
% 罰則をかける項を表す論理ベクトル（項数×1）
switch scope
    case 'All'
        isPenalized = true(size(design.exponents, 1), 1);
    case 'HO'
        isPenalized = design.isHighOrder;
    otherwise
        isPenalized = false(size(design.exponents, 1), 1);
end
end

function geometry = build_geometry(cfg)
% マーク位置と評価グリッドを、Shot中心を原点とする座標 [mm] で作る
[markXgrid, markYgrid] = meshgrid(cfg.markXmm, cfg.markYmm);
geometry.markX = markXgrid(:);
geometry.markY = markYgrid(:);
evalAxisX = linspace(-cfg.shotHalfWidthMm, cfg.shotHalfWidthMm, cfg.evalGridCount);
evalAxisY = linspace(-cfg.shotHalfHeightMm, cfg.shotHalfHeightMm, cfg.evalGridCount);
[evalXgrid, evalYgrid] = meshgrid(evalAxisX, evalAxisY);
geometry.evalX = evalXgrid(:);
geometry.evalY = evalYgrid(:);
end

function [variants, vifTable] = build_model_variants(cfg, geometry)
% ASML式・ニコン式と、それぞれのVIF削減版を作る。
% VIFで削除する項がない式は、削減版が元の式と同じになるので削減版を作らない。
baseNames = {'ASML', 'Nikon'};
axisNames = {'x', 'y'};
variants = struct('name', {}, 'x', {}, 'y', {});
vifRows = {};
for b = 1:numel(baseNames)
    base = model_terms(baseNames{b});
    fullDesign = struct();
    reducedDesign = struct();
    isAnyRemoved = false;
    for a = 1:numel(axisNames)
        axisName = axisNames{a};
        exponents = base.(axisName).exp;
        names = base.(axisName).name;
        Araw = design_matrix(geometry.markX, geometry.markY, exponents);
        vifFull = compute_vif(Araw, exponents);
        [keepMask, removalLog] = reduce_terms_by_vif(Araw, exponents, cfg.vifThreshold);
        isAnyRemoved = isAnyRemoved || ~all(keepMask);
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
    if isAnyRemoved
        variants(end + 1) = struct('name', [baseNames{b} '-VIF'], 'x', reducedDesign.x, 'y', reducedDesign.y); %#ok<AGROW>
    end
end
vifTable = cell2table(vifRows, 'VariableNames', ...
    {'Model', 'Axis', 'Term', 'Monomial', 'Order', 'VIF', 'Kept', 'RemovalOrder', 'VIFAfter'});
end

function design = make_design(exponents, names, geometry)
% 回帰用に列を標準化する（切片以外を平均0・標準偏差1）。評価グリッドにも同じ変換をかける。
% 標準化により、項ごとの座標の桁の違いによらず、罰則が各項に同じ強さで効く。
Araw = design_matrix(geometry.markX, geometry.markY, exponents);
isIntercept = all(exponents == 0, 2)';
center = mean(Araw, 1);
scale = std(Araw, 1, 1);
center(isIntercept) = 0;
scale(isIntercept) = 1;
design.exponents   = exponents;
design.names       = names;
design.Astd        = (Araw - center) ./ scale;
design.Aeval       = (design_matrix(geometry.evalX, geometry.evalY, exponents) - center) ./ scale;
design.isHighOrder = sum(exponents, 2) >= 2;   % 2次以上の項（線形成分以外）
end

function combinedTable = summarize_combined(rmsTrue, rmsApparent, variants, methods, scenarioNames)
% dx と dy の残差RMSを Shot ごとに二乗平均して1つにまとめ、分布を要約する
rows = {};
for sc = 1:numel(scenarioNames)
    for v = 1:numel(variants)
        for m = 1:numel(methods)
            rmsCombined = combined_rms(rmsTrue, sc, v, m);
            apparentCombined = combined_rms(rmsApparent, sc, v, m);
            rows(end + 1, :) = {scenarioNames{sc}, variants(v).name, methods(m).name, methods(m).loss, ...
                methods(m).penalty, methods(m).scope, median(rmsCombined), prctile(rmsCombined, 95), ...
                median(apparentCombined)}; %#ok<AGROW>
        end
    end
end
combinedTable = cell2table(rows, 'VariableNames', {'Scenario', 'Model', 'Method', 'Loss', 'Penalty', ...
    'Scope', 'MedianRmsNm', 'P95RmsNm', 'MedianApparentRmsNm'});
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
        Aeval = variants(v).(axisName).Aeval;
        bestResidual = trueAtEval - Aeval * (Aeval \ trueAtEval);
        floorRms{v, a} = sqrt(mean(bestResidual.^2, 1));
        rows(end + 1, :) = {['Floor_' variants(v).name], upper(axisName), ...
            median(floorRms{v, a}), prctile(floorRms{v, a}, 95)}; %#ok<AGROW>
    end
end
referenceTable = cell2table(rows, 'VariableNames', {'Condition', 'Axis', 'MedianRmsNm', 'P95RmsNm'});
end

function comparisonTable = compare_conditions(rmsTrue, variants, methods, scenarioNames)
% 条件Aと条件Bを同じShotで比べる。X・Yの残差RMSは二乗平均で1つにまとめる。
% 列: 観点, 区分, 条件A（補正式・手法）, 条件B, Aの中央値, Bの中央値, 中央値の低減率, AがBより良いShotの割合
variantNames = {variants.name};
methodNames = {methods.name};
pairs = comparison_pairs(variantNames, methods);
rows = {};
for sc = 1:numel(scenarioNames)
    for p = 1:size(pairs, 1)
        rmsA = combined_rms(rmsTrue, sc, find(strcmp(variantNames, pairs{p, 3})), ...
            find(strcmp(methodNames, pairs{p, 4})));
        rmsB = combined_rms(rmsTrue, sc, find(strcmp(variantNames, pairs{p, 5})), ...
            find(strcmp(methodNames, pairs{p, 6})));
        rows(end + 1, :) = {scenarioNames{sc}, pairs{p, :}, median(rmsA), median(rmsB), ...
            100 * (1 - median(rmsA) / median(rmsB)), 100 * mean(rmsA < rmsB)}; %#ok<AGROW>
    end
end
comparisonTable = cell2table(rows, 'VariableNames', {'Scenario', 'Group', 'Block', 'ModelA', 'MethodA', ...
    'ModelB', 'MethodB', 'MedianA_nm', 'MedianB_nm', 'MedianReductionPercent', 'WinRateA_percent'});
end

function pairs = comparison_pairs(variantNames, methods)
% 比較する条件の組 {観点, 区分, 補正式A, 手法A, 補正式B, 手法B} を作る（Aが評価したい条件）。
%   Scope          : 罰則を高次項のみにかける場合 対 全項にかける場合
%   Loss           : 損失関数どうし（LAD対OLS、Huber対OLS、Huber対LAD）。罰則の条件ごと
%   Penalty        : L2 対 L1
%   Regularization : 罰則あり 対 罰則なし
%   VIF            : VIFで項を削除した補正式 対 削除前の補正式（削減版がある式だけ）
lossNames = unique({methods.loss}, 'stable');
penaltyNames = {'L2', 'L1'};
scopeNames = {'HO', 'All'};
settings = {'none', 'none'; 'L2', 'HO'; 'L1', 'HO'; 'L2', 'All'; 'L1', 'All'};
lossPairs = {'LAD', 'OLS'; 'Huber', 'OLS'; 'Huber', 'LAD'};
pairs = cell(0, 6);
for v = 1:numel(variantNames)
    model = variantNames{v};
    for l = 1:numel(lossNames)
        for p = 1:numel(penaltyNames)
            pairs(end + 1, :) = {'Scope', lossNames{l}, model, method_name(lossNames{l}, penaltyNames{p}, 'HO'), ...
                model, method_name(lossNames{l}, penaltyNames{p}, 'All')}; %#ok<AGROW>
        end
    end
end
for v = 1:numel(variantNames)
    model = variantNames{v};
    for s = 1:size(settings, 1)
        block = method_name('', settings{s, 1}, settings{s, 2});
        for q = 1:size(lossPairs, 1)
            pairs(end + 1, :) = {'Loss', block, model, method_name(lossPairs{q, 1}, settings{s, 1}, settings{s, 2}), ...
                model, method_name(lossPairs{q, 2}, settings{s, 1}, settings{s, 2})}; %#ok<AGROW>
        end
    end
end
for v = 1:numel(variantNames)
    model = variantNames{v};
    for l = 1:numel(lossNames)
        for c = 1:numel(scopeNames)
            pairs(end + 1, :) = {'Penalty', lossNames{l}, model, method_name(lossNames{l}, 'L2', scopeNames{c}), ...
                model, method_name(lossNames{l}, 'L1', scopeNames{c})}; %#ok<AGROW>
        end
    end
end
for v = 1:numel(variantNames)
    model = variantNames{v};
    for l = 1:numel(lossNames)
        for s = 2:size(settings, 1)
            pairs(end + 1, :) = {'Regularization', lossNames{l}, model, ...
                method_name(lossNames{l}, settings{s, 1}, settings{s, 2}), model, lossNames{l}}; %#ok<AGROW>
        end
    end
end
for v = 1:numel(variantNames)
    reducedName = [variantNames{v} '-VIF'];
    if ~any(strcmp(variantNames, reducedName))
        continue;
    end
    for m = 1:numel(methods)
        pairs(end + 1, :) = {'VIF', methods(m).loss, reducedName, methods(m).name, ...
            variantNames{v}, methods(m).name}; %#ok<AGROW>
    end
end
end

function rmsCombined = combined_rms(rmsCell, scenarioIndex, variantIndex, methodIndex)
rmsX = rmsCell{scenarioIndex, variantIndex, methodIndex, 1};
rmsY = rmsCell{scenarioIndex, variantIndex, methodIndex, 2};
rmsCombined = sqrt((rmsX.^2 + rmsY.^2) / 2);
end

function sweepTable = sweep_lambda(cfg, data, variants, axisNames)
% 異常値ありのシナリオで、Huber回帰に組み合わせる罰則の強さλを固定値の何倍かに変える
methods = build_methods();
methods = methods(strcmp({methods.loss}, 'Huber') & ~strcmp({methods.penalty}, 'none'));
rows = {};
for scaleValue = cfg.lambdaScaleList
    for v = 1:numel(variants)
        for m = 1:numel(methods)
            rmsByAxis = cell(1, numel(axisNames));
            for a = 1:numel(axisNames)
                axisName = axisNames{a};
                design = variants(v).(axisName);
                theta = fit_batch(design.Astd, data.(axisName).measuredOutlier, methods(m), ...
                    penalty_mask(design, methods(m).scope), cfg.lambdaL1 * scaleValue, cfg.lambdaL2 * scaleValue, cfg);
                residual = data.(axisName).trueAtEval - design.Aeval * theta;
                rmsByAxis{a} = sqrt(mean(residual.^2, 1));
            end
            rmsCombined = sqrt((rmsByAxis{1}.^2 + rmsByAxis{2}.^2) / 2);
            rows(end + 1, :) = {scaleValue, variants(v).name, methods(m).name, methods(m).penalty, ...
                methods(m).scope, median(rmsCombined), prctile(rmsCombined, 95)}; %#ok<AGROW>
            fprintf('  sweep x%-5.2f %-9s %-13s median = %.3f nm\n', scaleValue, variants(v).name, ...
                methods(m).name, median(rmsCombined));
        end
    end
end
sweepTable = cell2table(rows, 'VariableNames', {'LambdaScale', 'Model', 'Method', 'Penalty', 'Scope', ...
    'MedianRmsNm', 'P95RmsNm'});
end

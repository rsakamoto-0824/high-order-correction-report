function write_tex_tables()
%WRITE_TEX_TABLES シミュレーション結果からLaTeXの表の中身（行）を作り、results/tex/ に保存する。
%   先に run_all を実行して results/simulation_results.mat を作っておくこと。
%   出力した行を report/main.tex の tabular 環境へ貼り付けて使う（数値の転記ミスを防ぐため）。
%   残差の値と比較（低減率・勝率）は、dx と dy の残差RMSを Shot ごとに二乗平均した値で求めている。
cfg = simulation_config();
resultPath = fullfile(cfg.resultDir, 'simulation_results.mat');
if ~isfile(resultPath)
    error('write_tex_tables:noResult', '結果ファイルがありません。先に run_all を実行してください: %s', resultPath);
end
result = load(resultPath);
texDir = fullfile(cfg.resultDir, 'tex');
if ~exist(texDir, 'dir')
    mkdir(texDir);
end

write_summary_rows(result, 'outlier', 'P95', fullfile(texDir, 'rows_summary_outlier.tex'));
write_summary_rows(result, 'clean', 'Apparent', fullfile(texDir, 'rows_summary_clean.tex'));
scenarioNames = {'outlier', 'clean'};
withinModelGroups = {'Scope', 'Loss', 'Penalty', 'Regularization'};
for sc = 1:numel(scenarioNames)
    for g = 1:numel(withinModelGroups)
        fileName = sprintf('rows_%s_%s.tex', lower(withinModelGroups{g}), scenarioNames{sc});
        write_within_model_rows(result, scenarioNames{sc}, withinModelGroups{g}, fullfile(texDir, fileName));
    end
    write_vif_effect_rows(result, scenarioNames{sc}, ...
        fullfile(texDir, sprintf('rows_vif_effect_%s.tex', scenarioNames{sc})));
end
write_vif_rows(result, fullfile(texDir, 'rows_vif.tex'));
write_reference_rows(result, fullfile(texDir, 'rows_reference.tex'));
fprintf('LaTeXの表の行を保存しました: %s\n', texDir);
end

% ------------------------------------------------------------------------
function write_summary_rows(result, scenario, secondMetric, filePath)
% 列: 手法 | 罰則の対象 | 補正式ごとに（中央値 | 2列目）
%   2列目は 'P95'（95パーセンタイル）または 'Apparent'（測定点での見かけの残差の中央値）。
%   太字は補正式ごとに中央値が最小の手法。
combined = result.combinedTable(strcmp(result.combinedTable.Scenario, scenario), :);
variantNames = {result.variants.name};
methods = result.methods;
bestMedian = zeros(1, numel(variantNames));
for v = 1:numel(variantNames)
    bestMedian(v) = min(combined.MedianRmsNm(strcmp(combined.Model, variantNames{v})));
end
fileId = open_file(filePath);
for m = 1:numel(methods)
    if m > 1 && ~strcmp(methods(m).loss, methods(m - 1).loss)
        fprintf(fileId, '\\midrule\n');
    end
    cells = {method_label(methods(m).name), scope_label(methods(m).name)};
    for v = 1:numel(variantNames)
        row = combined(strcmp(combined.Model, variantNames{v}) & strcmp(combined.Method, methods(m).name), :);
        cells{end + 1} = bold_if_best(row.MedianRmsNm, bestMedian(v)); %#ok<AGROW>
        if strcmp(secondMetric, 'P95')
            cells{end + 1} = sprintf('%.3f', row.P95RmsNm); %#ok<AGROW>
        else
            cells{end + 1} = sprintf('%.3f', row.MedianApparentRmsNm); %#ok<AGROW>
        end
    end
    fprintf(fileId, '%s \\\\\n', strjoin(cells, ' & '));
end
fclose(fileId);
end

function write_within_model_rows(result, scenario, group, filePath)
% 同じ補正式の中での比較。列: 条件A | 条件B | 補正式ごとに（低減率 | 勝率）
%   区分（Block）が変わるところに罫線を入れる。低減率が正なら条件Aが優れる。
comparison = result.comparisonTable(strcmp(result.comparisonTable.Scenario, scenario) & ...
    strcmp(result.comparisonTable.Group, group), :);
variantNames = {result.variants.name};
pairKeys = unique(comparison(:, {'Block', 'MethodA', 'MethodB'}), 'stable');
fileId = open_file(filePath);
for k = 1:height(pairKeys)
    if k > 1 && ~strcmp(pairKeys.Block{k}, pairKeys.Block{k - 1})
        fprintf(fileId, '\\midrule\n');
    end
    cells = {pretty_method(pairKeys.MethodA{k}), pretty_method(pairKeys.MethodB{k})};
    for v = 1:numel(variantNames)
        row = comparison(strcmp(comparison.ModelA, variantNames{v}) & ...
            strcmp(comparison.MethodA, pairKeys.MethodA{k}) & strcmp(comparison.MethodB, pairKeys.MethodB{k}), :);
        cells(end + 1:end + 2) = {sprintf('%+.1f', row.MedianReductionPercent), sprintf('%.1f', row.WinRateA_percent)};
    end
    fprintf(fileId, '%s \\\\\n', strjoin(cells, ' & '));
end
fclose(fileId);
end

function write_vif_effect_rows(result, scenario, filePath)
% VIFで項を削除した効果。列: 補正式 | 手法 | 罰則の対象 | 削除前の中央値 | 削除後の中央値 | 低減率 | 勝率
comparison = result.comparisonTable(strcmp(result.comparisonTable.Scenario, scenario) & ...
    strcmp(result.comparisonTable.Group, 'VIF'), :);
fileId = open_file(filePath);
for r = 1:height(comparison)
    isNewModel = r == 1 || ~strcmp(comparison.ModelB{r}, comparison.ModelB{r - 1});
    if r > 1 && (isNewModel || ~strcmp(comparison.Block{r}, comparison.Block{r - 1}))
        fprintf(fileId, '\\midrule\n');
    end
    modelText = '';
    if isNewModel
        modelText = comparison.ModelB{r};
    end
    fprintf(fileId, '%s & %s & %s & %.3f & %.3f & %+.1f & %.1f \\\\\n', modelText, ...
        method_label(comparison.MethodA{r}), scope_label(comparison.MethodA{r}), comparison.MedianB_nm(r), ...
        comparison.MedianA_nm(r), comparison.MedianReductionPercent(r), comparison.WinRateA_percent(r));
end
fclose(fileId);
end

function write_vif_rows(result, filePath)
% 列: 補正式 | 軸 | 項数（削除前→後） | 削除した項（削除した順） | 最大VIF（削除前） | 最大VIF（削除後）
vifTable = result.vifTable;
combos = {'ASML', 'X'; 'ASML', 'Y'; 'Nikon', 'X'; 'Nikon', 'Y'};
fileId = open_file(filePath);
for c = 1:size(combos, 1)
    rows = vifTable(strcmp(vifTable.Model, combos{c, 1}) & strcmp(vifTable.Axis, combos{c, 2}), :);
    removed = rows(~rows.Kept, :);
    [~, order] = sort(removed.RemovalOrder);
    removed = removed(order, :);
    if isempty(removed)
        removedText = 'なし';
    else
        removedText = strjoin(arrayfun(@(i) sprintf('$%s$', removed.Monomial{i}), 1:height(removed), ...
            'UniformOutput', false), ', ');
    end
    fprintf(fileId, '%s & $d_%s$ & %d $\\to$ %d & %s & %.1f & %.1f \\\\\n', combos{c, 1}, lower(combos{c, 2}), ...
        height(rows), sum(rows.Kept), removedText, max(rows.VIF, [], 'omitnan'), max(rows.VIFAfter, [], 'omitnan'));
end
fclose(fileId);
end

function write_reference_rows(result, filePath)
% 列: 条件 | dx中央値 | dx 95% | dy中央値 | dy 95%
reference = result.referenceTable;
conditions = unique(reference.Condition, 'stable');
fileId = open_file(filePath);
for c = 1:numel(conditions)
    rowX = reference(strcmp(reference.Condition, conditions{c}) & strcmp(reference.Axis, 'X'), :);
    rowY = reference(strcmp(reference.Condition, conditions{c}) & strcmp(reference.Axis, 'Y'), :);
    label = strrep(strrep(conditions{c}, 'NoCorrection', '補正なし'), 'Floor_', '表現限界 ');
    fprintf(fileId, '%s & %.3f & %.3f & %.3f & %.3f \\\\\n', label, rowX.MedianRmsNm, rowX.P95RmsNm, ...
        rowY.MedianRmsNm, rowY.P95RmsNm);
end
fclose(fileId);
end

function text = method_label(methodName)
% 手法名から罰則の対象を除いた表示名（例: 'Huber+L2(HO)' -> 'Huber+L2'）
text = regexprep(methodName, '\((HO|All)\)$', '');
end

function text = scope_label(methodName)
% 罰則の対象の表示名
if endsWith(methodName, '(HO)')
    text = '高次項';
elseif endsWith(methodName, '(All)')
    text = '全項';
else
    text = '---';
end
end

function text = pretty_method(methodName)
% 比較の表で使う表示名（例: 'Huber+L2(HO)' -> 'Huber+L2（高次項）'）
text = method_label(methodName);
if ~strcmp(scope_label(methodName), '---')
    text = sprintf('%s（%s）', text, scope_label(methodName));
end
end

function text = bold_if_best(value, bestValue)
if strcmp(sprintf('%.3f', value), sprintf('%.3f', bestValue))
    text = sprintf('\\textbf{%.3f}', value);
else
    text = sprintf('%.3f', value);
end
end

function fileId = open_file(filePath)
fileId = fopen(filePath, 'w', 'n', 'UTF-8');
if fileId < 0
    error('write_tex_tables:cannotOpen', 'ファイルを開けません: %s', filePath);
end
end

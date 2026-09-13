function write_tex_tables()
%WRITE_TEX_TABLES シミュレーション結果からLaTeXの表の中身（行）を作り、results/tex/ に保存する。
%   先に run_all を実行して results/simulation_results.mat を作っておくこと。
%   出力した行を report/main.tex の tabular 環境へ貼り付けて使う（数値の転記ミスを防ぐため）。
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

write_summary_rows(result, 'outlier', fullfile(texDir, 'rows_summary_outlier.tex'));
write_clean_rows(result, fullfile(texDir, 'rows_summary_clean.tex'));
write_comparison_rows(result, 'outlier', fullfile(texDir, 'rows_comparisons_outlier.tex'));
write_comparison_rows(result, 'clean', fullfile(texDir, 'rows_comparisons_clean.tex'));
write_vif_rows(result, fullfile(texDir, 'rows_vif.tex'));
write_reference_rows(result, fullfile(texDir, 'rows_reference.tex'));
fprintf('LaTeXの表の行を保存しました: %s\n', texDir);
end

% ------------------------------------------------------------------------
function write_summary_rows(result, scenario, filePath)
% 列: 補正式 | 手法 | dx中央値 | dx 95% | dx |平均|+3σ | dy中央値 | dy 95% | dy |平均|+3σ
summary = result.summaryTable(strcmp(result.summaryTable.Scenario, scenario), :);
variantNames = {result.variants.name};
methodNames = result.methodNames;
fileId = open_file(filePath);
for v = 1:numel(variantNames)
    if v > 1
        fprintf(fileId, '\\midrule\n');
    end
    block = summary(strcmp(summary.Model, variantNames{v}), :);
    bestX = min(block.MedianRmsNm(strcmp(block.Axis, 'X')));
    bestY = min(block.MedianRmsNm(strcmp(block.Axis, 'Y')));
    for m = 1:numel(methodNames)
        rowX = block(strcmp(block.Method, methodNames{m}) & strcmp(block.Axis, 'X'), :);
        rowY = block(strcmp(block.Method, methodNames{m}) & strcmp(block.Axis, 'Y'), :);
        fprintf(fileId, '%s & %s & %s & %.3f & %.3f & %s & %.3f & %.3f \\\\\n', ...
            model_cell(variantNames{v}, rowX.TermCount, rowY.TermCount, m), methodNames{m}, ...
            bold_if_best(rowX.MedianRmsNm, bestX), rowX.P95RmsNm, rowX.MeanPlus3SigmaNm, ...
            bold_if_best(rowY.MedianRmsNm, bestY), rowY.P95RmsNm, rowY.MeanPlus3SigmaNm);
    end
end
fclose(fileId);
end

function write_clean_rows(result, filePath)
% 列: 補正式 | 手法 | 異常値なし dx中央値 | dy中央値 | 測定点での見かけの残差 dx | dy
summary = result.summaryTable(strcmp(result.summaryTable.Scenario, 'clean'), :);
variantNames = {result.variants.name};
methodNames = result.methodNames;
fileId = open_file(filePath);
for v = 1:numel(variantNames)
    if v > 1
        fprintf(fileId, '\\midrule\n');
    end
    block = summary(strcmp(summary.Model, variantNames{v}), :);
    bestX = min(block.MedianRmsNm(strcmp(block.Axis, 'X')));
    bestY = min(block.MedianRmsNm(strcmp(block.Axis, 'Y')));
    for m = 1:numel(methodNames)
        rowX = block(strcmp(block.Method, methodNames{m}) & strcmp(block.Axis, 'X'), :);
        rowY = block(strcmp(block.Method, methodNames{m}) & strcmp(block.Axis, 'Y'), :);
        fprintf(fileId, '%s & %s & %s & %s & %.3f & %.3f \\\\\n', ...
            model_cell(variantNames{v}, rowX.TermCount, rowY.TermCount, m), methodNames{m}, ...
            bold_if_best(rowX.MedianRmsNm, bestX), bold_if_best(rowY.MedianRmsNm, bestY), ...
            rowX.MedianApparentRmsNm, rowY.MedianApparentRmsNm);
    end
end
fclose(fileId);
end

function write_comparison_rows(result, scenario, filePath)
% 列: 比較の観点 | 条件A | 条件B | A中央値 | B中央値 | 中央値の低減率 | Aの勝率
comparison = result.comparisonTable(strcmp(result.comparisonTable.Scenario, scenario), :);
groupLabel = containers.Map( ...
    {'C1_VIFbeforeL1', 'C2_L2vsL1_VIF', 'C2_L2vsL1_Full', 'C3_Robust', 'C4_VIFonly'}, ...
    {'事前削除（L1併用）', 'L2とL1（削除後）', 'L2とL1（削除前）', 'ロバスト回帰', '事前削除（正則化なし）'});
fileId = open_file(filePath);
previousGroup = '';
for r = 1:height(comparison)
    group = comparison.Group{r};
    label = '';
    if ~strcmp(group, previousGroup)
        if ~isempty(previousGroup)
            fprintf(fileId, '\\midrule\n');
        end
        label = groupLabel(group);
        previousGroup = group;
    end
    fprintf(fileId, '%s & %s & %s & %.3f & %.3f & %+.1f & %.1f \\\\\n', label, ...
        comparison.ConditionA{r}, comparison.ConditionB{r}, comparison.MedianA_nm(r), ...
        comparison.MedianB_nm(r), comparison.MedianReductionPercent(r), comparison.WinRateA_percent(r));
end
fclose(fileId);
end

function write_vif_rows(result, filePath)
% 列: 補正式 | 軸 | 項数（削除前→後） | 削除した項（削除時のVIF） | 最大VIF（削除前） | 最大VIF（削除後）
vifTable = result.vifTable;
combos = {'ASML', 'X'; 'ASML', 'Y'; 'Nikon', 'X'; 'Nikon', 'Y'};
fileId = open_file(filePath);
for c = 1:size(combos, 1)
    rows = vifTable(strcmp(vifTable.Model, combos{c, 1}) & strcmp(vifTable.Axis, combos{c, 2}), :);
    removed = rows(~rows.Kept, :);
    [~, order] = sort(removed.RemovalOrder);
    removed = removed(order, :);
    removedText = strjoin(arrayfun(@(i) sprintf('$%s$', removed.Monomial{i}), 1:height(removed), ...
        'UniformOutput', false), ', ');
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

function text = model_cell(variantName, termCountX, termCountY, rowIndex)
% ブロックの1行目に補正式名、2行目に項数（dx/dy）を書く
switch rowIndex
    case 1
        text = variantName;
    case 2
        text = sprintf('(%d/%d)', termCountX, termCountY);
    otherwise
        text = '';
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

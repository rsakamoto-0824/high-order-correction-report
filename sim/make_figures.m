function make_figures()
%MAKE_FIGURES シミュレーション結果からレポート用の図（PNG, 300 dpi）を作る。
%   先に run_all を実行して results/simulation_results.mat を作っておくこと。
%   図は report/figures/ に保存する。文字化けを避けるため図中の文字は英語にしている。
%   MATLABのベクターPDFはフォントが埋め込まれず、Overleaf等で文字が消えることがあるためPNGで出力する。
cfg = simulation_config();
resultPath = fullfile(cfg.resultDir, 'simulation_results.mat');
if ~isfile(resultPath)
    error('make_figures:noResult', '結果ファイルがありません。先に run_all を実行してください: %s', resultPath);
end
result = load(resultPath);
if ~exist(cfg.figureDir, 'dir')
    mkdir(cfg.figureDir);
end

style.fontSize       = 9;
style.fontName       = 'Helvetica';
style.resolutionDpi  = 300;
style.colorKept      = [0.00 0.35 0.60];
style.colorRemoved   = [0.90 0.55 0.10];
style.colorOutlier   = [0.80 0.10 0.10];
style.colorL1        = [0.85 0.33 0.10];
style.colorL2        = [0.00 0.45 0.74];
style.colorNoPenalty = [0.55 0.55 0.55];   % 残差分布図: 罰則なし
style.colorHighOrder = [0.20 0.60 0.20];   % 残差分布図: 高次項のみに罰則
style.colorAllTerms  = [0.55 0.25 0.65];   % 残差分布図: 全項に罰則
style.vectorMmPerNm  = 0.6;          % ベクトル図: 1 nm を何 mm の矢印で描くか
style.scaleBarNm     = 2;            % ベクトル図の目盛りの長さ [nm]
style.boxYLimitNm    = [0.05 3];     % 残差分布図の縦軸範囲 [nm]

plot_shot_example(result, cfg, style);
plot_vif(result, cfg, style);
plot_residual_box(result, cfg, style);
plot_lambda_sweep(result, cfg, style);
fprintf('図を保存しました: %s\n', cfg.figureDir);
end

% ------------------------------------------------------------------------
function plot_shot_example(result, cfg, style)
data = result.data;
geometry = result.geometry;
exampleShot = find(any(data.isOutlierMark, 1), 1);
if isempty(exampleShot)
    exampleShot = 1;
end
evalXmm = geometry.evalX;
evalYmm = geometry.evalY;
markXmm = geometry.markX;
markYmm = geometry.markY;
gridIndexX = round((evalXmm / cfg.shotHalfWidthMm + 1) * (cfg.evalGridCount - 1) / 2);
gridIndexY = round((evalYmm / cfg.shotHalfHeightMm + 1) * (cfg.evalGridCount - 1) / 2);
isShown = mod(gridIndexX, 2) == 0 & mod(gridIndexY, 2) == 0;   % 見やすさのため1点おきに表示
k = style.vectorMmPerNm;

fig = new_figure(style, 16, 8.5);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
prepare_shot_axes(ax, cfg, style, k);
quiver(ax, evalXmm(isShown), evalYmm(isShown), k * data.x.trueAtEval(isShown, exampleShot), ...
    k * data.y.trueAtEval(isShown, exampleShot), 0, 'Color', style.colorKept, 'MaxHeadSize', 0.4);
plot(ax, markXmm, markYmm, 'ks', 'MarkerSize', 4, 'MarkerFaceColor', [0.3 0.3 0.3]);
title(ax, '(a) True distortion (arrows) and marks (squares)', 'FontWeight', 'normal');

ax = nexttile(layout);
prepare_shot_axes(ax, cfg, style, k);
isOutlier = data.isOutlierMark(:, exampleShot);
measuredX = data.x.measuredOutlier(:, exampleShot);
measuredY = data.y.measuredOutlier(:, exampleShot);
quiver(ax, markXmm(~isOutlier), markYmm(~isOutlier), k * measuredX(~isOutlier), k * measuredY(~isOutlier), 0, ...
    'Color', style.colorKept, 'MaxHeadSize', 0.4);
quiver(ax, markXmm(isOutlier), markYmm(isOutlier), k * measuredX(isOutlier), k * measuredY(isOutlier), 0, ...
    'Color', style.colorOutlier, 'LineWidth', 1.2, 'MaxHeadSize', 0.4);
title(ax, '(b) Measured values (red: outlier marks)', 'FontWeight', 'normal');

save_figure(fig, cfg, style, 'fig_shot_example.png');
end

function prepare_shot_axes(ax, cfg, style, vectorMmPerNm)
hold(ax, 'on');
box(ax, 'on');
rectangle(ax, 'Position', [-cfg.shotHalfWidthMm, -cfg.shotHalfHeightMm, ...
    2 * cfg.shotHalfWidthMm, 2 * cfg.shotHalfHeightMm], 'EdgeColor', [0.6 0.6 0.6]);
barStartX = -cfg.shotHalfWidthMm;
barY = -cfg.shotHalfHeightMm - 3;
plot(ax, barStartX + [0, style.scaleBarNm * vectorMmPerNm], [barY barY], 'k-', 'LineWidth', 1.5);
text(ax, barStartX + style.scaleBarNm * vectorMmPerNm + 0.8, barY, sprintf('%g nm', style.scaleBarNm), ...
    'FontSize', style.fontSize - 1, 'VerticalAlignment', 'middle');
axis(ax, 'equal');
xlim(ax, [-cfg.shotHalfWidthMm - 4, cfg.shotHalfWidthMm + 4]);
ylim(ax, [-cfg.shotHalfHeightMm - 5, cfg.shotHalfHeightMm + 4]);
xlabel(ax, 'x [mm] (slit)');
ylabel(ax, 'y [mm] (scan)');
end

function plot_vif(result, cfg, style)
vifTable = result.vifTable;
combos = {'ASML', 'X'; 'ASML', 'Y'; 'Nikon', 'X'; 'Nikon', 'Y'};
fig = new_figure(style, 16, 11);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:size(combos, 1)
    isRow = strcmp(vifTable.Model, combos{c, 1}) & strcmp(vifTable.Axis, combos{c, 2}) & ~isnan(vifTable.VIF);
    subTable = vifTable(isRow, :);
    termCount = height(subTable);
    ax = nexttile(layout);
    hold(ax, 'on');
    box(ax, 'on');
    barHandle = bar(ax, 1:termCount, subTable.VIF, 'FaceColor', 'flat', 'BaseValue', 0.5);
    barColor = repmat(style.colorKept, termCount, 1);
    barColor(~subTable.Kept, :) = repmat(style.colorRemoved, sum(~subTable.Kept), 1);
    barHandle.CData = barColor;
    set(ax, 'YScale', 'log');
    yline(ax, cfg.vifThreshold, '--k');
    for t = find(~subTable.Kept)'
        text(ax, t, subTable.VIF(t) * 1.6, sprintf('#%d', subTable.RemovalOrder(t)), ...
            'HorizontalAlignment', 'center', 'FontSize', style.fontSize - 2, 'Color', style.colorRemoved);
    end
    % 対数軸の目盛りを10のべき乗ごとに明示する（自動では上位の目盛りラベルが省略されるため）
    upperExponent = ceil(log10(max(subTable.VIF) * 3));
    ylim(ax, [0.5, 10 ^ upperExponent]);
    yticks(ax, 10 .^ (0:upperExponent));
    xlim(ax, [0.4, termCount + 0.6]);
    xticks(ax, 1:termCount);
    xticklabels(ax, subTable.Monomial);
    ax.TickLabelInterpreter = 'tex';
    ax.XTickLabelRotation = 0;
    ax.FontSize = style.fontSize - 1;
    ylabel(ax, 'VIF');
    title(ax, sprintf('%s model, d%s', combos{c, 1}, lower(combos{c, 2})), 'FontWeight', 'normal', ...
        'FontSize', style.fontSize);
end
title(layout, sprintf('Blue: kept, orange: removed (#: removal order), dashed line: VIF = %g', ...
    cfg.vifThreshold), 'FontSize', style.fontSize);
save_figure(fig, cfg, style, 'fig_vif.png');
end

function plot_residual_box(result, cfg, style)
% 異常値ありのシナリオで、dx と dy をまとめた補正残差RMSの分布を、補正式ごとに1段で描く。
% 箱の色で罰則の対象（なし / 高次項のみ / 全項）を区別する。
scenarioIndex = find(strcmp(result.scenarioNames, 'outlier'));
variantNames = {result.variants.name};
methods = result.methods;
scopeNames  = {'none', 'HO', 'All'};
scopeColors = {style.colorNoPenalty, style.colorHighOrder, style.colorAllTerms};
scopeLabels = {'No penalty', 'Penalty on higher-order terms only', 'Penalty on all terms'};
lossBoundary = find(~strcmp({methods(2:end).loss}, {methods(1:end - 1).loss})) + 0.5;

fig = new_figure(style, 16, 19);
layout = tiledlayout(fig, numel(variantNames), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for v = 1:numel(variantNames)
    ax = nexttile(layout);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');
    for s = 1:numel(scopeNames)
        positions = [];
        values = [];
        for m = find(strcmp({methods.scope}, scopeNames{s}))
            rmsCombined = combined_rms(result.rmsTrue, scenarioIndex, v, m);
            positions = [positions; repmat(m, numel(rmsCombined), 1)]; %#ok<AGROW>
            values = [values; rmsCombined(:)]; %#ok<AGROW>
        end
        boxchart(ax, positions, values, 'BoxFaceColor', scopeColors{s}, 'MarkerStyle', 'none', ...
            'DisplayName', scopeLabels{s});
    end
    for boundary = lossBoundary
        xline(ax, boundary, '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    end
    set(ax, 'YScale', 'log');
    ylim(ax, style.boxYLimitNm);
    yticks(ax, [0.05 0.1 0.2 0.5 1 2]);
    xlim(ax, [0.4, numel(methods) + 0.6]);
    xticks(ax, 1:numel(methods));
    ax.TickLabelInterpreter = 'none';
    if v == numel(variantNames)
        xticklabels(ax, {methods.name});
        ax.XTickLabelRotation = 45;
    else
        xticklabels(ax, {});
    end
    ylabel(ax, 'Residual RMS [nm]');
    title(ax, sprintf('%s model', variantNames{v}), 'FontWeight', 'normal');
    if v == 1
        legend(ax, 'Location', 'northoutside', 'NumColumns', numel(scopeNames), 'Box', 'off');
    end
end
save_figure(fig, cfg, style, 'fig_residual_box.png');
end

function plot_lambda_sweep(result, cfg, style)
% Huber回帰に組み合わせる罰則の強さλを変えたときの残差の中央値（補正式ごとに1枚）
sweepTable = result.sweepTable;
scenarioIndex = find(strcmp(result.scenarioNames, 'outlier'));
variantNames = {result.variants.name};
huberIndex = find(strcmp({result.methods.name}, 'Huber'));
seriesList = {   % {罰則, 対象, 線種, 色, 凡例}
    'L1', 'HO',  '-o',  style.colorL1, 'Huber+L1, higher-order terms only'
    'L2', 'HO',  '-s',  style.colorL2, 'Huber+L2, higher-order terms only'
    'L1', 'All', '--o', style.colorL1, 'Huber+L1, all terms'
    'L2', 'All', '--s', style.colorL2, 'Huber+L2, all terms'
    };
fig = new_figure(style, 16, 8.5);
layout = tiledlayout(fig, 1, numel(variantNames), 'TileSpacing', 'compact', 'Padding', 'compact');
axesList = gobjects(1, numel(variantNames));
for v = 1:numel(variantNames)
    ax = nexttile(layout);
    axesList(v) = ax;
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');
    for s = 1:size(seriesList, 1)
        isRow = strcmp(sweepTable.Model, variantNames{v}) & strcmp(sweepTable.Penalty, seriesList{s, 1}) & ...
            strcmp(sweepTable.Scope, seriesList{s, 2});
        plot(ax, sweepTable.LambdaScale(isRow), sweepTable.MedianRmsNm(isRow), seriesList{s, 3}, ...
            'Color', seriesList{s, 4}, 'MarkerSize', 4, 'DisplayName', seriesList{s, 5});
    end
    rmsCombined = combined_rms(result.rmsTrue, scenarioIndex, v, huberIndex);
    yline(ax, median(rmsCombined), ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, 'DisplayName', 'Huber (no penalty)');
    xline(ax, 1, ':k', 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log');
    xticks(ax, cfg.lambdaScaleList);
    xlabel(ax, '\lambda / \lambda_{fixed}');
    if v == 1
        ylabel(ax, 'Median residual RMS [nm]');
    end
    title(ax, sprintf('%s model', variantNames{v}), 'FontWeight', 'normal');
end
linkaxes(axesList, 'y');
ylim(axesList(1), [0, max(sweepTable.MedianRmsNm) * 1.05]);
legendHandle = legend(axesList(end), 'NumColumns', 3, 'FontSize', style.fontSize - 1, 'Box', 'off');
legendHandle.Layout.Tile = 'south';
save_figure(fig, cfg, style, 'fig_lambda_sweep.png');
end

function rmsCombined = combined_rms(rmsCell, scenarioIndex, variantIndex, methodIndex)
% dx と dy の残差RMSを Shot ごとに二乗平均して1つにまとめる
rmsX = rmsCell{scenarioIndex, variantIndex, methodIndex, 1};
rmsY = rmsCell{scenarioIndex, variantIndex, methodIndex, 2};
rmsCombined = sqrt((rmsX.^2 + rmsY.^2) / 2);
end

function fig = new_figure(style, widthCm, heightCm)
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [1 1 widthCm heightCm], 'Color', 'w');
set(fig, 'DefaultAxesFontSize', style.fontSize, 'DefaultAxesFontName', style.fontName, ...
    'DefaultTextFontSize', style.fontSize, 'DefaultTextFontName', style.fontName);
end

function save_figure(fig, cfg, style, fileName)
exportgraphics(fig, fullfile(cfg.figureDir, fileName), 'Resolution', style.resolutionDpi);
close(fig);
end

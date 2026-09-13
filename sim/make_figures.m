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

style.fontSize      = 9;
style.fontName      = 'Helvetica';
style.resolutionDpi = 300;
style.colorKept     = [0.00 0.35 0.60];
style.colorRemoved  = [0.90 0.55 0.10];
style.colorOutlier  = [0.80 0.10 0.10];
style.colorL1       = [0.85 0.33 0.10];
style.colorL2       = [0.00 0.45 0.74];
style.vectorMmPerNm = 0.6;          % ベクトル図: 1 nm を何 mm の矢印で描くか
style.scaleBarNm    = 2;            % ベクトル図の目盛りの長さ [nm]
style.boxYLimitNm   = [0.03 3];     % 残差分布図の縦軸範囲 [nm]

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
evalXmm = geometry.evalX * cfg.shotHalfWidthMm;
evalYmm = geometry.evalY * cfg.shotHalfHeightMm;
markXmm = geometry.markX * cfg.shotHalfWidthMm;
markYmm = geometry.markY * cfg.shotHalfHeightMm;
gridIndexX = round((geometry.evalX + 1) * (cfg.evalGridCount - 1) / 2);
gridIndexY = round((geometry.evalY + 1) * (cfg.evalGridCount - 1) / 2);
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
scenarioIndex = find(strcmp(result.scenarioNames, 'outlier'));
variantNames = {result.variants.name};
methodNames = result.methodNames;
axisLabels = {'dx', 'dy'};
fig = new_figure(style, 16, 13);
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for a = 1:numel(axisLabels)
    values = [];
    modelGroup = {};
    methodGroup = {};
    for v = 1:numel(variantNames)
        for m = 1:numel(methodNames)
            rms = result.rmsTrue{scenarioIndex, v, m, a};
            values = [values; rms(:)]; %#ok<AGROW>
            modelGroup = [modelGroup; repmat(variantNames(v), numel(rms), 1)]; %#ok<AGROW>
            methodGroup = [methodGroup; repmat(methodNames(m), numel(rms), 1)]; %#ok<AGROW>
        end
    end
    ax = nexttile(layout);
    boxchart(ax, categorical(modelGroup, variantNames), values, ...
        'GroupByColor', categorical(methodGroup, methodNames), 'MarkerStyle', 'none');
    set(ax, 'YScale', 'log');
    ylim(ax, style.boxYLimitNm);
    yticks(ax, [0.05 0.1 0.2 0.5 1 2]);
    grid(ax, 'on');
    ylabel(ax, sprintf('Residual RMS of %s [nm]', axisLabels{a}));
    if a == 1
        legend(ax, 'Location', 'northoutside', 'NumColumns', numel(methodNames), 'Box', 'off');
    end
end
save_figure(fig, cfg, style, 'fig_residual_box.png');
end

function plot_lambda_sweep(result, cfg, style)
sweepTable = result.sweepTable;
scenarioIndex = find(strcmp(result.scenarioNames, 'outlier'));
variantNames = {result.variants.name};
huberIndex = find(strcmp(result.methodNames, 'Huber'));
families = {'ASML', 'Nikon'};
fig = new_figure(style, 16, 10);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for f = 1:numel(families)
    ax = nexttile(layout);
    hold(ax, 'on');
    box(ax, 'on');
    grid(ax, 'on');
    series = {
        families{f},           'Huber+L1', '-o',  style.colorL1
        families{f},           'Huber+L2', '-s',  style.colorL2
        [families{f} '-VIF'],  'Huber+L1', '--o', style.colorL1
        [families{f} '-VIF'],  'Huber+L2', '--s', style.colorL2
        };
    for s = 1:size(series, 1)
        isRow = strcmp(sweepTable.Model, series{s, 1}) & strcmp(sweepTable.Method, series{s, 2});
        plot(ax, sweepTable.LambdaScale(isRow), sweepTable.MedianRmsNm(isRow), series{s, 3}, ...
            'Color', series{s, 4}, 'MarkerSize', 4, 'DisplayName', [series{s, 1} ' ' series{s, 2}]);
    end
    referenceStyle = {':', '-.'};
    referenceModels = {families{f}, [families{f} '-VIF']};
    for r = 1:numel(referenceModels)
        v = find(strcmp(variantNames, referenceModels{r}));
        rmsX = result.rmsTrue{scenarioIndex, v, huberIndex, 1};
        rmsY = result.rmsTrue{scenarioIndex, v, huberIndex, 2};
        yline(ax, median(sqrt((rmsX.^2 + rmsY.^2) / 2)), referenceStyle{r}, 'Color', [0.4 0.4 0.4], ...
            'LineWidth', 1, 'DisplayName', [referenceModels{r} ' Huber (no penalty)']);
    end
    xline(ax, 1, ':k', 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log');
    xticks(ax, cfg.lambdaScaleList);
    xlabel(ax, '\lambda / \lambda_{fixed}');
    ylabel(ax, 'Median residual RMS [nm]');
    title(ax, sprintf('%s-based models', families{f}), 'FontWeight', 'normal');
    legend(ax, 'Location', 'southoutside', 'NumColumns', 2, 'FontSize', style.fontSize - 2, 'Box', 'off');
end
save_figure(fig, cfg, style, 'fig_lambda_sweep.png');
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

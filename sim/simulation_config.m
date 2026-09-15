function cfg = simulation_config()
%SIMULATION_CONFIG 実験条件を1か所にまとめて返す。
%   条件を変えるときは、この関数の値だけを編集する。

% ---- 再現性 ----
cfg.randomSeed = 20260913;

% ---- Shotとマーク配置（単位: mm、Shot中心基準）----
% x: スリット方向、y: スキャン方向。補正式の x, y にはこの座標をそのまま使う。
cfg.shotHalfWidthMm  = 13.0;              % Shot幅 26 mm の半分
cfg.shotHalfHeightMm = 16.5;              % Shot高さ 33 mm の半分
cfg.markXmm = linspace(-12.0, 12.0, 5);   % 5列
cfg.markYmm = linspace(-15.5, 15.5, 7);   % 7行（ニコン式の y^6 項を解くには7行以上が必要）
cfg.evalGridCount = 21;                   % 補正残差を評価する密グリッド（Shot全面 21×21点）

% ---- 生成するShot ----
cfg.shotCount = 1000;
% 真の歪み係数: 次数0〜6ごとの標準偏差 [nm] と、その項が現れる確率。
% 標準偏差は「Shot端（|x| = 13 mm、|y| = 16.5 mm）での変位量」で表し、
% generate_shots で mm座標の係数（nm/mm^次数）に換算する。
cfg.truthSigmaNmByOrder    = [2.0 2.0 0.8 0.6 0.4 0.4 0.4];
cfg.truthActiveProbByOrder = [1.0 1.0 0.5 0.5 0.3 0.3 0.3];

% ---- 測定誤差 ----
cfg.noiseSigmaNm       = 0.3;    % 計測の偶然誤差（1σ）
cfg.outlierProbability = 0.05;   % マークが異常値（フライヤー）になる確率
cfg.outlierMinNm       = 3.0;    % 異常値の大きさの範囲（符号はランダム）
cfg.outlierMaxNm       = 8.0;

% ---- VIFによる項削減 ----
cfg.vifThreshold = 10;

% ---- ロバスト回帰（LAD・Huber M推定、IRLS）----
cfg.huberTuning = 1.345;    % 正規誤差のもとで漸近効率95%となる定数
cfg.madToSigma  = 1.4826;   % MADを標準偏差に換算する係数
cfg.ladMinResidualRatio = 0.1;    % LADの重み s/|r| で、|r| を「この比 × s」未満として扱わない（重みの発散防止）
cfg.irlsScaleUpdateIter = 20;     % 尺度 s を更新する反復回数。以後は固定して収束させる
cfg.irlsMaxIter = 1000;
cfg.irlsTolNm   = 1e-5;           % 係数の変化がこれ未満になったShotから反復を止める

% ---- 正則化（ハイパーパラメータは固定）----
% 罰則の対象は「高次項のみ」と「切片・線形成分（並進・倍率・回転）を含む全項」の2通りで評価する。
markCount = numel(cfg.markXmm) * numel(cfg.markYmm);
cfg.lambdaL1 = cfg.noiseSigmaNm / sqrt(markCount);  % 係数の標準誤差1個分を打ち切り幅にする
cfg.lambdaL2 = 0.1;                                  % 標準化した直交列を 1/(1+0.1) に縮小する強さ
cfg.cdMaxSweep = 20000;
cfg.cdTolNm    = 1e-7;

% ---- λ感度解析（考察用）----
cfg.lambdaScaleList = [0.25 0.5 1 2 4 8];

% ---- 出力先 ----
simDir = fileparts(mfilename('fullpath'));
cfg.resultDir = fullfile(simDir, '..', 'results');
cfg.figureDir = fullfile(simDir, '..', 'report', 'figures');
end

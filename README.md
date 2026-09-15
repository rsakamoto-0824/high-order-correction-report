# high-order-correction-report

リソグラフィー工程における **Shot高次補正の回帰手法** を検討した技術レポート（TeX）と、
その実験に使ったシミュレーション（MATLAB）一式。

- VIF（分散拡大係数）によるkパラメータの事前選定
- ロバスト回帰（損失関数: OLS / LAD / Huber）
- L1 / L2 正則化（罰則を高次項のみにかける場合と、線形成分を含む全項にかける場合）

を組み合わせた15種類の回帰手法を、露光機メーカーごとのShot高次補正演算式
（ASML式: dx 7項・dy 8項、ニコン式: dx 12項・dy 15項）に適用し、乱数で作った1000通りのShotで比較する。

## 利用手順書

シミュレーションの実行手順と、Overleafでの `main.tex` の編集手順は
[docs/usage-guide.html](docs/usage-guide.html) にまとめている。
GitHub上ではHTMLのソースが表示されるため、リポジトリを取得（clone）してからブラウザで開く。

## 処理フロー資料

評価に使ったスクリプト（`sim/`）の処理の流れ、計算式、出力ファイルの列は
[docs/processing-flow.html](docs/processing-flow.html) にまとめている（同じくcloneしてからブラウザで開く）。

## フォルダ構成

```
high-order-correction-report/
├── docs/
│   ├── usage-guide.html      # 利用手順書（シミュレーション・Overleaf編集）
│   └── processing-flow.html  # 処理フロー資料（評価スクリプトの詳細）
├── report/
│   ├── main.tex          # レポート本体（Overleafにアップロードする）
│   └── figures/          # 図（make_figures.m で生成、PNG 300 dpi）
├── sim/                  # MATLABシミュレーション
│   ├── run_all.m             # 実行の入口（回帰・評価・CSV出力）
│   ├── make_figures.m        # 図を作る
│   ├── write_tex_tables.m    # LaTeXの表の行を作る
│   ├── simulation_config.m   # 実験条件（ここだけ編集すれば条件を変えられる）
│   ├── model_terms.m         # ASML式・ニコン式の項の定義
│   ├── generate_shots.m      # 高次歪みShotと測定誤差の生成
│   ├── compute_vif.m / reduce_terms_by_vif.m
│   ├── fit_batch.m           # 回帰（OLS / LAD / Huber × 罰則なし / L2 / L1 × 罰則の対象）
│   └── design_matrix.m / monomial_label.m
└── results/              # 計算結果（CSV、.mat、表の行）
```

## 実行方法（シミュレーション）

必要な環境: MATLAB R2022b 以降、Statistics and Machine Learning Toolbox（`prctile` を使用）

```matlab
cd sim
run_all          % 1000 Shot（実測で約9分。λ感度解析を含む）。動作確認だけなら run_all(50)（約46秒）
make_figures     % report/figures/ に図を出力
write_tex_tables % results/tex/ に表の行を出力
```

乱数シードは `simulation_config.m` で固定しているため、同じ結果が再現される。

## Overleafでのコンパイル

1. `report/` フォルダ（`main.tex` と `figures/`）をZIPにしてOverleafの「New Project → Upload Project」でアップロードする
2. メニュー → **Compiler を「LuaLaTeX」** にする（日本語クラス `ltjsarticle` を使うため）
3. Recompile

表紙の著者名と謝辞は `main.tex` 内の「← 記入する」のコメント箇所を書き換える。

## 注意事項

- 計測データはすべて乱数で生成したもので、実データは含まない。
- **補正式は露光機メーカーごとに定められた式を使う。** `sim/model_terms.m` の ASML式・ニコン式は、
  前提条件の式から項・係数番号・変数・符号を変えずに写したもの。メーカー間で共通化したり、項を補ったりしない。
- 補正式の `x`, `y` はShot中心を原点とするマーク座標 [mm]、測定値は [nm]。
- 正則化（L1 / L2）は、罰則を高次項（2次以上）のみにかける場合と、切片・線形成分（並進・倍率・回転）を
  含む全項にかける場合の両方を計算する。
- ハイパーパラメータ（λ1, λ2）は固定値であり、結論はその値に依存する部分がある（レポートの考察を参照）。
- 認証情報は使用していない。環境変数の設定は不要。

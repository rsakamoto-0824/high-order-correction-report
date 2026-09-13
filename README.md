# high-order-correction-report

リソグラフィー工程における **Shot高次補正の回帰手法** を検討した技術レポート（TeX）と、
その実験に使ったシミュレーション（MATLAB）一式。

- VIF（分散拡大係数）によるkパラメータの事前選定
- Huber損失によるロバスト回帰
- L1 / L2 正則化

を、ASML式（k1〜k20）とニコン式（Shot線形＋高次係数）に適用し、乱数で作った1000通りのShotで比較する。

## フォルダ構成

```
high-order-correction-report/
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
│   ├── fit_batch.m           # OLS / Huber / L2 / L1 / Huber+L2 / Huber+L1
│   └── design_matrix.m / monomial_label.m
└── results/              # 計算結果（CSV、.mat、表の行）
```

## 実行方法（シミュレーション）

必要な環境: MATLAB R2022b 以降、Statistics and Machine Learning Toolbox（`prctile` を使用）

```matlab
cd sim
run_all          % 1000 Shot（実測で約86分。λ感度解析を含む）。動作確認だけなら run_all(50)（約70秒）
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
- **ニコン式の項構成は仮定を含む。** 装置のShot高次補正係数名の添字 `mn` を `x^m y^n` の次数と解釈している。
  装置仕様で確認し、違っていれば `sim/model_terms.m` の `'Nikon'` の指数を修正して再実行する。
- ハイパーパラメータ（λ1, λ2）は固定値であり、結論はその値に依存する部分がある（レポートの考察を参照）。
- 認証情報は使用していない。環境変数の設定は不要。

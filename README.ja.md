### 🇯🇵 Japanese (`README.ja.md`)

# オイラーの弾性曲線 & 非線形梁アナライザー
[English](README.md) | [Français](README.fr.md) |  [Español](README.es.md) |[繁體中文](README.zh-TW.md) | [简体中文](README.zh-CN.md) | [Deutsch](README.de.md) | [日本語](README.ja.md) | [한국어](README.ko.md)
<div align="center">
  <img src="results/batch_runs/3d_renders/3D_Stress_FEM_Case8_F-20.0N.jpg" width="800" alt="非線形梁の3Dアイソメトリックレンダー">
  <p><em>極端な中央点荷重下での単純支持梁の大変形（Corotational FEM）。</em></p>
</div>

**大変形幾何学的非線形梁（オイラーの弾性曲線）** を解析、シミュレーション、および検証するために設計された、堅牢性の高いPythonベースの計算力学フレームワークです。

3つの異なる数学的側面（解析的、ルンゲ＝クッタ・シューティング法、およびCorotational FEM）を相互検証することにより、このツールは線形仮定と非線形現実との間の正確な数値的境界を明示的にマッピングします。

## 🧮 理論的基礎

従来の線形力学（オイラー・ベルヌーイの梁理論）は、微小回転（$w' \approx 0$）を仮定して正確な曲率方程式を単純化するため、極端な荷重下でのたわみを大幅に過大評価する原因となります。このフレームワークは、完全な幾何学的非線形性を維持することで、**オイラーの弾性曲線**問題を根本的に解決します：

$$ \kappa = \frac{M(x)}{EI} = \frac{w''}{(1 + (w')^2)^{3/2}} $$

この高度に連成された非線形微分方程式は、ほとんどの複雑な荷重ケースに対して閉形式の解析解が存在しないため、本プロジェクトで実装された高度な数値ソルバーが必要となります。

## ⚙️ 技術詳細とアルゴリズム

このフレームワークは、3つの独立したソルバーエンジンと高度な空間マッピングを通じて、高忠実度の検証を実現します：

### 1. 三重の数学的相互検証とラグランジュ座標系での整合
* **線形ベースライン：** 古典的なオイラー・ベルヌーイの微小たわみ理論に基づく閉形式解。これは非線形発散のしきい値を定量化するための比較基準として機能します。
* **厳密な非線形エンジン：** オイラーの弾性曲線の厳密な微分方程式のルンゲ＝クッタ積分。
* **FEMのゴールドスタンダード：** `Corotational`幾何学的非線形定式化を利用したOpenSeesPy。
* **原子レベルでのラグランジュ整合：** 極端な曲げにより、梁の水平投影は大幅に縮小します。すべての誤差評価と可視化は、厳密に**ラグランジュ座標系**（初期弧長 $s$ に基づき、`scipy.interpolate.interp1d` を介して）で実行されます。これにより、オイラー座標系での極端な幾何学的変形に固有の空間的な不整合問題が効果的に解決されます。

### 2. RKシューティング法エンジン
* **状態空間定式化：** 高階微分方程式を1階常微分方程式の系に変換し、状態ベクトル $\mathbf{y} =[w, \theta, M, V, N]^T$ を確立します。
* **境界値問題(BVP)から初期値問題(IVP)への変換：** シューティング法を介して境界値問題(BVP)を解きます。**レーベンバーグ・マルカート法 (`lm`)** (`scipy.optimize.least_squares`) を使用して初期推定値を反復的に精密化し、深い曲げ領域でのヤコビ行列の特異性を効果的に防止します。
* **不連続性の処理：** 区分的に連続な境界適合条件を実装し、内部の力/モーメントのジャンプをスムーズに解決します。

### 3. 自動境界探索エンジン
* **求根アルゴリズム：** **ブレント法** (`scipy.optimize.brentq`) を統合し、線形解析モデルと非線形FEMモデルの間の相対誤差が厳密な**5%のしきい値**に達する正確な適用荷重（または荷重位置）を動的に探索します。

### 4. ゼロ依存の高忠実度3Dレンダリング
* 重い3D科学可視化ライブラリ（VTK、Mayavi、ParaViewなど）を回避します。1D梁要素を数学的に3D物理ソリッドに再構築し、等価応力テンソルを表面にマッピングすることで、**純粋なMatplotlib**を使用して出版品質の3Dアイソメトリック応力コンターを生成します。

## 📦 環境設定

このプロジェクトは、計算用のPythonスクリプトのコレクションです。依存関係が満たされれば、直接実行できます。

### オプションA：標準セットアップ (Windows / macOS / Debian系 / Red Hat系)
標準的なOS環境では、仮想環境の使用は任意ですが、依存関係の競合を避けるために推奨されます。
```bash
# 1. リポジトリをクローンする
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. (任意) 仮想環境を作成して有効化する
python3 -m venv venv
source venv/bin/activate  # Windowsでは次を使用: venv\Scripts\activate

# 3. 主要な依存関係をインストールする
pip install --upgrade pip
pip install -r requirements.txt
```
> **⚠️ macOS Apple Silicon (M1/M2/M3) ユーザーへの注意：** `openseespy` はC++でラップされたフレームワークです。`pip` が互換性のあるARM64ホイールを見つけられない場合、Rosetta 2を使用してターミナルを実行するか、[OpenSeesPyの公式インストール手順](https://openseespydoc.readthedocs.io/en/latest/src/installation.html)を参照する必要があるかもしれません。

### オプションB：Arch系Linux (AUR)
Arch系Linuxを使用している場合、グローバルな `pip` インストールは外部で管理されています (PEP 668)。仮想環境を安全にスキップし、システムのパッケージマネージャーを介して直接依存関係をインストールできます。
*(注意：AURパッケージ `python-openseespy` は、このリポジトリの作者 [@LwhJesse](https://aur.archlinux.org/packages/python-openseespy) によって公式にメンテナンスされています)。*
```bash
# 1. リポジトリをクローンする
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. pacmanとAURヘルパー (例: yay) を介して依存関係をインストールする
sudo pacman -S python-numpy python-pandas python-scipy python-matplotlib python-rich
yay -S python-openseespy
```

## 🚀 使用ガイド

プロジェクトは `core/config.py` にある単一のパラメータ設定ファイルに依存しています。

### 1. 単一ケース分析
単一のケースを実行して、2Dラグランジュ比較プロットを生成します：
```bash
python run_single.py
```
現在の設定で3Dアイソメトリック応力コンターを生成します：
```bash
python run_3d_render.py
```

### 2. バッチ実行とベンチマーク
10個の事前定義されたベンチマークケースをすべて実行して、包括的な検証スイートを生成します：
```bash
python run_batch.py
python run_batch_3d.py
```
<details>
<summary><b>クリックして表示：厳密RK vs FEM多物理整合 (ケース8)</b></summary>
<br>
<div align="center">
  <img src="results/batch_runs/2d_plots/case8-1_comparison.jpg" width="800" alt="2D比較プロット">
  <p><em>線形解析解（緑）が大幅に発散しているのに対し、RK（赤）とFEM（黒）が完全に一致していることに注目してください。</em></p>
</div>
</details>

### 3. 臨界境界探索 (ライブクラスター)
`rich` を利用したライブターミナルダッシュボードを使用して、すべてのCPUコアで同時に5%の非線形誤差位相境界をマッピングします：
```bash
python run_multiprocess.py
```
<details>
<summary><b>クリックして表示：二変数安全/危険ゾーンマップ (ケース3)</b></summary>
<br>
<div align="center">
  <img src="results/boundary_analysis/pure_critical_boundary_fem_case3.jpg" width="600" alt="境界包絡線図">
  <p><em>赤い臨界境界は、幾何学的非線形モデルを採用すべき時期を明示的に示しています。</em></p>
</div>
</details>

## 🛠️ サポートされているベンチマーク荷重ケース

| ケースID | 境界条件 | 荷重タイプ | アクティブ変数 |
| :---: | :--- | :--- | :--- |
| **1** | 片持ち梁 | 端部曲げモーメント | $M_e$ |
| **2** | 片持ち梁 | 端部集中荷重 | $F$ |
| **3** | 片持ち梁 | 中間集中荷重 | $F, a$ |
| **4** | 片持ち梁 | 等分布荷重 | $q$ |
| **5** | 単純支持 | 左端モーメント | $M_e$ |
| **6** | 単純支持 | 右端モーメント | $M_e$ |
| **7** | 単純支持 | 中間モーメント | $M_e, a$ |
| **8** | 単純支持 | 中点集中荷重 | $F$ |
| **9** | 単純支持 | 中間集中荷重 | $F, a$ |
| **10** | 単純支持 | 等分布荷重 | $q$ |

## 📁 出力ディレクトリ構造
結果はソースコードから完全に分離され、自動的に整理されます：
- `results/single_runs/` : 手動設定のスナップショットとレンダー。
- `results/batch_runs/` : 10個のベンチマークケースの完全な検証スイート。
- `results/boundary_analysis/` : 動的な誤差の進化と臨界しきい値曲線。
- `results/mesh_convergence/` : メッシュ独立性の検証プロット。

## 📄 ライセンス
このプロジェクトはMITライセンスの下でライセンスされています - 詳細は[LICENSE](LICENSE)ファイルをご覧ください。


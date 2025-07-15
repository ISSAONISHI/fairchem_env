# fairchem_env 仮想環境構築マニュアル（日本語）

このドキュメントは [FAIRChem (facebookresearch/fairchem)](https://github.com/facebookresearch/fairchem) のコードを Python 3.12 環境で実行するための仮想環境 `fairchem_env` の構築手順をまとめたものです。

---

## ✓ 前提条件

* Anaconda または Miniconda がインストールされていること
* Windows 環境を想定

---

## 1. 仮想環境の作成

GitHub リポジトリ内の `fairchem_env.yml` を使用して仮想環境を作成します。

```bash
conda env create -f fairchem_env.yml
conda activate fairchem_env
```

---

## 2. HuggingFace 認証 (モデル利用時)

FAIRChem の UMA モデル (例: `uma-s-1p1`) を使用するには、HuggingFace のアカウントにログインしてアクセストークンを登録する必要があります。

```bash
huggingface-cli login
```

以下のリンクから UMA モデルへのアクセス申請を行ってください:
[https://huggingface.co/facebook/UMA](https://huggingface.co/facebook/UMA)

---

## 3. UMA モデルの実行例

Cu 表面上の CO 吸着構造を伸縮する例:

```python
from ase.build import fcc100, add_adsorbate, molecule
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator

slab = fcc100("Cu", (3, 3, 3), vacuum=8)
add_adsorbate(slab, molecule("CO"), height=2.0, position="bridge")

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
slab.calc = FAIRChemCalculator(predictor, task_name="oc20")

opt = LBFGS(slab)
opt.run(fmax=0.05, steps=100)
```

---

## 4. Spyder (任意)

GUI で実行する場合、Spyder を conda-forge から最新バージョンで設置します:

```bash
conda install -c conda-forge spyder
```

---

## 5. Visual Studio Build Tools (必要な場合)

パッケージによっては C++ コンパイラが必要な場合があります。

1. [Visual Studio Installer](https://visualstudio.microsoft.com/ja/visual-cpp-build-tools/) をダウンロード
2. 次のワークロードを選択:

   * ☑ C++ によるデスクトップ開発
3. 完了後、Anaconda Prompt を再起動

---

## 6. 備考

* 環境が壊れた場合は前て削除して作り直すこと

```bash
conda remove --name fairchem_env --all
```

* 環境再現性を保つため、下記のようにエクスポートすること

```bash
conda list --explicit > fairchem_env.txt
```




# fairchem\_env Virtual Environment Setup Manual (English)

This document summarizes the setup procedure for the `fairchem_env` virtual environment to run code from [FAIRChem (facebookresearch/fairchem)](https://github.com/facebookresearch/fairchem) using Python 3.12.

---

## ✓ Prerequisites

* Anaconda or Miniconda must be installed
* Assumes Windows OS

---

## 1. Create the Environment

Use the `fairchem_env.yml` file in the GitHub repository to create the environment.

```bash
conda env create -f fairchem_env.yml
conda activate fairchem_env
```

---

## 2. HuggingFace Authentication (for UMA model access)

To use FAIRChem's UMA models (e.g., `uma-s-1p1`), you need to log in to your HuggingFace account and register an access token.

```bash
huggingface-cli login
```

Apply for access to the UMA models at the following URL:
[https://huggingface.co/facebook/UMA](https://huggingface.co/facebook/UMA)

---

## 3. Example: Using UMA Model

Example of relaxing CO adsorption structure on Cu surface:

```python
from ase.build import fcc100, add_adsorbate, molecule
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator

slab = fcc100("Cu", (3, 3, 3), vacuum=8)
add_adsorbate(slab, molecule("CO"), height=2.0, position="bridge")

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
slab.calc = FAIRChemCalculator(predictor, task_name="oc20")

opt = LBFGS(slab)
opt.run(fmax=0.05, steps=100)
```

---

## 4. Spyder (Optional)

If you wish to use a GUI, install Spyder from conda-forge with the latest version:

```bash
conda install -c conda-forge spyder
```

---

## 5. Visual Studio Build Tools (if required)

Some packages may require a C++ compiler. Follow the steps below:

1. Download the [Visual Studio Installer](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
2. Select the following workload:

   * ✓ Desktop development with C++
3. After installation, restart Anaconda Prompt

---

## 6. Notes

* If the environment breaks, remove and recreate it:

```bash
conda remove --name fairchem_env --all
```

* To preserve reproducibility, export the environment as follows:

```bash
conda list --explicit > fairchem_env.txt
```

To restore the environment:

```bash
conda create --name fairchem_env --file fairchem_env.txt
```


復元方法:

```bash
conda create --name fairchem_env --file fairchem_env.txt
```

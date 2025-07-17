# fairchem_env 仮想環境構築マニュアル（日本語）

このマニュアルでは、[FAIRChem (facebookresearch/fairchem)](https://github.com/facebookresearch/fairchem) の UMA モデルを Windows 環境で利用するための仮想環境 `fairchem_env` の構築手順を説明します。

---

## ✓ 前提条件

- Windows 10 または 11（64-bit）
- Anaconda または Miniconda がインストール済み
- Git がインストールされており `git --version` が通ること（URL：https://git-scm.com/downloads/win）
- GPU 環境推奨（CUDA 12.1）

---

## 1. 仮想環境の作成とパッケージ導入

```bash
conda create -n fairchem_env python=3.12
conda activate fairchem_env

# PyTorch + CUDA（適宜バージョン調整）
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# 依存パッケージ
pip install git+https://github.com/facebookresearch/e3nn.git
pip install omegaconf hydra-core einops torchtnt pymatgen huggingface_hub tqdm
```

---

## 2. fairchem のクローンとパス設定

```bash
git clone https://github.com/facebookresearch/fairchem.git
```

Python スクリプトの冒頭に以下を追加：

```python
import sys
sys.path.append("C:/Users/ユーザー名/fairchem/src")  # ← 自分の環境に合わせて修正
```

---

## 3. HuggingFace 認証と UMA モデルのアクセス

```bash
huggingface-cli login
```

以下ページで UMA モデルの使用申請を行う：
[https://huggingface.co/facebook/UMA](https://huggingface.co/facebook/UMA)

承認後、以下のようにスクリプトにトークンを設定：

```python
from huggingface_hub import login
login("hf_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
```

---

## 4. UMA モデルの使用例（Cu 表面への CO 吸着構造最適化）

```python
from ase.build import fcc100, add_adsorbate, molecule
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator

slab = fcc100("Cu", (3, 3, 3), vacuum=8)
add_adsorbate(slab, molecule("CO"), height=2.0, position="bridge")

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda") # CPUの場合は "cpu"
slab.calc = FAIRChemCalculator(predictor, task_name="oc20")

opt = LBFGS(slab)
opt.run(fmax=0.05, steps=100)
```

---

## 5. オプション（GUI 実行用）

Spyder を導入する場合は以下を実行：

```bash
conda install -c conda-forge spyder
```

---

## 6. 環境のエクスポート／復元

### エクスポート：

```bash
conda list --explicit > fairchem_env.txt
```

### 復元：

```bash
conda create --name fairchem_env --file fairchem_env.txt
```

---

## 7. トラブル時の削除：

```bash
conda remove --name fairchem_env --all
```

---

# fairchem_env Virtual Environment Setup Manual (English)

This manual describes how to set up the `fairchem_env` virtual environment to use [FAIRChem (facebookresearch/fairchem)](https://github.com/facebookresearch/fairchem) UMA models on Windows.

---

## ✓ Prerequisites

- Windows 10 or 11 (64-bit)
- Anaconda or Miniconda installed
- Git installed and `git --version` works (URL：https://git-scm.com/downloads/win)
- GPU environment recommended (CUDA 12.1)

---

## 1. Create Environment and Install Dependencies

```bash
conda create -n fairchem_env python=3.12
conda activate fairchem_env

# PyTorch with CUDA support (adjust version if needed)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Additional dependencies
pip install git+https://github.com/facebookresearch/e3nn.git
pip install omegaconf hydra-core einops torchtnt pymatgen huggingface_hub tqdm
```

---

## 2. Clone fairchem and Set Python Path

```bash
git clone https://github.com/facebookresearch/fairchem.git
```

Add the following to your script:

```python
import sys
sys.path.append("C:/Users/YourUsername/fairchem/src")  # ← Adjust to your environment
```

---

## 3. HuggingFace Authentication and UMA Access

```bash
huggingface-cli login
```

Apply for UMA model access at:
[https://huggingface.co/facebook/UMA](https://huggingface.co/facebook/UMA)

After approval, use your token as:

```python
from huggingface_hub import login
login("hf_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
```

---

## 4. Example: CO Adsorption Relaxation on Cu Surface

```python
from ase.build import fcc100, add_adsorbate, molecule
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip, FAIRChemCalculator

slab = fcc100("Cu", (3, 3, 3), vacuum=8)
add_adsorbate(slab, molecule("CO"), height=2.0, position="bridge")

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda") # Use "cpu" if no GPU
slab.calc = FAIRChemCalculator(predictor, task_name="oc20")

opt = LBFGS(slab)
opt.run(fmax=0.05, steps=100)
```

---

## 5. (Optional) Spyder GUI

```bash
conda install -c conda-forge spyder
```

---

## 6. Export/Restore Environment

### Export:

```bash
conda list --explicit > fairchem_env.txt
```

### Restore:

```bash
conda create --name fairchem_env --file fairchem_env.txt
```

---

## 7. Remove Environment

```bash
conda remove --name fairchem_env --all
```

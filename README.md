![Banner GIF](img/banner.gif)

## Environment Setup

```bash
git clone git@github.com:sujitdebnath/fau-bionets-project-ws23.git
cd fau-bionets-project-ws23

python3 -m venv <env_name>
source <env_name>/bin/activate

pip install --upgrade pip

pip install pandas numpy scipy scikit-learn seaborn matplotlib jupyter openpyxl scanpy anndata leidenalg louvain

#-------------
# CPU only
pip install torch===2.0.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install torch_geometric
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cpu.html

pip install -U omicverse
#-------------

deactivate
rm -rf <env_name>
```

MacOS Silicon Based

```bash
conda create -n bionets python=3.10
conda activate bionets

conda install -c conda-forge pandas numpy scipy scikit-learn seaborn matplotlib jupyterlab scanpy anndata pymde python-igraph leidenalg
pip install louvain

# https://omicverse.readthedocs.io/en/latest/Installation_guild/
conda install pytorch torchvision torchaudio cpuonly -c pytorch
pip install torch_geometric
conda install s_gd2 -c conda-forge
pip install -U omicverse

conda deactivate
conda remove -n bionets --all
```
> _**Disclaimer:** Necessary data, documents, and implemented pipeline for the BioNets Project are the intellectual property of [Prof. Dr. David B. Blumenthal](https://www.bionets.tf.fau.de/person/david-b-blumenthal/), and [Dr. Anne Hartebrodt](https://www.bionets.tf.fau.de/person/anne-hartebrodt/) at [FAU Erlangen-Nürnberg](https://www.fau.eu/). Please be aware that copying content from here holds you accountable._

![Banner GIF](img/banner.gif)

Welcome to the Biomedical Network Science (BioNets) Project repository for the Winter'23/24 semester at [Friedrich-Alexander University Erlangen-Nürnberg](https://www.fau.eu/). This repository contains necessary data, documents, and implemented pipeline for the BioNets Project, proposed by Biomedical Network Science ([BIONETS](https://www.bionets.tf.fau.de/)) lab, supervised by [Prof. Dr. David B. Blumenthal](https://www.bionets.tf.fau.de/person/david-b-blumenthal/), and [Dr. Anne Hartebrodt](https://www.bionets.tf.fau.de/person/anne-hartebrodt/) at FAU Erlangen-Nürnberg.

## Project Goals
The core goals of the project are as follows:

- Retrieve scRNA-seq data from case-control studies for one fixed disease.
- Implemented and run pipelines for
    - Clustering and automatic cell type annotation.
    - Identification of DEGs, comparing cells from the same cell type between case and control.
- Make a dashboard to interactively visualize the results.
- (Possibly) Extend to > 1 diseases.

This project consists of _**two**_ parts:

1. scRNA-seq Analysis: Which consists of disease research, relevant data collection, automatic cell type annotation, and differentially expressed genes (DEGs) analysis.
2. Dashboard: Preparing a web-based template to demonstrate the finding which should be interactive and easy to use.

## Project Structure

```bash
fau-bionet-project-ws23/
├── dataset/                  # Data directory
│   ├── disease_id1/          # Disease Name 1, e.g. Diabetes II
│   │   ├── dataset_id1/      # Dataset 1 for disease 1
│   │   │   ├── case/         # Case samples
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── control/      # Control samples
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── dataset_id2/
│   │   └── ...
│   ├── disease_id2/          # Disease Name 2, e.g. MPN
│   └── ...
├── pipelines/                # Main pipeline modules
│   ├── services/             # Pipeline services
│   │   ├── adata_handler.py
│   │   ├── adata_preprocessor.py
│   │   ├── cell_type_annotation.py
│   │   └── diff_gene_exp_analysis.py
│   ├── run_pipelines.sh      # Script to run pipelines
│   └── cleaner.sh            # Script to clean temporary files
├── results/                  # Results directory
│   ├── disease_id1/          # Store results for Disease 1
│   │   ├── dataset_id1/      # Store results of dataset1
│   │   ├── dataset_id2/      # Store results of dataset2
│   │   └── ...
│   ├── disease_id2/          # Store results for Disease 2
│   └── ...
├── dashboard/                # Dashboard directory
└── README.md                 # Project documentation
```

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
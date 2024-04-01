> _**Disclaimer:** Necessary data, documents, and implemented pipeline for the BioNets Project are the intellectual property of [Prof. Dr. David B. Blumenthal](https://www.bionets.tf.fau.de/person/david-b-blumenthal/), and [Dr. Anne Hartebrodt](https://www.bionets.tf.fau.de/person/anne-hartebrodt/) at [FAU Erlangen-Nürnberg](https://www.fau.eu/). Please be aware that copying content from here holds you accountable._

![Banner GIF](img/banner.gif)

Welcome to the Biomedical Network Science (BioNets) Project repository for the Winter'23/24 semester at [Friedrich-Alexander University Erlangen-Nürnberg](https://www.fau.eu/). This repository contains necessary files and documents for the BioNets project called _**"Large-scale Differential Gene Expression Analysis in scRNA-seq Data"**_, proposed by Biomedical Network Science ([BIONETS](https://www.bionets.tf.fau.de/)) lab, supervised by [Prof. Dr. David B. Blumenthal](https://www.bionets.tf.fau.de/person/david-b-blumenthal/), and [Dr. Anne Hartebrodt](https://www.bionets.tf.fau.de/person/anne-hartebrodt/) at FAU Erlangen-Nürnberg.

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
├── dataset/                            # Data directory
│   ├── disease_id1/                    # Disease Name 1, e.g. Diabetes II
│   │   ├── dataset_id1/                # Dataset 1 for disease 1
│   │   │   ├── case/                   # Case samples
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── control/                # Control samples
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── dataset_id2/
│   │   └── ...
│   ├── disease_id2/                    # Disease Name 2, e.g. MPN
│   └── ...
├── pipelines/                          # Main pipeline modules
│   ├── services/                       # Pipeline services
│   │   ├── adata_handler.py            # Script for data loading, saving, etc
│   │   ├── adata_preprocessor.py       # Script for preprocessing
│   │   ├── cell_type_annotation.py     # Script for automatic cell annotation
│   │   └── diff_gene_exp_analysis.py   # Script for DGE analysis
│   ├── run_pipelines.sh                # Script to run pipelines
│   └── cleaner.sh                      # Script to clean temporary files
├── results/                            # Results directory
│   ├── disease_id1/                    # Store results for Disease 1
│   │   ├── dataset_id1/                # Store results of dataset1
│   │   ├── dataset_id2/                # Store results of dataset2
│   │   └── ...
│   ├── disease_id2/                    # Store results for Disease 2
│   └── ...
├── dashboard/                          # Dashboard directory
└── README.md                           # Project documentation
```

## Project Details

### 1. Dataset

### 2. Pipelines

#### 2.1. Preprocessing adata:
#### 2.2. Automatic Cell-type Annotation:
#### 2.3. DGE Analysis:

### 3. Dashboard

## Environment Setup

This section provides step-by-step instructions for setting up the required environment on Linux, and MacOS systems. Please note that the setup process for MacOS systems with Silicon-based processors may vary slightly.

### Used Technology
1. [Python3.x](https://www.python.org), and [Anaconda](https://anaconda.org) distribution (for Silicon Based MacOS)
2. [Scanpy](https://scanpy.readthedocs.io/en/stable/) - it is a Python package and a scalable toolkit for analyzing single-cell gene expression data built jointly with [anndata](https://anndata.readthedocs.io/en/latest/).
3. [Omicverse](https://omicverse.readthedocs.io/en/latest/) - Omicverse is the fundamental package for multi omics included bulk and single cell RNA-seq analysis with Python.

### Linux

```bash
# Clone the repository
git clone git@github.com:sujitdebnath/fau-bionet-project-ws23.git
cd fau-bionets-project-ws23

# Create a virtual environment and activate
python3 -m venv <env_name>
source <env_name>/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install required Python packages
pip install pandas numpy scipy scikit-learn seaborn matplotlib jupyter openpyxl scanpy anndata leidenalg louvain

# Install PyTorch, PyTorch Geometric and additional packages for CPU-only operations
pip install torch===2.0.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install torch_geometric
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cpu.html

# Install Omicverse package
pip install -U omicverse

# [if needed] Deactivate and remove virtual environment
deactivate
rm -rf <env_name>
```

### MacOS (Silicon Based)

```bash
# Clone the repository
git clone git@github.com:sujitdebnath/fau-bionet-project-ws23.git
cd fau-bionets-project-ws23

# Create a conda environment and activate
conda create -n <conda_env_name> python=<python_version>
conda activate <conda_env_name>

# Install required Python packages
conda install -c conda-forge pandas numpy scipy scikit-learn seaborn matplotlib jupyterlab scanpy anndata pymde python-igraph leidenalg
pip install louvain

# Install PyTorch, PyTorch Geometric and additional packages for CPU-only operations
conda install pytorch torchvision torchaudio cpuonly -c pytorch
pip install torch_geometric
conda install s_gd2 -c conda-forge

# Install Omicverse package
pip install -U omicverse

# [if needed] Deactivate and remove conda environment
conda deactivate
conda remove -n bionets --all
```
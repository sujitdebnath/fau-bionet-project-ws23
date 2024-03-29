import os
import sys
import argparse
import anndata as ad
import scanpy as sc
import omicverse as ov
from typing import List, Callable, Any, Optional


BASE_DIR      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE_DATA_DIR = os.path.join(BASE_DIR, 'dataset')
BASE_FIG_DIR  = os.path.join(BASE_DIR, 'results', 'figures')
BASE_DEG_DIR  = os.path.join(BASE_DIR, 'results', 'deg_results')


def parse_add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument('--disease_id', type=str, help='Unique identifier for disease, e.g., diabetesII, mpn')
    parser.add_argument('--dataset_id', type=str, help='Unique identifier for the dataset')

    return parser

def create_directory(dir_path: str) -> str:
    try:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(f"Succeed: Directory created at {dir_path}")
        else:
            print(f"Succeed: Directory already exists at {dir_path}")
    except Exception as e:
        print(f"Failed: Error creating directory at {dir_path}. {e}")
        sys.exit(1)
    
    return dir_path

def load_data_case_and_control(dataset_path: str) -> sc.AnnData:
    case_dir_path    = None
    control_dir_path = None

    for dir_name in os.listdir(dataset_path):
        if dir_name.startswith('case'):
            case_dir_path = os.path.join(dataset_path, dir_name)
        elif dir_name.startswith('control'):
            control_dir_path = os.path.join(dataset_path, dir_name)

    if case_dir_path is not None:
        adata_case              = sc.read_10x_mtx(case_dir_path, var_names='gene_symbols', cache=True)
        adata_case.obs['donor'] = 'case'
    
    if control_dir_path is not None:
        adata_control              = sc.read_10x_mtx(control_dir_path, var_names='gene_symbols', cache=True)
        adata_control.obs['donor'] = 'control'

    if case_dir_path is not None and control_dir_path is not None:
        adata = ad.concat([adata_case, adata_control])
    elif case_dir_path is not None:
        adata = adata_case
    elif control_dir_path is not None:
        adata = adata_control
    
    # adata.obs_names_make_unique()
    
    return adata

def run_data_handler() -> None:
    parser = argparse.ArgumentParser(description='Data Handler Script')
    parser = parse_add_args(parser)
    args   = parser.parse_args()
    
    disease_id = args.disease_id
    dataset_id = args.dataset_id

    dataset_dir     = os.path.join(BASE_DATA_DIR, disease_id, dataset_id)
    disease_fig_dir = create_directory(dir_path=os.path.join(BASE_FIG_DIR, disease_id))
    dataset_fig_dir = create_directory(dir_path=os.path.join(disease_fig_dir, dataset_id))
    disease_deg_dir = create_directory(dir_path=os.path.join(BASE_DEG_DIR, disease_id))
    dataset_deg_dir = create_directory(dir_path=os.path.join(disease_deg_dir, dataset_id))

    adata = load_data_case_and_control(dataset_path=dataset_dir)
    print(adata)


if __name__ == '__main__':
    run_data_handler()

import os
import sys
import argparse
import anndata as ad
import scanpy as sc
import omicverse as ov
from typing import List, Callable, Any, Optional


BASE_DIR      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE_DATA_DIR = os.path.join(BASE_DIR, 'dataset')
BASE_RES_DIR  = os.path.join(BASE_DIR, 'results')


def parse_add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument('--disease_id', type=str, help='Unique identifier for disease, e.g., diabetesII, mpn')
    parser.add_argument('--dataset_id', type=str, help='Unique identifier of the dataset for that specific disease')

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

def prepare_all_dirs(disease_id: str, dataset_id: str) -> None:
    # create fig and deg_res dir
    fig_dir     = create_directory(dir_path=os.path.join(BASE_RES_DIR, 'figures'))
    deg_res_dir = create_directory(dir_path=os.path.join(BASE_RES_DIR, 'deg_results'))

    # create disease dir
    disease_fig_dir = create_directory(dir_path=os.path.join(fig_dir, disease_id))
    disease_deg_dir = create_directory(dir_path=os.path.join(deg_res_dir, disease_id))

    # create fig and deg dirs for specific dataset of specific disease
    dataset_fig_dir = create_directory(dir_path=os.path.join(disease_fig_dir, dataset_id))
    dataset_deg_dir = create_directory(dir_path=os.path.join(disease_deg_dir, dataset_id))

    # create temp dirs to store temp adata in pipelines
    temp_adata_dir = create_directory(dir_path=os.path.join(BASE_DIR, 'pipelines', 'temp_adata'))

def load_data_case_and_control(disease_id: str, dataset_id: str) -> sc.AnnData:
    dataset_path = os.path.join(BASE_DIR, 'pipelines', 'temp_adata', f'{disease_id}_{dataset_id}.h5ad')

    if os.path.exists(dataset_path):
        adata = sc.read_h5ad(dataset_path) 
    else:
        dataset_path     = os.path.join(BASE_DATA_DIR, disease_id, dataset_id)
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

def save_adata(adata: sc.AnnData, disease_id: str, dataset_id: str) -> None:
    try:
        temp_adata_dir = os.path.join(BASE_DIR, 'pipelines', 'temp_adata')
        res_fpath      = os.path.join(temp_adata_dir, f'{disease_id}_{dataset_id}.h5ad')
        
        adata.write(res_fpath)
        print(f"Succeed: Successfully saved AnnData object as {os.path.basename(res_fpath)} in the temporary adata dir.")
    except Exception as e:
        print(f"Failed: Error while saving AnnData object. {str(e)}")
        sys.exit(1)

def run_adata_handler() -> None:
    parser = argparse.ArgumentParser(description='Data Handler Script')
    parser = parse_add_args(parser)
    args   = parser.parse_args()
    
    disease_id = args.disease_id
    dataset_id = args.dataset_id

    prepare_all_dirs(disease_id=disease_id, dataset_id=dataset_id)

    adata = load_data_case_and_control(disease_id=disease_id, dataset_id=dataset_id)
    print(adata)

    save_adata(adata=adata, disease_id=disease_id, dataset_id=dataset_id)


if __name__ == '__main__':
    run_adata_handler()

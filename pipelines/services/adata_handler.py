# Python imports
import os
import sys
import argparse

# Third party imports
import anndata as ad
import scanpy as sc
import omicverse as ov

# Self imports


# Define the necessary directories of the project
BASE_DIR      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE_DATA_DIR = os.path.join(BASE_DIR, 'dataset')
BASE_RES_DIR  = os.path.join(BASE_DIR, 'results')


def parse_add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add arguments related to disease and dataset identifiers to the provided argument parser.

    Parameters:
        parser (argparse.ArgumentParser): The argument parser to which the arguments will be added.

    Returns:
        argparse.ArgumentParser: The modified argument parser.
    """
    parser.add_argument('--disease_id', type=str, help='Unique identifier for disease, e.g., diabetesII, mpn')
    parser.add_argument('--dataset_id', type=str, help='Unique identifier of the dataset for that specific disease')

    return parser

def create_directory(dir_path: str) -> str:
    """
    Create a directory at the specified path if it does not exist already.

    Parameters:
        dir_path (str): The path of the directory to be created.

    Returns:
        str: The path of the directory that was created or already existed.
    """
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
    """
    Prepare directories for storing results and temporary data.

    Parameters:
        disease_id (str): Unique identifier for the disease.
        dataset_id (str): Unique identifier of the dataset for the specific disease.

    Returns:
        None
    """
    # create disease and dataset dir in results
    disease_res_dir = create_directory(dir_path=os.path.join(BASE_RES_DIR, disease_id))
    dataset_res_dir = create_directory(dir_path=os.path.join(disease_res_dir, dataset_id))

    # create temp dirs to store temp adata in pipelines
    temp_adata_dir = create_directory(dir_path=os.path.join(BASE_DIR, 'pipelines', 'temp_adata'))

def load_data_case_and_control(disease_id: str, dataset_id: str) -> sc.AnnData:
    """
    Load case and control data for a specific disease and dataset.

    Parameters:
        disease_id (str): Unique identifier for the disease.
        dataset_id (str): Unique identifier of the dataset for the specific disease.

    Returns:
        sc.AnnData: Anndata object containing the loaded data.
    """
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
        
        # make indices unique for the whole adata object
        # as the case vs control dataset is concatenated, the duplicate indices can occur
        adata.obs_names_make_unique()
    
    return adata

def save_adata(adata: sc.AnnData, disease_id: str, dataset_id: str) -> None:
    """
    Save an AnnData object to a file.

    Parameters:
        adata (sc.AnnData): Anndata object to be saved.
        disease_id (str): Unique identifier for the disease.
        dataset_id (str): Unique identifier of the dataset for the specific disease.

    Returns:
        None
    """
    try:
        temp_adata_dir = os.path.join(BASE_DIR, 'pipelines', 'temp_adata')
        res_fpath      = os.path.join(temp_adata_dir, f'{disease_id}_{dataset_id}.h5ad')
        
        adata.write(res_fpath)
        print(f"Succeed: Successfully saved AnnData object as {os.path.basename(res_fpath)} in the temporary adata dir.")
    except Exception as e:
        print(f"Failed: Error while saving AnnData object. {str(e)}")
        sys.exit(1)

def main() -> None:
    """
    Main function to handle data processing tasks.

    This function parses command-line arguments, prepares directories, loads data,
    saves data, and performs other necessary tasks.

    Parameters:
        None

    Returns:
        None
    """
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
    main()

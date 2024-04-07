import os
import pandas as pd
import scanpy as sc


BASE_DIR      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_DATA_DIR = os.path.join(BASE_DIR, 'dataset')
BASE_RES_DIR  = os.path.join(BASE_DIR, 'results')


def load_adata(disease_id: str, dataset_id: str) -> sc.AnnData:
    dataset_path = os.path.join(BASE_DIR, 'pipelines', 'temp_adata', f'{disease_id}_{dataset_id}.h5ad')
    adata        = sc.read_h5ad(dataset_path)
    
    return adata

def combine_all_adata_dfs() -> pd.DataFrame:
    disease_ids   = sorted([dir for dir in os.listdir(BASE_RES_DIR) if os.path.isdir(os.path.join(BASE_RES_DIR, dir))])
    all_adata_dfs = []

    for disease_id in disease_ids:
        dataset_ids = sorted(os.listdir(os.path.join(BASE_RES_DIR, disease_id)))

        for dataset_id in dataset_ids:
            adata = load_adata(disease_id=disease_id, dataset_id=dataset_id)

            adata.obs['disease_id'] = disease_id
            adata.obs['dataset_id'] = dataset_id

            adata_df = pd.DataFrame(adata.obs[['disease_id', 'dataset_id', 'donor', 'scsa_celltype_cellmarker', 'scsa_celltype_panglaodb', 'MetaTiME', 'Major_MetaTiME']])
            print(f"({disease_id}, {dataset_id}) => {adata_df.shape}")
            all_adata_dfs.append(adata_df)
    
    combined_df = pd.concat(all_adata_dfs, ignore_index=True)

    return combined_df

def standardize_cell_type_names(df: pd.DataFrame) -> pd.DataFrame:
    # define mapping for standardizing cell type names
    cell_type_mapping = {
        'T-cell': ['T cell', 'T Cells', 'T', 'T follicular helper(Tfh) cell', 'T Helper Cells', 'CD4+ T cell', 'Naive CD8+ T cell', 'Natural killer T (NKT) cell', 'Gamma Delta T Cells', 'Regulatory T(Treg) cell', 'T Regulatory Cells', 'T Memory Cells'],
        'B-cell': ['B cell', 'B Cells', 'B', 'Plasma Cells'],
        'Monocyte': ['Monocyte', 'Monocytes', 'M', 'Myeloid'],
        'Megakaryocyte': ['Megakaryocyte', 'Platelets'],
        'Natural killer cell': ['Natural killer cell', 'NK Cells'],
        'Endothelial cell': ['Endothelial cell', 'Stroma', 'Decidual Cells'],
        'Hematopoietic stem cell': ['Hematopoietic stem cell', 'Hematopoietic Stem Cells'],
        'Dendritic cell': ['Dendritic cell', 'Dendritic Cells', 'DC', 'Plasmacytoid dendritic cell(pDC)', 'Plasmacytoid Dendritic Cells'],
        'Pancreatic cell': ['Pan'],
        'Unknown': ['Unknown']
    }

    # apply the mapping to all columns containing cell type names
    for col in df.columns:
        for cell_type, aliases in cell_type_mapping.items():
            df[col] = df[col].replace(aliases, cell_type)

    return df

def main():
    adata_dfs             = combine_all_adata_dfs()
    adata_standardize_dfs = standardize_cell_type_names(adata_dfs)
    print(f"Combined df shape: {adata_dfs.shape}")
    adata_standardize_dfs.to_csv(os.path.join(BASE_DIR, 'dashboard', 'cell_anno_res.csv'))


if __name__ == '__main__':
    main()

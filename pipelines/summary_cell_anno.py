# Python imports
import os
import sys
from typing import List

# Third-party imports
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# Self imports


# Define the necessary directories of the project
BASE_DIR     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_RES_DIR = os.path.join(BASE_DIR, 'results')


def create_directory(dir_path: str) -> str:
    """
    Create a directory if it does not exist.

    Parameters:
        dir_path (str): The path of the directory to be created.

    Returns:
        str: The path of the created directory.
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

def load_adata(disease_id: str, dataset_id: str) -> sc.AnnData:
    """
    Load an AnnData object from a file.

    Parameters:
        disease_id (str): Unique identifier for the disease.
        dataset_id (str): Unique identifier for the dataset.

    Returns:
        sc.AnnData: The loaded AnnData object.
    """
    dataset_path = os.path.join(BASE_DIR, 'pipelines', 'temp_adata', f'{disease_id}_{dataset_id}.h5ad')
    adata        = sc.read_h5ad(dataset_path)
    
    return adata

def combine_all_adata_dfs() -> pd.DataFrame:
    """
    Combine all AnnData objects into a single DataFrame.

    Parameters:
        None
    
    Returns:
        pd.DataFrame: The combined DataFrame.
    """
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

def standardize_cell_type_names(df: pd.DataFrame, cell_type_mapping: dict) -> pd.DataFrame:
    """
    Standardize cell type names in a DataFrame using a given mapping.

    Parameters:
        df (pd.DataFrame): DataFrame containing cell type names to be standardized.
        cell_type_mapping (dict): Dictionary mapping original cell type names to standardized names.

    Returns:
        pd.DataFrame: DataFrame with standardized cell type names.
    """
    # apply the mapping to all columns containing cell type names
    for col in df.columns:
        for cell_type, aliases in cell_type_mapping.items():
            df[col] = df[col].replace(aliases, cell_type)

    return df

def map_cell_anno_method_to_col_name(cell_anno_method: str) -> str:
    """
    Map the cell type annotation method to the corresponding column name.

    Parameters:
        cell_anno_method (str): The name of the cell type annotation method.

    Returns:
        str: The corresponding column name in the AnnData object.
    """
    if cell_anno_method == "SCSA - cellmarker":
        cell_anno_col_name = "scsa_celltype_cellmarker"
    elif cell_anno_method == "SCSA - panglaodb":
        cell_anno_col_name = "scsa_celltype_panglaodb"
    elif cell_anno_method == "MetaTiME":
        cell_anno_col_name = "Major_MetaTiME"
    
    return cell_anno_col_name

def map_col_name_to_cell_anno_method(cell_anno_col: str) -> str:
    """
    Map the column name in the AnnData object to the corresponding cell type annotation method.

    Parameters:
        cell_anno_col (str): The name of the column in the AnnData object.

    Returns:
        str: The corresponding cell type annotation method.
    """
    if cell_anno_col == "scsa_celltype_cellmarker":
        cell_anno_method = "SCSA - cellmarker"
    elif cell_anno_col == "scsa_celltype_panglaodb":
        cell_anno_method = "SCSA - panglaodb"
    elif cell_anno_col == "Major_MetaTiME":
        cell_anno_method = "MetaTiME"
    
    return cell_anno_method

def cell_anno_box_plot_case_control_subplots(
        df: pd.DataFrame,
        disease_ids: List,
        cell_anno_methods: List,
        figure_path: str
    ) -> None:
    """
    Generate subplots of box plots comparing the distribution of cell types between case and control donors.

    Parameters:
        df (pd.DataFrame): DataFrame containing the data.
        disease_ids (List): List of disease IDs.
        cell_anno_methods (List): List of cell annotation methods.
        figure_path (str): Path to save the generated figure.

    Returns:
        None
    """
    num_diseases = len(disease_ids)
    num_methods  = len(cell_anno_methods)
    
    fig, axes = plt.subplots(num_diseases, num_methods, figsize=(32, 16))
    cap_alpha = ord('a') - 1

    for i, disease_id in enumerate(disease_ids):
        for j, cell_anno_method in enumerate(cell_anno_methods):
            cell_anno_col_name = map_cell_anno_method_to_col_name(cell_anno_method)

            # Aggregate counts of cell types for each disease, dataset, and donor
            agg_df = df.groupby(['disease_id', 'dataset_id', 'donor', cell_anno_col_name]).size().reset_index(name='count')

            # Filter rows with given disease
            agg_df = agg_df[agg_df['disease_id'] == disease_id]

            # Filter out rows with 'Unknown' cell types
            agg_df = agg_df[agg_df[cell_anno_col_name] != 'Unknown']

            # Calculate total count of each cell type across all datasets and diseases
            total_counts = agg_df.groupby(cell_anno_col_name)['count'].sum().reset_index()
            total_counts = total_counts.sort_values(by='count', ascending=False)  # Sort by count

            # Sort agg_df based on total count order
            agg_df[cell_anno_col_name] = pd.Categorical(agg_df[cell_anno_col_name], categories=total_counts[cell_anno_col_name], ordered=True)
            agg_df = agg_df.sort_values(by=cell_anno_col_name)

            # Create a Seaborn box plot
            ax = sns.boxplot(
                data=agg_df,
                x='count',
                y=cell_anno_col_name,
                hue='donor',
                hue_order=['control', 'case'],
                gap=0.2,
                palette='Set2',
                ax=axes[i, j]
            )
            ax.set_title(f'({chr(cap_alpha + 1)}) Disease: {disease_id.upper()} and Method: {cell_anno_method}')
            ax.set_xlabel('Count')
            ax.set_ylabel('Cell Types')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.legend(title='Donor')
            cap_alpha += 1

    plt.tight_layout()
    fig.subplots_adjust(top=0.94)
    plt.suptitle(f"Distribution of Cell Types for Case vs Control", fontsize=16, fontweight='bold')
    plt.savefig(figure_path)
    plt.close()

def cell_anno_box_plot_for_diseases_subplots(
        df: pd.DataFrame,
        disease_ids: List,
        cell_anno_methods: List,
        figure_path: str
    ) -> None:
    """
    Generate subplots of box plots comparing the distribution of cell types across different annotation methods,
    for each disease.

    Parameters:
        df (pd.DataFrame): DataFrame containing the data.
        disease_ids (List): List of disease IDs.
        cell_anno_methods (List): List of cell annotation methods.
        figure_path (str): Path to save the generated figure.

    Returns:
        None
    """
    fig, axes = plt.subplots(1, len(disease_ids), figsize=(20, 6))
    cap_alpha = ord('a')

    for i, disease_id in enumerate(disease_ids):
        merged_df = None

        for cell_anno_method in cell_anno_methods:
            cell_anno_col_name = map_cell_anno_method_to_col_name(cell_anno_method)

            agg_df = df.groupby(['disease_id', 'dataset_id', cell_anno_col_name]).size().reset_index(name='count')
            agg_df = agg_df.rename(columns={cell_anno_col_name: 'cell_type'})
            agg_df['method'] = cell_anno_method

            if merged_df is None:
                merged_df = agg_df
            else:
                merged_df = pd.concat([merged_df, agg_df])

        # Filter rows with given disease
        merged_df = merged_df[merged_df['disease_id'] == disease_id]

        # Calculate total count of each cell type across all annotation methods
        total_counts = merged_df.groupby('cell_type')['count'].sum().reset_index()
        total_counts = total_counts.sort_values(by='count', ascending=False)

        # Sort merged_df based on method order and total count order
        merged_df['cell_type'] = pd.Categorical(merged_df['cell_type'], categories=total_counts['cell_type'], ordered=True)
        merged_df['method'] = pd.Categorical(merged_df['method'], categories=cell_anno_methods, ordered=True)
        merged_df = merged_df.sort_values(by=['method', 'cell_type'])

        # Create a Seaborn box plot
        sns.boxplot(
            data=merged_df,
            x='count',
            y='cell_type',
            hue='method',
            gap=0.25,
            palette='tab20',
            ax=axes[i]
        )

        axes[i].set_title(f'({chr(cap_alpha+i)}) {disease_id.upper()}')
        axes[i].set_xlabel('Count')
        axes[i].set_ylabel('Cell Types')
        axes[i].legend()

    plt.tight_layout()
    fig.subplots_adjust(top=0.9)
    plt.suptitle(f"Distribution of Cell Types According to All Annotation Methods", fontsize=16, fontweight='bold')
    plt.savefig(figure_path)
    plt.close()

def cell_anno_heatmap_subplots(
        df: pd.DataFrame,
        disease_ids: List,
        cell_anno_methods: List,
        figure_path: str
    ) -> None:
    """
    Generate subplots of heatmaps showing the co-occurrence of cell types across different annotation methods,
    for each disease.

    Parameters:
        df (pd.DataFrame): DataFrame containing the data.
        disease_ids (List): List of disease IDs.
        cell_anno_methods (List): List of cell annotation methods.
        figure_path (str): Path to save the generated figure.

    Returns:
        None
    """
    fig, axes = plt.subplots(len(disease_ids), len(cell_anno_methods), figsize=(35, 16))
    cap_alpha = ord('a') - 1

    for xax, disease_id in enumerate(disease_ids):
        yax = 0

        for i, cell_anno_method1 in enumerate(cell_anno_methods):
            for cell_anno_method2 in cell_anno_methods[i+1:]:
                cell_anno_col1 = map_cell_anno_method_to_col_name(cell_anno_method1)
                cell_anno_col2 = map_cell_anno_method_to_col_name(cell_anno_method2)

                # Filter data by disease_id
                df_filtered = df[df['disease_id'] == disease_id]

                # Compute the co-occurrence matrix
                co_occurrence_matrix = df_filtered.groupby([cell_anno_col1, cell_anno_col2]).size().unstack(fill_value=0)

                # Filter out rows and columns with all zero values
                co_occurrence_matrix = co_occurrence_matrix.loc[(co_occurrence_matrix != 0).any(axis=1), (co_occurrence_matrix != 0).any(axis=0)]

                # Create the heatmap
                sns.heatmap(
                    co_occurrence_matrix,
                    annot=True,
                    fmt='d',  # Display integer values
                    linewidths=0.4,
                    linecolor='grey',
                    cmap='inferno_r',
                    annot_kws={"size": 14},
                    ax=axes[xax, yax]
                )

                # Customize the plot
                axes[xax, yax].set_title(f"({chr(cap_alpha + 1)})")
                axes[xax, yax].set_xlabel(cell_anno_method2)
                axes[xax, yax].set_ylabel(cell_anno_method1)
                yax += 1
                cap_alpha += 1
    
    plt.tight_layout()
    fig.subplots_adjust(top=0.94)
    plt.suptitle(f"Cell Type Co-occurrence Heatmap for {', '.join([d.upper() for d in disease_ids])} Diseases", fontsize=16, fontweight='bold')
    plt.savefig(figure_path)
    plt.close()

def main():
    """
    Main function for processing and standardizing cell type annotations.

    Combines cell type annotations from different datasets, standardizes cell type names,
    and saves the results to a CSV file. Additionally, generates and saves visualizations
    depicting the distribution of cell types for case vs control, distribution of cell types
    according to all annotation methods, and a heatmap of cell type co-occurrence.

    Parameters:
        None

    Returns:
        None
    """
    adata_dfs = combine_all_adata_dfs()

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
    adata_standardize_dfs = standardize_cell_type_names(df=adata_dfs, cell_type_mapping=cell_type_mapping)

    summary_dir = create_directory(dir_path=os.path.join(BASE_DIR, 'results', 'summary'))
    adata_standardize_dfs.to_csv(os.path.join(summary_dir, 'cell_anno_res.csv'))
    
    print(f"Combined df shape: {adata_dfs.shape}")
    print(adata_standardize_dfs['scsa_celltype_cellmarker'].unique())
    print(adata_standardize_dfs['scsa_celltype_panglaodb'].unique())
    print(adata_standardize_dfs['Major_MetaTiME'].unique())

    disease_ids       = sorted([dir for dir in os.listdir(BASE_RES_DIR) if os.path.isdir(os.path.join(BASE_RES_DIR, dir)) and dir != 'summary'])
    cell_anno_methods = ['SCSA - cellmarker', 'SCSA - panglaodb', 'MetaTiME']

    # distribution of cell types for case vs control (excluded unknown labels)
    cell_anno_box_plot_case_control_subplots(
        df=adata_standardize_dfs,
        disease_ids=disease_ids,
        cell_anno_methods=cell_anno_methods,
        figure_path=os.path.join(BASE_RES_DIR, 'summary', 'dist_case_control.jpeg')
    )
    
    # distribution of cell types according to all annotation methods
    cell_anno_box_plot_for_diseases_subplots(
        df=adata_standardize_dfs,
        disease_ids=disease_ids,
        cell_anno_methods=cell_anno_methods,
        figure_path=os.path.join(BASE_RES_DIR, 'summary', 'dist_anno_methods.jpeg')
    )

    # heatmap of cell type co-occurrence
    cell_anno_heatmap_subplots(
        df=adata_standardize_dfs,
        disease_ids=disease_ids,
        cell_anno_methods=cell_anno_methods,
        figure_path=os.path.join(BASE_RES_DIR, 'summary', 'heatmap.jpeg')
    )


if __name__ == '__main__':
    main()

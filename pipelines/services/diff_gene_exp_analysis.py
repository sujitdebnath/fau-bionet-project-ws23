# Python imports
import os
import argparse
from typing import List

# Third-party imports
import numpy as np
import pandas as pd
import scanpy as sc
import omicverse as ov

# Self imports
import adata_handler


def parse_add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Parse and add command-line arguments related to differential gene expression analysis.

    Parameters:
        parser (argparse.ArgumentParser): Argument parser object.

    Returns:
        argparse.ArgumentParser: Updated argument parser object.
    """
    parser.add_argument('--disease_id', type=str, help='Unique identifier for disease, e.g., diabetesII, mpn')
    parser.add_argument('--dataset_id', type=str, help='Unique identifier of the dataset for that specific disease')
    parser.add_argument('--dge_methods', nargs='+', type=str,
                        default=['t-test', 'wilcoxon', 'logreg', 't-test_overestim_var'],
                        help='Possible method list of differential gene expression analysis')
    return parser

def perform_diff_gene_exp_analysis(
        adata: sc.AnnData,
        disease_id: str,
        dataset_id: str,
        dge_methods: List,
        specific_celltype: str
    ) -> None:
    """
    Perform differential gene expression analysis.

    Parameters:
        adata (sc.AnnData): Annotated data object.
        disease_id (str): Unique identifier for disease.
        dataset_id (str): Unique identifier of the dataset for that specific disease.
        dge_methods (List): List of differential gene expression analysis methods.
        specific_celltype (str): Specific cell type for analysis.

    Returns:
        None
    """
    all_results = pd.DataFrame()
    
    for method in dge_methods:
        print(f"{'*'*10} Performing {method} {'='*10}")
        
        sub_adata = adata.copy()
    
        if specific_celltype != 'All':
            sub_adata = sub_adata[sub_adata.obs['scsa_celltype_cellmarker'] == specific_celltype]        
    
        sc.tl.rank_genes_groups(sub_adata, groupby='donor', method=method)
    
        result = sub_adata.uns['rank_genes_groups']
        donors = result['names'].dtype.names
        out    = np.zeros((1, 7))
        print(donors)
        
        for donor in donors:
            if 'pvals_adj' in result.keys():
                out = np.vstack((out, np.vstack((
                        result['names'][donor],
                        np.array([specific_celltype] * len(result['names'][donor])).astype('object'),
                        np.array([donor] * len(result['names'][donor])).astype('object'),
                        np.array([result['params']['method']] * len(result['names'][donor])).astype('object'),
                        result['scores'][donor],
                        result['pvals_adj'][donor],
                        result['logfoldchanges'][donor]
                    )).T))
            else:
                out = np.vstack((out, np.vstack((
                        result['names'][donor],
                        np.array([specific_celltype] * len(result['names'][donor])).astype('object'),
                        np.array([donor] * len(result['names'][donor])).astype('object'),
                        np.array([result['params']['method']] * len(result['names'][donor])).astype('object'),
                        result['scores'][donor],
                        np.array([np.NaN] * len(result['names'][donor])).astype('object'),
                        np.array([np.NaN] * len(result['names'][donor])).astype('object'),
                    )).T))    
        
        results_df  = pd.DataFrame(out[1:], columns=['Gene', 'Target Cell Type', 'Donor', 'Method', 'Score', 'pval_adj', 'LFC'])
        all_results = pd.concat([all_results, results_df])
    
    all_results.to_csv(os.path.join(adata_handler.BASE_RES_DIR, disease_id, dataset_id, f'dge_analysis_result.csv'), index=False)

def main() -> None:
    """
    Main function for performing automatic differential gene expression analysis.

    Parses command-line arguments, loads annotated data, performs differential gene expression analysis,
    and saves the results to a CSV file.

    Parameters:
        None

    Returns:
        None
    """
    parser = argparse.ArgumentParser(description='Automatic Cell-type Annotation Script')
    parser = parse_add_args(parser)
    args   = parser.parse_args()

    disease_id  = args.disease_id
    dataset_id  = args.dataset_id
    dge_methods = args.dge_methods

    # load existing anndata with case and control
    adata = adata_handler.load_data_case_and_control(disease_id=disease_id, dataset_id=dataset_id)

    # perform differential gene expression analysis and store the data in a csv
    perform_diff_gene_exp_analysis(adata=adata, disease_id=disease_id, dataset_id=dataset_id,
                                   dge_methods=dge_methods, specific_celltype='All')

    # Save adata
    adata_handler.save_adata(adata=adata, disease_id=disease_id, dataset_id=dataset_id)


if __name__ == '__main__':
    main()

import os
import argparse
import pandas as pd
import scanpy as sc
import omicverse as ov
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Callable, Any, Tuple
import adata_handler


def parse_add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument('--disease_id', type=str, help='Unique identifier for disease, e.g., diabetesII, mpn')
    parser.add_argument('--dataset_id', type=str, help='Unique identifier of the dataset for that specific disease')
    parser.add_argument('--cellmarker', action='store_true', help='Use cellmarker as the reference dataset for annotation')
    parser.add_argument('--panglaodb', action='store_true', help='Use panglaodb as the reference dataset for annotation')
    parser.add_argument('--foldchange', type=float, default=1.5, help='Fold change threshold')
    parser.add_argument('--pvalue', type=float, default=0.01, help='P-value threshold')
    parser.add_argument('--celltype', type=str, default='normal', help='Cell type for annotation')
    parser.add_argument('--clustertype', type=str, default='leiden', help='Clustering name used in scanpy')
    parser.add_argument('--rank_rep', action='store_true', help='Rank representation')

    return parser

def apply_cell_type_anno_scsa(adata: sc.AnnData, cellmarker: bool, panglaodb: bool,
                              foldchange: float, pvalue: float, celltype: str,
                              clustertype: str, rank_rep: bool) -> None:
    if cellmarker:
        scsa = ov.single.pySCSA(adata=adata, foldchange=foldchange, pvalue=pvalue,
                                celltype=celltype, target='cellmarker', tissue='All')
        
        scsa.cell_anno(clustertype=clustertype, cluster='all', rank_rep=rank_rep)
        scsa.cell_auto_anno(adata, key='scsa_celltype_cellmarker')

    if panglaodb:
        scsa = ov.single.pySCSA(adata=adata, foldchange=foldchange, pvalue=pvalue,
                                celltype=celltype, target='panglaodb', tissue='All')
        
        scsa.cell_anno(clustertype=clustertype, cluster='all', rank_rep=rank_rep)
        scsa.cell_auto_anno(adata, key='scsa_celltype_panglaodb')

def plot_embeddings(
        adata: sc.AnnData,
        disease_id: str,
        dataset_id: str,
        plot_function: Callable[[sc.AnnData, Any], Any],
        fig_size: Tuple[int, int],
        plot_title: str = None,
        xlabel: str = None,
        ylabel: str = None,
        plot_fname: str = None,
        *args, **kwargs
    ) -> None:
    figures_dir = os.path.join(adata_handler.BASE_RES_DIR, disease_id, dataset_id)

    try:
        fig, ax = plt.subplots(figsize=fig_size)
        plot_function(adata, *args, **kwargs, ax=ax, show=False)
        
        if xlabel: ax.set_xlabel(xlabel)
        if ylabel: ax.set_ylabel(ylabel)
        if plot_title: ax.set_title(plot_title)
        if plot_fname is None: plot_fname = f"{disease_id}_{dataset_id}_{str(len(os.listdir(figures_dir))+1)}.png"
        
        plt.savefig(os.path.join(figures_dir, plot_fname), bbox_inches='tight')
        plt.close()
        print(f"Succeed: Plot saved as {plot_fname} in the corresponding figure dir.")
    except Exception as e:
        print(f"Failed: Error during plot generation, {str(e)}")

def calculate_roe(adata: sc.AnnData, sample_key: str, cell_type_key: str) -> pd.DataFrame:
    roe = None

    try:
        roe = ov.utils.roe(adata, sample_key=sample_key, cell_type_key=cell_type_key)
    except Exception as e:
        print(f"Failed: Error during Ro/e calculation, {str(e)}")

    return roe

def plot_roe(roe: pd.DataFrame, disease_id: str, dataset_id: str,):
    figures_dir = os.path.join(adata_handler.BASE_RES_DIR, disease_id, dataset_id)
    fig, ax     = plt.subplots(figsize=(2, 4))

    transformed_roe = roe.copy()
    transformed_roe = transformed_roe.applymap(lambda x: '+++' if x >= 2 else ('++' if x >= 1.5 else ('+' if x >= 1 else '+/-')))
    
    sns.heatmap(roe, annot=transformed_roe, cmap='RdBu_r', fmt='', cbar=True, ax=ax, vmin=0.5, vmax=1.5, cbar_kws={'shrink':0.5})
    
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.xlabel('Donor', fontsize=12)
    plt.ylabel('Cell type', fontsize=12)
    plt.title('Ro/e', fontsize=12)
    plt.savefig(os.path.join(figures_dir, 'roe.png'), bbox_inches='tight')

def run_cell_type_annotation() -> None:
    parser = argparse.ArgumentParser(description='Automatic Cell-type Annotation Script')
    parser = parse_add_args(parser)
    args   = parser.parse_args()

    disease_id  = args.disease_id
    dataset_id  = args.dataset_id
    cellmarker  = args.cellmarker
    panglaodb   = args.panglaodb
    foldchange  = args.foldchange
    pvalue      = args.pvalue
    celltype    = args.celltype
    clustertype = args.clustertype
    rank_rep    = args.rank_rep

    # load existing anndata with case and control
    adata = adata_handler.load_data_case_and_control(disease_id=disease_id, dataset_id=dataset_id)

    # apply SCSA cell type annotation
    apply_cell_type_anno_scsa(adata=adata, cellmarker=cellmarker, panglaodb=panglaodb,
                              foldchange=foldchange, pvalue=pvalue, celltype=celltype,
                              clustertype=clustertype, rank_rep=rank_rep)
    
    # plot embeddings
    plot_embeddings(adata=adata, disease_id=disease_id, dataset_id=dataset_id,
                 plot_function=ov.utils.embedding, fig_size=(4,4), plot_fname='leiden.png',
                 basis='X_mde', color=['leiden'], palette=ov.utils.palette(),
                 plot_title='Leiden Cluster', xlabel='X_mde1', ylabel='X_mde2',
                 frameon='small', legend_loc='on data', legend_fontoutline=0.01)
    plot_embeddings(adata=adata, disease_id=disease_id, dataset_id=dataset_id,
                 plot_function=ov.utils.embedding, fig_size=(4,4), plot_fname='scsa_cellmarker.png',
                 plot_title='SCSA Celltype Annotation (cellmarker)', xlabel='X_mde1', ylabel='X_mde2',
                 basis='X_mde', color=['scsa_celltype_cellmarker'], palette=ov.utils.palette(),
                 frameon='small', legend_fontoutline=0.01)
    plot_embeddings(adata=adata, disease_id=disease_id, dataset_id=dataset_id,
                 plot_function=ov.utils.embedding, fig_size=(4,4), plot_fname='scsa_panglaodb.png',
                 plot_title='SCSA Celltype Annotation (panglaodb)', xlabel='X_mde1', ylabel='X_mde2',
                 basis='X_mde', color=['scsa_celltype_panglaodb'], palette=ov.utils.palette(),
                 frameon='small', legend_fontoutline=0.01)
    plot_embeddings(adata=adata, disease_id=disease_id, dataset_id=dataset_id,
                 plot_function=ov.utils.embedding, fig_size=(4,4), plot_fname='donor_cells.png',
                 plot_title='Donors (Case vs Control)', xlabel='X_mde1', ylabel='X_mde2',
                 basis='X_mde', color=['donor'], palette=ov.utils.red_palette(),
                 frameon='small', legend_fontoutline=2)
    
    # calculate Ro/e and plot if Ro/e calculated successfully
    roe = calculate_roe(adata=adata, sample_key='donor', cell_type_key='scsa_celltype_cellmarker')
    if roe is not None:
        plot_roe(roe=roe, disease_id=disease_id, dataset_id=dataset_id)

    # Save adata
    adata_handler.save_adata(adata=adata, disease_id=disease_id, dataset_id=dataset_id)


if __name__ == '__main__':
    run_cell_type_annotation()

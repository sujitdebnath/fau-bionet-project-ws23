import argparse
import scanpy as sc
import omicverse as ov
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
    print(adata)

    # apply SCSA cell type annotation
    apply_cell_type_anno_scsa(adata=adata, cellmarker=cellmarker, panglaodb=panglaodb,
                              foldchange=foldchange, pvalue=pvalue, celltype=celltype,
                              clustertype=clustertype, rank_rep=rank_rep)

    # Save adata
    adata_handler.save_adata(adata=adata, disease_id=disease_id, dataset_id=dataset_id)
    print(adata)


if __name__ == '__main__':
    run_cell_type_annotation()

# Python imports
import argparse

# Third party imports
import scanpy as sc
import omicverse as ov

# Self imports
import adata_handler


def parse_add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add arguments to the argument parser for preprocessing AnnData objects.

    Parameters:
        parser (argparse.ArgumentParser): Argument parser object.

    Returns:
        argparse.ArgumentParser: Updated argument parser object.
    """
    parser.add_argument('--disease_id', type=str, help='Unique identifier for disease, e.g., diabetesII, mpn')
    parser.add_argument('--dataset_id', type=str, help='Unique identifier of the dataset for that specific disease')
    parser.add_argument('--mito_perc', type=float, default=0.05, help='Percentage of mitochondrial genes threshold for quality control')
    parser.add_argument('--n_umis', type=int, default=500, help='Number of UMIs threshold for quality control')
    parser.add_argument('--detected_genes', type=int, default=250, help='Number of detected genes threshold for quality control')
    parser.add_argument('--preprocess_mode', type=str, default='shiftlog|pearson', help='Preprocessing mode')
    parser.add_argument('--n_hvgs', type=int, default=2000, help='Number of highly variable genes to select')
    parser.add_argument('--n_pcs', type=int, default=50, help='Number of principal components for PCA')
    parser.add_argument('--n_neighbors', type=int, default=15, help='Number of neighbors for constructing neighborhood graph')
    parser.add_argument('--cluster_mode', type=str, default='leiden', choices=['leiden'], help='Clustering method')

    return parser

def apply_quantity_control(adata: sc.AnnData, mito_perc: float, n_umis: int, detected_genes: int) -> sc.AnnData:
    """
    Apply quality control to the AnnData object.

    Parameters:
        adata (sc.AnnData): AnnData object.
        mito_perc (float): Percentage of mitochondrial genes threshold for quality control.
        n_umis (int): Number of UMIs threshold for quality control.
        detected_genes (int): Number of detected genes threshold for quality control.

    Returns:
        sc.AnnData: Updated AnnData object after quality control.
    """
    adata = ov.pp.qc(adata, tresh = {'mito_perc': mito_perc, 'nUMIs': n_umis, 'detected_genes': detected_genes})

    return adata

def preprocess_adata(adata: sc.AnnData, preprocess_mode: str, n_hvgs: int) -> sc.AnnData:
    """
    Preprocess the AnnData object.

    Parameters:
        adata (sc.AnnData): AnnData object.
        preprocess_mode (str): Preprocessing mode.
        n_hvgs (int): Number of highly variable genes to select.

    Returns:
        sc.AnnData: Updated AnnData object after preprocessing.
    """
    adata     = ov.pp.preprocess(adata, mode=preprocess_mode, n_HVGs=n_hvgs)
    adata.raw = adata
    adata     = adata[:, adata.var.highly_variable_features]

    return adata

def apply_pca(adata: sc.AnnData, n_pcs: int) -> sc.AnnData:
    """
    Apply Principal Component Analysis (PCA) to the AnnData object.

    Parameters:
        adata (sc.AnnData): AnnData object.
        n_pcs (int): Number of principal components for PCA.

    Returns:
        sc.AnnData: Updated AnnData object after PCA.
    """
    ov.pp.scale(adata)
    ov.pp.pca(adata, layer='scaled', n_pcs=n_pcs)

    return adata

def construct_neighbourhood_graph(adata: sc.AnnData, n_neighbors: int, n_pcs: int) -> sc.AnnData:
    """
    Construct a neighborhood graph for the AnnData object.

    Parameters:
        adata (sc.AnnData): AnnData object.
        n_neighbors (int): Number of neighbors for constructing neighborhood graph.
        n_pcs (int): Number of principal components for constructing neighborhood graph.

    Returns:
        sc.AnnData: Updated AnnData object with neighborhood graph.
    """
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep='scaled|original|X_pca')

    # apply dimensionality reduction for visualization (X_mde=X_umap+GPU)
    adata.obsm["X_mde"] = ov.utils.mde(adata.obsm["scaled|original|X_pca"])

    return adata

def apply_clustering(adata: sc.AnnData, cluster_mode: str) -> sc.AnnData:
    """
    Apply clustering to the AnnData object.

    Parameters:
        adata (sc.AnnData): AnnData object.
        cluster_mode (str): Clustering method.

    Returns:
        sc.AnnData: Updated AnnData object after clustering.
    """
    if cluster_mode == 'leiden':
        sc.tl.leiden(adata)
    
    return adata

def main() -> None:
    """
    Main function for data preprocessing.

    Parses command-line arguments, loads data, applies preprocessing steps,
    and saves the preprocessed data.

    Parameters:
        None

    Returns:
        None
    """
    parser = argparse.ArgumentParser(description='Data Preprocessor Script')
    parser = parse_add_args(parser)
    args   = parser.parse_args()

    disease_id      = args.disease_id
    dataset_id      = args.dataset_id
    mito_perc       = args.mito_perc
    n_umis          = args.n_umis
    detected_genes  = args.detected_genes
    preprocess_mode = args.preprocess_mode
    n_hvgs          = args.n_hvgs
    n_pcs           = args.n_pcs
    n_neighbors     = args.n_neighbors
    cluster_mode    = args.cluster_mode

    # load existing anndata with case and control
    adata = adata_handler.load_data_case_and_control(disease_id=disease_id, dataset_id=dataset_id)

    # apply quantity control
    adata = apply_quantity_control(adata=adata, mito_perc=mito_perc, n_umis=n_umis, detected_genes=detected_genes)
    
    # calculate high variable genes (HVGs), save the whole genes and filter the non-HVGs
    adata = preprocess_adata(adata=adata, preprocess_mode=preprocess_mode, n_hvgs=n_hvgs)
    
    # scale and apply dimensionality reduction
    adata = apply_pca(adata=adata, n_pcs=n_pcs)
    
    # construct neighbourhood graph and apply dimensionality reduction for visualization (X_mde=X_umap+GPU)
    adata = construct_neighbourhood_graph(adata=adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # apply clusters
    adata = apply_clustering(adata=adata, cluster_mode=cluster_mode)

    # save preprocessed data
    adata_handler.save_adata(adata=adata, disease_id=disease_id, dataset_id=dataset_id)


if __name__ == '__main__':
    main()

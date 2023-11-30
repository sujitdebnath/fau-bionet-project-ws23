import numpy as np
import pandas as pd
import scanpy as sc


class Pipeline:
    def __init__(self, verbosity_lv, source_file_path, results_file_path):
        self.verbosity_lv = verbosity_lv
        self.source_file_path = source_file_path
        self.result_file_path = results_file_path
        self._set_settings()
        self._load_data()



    def _set_settings(self):
        sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=80, facecolor='white')

    def _load_data(self):
        adata = sc.read_10x_mtx(
            self.source_file_path,  # the directory with the `.mtx` file
            var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
            cache=True)

        self.adata = adata

    def plot_highest_expr_genes(self):
        sc.pl.highest_expr_genes(self.adata, n_top=20)

    def plot_highly_variable_genes(self):
        sc.pp.highly_variable_genes(self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(self.adata)

    def preprocessing(self):
        sc.pp.filter_cells(self.adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)

        self.adata.var['mt'] = self.adata.var_names.str.startswith(
            'MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    def plot_scatter_adata(self):
        sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_mt')
        sc.pl.scatter(self.adata, x='total_counts', y='n_genes_by_counts')

    def filter_data(self):
        self.adata = self.adata[self.adata.obs.n_genes_by_counts < 2500, :]
        self.adata = self.adata[self.adata.obs.pct_counts_mt < 5, :]

    def normalize_data(self):
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)

    def select_highly_variable_data(self):
        self.adata.raw = self.adata
        self.adata = self.adata[:, self.adata.var.highly_variable]

    def scale_data(self):
        sc.pp.scale(self.adata, max_value=10)

    def dimension_reduction(self):
        sc.tl.pca(self.adata, svd_solver='arpack')

    def plot_pca(self):
        sc.pl.pca(self.adata, color='CST3')

    def plot_pca_variance_ration(self):
        sc.pl.pca_variance_ratio(self.adata, log=True)

    def write_result_file(self):
        self.adata.write(self.result_file_path)

    def find_neighbours(self, n_neighbors, n_pcs):
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    def cluster(self):
        sc.tl.leiden(self.adata)

    def compute_UMAP(self):
        sc.tl.umap(self.adata)

    def plot_UMAP(self, use_raw, colors):
        sc.pl.umap(self.adata, color=colors, use_raw=use_raw)

    def rank_gene_groups(self, n_genes, group_by='leiden', method='t-test'):
        sc.tl.rank_genes_groups(self.adata, groupby=group_by, method=method)
        sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=False)

    def plot_violin_data(self):
        sc.pl.rank_genes_groups_violin(self.adata, groups='0', n_genes=8)

    def cell_type_annotation(self):
        pass

    def differential_gene_expression(self):
        pass

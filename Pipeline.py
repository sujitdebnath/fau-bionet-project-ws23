import os.path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


class Pipeline:
    def __init__(self, verbosity_lv, source_file_path, results_file_path, name):
        self.verbosity_lv = verbosity_lv
        self.source_file_path = source_file_path
        self.result_file_path = results_file_path
        self.name = name
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
        plt.clf()
        sc.pl.highest_expr_genes(self.adata, n_top=20, show=False)
        # plt.title(f'Highest Expressed Genes in {self.name}')
        file_name = f'highest_expr_genes-{self.name}.png'
        self.highest_expr_genes_url = f'www/{file_name}'
        plt.savefig(self.highest_expr_genes_url)

    def plot_highly_variable_genes(self):
        plt.clf()
        sc.pp.highly_variable_genes(self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(self.adata, show=False)
        # plt.title(f'Highly Variable Genes in {self.name}')
        file_name = f'highly_variable_genes-{self.name}.png'
        self.highly_variable_genes_url = f'www/{file_name}'
        plt.savefig(self.highly_variable_genes_url)

    def preprocessing(self):
        sc.pp.filter_cells(self.adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)

        self.adata.var['mt'] = self.adata.var_names.str.startswith(
            'MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    def plot_scatter_adata(self):
        plt.clf()
        sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_mt', show=False)
        # plt.title(f'PCT Counts (MT) in {self.name}')
        file_name = f'pct_counts_mt-{self.name}.png'
        self.pct_counts_mt_url = f'www/{file_name}'
        plt.savefig(self.pct_counts_mt_url)

        plt.clf()
        sc.pl.scatter(self.adata, x='total_counts', y='n_genes_by_counts', show=False)
        # plt.title(f'Number of Genes by Count in {self.name}')
        file_name = f'n_genes_by_counts-{self.name}.png'
        self.n_genes_by_counts_url = f'www/{file_name}'
        plt.savefig(self.n_genes_by_counts_url)

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
        plt.clf()
        sc.pl.pca(self.adata, color='CST3', show=False)
        # plt.title(f'PCA in {self.name}')
        file_name = f'pca-{self.name}.png'
        self.pca_url = f'www/{file_name}'
        plt.savefig(self.pca_url)

    def plot_pca_variance_ration(self):
        plt.clf()
        sc.pl.pca_variance_ratio(self.adata, log=True, show=False)
        # plt.title(f'PCA Variance in {self.name}')
        file_name = f'pca_variance-{self.name}.png'
        self.pca_variance_url = f'www/{file_name}'
        plt.savefig(self.pca_variance_url)

    def write_result_file(self):
        if not os.path.exists('write'):
            os.mkdir('write')

        self.adata.write(self.result_file_path)

    def find_neighbours(self, n_neighbors, n_pcs):
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    def cluster(self):
        sc.tl.leiden(self.adata)

    def compute_UMAP(self):
        sc.tl.umap(self.adata)

    def plot_UMAP(self, use_raw, colors):
        plt.clf()
        sc.pl.umap(self.adata, color=colors, use_raw=use_raw, show=False)
        # plt.title(f'UMAP in {self.name}')
        file_name = f'umap-{colors}-{self.name}.png'
        self.umap_url = f'www/{file_name}'
        plt.savefig(self.umap_url)

    def rank_gene_groups(self, n_genes, group_by='leiden', method='t-test'):
        plt.clf()
        sc.tl.rank_genes_groups(self.adata, groupby=group_by, method=method)
        sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=False, show=False)
        # plt.title(f'Rank Genes Group in {self.name}')
        file_name = f'rank_genes_group_by-{group_by}-{method}-{self.name}.png'
        self.rank_genes_groups_url = f'www/{file_name}'
        plt.savefig(self.rank_genes_groups_url)

    def plot_violin_data(self):
        plt.clf()
        sc.pl.rank_genes_groups_violin(self.adata, groups='0', n_genes=8, show=False)
        # plt.title(f'Rank Genes Group Violin in {self.name}')
        file_name = f'rank_genes_groups_violin-{self.name}.png'
        self.rank_genes_groups_violin_url = f'www/{file_name}'
        plt.savefig(self.rank_genes_groups_violin_url)

    def cell_type_annotation(self):
        pass

    def differential_gene_expression(self):
        pass

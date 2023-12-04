import os.path
import scanpy as sc
import matplotlib.pyplot as plt
from PIL import Image


def merge_umap_plots_vertically(image_paths, output_path):
    images = [Image.open(path) for path in image_paths]

    # Get the maximum width among all images
    max_width = max(img.width for img in images)

    # Calculate the total height by summing individual image heights
    total_height = sum(img.height for img in images)

    # Create a new image with the maximum width and total height
    result = Image.new("RGB", (max_width, total_height), (255, 255, 255))

    # Paste each image onto the result image
    current_height = 0
    for img in images:
        result.paste(img, (0, current_height))
        current_height += img.height

    # Save the result image
    result.save(output_path)


class Pipeline:
    def __init__(self, verbosity_lv, source_file_path, name):
        self.verbosity_lv = verbosity_lv
        self.source_file_path = source_file_path
        self.name = name
        self.result_file_path = f'write/{self.name}/data.h5ad'
        self._set_settings()
        self._load_data()

    def _set_settings(self):
        sc.settings.verbosity = self.verbosity_lv  # verbosity: errors (0), warnings (1), info (2), hints (3)
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=80, facecolor='white')

        if not os.path.exists('figures'):
            os.mkdir('figures')

        if not os.path.exists(f'figures/{self.name}'):
            os.mkdir(f'figures/{self.name}')

        if not os.path.exists('write'):
            os.mkdir('write')

        if not os.path.exists(f'write/{self.name}'):
            os.mkdir(f'write/{self.name}')

    def _load_data(self):
        adata = sc.read_10x_mtx(
            self.source_file_path,  # the directory with the `.mtx` file
            var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
            cache=True
        )

        self.adata = adata

    def _plot_highest_expr_genes(self):
        plt.clf()
        sc.pl.highest_expr_genes(self.adata, n_top=20, show=False)
        # plt.title(f'Highest Expressed Genes in {self.name}')
        file_name = f'highest_expr_genes.png'
        self.highest_expr_genes_url = f'figures/{self.name}/{file_name}'
        plt.savefig(self.highest_expr_genes_url)

    def _plot_highly_variable_genes(self):
        plt.clf()
        sc.pp.highly_variable_genes(self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(self.adata, show=False)
        file_name = f'highly_variable_genes.png'
        self.highly_variable_genes_url = f'figures/{self.name}/{file_name}'
        plt.savefig(self.highly_variable_genes_url)

    def _preprocessing(self):
        sc.pp.filter_cells(self.adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)

        self.adata.var['mt'] = self.adata.var_names.str.startswith(
            'MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    def _plot_scatter_adata(self):
        plt.clf()
        sc.pl.scatter(self.adata, x='total_counts', y='pct_counts_mt', show=False)
        file_name = f'pct_counts_mt.png'
        self.pct_counts_mt_url = f'figures/{self.name}/{file_name}'
        plt.savefig(self.pct_counts_mt_url)

        plt.clf()
        sc.pl.scatter(self.adata, x='total_counts', y='n_genes_by_counts', show=False)
        file_name = f'n_genes_by_counts.png'
        self.n_genes_by_counts_url = f'figures/{self.name}/{file_name}'
        plt.savefig(self.n_genes_by_counts_url)

    def _filter_data(self):
        self.adata = self.adata[self.adata.obs.n_genes_by_counts < 2500, :]
        self.adata = self.adata[self.adata.obs.pct_counts_mt < 5, :]

    def _normalize_data(self):
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)

    def _select_highly_variable_data(self):
        self.adata.raw = self.adata
        self.adata = self.adata[:, self.adata.var.highly_variable]

    def _scale_data(self):
        sc.pp.scale(self.adata, max_value=10)

    def _dimension_reduction(self):
        sc.tl.pca(self.adata, svd_solver='arpack')

    def _plot_pca(self):
        plt.clf()
        sc.pl.pca(self.adata, color='CST3', show=False)
        file_name = f'pca.png'
        self.pca_url = f'figures/{self.name}/{file_name}'
        plt.savefig(self.pca_url)

    def _plot_pca_variance_ration(self):
        plt.clf()
        sc.pl.pca_variance_ratio(self.adata, log=True, show=False)
        file_name = f'pca_variance.png'
        self.pca_variance_url = f'figures/{self.name}/{file_name}'
        plt.savefig(self.pca_variance_url)

    def _write_result_file(self):
        self.adata.write(self.result_file_path)

    def _find_neighbours(self, n_neighbors, n_pcs):
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    def _cluster(self):
        sc.tl.leiden(self.adata)

    def _compute_UMAP(self):
        sc.tl.umap(self.adata)

    def plot_UMAP(self, use_raw, colors, is_after_clustering=False):
        umap_plots = []
        umap_plots_paths = []

        for c in colors:
            plt.clf()
            sc.pl.umap(self.adata, color=[c], use_raw=use_raw, show=False)
            umapt = plt.gcf()
            umap_plots.append(umapt)
            file_name = f'umap-{c}.png'
            plt.savefig(f'figures/{self.name}/{file_name}')
            umap_plots_paths.append(f'figures/{self.name}/{file_name}')

        output_umap_plots_path = f'figures/{self.name}/umap.jpg'

        self.umap_url = output_umap_plots_path

        merge_umap_plots_vertically(image_paths=umap_plots_paths, output_path=output_umap_plots_path)

    def _rank_gene_groups(self, n_genes, group_by='leiden', method='t-test'):
        plt.clf()
        sc.tl.rank_genes_groups(self.adata, groupby=group_by, method=method)
        sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=False, show=False)
        file_name = f'rank_genes_group_by-{group_by}-{method}.png'
        self.rank_genes_groups_url = f'figures/{self.name}/{file_name}'
        plt.savefig(self.rank_genes_groups_url)

    def _plot_violin_data(self):
        plt.clf()
        sc.pl.rank_genes_groups_violin(self.adata, groups='0', n_genes=8, show=False)
        # plt.title(f'Rank Genes Group Violin in {self.name}')
        file_name = f'rank_genes_groups_violin.png'
        self.rank_genes_groups_violin_url = f'figures/{self.name}/{file_name}'
        plt.savefig(self.rank_genes_groups_violin_url)

    def run(self):
        self._plot_highest_expr_genes()

        self._preprocessing()

        self._plot_scatter_adata()

        self._filter_data()

        self._normalize_data()

        self._plot_highly_variable_genes()

        self._select_highly_variable_data()

        self._scale_data()

        self._dimension_reduction()

        self._plot_pca()

        self._plot_pca_variance_ration()

        self._write_result_file()

        self._find_neighbours(n_neighbors=10, n_pcs=40)

        self._compute_UMAP()

        self._cluster()

        self._write_result_file()

        self._rank_gene_groups(n_genes=25, group_by='leiden', method='t-test')

        self._rank_gene_groups(n_genes=25, group_by='leiden', method='wilcoxon')

        self._rank_gene_groups(n_genes=25, group_by='leiden', method='logreg')

        self._plot_violin_data()

    def cell_type_annotation(self):
        pass

    def differential_gene_expression(self):
        pass

from Pipeline import Pipeline

if __name__ == '__main__':
    pipeline = Pipeline(verbosity_lv=3,
                        source_file_path='data/filtered_gene_bc_matrices/hg19/',
                        results_file_path='write/SC3.h5ad')

    pipeline.plot_highest_expr_genes()

    pipeline.preprocessing()

    pipeline.plot_scatter_adata()

    pipeline.filter_data()

    pipeline.normalize_data()

    pipeline.plot_highly_variable_genes()

    pipeline.select_highly_variable_data()

    pipeline.scale_data()

    pipeline.dimension_reduction()

    pipeline.plot_pca()

    pipeline.plot_pca_variance_ration()

    pipeline.write_result_file()

    pipeline.find_neighbours(n_neighbors=10, n_pcs=40)

    pipeline.compute_UMAP()

    pipeline.plot_UMAP(use_raw=False, colors=['CST3', 'NKG7', 'PPBP'])

    pipeline.cluster()

    pipeline.plot_UMAP(use_raw=True, colors=['leiden', 'CST3', 'NKG7'])

    pipeline.write_result_file()

    pipeline.rank_gene_groups(n_genes=25, group_by='leiden', method='t-test')

    pipeline.rank_gene_groups(n_genes=25, group_by='leiden', method='wilcoxon')

    pipeline.rank_gene_groups(n_genes=25, group_by='leiden', method='logreg')

    pipeline.plot_violin_data()
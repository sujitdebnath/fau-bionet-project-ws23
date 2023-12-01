import pandas as pd
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo
from Pipeline import Pipeline
import os


app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.panel_sidebar(

            ui.row(
                ui.input_text("txt", "Enter the text to display below:", "delete me"),
                ui.output_text('text'),
            ),

            ui.row(
                ui.input_select(
                    "dataset",
                    "Choose a dataset:",
                    ['WB Lysis Granulocytes 5p Introns 8kCells', 'PBMC3k'],
                ),
            ),
            ui.input_action_button("run", "Begin Analysis!", class_="btn-success"),
        ),
        ui.output_plot('run_pipeline')
    ),
)


def server(input: Inputs, output: Outputs, session: Session):
    @output
    @render.text
    def text():
        return input.txt()

    @output
    @render.plot
    @reactive.event(input.run, ignore_none=True)
    def run_pipeline():
        dataset_name = str(input.dataset())
        if dataset_name == 'PBMC3k':
            pipeline = Pipeline(verbosity_lv=3,
                                source_file_path='data/filtered_gene_bc_matrices/hg19',
                                results_file_path='write/filtered_gene_bc_matrices.h5ad',
                                name='PBMC3k')

        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = Pipeline(verbosity_lv=3,
                                source_file_path='data/WB_Lysis_Granulocytes_5p_Introns_8kCells_filtered_feature_bc_matrix/filtered_feature_bc_matrix',
                                results_file_path='write/WB_Lysis_Granulocytes_5p_Introns_8kCells.h5ad',
                                name='WB-Lysis')

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
        return


app = App(app_ui, server)
app.run()

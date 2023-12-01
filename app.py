import pandas as pd
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo, ImgData
from Pipeline import Pipeline
import os
from PIL import Image
from datetime import datetime

app_ui = ui.page_fluid(
    # ui.head_content(ui.include_js("js/main.js", method="inline")),

    ui.h1("BioNet Project"),
    ui.p("Created by Farzam"),
    ui.layout_sidebar(
        ui.panel_sidebar(

            ui.row(
                ui.input_select(
                    "dataset",
                    "Choose a dataset:",
                    ['WB Lysis Granulocytes 5p Introns 8kCells', 'PBMC3k'],
                ),
            ),
            ui.input_action_button("run", "Begin Analysis", class_="btn-success"),
            ui.input_action_button("plot_figures", "Plot Figures", class_="btn-primary"),
        ),
        ui.navset_tab(
            ui.nav("Highest Expressed Genes", ui.output_image('plot_highest_expressed_genes')),
            ui.nav("Highly Variable Genes", ui.output_image('plot_highly_variable_genes')),
            ui.nav("PCT Counts (MT)", ui.output_image('plot_pct_count')),
            ui.nav("Number of Genes by Count", ui.output_image('plot_num_of_genes_by_count')),
            ui.nav("PCA", ui.output_image('plot_pca')),
            ui.nav("PCA Variance", ui.output_image('plot_pca_variance')),
            ui.nav("UMAP", ui.output_image('plot_umap')),
            ui.nav("Rank Genes Group", ui.output_image('plot_rank_genes_group')),
            ui.nav("Rank Genes Group Violin", ui.output_image('plot_rank_genes_group_violin')),
        ),
        ui.output_ui('run_pipeline'),
    ),
)

pipelines_info = {
    'PBMC3k': None,
    'WB_Lysis_Granulocytes_5p_Introns_8kCells': None
}


def server(input: Inputs, output: Outputs, session: Session):
    @output
    @render.text
    def text():
        return input.txt()

    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_highest_expressed_genes():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            if pipeline.highest_expr_genes_url is None:
                ui.notification_show('Figure has not been yet plotted.', type='error')
                return
            else:
                return ImgData(src=pipeline.highest_expr_genes_url, height='100%')

    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_highly_variable_genes():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            if pipeline.highest_expr_genes_url is None:
                ui.notification_show('Figure has not been yet plotted.', type='error')
                return
            else:
                return ImgData(src=pipeline.highly_variable_genes_url, height='100%')

    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_pct_count():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            if pipeline.highest_expr_genes_url is None:
                ui.notification_show('Figure has not been yet plotted.', type='error')
                return
            else:
                return ImgData(src=pipeline.pct_counts_mt_url, height='100%')


    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_num_of_genes_by_count():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            if pipeline.highest_expr_genes_url is None:
                ui.notification_show('Figure has not been yet plotted.', type='error')
                return
            else:
                return ImgData(src=pipeline.n_genes_by_counts_url, height='100%')


    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_pca():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            if pipeline.highest_expr_genes_url is None:
                ui.notification_show('Figure has not been yet plotted.', type='error')
                return
            else:
                return ImgData(src=pipeline.pca_url, height='100%')


    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_pca_variance():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            if pipeline.highest_expr_genes_url is None:
                ui.notification_show('Figure has not been yet plotted.', type='error')
                return
            else:
                return ImgData(src=pipeline.pca_variance_url, height='100%')

    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_umap():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            if pipeline.highest_expr_genes_url is None:
                ui.notification_show('Figure has not been yet plotted.', type='error')
                return
            else:
                return ImgData(src=pipeline.umap_url, height='100%')


    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_rank_genes_group():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            if pipeline.highest_expr_genes_url is None:
                ui.notification_show('Figure has not been yet plotted.', type='error')
                return
            else:
                return ImgData(src=pipeline.rank_genes_groups_url, height='100%')

    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_rank_genes_group_violin():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            if pipeline.highest_expr_genes_url is None:
                ui.notification_show('Figure has not been yet plotted.', type='error')
                return
            else:
                return ImgData(src=pipeline.rank_genes_groups_violin_url, height='100%')

    @output
    @render.ui
    @reactive.event(input.run, ignore_none=True)
    def run_pipeline():
        dataset_name = str(input.dataset())

        if dataset_name == 'PBMC3k':
            pipeline = pipelines_info['PBMC3k']
        elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
            pipeline = pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells']

        if pipeline is not None:
            ui.notification_show('Analysis is already done previously. Try to plot the figures.', type='warning')

        else:
            ui.notification_show(f'Analysis has started!', type='message')

            if dataset_name == 'PBMC3k':
                pipeline = Pipeline(verbosity_lv=1,
                                    source_file_path='data/filtered_gene_bc_matrices/hg19',
                                    results_file_path='write/filtered_gene_bc_matrices.h5ad',
                                    name='PBMC3k')

                pipelines_info['PBMC3k'] = pipeline

            elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
                pipeline = Pipeline(verbosity_lv=1,
                                    source_file_path='data/WB_Lysis_Granulocytes_5p_Introns_8kCells_filtered_feature_bc_matrix/filtered_feature_bc_matrix',
                                    results_file_path='write/WB_Lysis_Granulocytes_5p_Introns_8kCells.h5ad',
                                    name='WB-Lysis')

                pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells'] = pipeline

            t0 = datetime.now()

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

            t1 = datetime.now()

            elapsed_time = (t1 - t0).total_seconds()

            ui.notification_show(f'Analysis is done in {round(elapsed_time, 3)} seconds!', type='message')


app = App(app_ui, server)
# app.run()

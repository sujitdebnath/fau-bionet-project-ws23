import matplotlib.pyplot as plt
from shiny import App, Inputs, Outputs, Session, reactive, render, ui, req
from shiny.types import FileInfo, ImgData
from Pipeline import Pipeline
import os
from PIL import Image
from datetime import datetime

app_ui = ui.page_fluid(
    ui.head_content(ui.include_js("js/main.js", method="inline")),

    ui.output_ui('show_welcome_modal'),
    ui.h1("BioNet Project"),
    ui.p("Created by Farzam"),

    ui.layout_sidebar(
        ui.panel_sidebar(
            ui.h5('Configurations for Analysis'),

            ui.row(
                ui.input_select(
                    "dataset",
                    "Choose a dataset:",
                    ['WB Lysis Granulocytes 5p Introns 8kCells', 'PBMC3k'],
                )
            ),

            ui.input_slider("min_genes", "Minimum Genes for Filtering", min=1, max=100, value=3),

            ui.input_slider("min_cells", "Minimum Cells for Filtering", min=1, max=1000, value=200),

            ui.br(),

            ui.input_action_button("run", "Begin Analysis", class_="btn-success", onclick='freeze_buttons();'),

            ui.HTML('<hr width="100%;" color="blue" size="8">'),

            ui.h5('Configurations for Plotting'),

            ui.input_slider("most_expr_genes", "Most expressed genes to plot", min=5, max=100, value=20),

            ui.input_checkbox_group(
                "umap_colors",
                "Choose color(s) for UMAP plots:",
                {
                    "CST3": ui.span("CST3"),
                    "NKG7": ui.span("NKG7"),
                    "PPBP": ui.span("PPBP"),
                    "leiden": ui.span("Leiden")
                }, selected=['CST3', 'NKG7', 'PPBP', 'leiden']
            ),

            # ui.input_slider("num_neighbours", "Number of neighbours to compute neighborhood graph",
            #                 min=5, max=50, value=10),
            # ui.input_slider("num_pcs", "Number of PC(s) to compute neighborhood graph",
            #                 min=5, max=50, value=40),

            ui.input_action_button("plot_figures", "Plot Figures", class_="btn-primary"),

            width=2,
        ),
        {'style': 'height: 1000px;'},

        ui.navset_tab(
            ui.nav(
                "Highest Expressed Genes",
                ui.row(
                    ui.column(6, ui.output_image('plot_highest_expressed_genes')),
                    ui.column(6, ui.output_ui('desc_highest_expressed_genes'))),
                ),

            ui.nav("Highly Variable Genes",
                   ui.row(
                       ui.column(6, ui.output_image('plot_highly_variable_genes')),
                       ui.column(6, ui.output_ui('desc_highly_variable_genes'))),
                   ),

            ui.nav("PCT Counts (MT)",
                   ui.row(
                       ui.column(6, ui.output_image('plot_pct_count')),
                       ui.column(6, ui.output_ui('desc_pct_count'))),
                   ),

            ui.nav("Number of Genes by Count",
                   ui.row(
                       ui.column(6, ui.output_image('plot_num_of_genes_by_count')),
                       ui.column(6, ui.output_ui('desc_num_of_genes_by_count'))),
                   ),

            ui.nav("PCA",
                   ui.row(
                       ui.column(6, ui.output_image('plot_pca')),
                       ui.column(6, ui.output_ui('desc_pca'))),
                   ),

            ui.nav("PCA Variance",
                   ui.row(
                       ui.column(6, ui.output_image('plot_pca_variance')),
                       ui.column(6, ui.output_ui('desc_pca_variance'))),
                   ),

            ui.nav("UMAP",
                   ui.row(
                       ui.column(6, ui.output_image('plot_umap')),
                       ui.column(6, ui.output_ui('desc_umap'))),
                   ),

            ui.nav("Rank Genes Group",
                   ui.row(
                       ui.column(6, ui.output_image('plot_rank_genes_group')),
                       ui.column(6, ui.output_ui('desc_rank_genes_group'))),
                   ),

            ui.nav("Rank Genes Group Violin",
                   ui.row(
                       ui.column(6, ui.output_image('plot_rank_genes_group_violin')),
                       ui.column(6, ui.output_ui('desc_rank_genes_group_violin'))),
                   ),
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
    @render.ui
    def show_welcome_modal():
        m = ui.modal(
            "In this project, I tried to perform sc-RNA data analysis of the selected datasets. This is a under-development project and thus, it's done fully done. Thank you.",
            title="Welcome to the Project Biomedical Network Science!",
            easy_close=True,
            footer='Farzam',
        )
        return ui.modal_show(m)

    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_highest_expressed_genes():
        return ui.markdown("### Hi, this is a sample description for highest_expressed_genes!")

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
                return ImgData(src=pipeline.highest_expr_genes_url, height='auto', width='100%')

    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_highly_variable_genes():
        return ui.markdown("### Hi, this is a sample description for highly_variable_genes!")

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
                return ImgData(src=pipeline.highly_variable_genes_url, height='auto', width='100%')


    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_pct_count():
        return ui.markdown("### Hi, this is a sample description for pct_count!")


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
                return ImgData(src=pipeline.pct_counts_mt_url, height='auto', width='100%')


    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_num_of_genes_by_count():
        return ui.markdown("### Hi, this is a sample description for num_of_genes_by_count!")


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
                return ImgData(src=pipeline.n_genes_by_counts_url, height='auto', width='100%')


    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_pca():
        return ui.markdown("### Hi, this is a sample description for pca!")


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
                return ImgData(src=pipeline.pca_url, height='auto', width='100%')


    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_pca_variance():
        return ui.markdown("### Hi, this is a sample description for pca_variance!")


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
                return ImgData(src=pipeline.pca_variance_url, height='auto', width='100%')


    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_umap():
        return ui.markdown("### Hi, this is a sample description for umap!")


    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_umap():
        dataset_name = str(input.dataset())
        umap_colors = input.umap_colors()

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
                pipeline.plot_UMAP(use_raw=False, colors=umap_colors, is_after_clustering=False)
                return ImgData(src=pipeline.umap_url, height='auto', width='60%')


    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_rank_genes_group():
        return ui.markdown("### Hi, this is a sample description rank genes group!")


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
                return ImgData(src=pipeline.rank_genes_groups_url, height='auto', width='100%')


    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_rank_genes_group_violin():
        return ui.markdown("### Hi, this is a sample description for rank genes group violin!")


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
                return ImgData(src=pipeline.rank_genes_groups_violin_url, height='auto', width='100%')

    @output
    @render.ui
    @reactive.event(input.run, ignore_none=True)
    def run_pipeline():
        dataset_name = str(input.dataset())
        n_most_expr_genes = input.most_expr_genes()
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()

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
                                    name='PBMC3k', n_most_expr_genes=n_most_expr_genes,
                                    min_num_genes_for_filtering=min_num_genes_for_filtering,
                                    min_num_cells_for_filtering=min_num_cells_for_filtering)

                pipelines_info['PBMC3k'] = pipeline

            elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
                pipeline = Pipeline(verbosity_lv=1,
                                    source_file_path='data/WB_Lysis_Granulocytes_5p_Introns_8kCells_filtered_feature_bc_matrix/filtered_feature_bc_matrix',
                                    name='WB-Lysis', n_most_expr_genes=n_most_expr_genes,
                                    min_num_genes_for_filtering=min_num_genes_for_filtering,
                                    min_num_cells_for_filtering=min_num_cells_for_filtering)

                pipelines_info['WB_Lysis_Granulocytes_5p_Introns_8kCells'] = pipeline

            t0 = datetime.now()

            pipeline.run()

            t1 = datetime.now()

            elapsed_time = (t1 - t0).total_seconds()

            ui.notification_show(f'Analysis is done in {round(elapsed_time, 3)} seconds!',
                                 type='message', id='success_message',
                                 action=ui.HTML('<script>release_buttons();</script>'))


app = App(app_ui, server)
# app.run()

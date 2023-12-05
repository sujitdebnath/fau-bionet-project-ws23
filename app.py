from shiny import App, Inputs, Outputs, Session, reactive, render, ui, req
from shiny.types import FileInfo, ImgData
from Pipeline import Pipeline
from datetime import datetime

app_ui = ui.page_fluid(
    ui.head_content(ui.include_js("assets/js/main.js", method="inline"),
                    ui.include_js("assets/js/bootstrap.bundle.min.js", method="inline"),
                    ui.include_css("assets/css/main.css", method="inline"),
                    ui.include_css("assets/css/bootstrap.min.css", method="inline"), ),

    ui.output_ui('show_welcome_modal'),

    ui.layout_sidebar(
        ui.panel_sidebar(
            ui.h5('Configurations for Analysis'),

            ui.row(
                ui.input_select(
                    "dataset",
                    "Choose a dataset:",
                    ['PBMC3k', 'WB Lysis Granulocytes 5p Introns 8kCells'],
                )
            ),

            ui.input_slider("min_genes", "Minimum Genes for Filtering", min=1, max=100, value=3),

            ui.input_slider("min_cells", "Minimum Cells for Filtering", min=1, max=1000, value=200),

            ui.input_slider("num_neighbours", "Number of neighbours to compute neighborhood graph",
                            min=5, max=50, value=10),
            ui.input_slider("num_pcs", "Number of PC(s) to compute neighborhood graph",
                            min=5, max=50, value=40),

            ui.br(),

            ui.input_action_button("run", "Begin Analysis", class_="btn-success"),

            ui.HTML('<hr width="100%;" color="black" size="8">'),

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

            ui.input_action_button("plot_figures", "Plot Figures", class_="btn-primary"),

            width=2,
        ),
        {'style': 'height: 1000px;'},

        ui.navset_tab(
            ui.nav(
                "Introduction",
                ui.output_ui('desc_introduction'),
            ),

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

pipelines = []


def check_if_pipeline_exists(pipeline_name, min_genes, min_cells, n_neighbours, n_pcs):
    for pipeline in pipelines:
        if pipeline.name == pipeline_name and pipeline.min_num_genes_for_filtering == min_genes and pipeline.min_num_cells_for_filtering == min_cells and pipeline.num_neighbours == n_neighbours and pipeline.num_pcs == n_pcs:
            return pipeline
    return None


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
    def desc_introduction():
        return ui.HTML("""
        <br>
        <div class="container">
            <div class='jumbotron'>
                <h2> Welcome to BioNet </h1>
                <p>This is a introduction to how to use this platform.
                <br>
                Basically, this platform's functions can be divided into 
                three groups:
                <br>
                <ol>
                    <li> Setting the hyperparameters
                    <li> Creating pipeline for Analysis
                    <li> Visualizing the plots
                </ol>
                
                At the moment, this platform supports <b>only two</b> datasets. However, any dataset with a 10x format that can
                be read with <code>sc.read_10x_mtx()</code> can be added.</p>
                
                <h3>First Step</h3>
                <p>To start with, please <b>choose a dataset</b> from the top-left dropdown.
                After selecting the dataset, please feel free to customize the hyperparameters as you wish.
                When you are done with the parameters, click on the <span style='color: green'><b>Begin Analysis</b></span> button. In this step, two occurrences are possible:</p>
                 
                 <ul>
                    <li> A pipeline will be created and will begin analyzing the data with the given parameters. (takes ~ 20 - 40 seconds)</li>
                    <li> Already a pipeline is created, and the results instantly will be shown to you.</li>
                 </ul>
                 
                 <p>In both cases, a pipeline has already been created and has analyzed the data. Thus, the results are now ready to be plotted. After the analysis are done, all plots will be shown automatically.</p>
                 
                 <h3>Second Step</h3>
                 <p>On the top of this section you can find the navigation bar, in which different plots are grouped and shown. In every window, you will find the plot related to that topic and a brief description about the plot and its analysis.
                 For example, you can learn about the <b>most expressed genes in the dataset</b> by clicking on the next window called <span style='color: red;'><b>Highest Expressed Genes</b></span>.</p>
                 
                 <h3>Third Step</h3>
                 <p>Some paramerts are only related to plotting data, for example <b>Most Expressed Genes</b> to plot, or the color set for the <b>UMap</b>. You can twist these parameters in the left column in the <b>Configuration for Plotting</b> section.
                 After any changes, click on the <span style='color: blue'><b>Plot Figures</b></span> to re-plot all the analysis in the desired format.
                 Please be aware that changing these settings doesn't have to do anything with the analysis and you don't need to analyze the data again.</p>
                 
                 <br><br>
                 <h5>Note:</h5>
                 <p>This is a provisional version of this project and is still under development.<br>
                 Latest update: 5/Dec/2023
                 <br>
                 Farzam</p>
            </div>
        </div>
        """)

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
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()
        n_most_expr_genes = input.most_expr_genes()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset with the given configuration.',
                                 type='error')
        else:
            pipeline.plot_highest_expr_genes(n_most_expr_genes=n_most_expr_genes)
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
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

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
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

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
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

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
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

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
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

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
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')

        else:
            pipeline.plot_umap(use_raw=False, colors=umap_colors, is_after_clustering=False)
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
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')
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
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is None:
            ui.notification_show('Please first run analysis on the selected dataset.', type='error')
        else:
            return ImgData(src=pipeline.rank_genes_groups_violin_url, height='auto', width='100%')

    @output
    @render.ui
    @reactive.event(input.run, ignore_none=True)
    def run_pipeline():
        dataset_name = str(input.dataset())
        min_num_genes_for_filtering = input.min_genes()
        min_num_cells_for_filtering = input.min_cells()
        num_neighbours = input.num_neighbours()
        num_pcs = input.num_pcs()

        pipeline = check_if_pipeline_exists(pipeline_name=dataset_name, min_genes=min_num_genes_for_filtering,
                                            min_cells=min_num_cells_for_filtering, n_neighbours=num_neighbours,
                                            n_pcs=num_pcs)

        if pipeline is not None:
            ui.notification_show('Analysis is already done previously. Try to plot the figures.', type='warning')

        else:
            ui.notification_show(f'Analysis has started!', type='message',
                                 action=ui.HTML('<script>freeze_buttons();</script>'))

            if dataset_name == 'PBMC3k':
                pipeline = Pipeline(verbosity_lv=1,
                                    source_file_path='data/filtered_gene_bc_matrices/hg19',
                                    name='PBMC3k',
                                    min_num_genes_for_filtering=min_num_genes_for_filtering,
                                    min_num_cells_for_filtering=min_num_cells_for_filtering,
                                    num_neighbours=num_neighbours, num_pcs=num_pcs)

            elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
                pipeline = Pipeline(verbosity_lv=1,
                                    source_file_path='data/WB_Lysis_Granulocytes_5p_Introns_8kCells_filtered_feature_bc_matrix/filtered_feature_bc_matrix',
                                    name='WB-Lysis',
                                    min_num_genes_for_filtering=min_num_genes_for_filtering,
                                    min_num_cells_for_filtering=min_num_cells_for_filtering,
                                    num_neighbours=num_neighbours, num_pcs=num_pcs)

            pipelines.append(pipeline)

            t0 = datetime.now()

            pipeline.run()

            t1 = datetime.now()

            elapsed_time = (t1 - t0).total_seconds()

            ui.notification_show(f'Analysis is done in {round(elapsed_time, 3)} seconds!',
                                 type='message', id='success_message',
                                 action=ui.HTML('<script>release_buttons(); auto_plot(); goto_first_tab();</script>'))


app = App(app_ui, server)
# app.run()

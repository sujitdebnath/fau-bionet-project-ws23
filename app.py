from shiny import App, Inputs, Outputs, Session, reactive, render, ui, req
from shiny.types import FileInfo, ImgData
from Pipeline import Pipeline
from datetime import datetime

app_ui = ui.page_fluid(
    ui.HTML('<html data-bs-theme="dark">'),
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

            ui.input_slider("most_expr_genes", "Most expressed genes to plot", min=5, max=100, value=10),

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

            ui.nav("Genes / PCT Count (MT)",
                   ui.row(
                       ui.column(6, ui.output_image('plot_genes_pct_count')),
                       ui.column(6, ui.output_ui('desc_genes_pct_count'))),
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
                <h2> Welcome to BioNet </h2>
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
        return ui.HTML(
            """
            <br>
            <div class="container">
                <div class='jumbotron'>
                    <h2> Preprocessing </h2>
                    
                    <p>
                    The code <code>sc.pl.highest_expr_genes()</code> is used to create a plot that visually represents the 
                    expression levels of the top genes in a single-cell RNA sequencing dataset. In biological terms, 
                    this plot is valuable for identifying and exploring the most highly expressed genes across individual 
                    cells within the dataset. Biologically, the plot provides insights into the genes that play a crucial 
                    role in cellular functions or responses. Genes with high expression levels are often associated with key 
                    biological processes, and their identification can aid in understanding the molecular mechanisms 
                    underlying cellular diversity and function. By visualizing the highest expressing genes, researchers can:
                    </p>
                    
                    <br>
                    
                    <ol>
                        <li><b>Identify Biomarkers:</b>High expression levels of specific genes may indicate cell types, states, or conditions. This information is vital for identifying potential biomarkers associated with certain biological phenomena.</li>
                        <li><b>Characterize Cell Types</b>: The plot helps in characterizing different cell types based on their distinctive gene expression profiles. This is particularly useful in deciphering the heterogeneity within a cell population.</li>
                        <li><b>Uncover Biological Processes</b>: Genes with elevated expression levels often participate in critical biological processes. Examining the top expressing genes allows researchers to uncover and prioritize pathways and functions relevant to the studied cells.</li>
                        <li><b>Refine Experimental Design</b>: The information obtained from the plot can guide the design of further experiments or analyses, helping researchers focus on specific genes or pathways of interest.</li>
                    </ol>
                </div>
            </div>
            """
        )

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
            return ImgData(src=pipeline.highest_expr_genes_url,
                           style='display: block;margin-left: auto;margin-right: auto;width: 80%;')

    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_highly_variable_genes():
        return ui.HTML(
            """
            <br>
            <div class="container">
                <div class='jumbotron'>
                    <h2> Preprocessing </h2>

                    <p>
                    Citing from “Simple Single Cell” workflows <a href='https://master.bioconductor.org/packages/release/workflows/html/simpleSingleCell.html#examining-gene-level-metrics'>(Lun, McCarthy & Marioni, 2017)</a>:
                    
                    <blockquote cite="http://www.worldwildlife.org/who/index.html">
                    High proportions are indicative of poor-quality cells (Islam et al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic RNA from perforated cells. The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane.
                    </blockquote>
                    
                    The code <code>sc.pp.highly_variable_genes()</code> followed by 
                    <code>sc.pl.highly_variable_genes()</code> is employed in single-cell RNA sequencing analysis to 
                    identify and visualize highly variable genes within the dataset. Biologically, this process 
                    serves the following purposes:
                    </p>

                    <br>

                    <ol>
                        <li><b>Identification of Variable Genes</b>:
                            <ul>
                                <li>The function <code>sc.pp.highly_variable_genes()</code> assesses the mean expression and dispersion of genes in the dataset. Genes with variable expression levels across cells are identified based on specified criteria such as minimum mean expression <code>(min_mean)</code>, maximum mean expression <code>(max_mean)</code>, and minimum dispersion <code>(min_disp)</code>.</li>
                                <li>Variable genes are crucial for understanding cellular diversity, as they often contribute to the distinct functional states of individual cells.</li>
                            </ul>
                        </li>
                        <li><b>Characterization of Cellular Heterogeneity</b>: Highly variable genes are indicative of differences in expression patterns among cells. Analyzing these genes helps researchers characterize and quantify cellular heterogeneity within a population of cells.</li>
                        <li><b>Selection of Discriminatory Features</b>: The identification of highly variable genes assists in selecting discriminatory features for downstream analyses, such as dimensionality reduction or clustering. Focusing on genes with substantial expression variation enhances the ability to capture biologically relevant differences.</li>
                        <li><b>Enhancement of Statistical Power</b>: By filtering out less informative genes, the analysis of highly variable genes enhances the statistical power of downstream analyses. This is particularly important in uncovering meaningful patterns in datasets with large cellular heterogeneity.</li>
                        <li><b>Visualization of Variable Genes</b>: The plotting function sc.pl.highly_variable_genes generates a visual representation of the identified highly variable genes. This plot provides a quick overview of the distribution of mean expression versus dispersion for each gene, aiding in the interpretation of gene selection criteria.</li>

                    </ol>
                </div>
            </div>
            """
        )

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
            return ImgData(src=pipeline.highly_variable_genes_url,
                           style='display: block;margin-left: auto;margin-right: auto;width: 80%;')

    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_genes_pct_count():
        return ui.HTML(
            """
            <br>
            <div class="container">
                <div class='jumbotron'>
                    <h2> Preprocessing </h2>

                    <p>The code <code>sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')</code> generates a 
                    scatter plot in the context of single-cell RNA sequencing analysis. This plot specifically 
                    illustrates the relationship between total counts (representing overall gene expression) and the 
                    percentage of counts attributed to mitochondrial genes (pct_counts_mt).

                    Biologically, this scatter plot serves the following purposes:
                    </p>

                    <br>

                    <ol>
                        <li><b>Quality Control Assessment</b>:The scatter plot aids in assessing the quality of individual cells based on their gene expression characteristics. It specifically focuses on the interplay between total counts and the percentage of counts associated with mitochondrial genes.</li>
                        <li><b>Mitochondrial Content as a Quality Indicator</b>: Cells with an unusually high percentage of counts derived from mitochondrial genes may indicate poor-quality cells, potential cell stress, or contamination during sample preparation. Monitoring this metric helps identify and filter out cells that might introduce noise into downstream analyses.</li>
                        <li><b>Normalization Considerations</b>: The scatter plot assists in evaluating the impact of normalization procedures on the dataset. Deviations from the expected relationship between total counts and mitochondrial content may prompt researchers to explore or adjust normalization strategies to ensure accurate biological interpretations.</li>
                        <li><b>Identification of Cell Subpopulations</b>: Patterns in the scatter plot can reveal distinct subpopulations of cells with varying mitochondrial content. This information is valuable for understanding heterogeneity within the cell population.</li>
                        <li><b>Data Filtering and Preprocessing</b>: Researchers can use this scatter plot to guide data filtering and preprocessing steps. Cells falling outside acceptable ranges for mitochondrial content may be flagged for further investigation or excluded from downstream analyses.</li>
                    </ol>
                </div>
            </div>
            """
        )

    @output
    @render.image
    @reactive.event(input.plot_figures, ignore_none=True)
    def plot_genes_pct_count():
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
            return ImgData(src=pipeline.gene_pct_count_url,
                           style='display: block;margin-left: auto;margin-right: auto;width: 60%;')

    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_pca():
        return ui.HTML(
            """
            <br>
            <div class="container">
                <div class='jumbotron'>
                    <h2> Principal Component Analysis </h2>

                    <p>
                    The code <code>sc.pl.pca(adata)</code> generates a visualization of principal component 
                    analysis (PCA) in the context of single-cell RNA sequencing analysis. PCA is a dimensionality 
                    reduction technique that transforms high-dimensional gene expression data into a 
                    lower-dimensional space, capturing the most significant sources of variation.

                    Biologically, the PCA plot serves the following purposes:
                    </p>

                    <br>

                    <ol>
                        <li><b>Dimensionality Reduction</b>:PCA condenses the complex gene expression data into a reduced set of principal components, each representing a linear combination of genes. This facilitates the visualization of the most significant sources of variation within the dataset.</li>
                        <li><b>Cellular Heterogeneity</b>: The PCA plot reveals patterns of cellular heterogeneity by representing individual cells in the reduced principal component space. Distinct clusters or groupings of cells may indicate different cell types, states, or biological conditions.</li>
                        <li><b>Identification of Cell Populations</b>: Clustering of cells in the PCA plot can suggest the presence of distinct cell populations with similar transcriptional profiles. This aids in the identification and characterization of cell types or subpopulations within the analyzed sample.</li>
                        <li><b>Quality Control Assessment</b>: Anomalies or outliers in the PCA plot may indicate issues such as batch effects, technical artifacts, or suboptimal sample quality. Visual inspection of the plot allows for the identification and potential correction of such issues.</li>
                        <li><b>Informing Downstream Analyses</b>: Results from PCA can inform downstream analyses, such as clustering or trajectory inference. Clusters identified in the PCA space can guide subsequent analyses focused on understanding the biological significance of the observed cell groupings.</li>
                    </ol>
                </div>
            </div>
            """
        )

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
            return ImgData(src=pipeline.pca_url,
                           style='display: block;margin-left: auto;margin-right: auto;width: 60%;')

    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_pca_variance():
        return ui.HTML(
            """
            <br>
            <div class="container">
                <div class='jumbotron'>
                    <h2> Principal Component Analysis </h2>

                    <p>
                    The code <code>sc.pl.pca_variance_ratio(adata)</code> generates a visualization of the variance 
                    explained by each principal component in the context of single-cell RNA sequencing analysis. This 
                    plot provides information about the contribution of each principal component to the overall 
                    variance in the dataset. Biologically, the PCA variance ratio plot serves the following purposes: 
                    </p>

                    <br>

                    <ol>
                        <li><b>Variance Distribution Across Principal Components</b>: The plot displays the proportion of variance explained by each principal component. This information helps in understanding the relative importance of individual components in capturing the overall variability within the gene expression data.</li>
                        <li><b>Dimensionality Reduction Assessment</b>: Researchers can assess how effectively the top principal components summarize the dataset's variance. A steep drop-off in variance explained suggests that a reduced number of principal components may be sufficient to capture the majority of the dataset's variability.</li>
                        <li><b>Selection of Informative Principal Components</b>: By examining the variance ratio plot, researchers can identify the principal components that contribute most significantly to the dataset's variability. This knowledge aids in the selection of the most informative components for downstream analyses.</li>
                        <li><b>Quality Control</b>: Anomalies in the variance ratio plot, such as unexpected spikes or plateaus, may indicate issues with the data, such as batch effects or technical artifacts. Detecting such patterns allows for the identification and potential correction of data quality issues.</li>
                        <li><b>Informing Downstream Analyses</b>: Understanding the distribution of variance across principal components informs subsequent analyses, such as clustering or trajectory inference. It helps researchers decide how many components to consider for downstream exploration and interpretation.</li>
                    </ol>
                </div>
            </div>
            """
        )

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
            return ImgData(src=pipeline.pca_variance_url,
                           style='display: block;margin-left: auto;margin-right: auto;width: 60%;')

    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_umap():
        return ui.HTML(
            """
            <br>
            <div class="container">
                <div class='jumbotron'>
                    <h2> Embedding the Neighborhood Graph </h2>

                    <p> The code <code>sc.tl.umap(adata)</code> generates a Uniform Manifold Approximation and 
                    Projection (UMAP) visualization in the context of single-cell RNA sequencing analysis. UMAP is a 
                    dimensionality reduction technique that preserves local and global relationships between data 
                    points, providing an effective tool for visualizing complex biological structures within 
                    single-cell datasets.
                    
                    Biologically, the UMAP plot serves the following purposes:
                    </p>

                    <br>

                    <ol>
                        <li><b>Cellular Heterogeneity Visualization</b>: UMAP projects high-dimensional gene expression data into a lower-dimensional space, revealing the relationships between individual cells. The resulting plot provides a visually interpretable representation of cellular heterogeneity within the dataset.</li>
                        <li><b>Identification of Cell Clusters</b>: Clusters of cells with similar transcriptional profiles appear as distinct groupings in the UMAP plot. This aids in identifying and characterizing different cell types or states present in the analyzed sample.</li>
                        <li><b>Preservation of Local and Global Structure</b>: UMAP is particularly effective at preserving both local and global structure in the data. Cells with similar gene expression profiles are positioned closer together, maintaining meaningful relationships observed in the original high-dimensional space.</li>
                        <li><b>Trajectory Inference and Developmental Pseudotime</b>: UMAP can be used to infer trajectories and developmental pseudotime in single-cell datasets. By capturing the progression of cells in the reduced space, researchers gain insights into potential cellular transitions and developmental processes.</li>
                        <li><b>Biological Insights from Spatial Arrangement</b>: The spatial arrangement of cells in the UMAP plot may provide biological insights. Proximity in the plot can indicate functional relationships or shared regulatory programs between cells, contributing to a deeper understanding of cellular biology.</li>
                        <li><b>Quality Control Assessment</b>: Anomalies or unexpected patterns in the UMAP plot may indicate issues with the data, such as batch effects or technical artifacts. Visual inspection of the plot allows for the identification and potential correction of data quality issues.</li>
                    </ol>
                </div>
            </div>
            """
        )

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
            return ImgData(src=pipeline.umap_url,
                           style='display: block;margin-left: auto;margin-right: auto;width: 60%;')

    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_rank_genes_group():
        return ui.HTML(
            """
            <br>
            <div class="container">
                <div class='jumbotron'>
                    <h2> Finding Marker Genes </h2>

                    <p>
                    The code <code>sc.tl.rank_genes_groups(adata)</code> followed by 
                    <code>sc.pl.rank_genes_groups(adata)</code> is used in single-cell RNA sequencing analysis to 
                    identify and visualize differentially expressed genes between groups of cells. This process is 
                    commonly employed to uncover genes that are significantly associated with specific cell types, 
                    conditions, or experimental factors.

                    Biologically, these steps serve the following purposes:
                    </p>

                    <br>

                    <ol>
                        <li><b>Identification of Differentially Expressed Genes</b>: <code>sc.tl.rank_genes_groups</code> performs statistical tests to identify genes that are differentially expressed between predefined groups of cells. These groups may represent different cell types, conditions, or experimental factors.</li>
                        <li><b>Ranking Genes Based on Statistical Significance</b>: Genes are ranked based on their statistical significance, providing a list of genes associated with the specified groups. The ranking considers factors such as fold change and p-values.</li>
                        <li><b>Visualization of Differentially Expressed Genes</b>: <code>sc.pl.rank_genes_groups</code> generates plots to visualize the results of the differential expression analysis. This includes dot plots or similar visualizations that display the top-ranked genes and their associated statistics.</li>
                        <li><b>Cell Type Marker Identification</b>: Differentially expressed genes often serve as markers for specific cell types. This analysis helps in identifying genes that are characteristic of certain cell populations, contributing to the molecular characterization of cell types.</li>
                        <li><b>Biological Interpretation of Results</b>: Visual inspection of the rank_genes_groups plot allows researchers to interpret the biological significance of differentially expressed genes. The visualization highlights the genes that contribute most significantly to the observed differences between groups.</li>
                        <li><b>Quality Control Assessment</b>: The analysis can be used for quality control and validation, ensuring that the identified differentially expressed genes align with biological expectations and experimental design.</li>
                    </ol>
                </div>
            </div>
            """
        )

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
            return ImgData(src=pipeline.rank_genes_groups_url,
                           style='display: block;margin-left: auto;margin-right: auto;width: 80%;')

    @output
    @render.ui
    @reactive.event(input.plot_figures, ignore_none=True)
    def desc_rank_genes_group_violin():
        return ui.HTML(
            """
            <br>
            <div class="container">
                <div class='jumbotron'>
                    <h2> Finding Marker Genes </h2>

                    <p>
                    The code <code>sc.pl.rank_genes_groups_violin(adata)</code> is used in single-cell RNA 
                    sequencing analysis to visualize the expression levels of top-ranked differentially expressed 
                    genes identified by the <code>sc.tl.rank_genes_groups</code> analysis. This type of 
                    visualization, often referred to as a violin plot, provides insights into the distribution of 
                    gene expression across different cell groups or conditions.

                    Biologically, these steps serve the following purposes:
                    </p>

                    <br>

                    <ol>
                        <li><b>Visualizing Gene Expression Patterns</b>: The plot visualizes the expression patterns of top-ranked genes that are differentially expressed between predefined groups of cells. Each violin plot represents the distribution of expression levels for a specific gene across the defined cell groups.</li>
                        <li><b>Identification of Cell Type Markers</b>: Differentially expressed genes identified by the <code>sc.tl.rank_genes_groups</code> analysis are often cell type markers. The violin plot allows researchers to observe how the expression of these markers varies across different cell types or conditions.</li>
                        <li><b>Understanding Expression Variability</b>: Violin plots provide information about the variability in gene expression within each group. The width of the violin plot at a particular expression level indicates the density of cells with that expression level.</li>
                        <li><b>Comparison of Expression Levels</b>: Researchers can compare the expression levels of specific genes across different cell types or conditions. This comparison helps in understanding the molecular distinctions between cell populations and their associated gene expression signatures.</li>
                        <li><b>Biological Interpretation of Results</b>: Visual inspection of the violin plots allows for the interpretation of the biological significance of differentially expressed genes. Researchers can observe how the expression of these genes contributes to the observed phenotypic differences between groups.</li>
                        <li><b>Quality Control and Validation</b>: The violin plots can be used for quality control and validation, ensuring that the identified differentially expressed genes exhibit the expected expression patterns across experimental conditions.</li>
                    </ol>
                </div>
            </div>
            """
        )

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
            return ImgData(src=pipeline.rank_genes_groups_violin_url,
                           style='display: block;margin-left: auto;margin-right: auto;width: 80%;')

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
                                    name=dataset_name,
                                    min_num_genes_for_filtering=min_num_genes_for_filtering,
                                    min_num_cells_for_filtering=min_num_cells_for_filtering,
                                    num_neighbours=num_neighbours, num_pcs=num_pcs)

            elif dataset_name == 'WB Lysis Granulocytes 5p Introns 8kCells':
                pipeline = Pipeline(verbosity_lv=1,
                                    source_file_path='data/WB_Lysis_Granulocytes_5p_Introns_8kCells_filtered_feature_bc_matrix/filtered_feature_bc_matrix',
                                    name=dataset_name,
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

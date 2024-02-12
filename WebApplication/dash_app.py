import os
import shutil

from dash import Dash, html, dcc, Input, Output, no_update, callback
import dash_bootstrap_components as dbc
import dash_ag_grid as dag
import pandas as pd

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], serve_locally=True)

if os.path.exists('DEG_results'):
    shutil.rmtree('DEG_results')
    shutil.copytree('./../Cell_Annotation_and_DGE/DEG_results', 'DEG_results')
else:
    shutil.copytree('./../Cell_Annotation_and_DGE/DEG_results', 'DEG_results')

if os.path.exists('assets/figures'):
    shutil.rmtree('assets/figures')
    shutil.copytree('./../Cell_Annotation_and_DGE/figures', 'assets/figures')
else:
    shutil.copytree('./../Cell_Annotation_and_DGE/figures', 'assets/figures')


diseases_list = ['DiabetesII', 'MPN']

datasets = [file_name.replace('.csv', '') for file_name in os.listdir(f'./DEG_results/{diseases_list[0]}')]

datasets.sort()

columnDefs = [
    {"field": "Gene", 'text-align': 'center'},
    {"field": "Target Cell Type", 'text-align': 'center'},
    {"field": "Donor", 'text-align': 'center'},
    {"field": "Method"},
    {"field": "Score", "filter": "agNumberColumnFilter", 'text-align': 'center'},
    {"field": "pval_adj", "filter": "agNumberColumnFilter", 'text-align': 'center'},
    {"field": "LFC", "filter": "agNumberColumnFilter", 'text-align': 'center'},
]

first_table_grid = dag.AgGrid(
    className="ag-theme-alpine-dark",
    id="first_table_grid",
    rowData=None,
    columnDefs=columnDefs,
    defaultColDef={"resizable": True, "sortable": True, "filter": True, "Width": 125},
    columnSize="autoSize",
    # dashGridOptions={"pagination": True, 'paginationPageSize': 12, "domLayout": "autoHeight"},
    dashGridOptions={"pagination": True, 'AutoPaginationPageSize': True},
    style={"height": 600}
)

specific_cell_DEG_grid = dag.AgGrid(
    className="ag-theme-alpine-dark",
    id="second_table_grid",
    rowData=None,
    columnDefs=columnDefs,
    defaultColDef={"resizable": True, "sortable": True, "filter": True, "Width": 125, },
    columnSize="autoSize",
    # dashGridOptions={"pagination": True, 'paginationPageSize': 12, "domLayout": "autoHeight"},
    dashGridOptions={"pagination": True, 'AutoPaginationPageSize': True},
    style={"height": 600},
)

app.layout = dbc.Container([
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H1("BioNet Project", style={'color': 'white', 'text-align': 'center'}),
                    html.Div([
                        html.Br(),
                        html.H3('Select the disease to load analysis', style={'color': 'white', 'margin': '30px'}),
                        dcc.Dropdown(diseases_list, value=diseases_list[0], id='disease_name',
                                     style={'margin': '30px', 'width': '50%'}),
                        html.H3(id='disease_name_banner',
                                style={'color': 'white', 'margin': '30px'}),
                    ], style={'border-style': 'solid', 'border-color': 'red'}),

                ],
                width=12,
            ),
        ],
        style={"margin": 30}
    ),

    dbc.Row([
        dbc.Col([
            html.H4("First Dataset:", style={'color': 'white'}),
            dcc.Dropdown(datasets, value=datasets[0], id='first_table_dropdown'),

            html.Br(),
            html.Br(),

            dbc.Col([
                html.H4('SCSA annotated cells vs. Panglaodb', style={'color': 'white'}),
                dbc.Row([html.Img(id='first_SCSA_cellType_annotation')]),
                html.Br(),
                html.H4('Cell Distribution between Case and Control', style={'color': 'white'}),
                dbc.Row([html.Img(id='first_donor_cells')]),
                html.Br(),
                html.H4('Cell Types based on SCSA', style={'color': 'white'}),
                dbc.Row([html.Img(id='first_embedding_cell_types')]),
                html.Br(),
                html.H4('Calculated Ro/e between Case and Control', style={'color': 'white'}),
                dbc.Row([html.Img(id='first_roe')]),
            ]),

            html.Br(),
            html.H6("Go to page:", style={'color': 'white'}),
            dcc.Input(id="first_table_page_input", type="number", style={'font-size': 'large'}, value=1),
            html.Br(),
            html.Br(),
            first_table_grid,
        ],
            width=6,
            style={"margin": 50, 'width': 550}),

        dbc.Col([
            html.H4("Second Dataset:", style={'color': 'white'}),
            dcc.Dropdown(datasets, value=datasets[1], id='second_table_dropdown'),

            html.Br(),
            html.Br(),

            dbc.Col([
                html.H4('SCSA annotated cells vs. Panglaodb', style={'color': 'white'}),
                dbc.Row([html.Img(id='second_SCSA_cellType_annotation')]),
                html.Br(),
                html.H4('Cell Distribution between Case and Control', style={'color': 'white'}),
                dbc.Row([html.Img(id='second_donor_cells')]),
                html.Br(),
                html.H4('Cell Types based on SCSA', style={'color': 'white'}),
                dbc.Row([html.Img(id='second_embedding_cell_types')]),
                html.Br(),
                html.H4('Calculated Ro/e between Case and Control', style={'color': 'white'}),
                dbc.Row([html.Img(id='second_roe')]),
            ]),

            html.Br(),
            html.H6("Go to page:", style={'color': 'white'}),
            dcc.Input(id="second_table_goto_page_input", type="number", style={'font-size': 'large'}, value=1),
            html.Br(),
            html.Br(),
            specific_cell_DEG_grid,
        ],
            width=6,
            style={"margin": 50, 'width': 550}),

    ]),
])


@callback(
    [Output('first_table_dropdown', 'options'),
     Output('second_table_dropdown', 'options'),
     Output('disease_name_banner', 'children')],
    Input('disease_name', 'value')
)
def update_disease(disease_name):
    datasets = [file_name.replace('.csv', '') for file_name in os.listdir(f'./DEG_results/{disease_name}')]
    datasets.sort()
    return datasets, datasets, f'Results for Differential Gene Analysis for {disease_name} Disease'


@callback(
    [Output('first_table_grid', 'rowData'),
     Output('first_SCSA_cellType_annotation', 'src'),
     Output('first_donor_cells', 'src'),
     Output('first_embedding_cell_types', 'src'),
     Output('first_roe', 'src')],
    [Input('first_table_dropdown', 'value'),
     Input('disease_name', 'value')]
)
def update_first_table(dataset_name, disease_name):
    dataset = pd.read_csv(f'./DEG_results/{disease_name}/{dataset_name}/DEG_All.csv')
    plot_base_url = f'assets/figures/X_mde/{disease_name}/{dataset_name}'
    return dataset.to_dict("records"), f'{plot_base_url}/SCSA_cellType_annotation.png', f'{plot_base_url}/donor_cells.png', f'{plot_base_url}/embedding_cell_types.png', f'{plot_base_url}/ROE.png'

@callback(
    [Output('second_table_grid', 'rowData'),
     Output('second_SCSA_cellType_annotation', 'src'),
     Output('second_donor_cells', 'src'),
     Output('second_embedding_cell_types', 'src'),
     Output('second_roe', 'src')],
    [Input('second_table_dropdown', 'value'),
     Input('disease_name', 'value')]
)
def update_second_table(dataset_name, disease_name):
    dataset = pd.read_csv(f'./DEG_results/{disease_name}/{dataset_name}/DEG_All.csv')
    plot_base_url = f'assets/figures/X_mde/{disease_name}/{dataset_name}'
    return dataset.to_dict("records"), f'{plot_base_url}/SCSA_cellType_annotation.png', f'{plot_base_url}/donor_cells.png', f'{plot_base_url}/embedding_cell_types.png', f'{plot_base_url}/ROE.png'


@callback(
    Output("first_table_grid", "paginationGoTo"),
    Input("first_table_page_input", "value"),
)
def update_page_size(goto):
    if goto is None:
        return no_update
    # grid page numbers start at zero
    return goto - 1


@callback(
    Output("second_table_grid", "paginationGoTo"),
    Input("second_table_goto_page_input", "value"),
)
def update_page_size(goto):
    if goto is None:
        return no_update
    # grid page numbers start at zero
    return goto - 1


if __name__ == "__main__":
    app.run_server(debug=True)

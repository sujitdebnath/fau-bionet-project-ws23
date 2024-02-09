import os

from dash import Dash, html, dcc, Input, Output, no_update, callback
import dash_bootstrap_components as dbc
import dash_ag_grid as dag
import pandas as pd
from dash.exceptions import PreventUpdate

# all_cell_DEG_df = pd.read_csv("DEG_All.csv")
# specific_cell_DEG_df = pd.read_csv("DEG_T cell.csv")

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], serve_locally=True)

dataset_files = [file_name.replace('.csv', '') for file_name in os.listdir('./resulted_DEG_files')]

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
                        html.H3("Results for Differential Gene Analysis for MPN Disease",
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
            dcc.Dropdown(dataset_files, value=dataset_files[0], id='first_table_dropdown'),
            html.Br(),
            html.H6("Go to page:", style={'color': 'white'}),
            dcc.Input(id="first_table_page_input", type="number", style={'font-size': 'large'}, value=1),
            html.Br(),
            html.Br(),
            first_table_grid,
        ],
            width=6,
            style={"margin": 50, 'width': 500}),

        dbc.Col([
            html.H4("Second Dataset:", style={'color': 'white'}),
            dcc.Dropdown(dataset_files, value=dataset_files[1], id='second_table_dropdown'),
            html.Br(),
            html.H6("Go to page:", style={'color': 'white'}),
            dcc.Input(id="second_table_goto_page_input", type="number", style={'font-size': 'large'}, value=1),
            html.Br(),
            html.Br(),
            specific_cell_DEG_grid,
        ],
            width=6,
            style={"margin": 50, 'width': 500}),

    ]),
])


@callback(
    Output('first_table_grid', 'rowData'),
    Input('first_table_dropdown', 'value')
)
def update_first_table(value):
    first_dataset = pd.read_csv(f'./resulted_DEG_files/{str(value)}.csv')
    return first_dataset.to_dict("records")


@callback(
    Output('second_table_grid', 'rowData'),
    Input('second_table_dropdown', 'value')
)
def update_first_table(value):
    first_dataset = pd.read_csv(f'./resulted_DEG_files/{str(value)}.csv')
    return first_dataset.to_dict("records")


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

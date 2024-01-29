from dash import Dash, html, dcc, Input, Output, no_update, callback
import dash_bootstrap_components as dbc
import dash_ag_grid as dag
import pandas as pd

all_cell_DEG_df = pd.read_csv("DEG_All.csv")
specific_cell_DEG_df = pd.read_csv("DEG_T cell.csv")

app = Dash(__name__, serve_locally=True)

columnDefs = [
    {"field": "Gene", 'text-align': 'center'},
    {"field": "Target Cell Type", 'text-align': 'center'},
    {"field": "Donor", 'text-align': 'center'},
    {"field": "Method"},
    {"field": "Score", "filter": "agNumberColumnFilter", 'text-align': 'center'},
    {"field": "pval_adj", "filter": "agNumberColumnFilter", 'text-align': 'center'},
    {"field": "LFC", "filter": "agNumberColumnFilter", 'text-align': 'center'},
]

all_DEG_grid = dag.AgGrid(
    className="ag-theme-alpine-dark",
    id="all_cell_DEG_grid",
    rowData=all_cell_DEG_df.to_dict("records"),
    columnDefs=columnDefs,
    defaultColDef={"resizable": True, "sortable": True, "filter": True, "Width": 125, },
    columnSize="sizeToFit",
    # dashGridOptions={"pagination": True, 'paginationPageSize': 12, "domLayout": "autoHeight"},
    dashGridOptions={"pagination": True, 'AutoPaginationPageSize': True},
    style={"height": 600}
)

specific_cell_DEG_grid = dag.AgGrid(
    className="ag-theme-alpine-dark",
    id="specific_cell_DEG_grid",
    rowData=specific_cell_DEG_df.to_dict("records"),
    columnDefs=columnDefs,
    defaultColDef={"resizable": True, "sortable": True, "filter": True, "Width": 125, },
    columnSize="sizeToFit",
    # dashGridOptions={"pagination": True, 'paginationPageSize': 12, "domLayout": "autoHeight"},
    dashGridOptions={"pagination": True, 'AutoPaginationPageSize': True},
    style={"height": 600}
)

app.layout = dbc.Container([
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H1("BioNet Project", style={'color': 'white', 'text-align': 'center'}),
                    html.Div([
                        html.H2("Results for Differential Gene Analysis", style={'color': 'white', 'margin': '30px'}),
                        html.P(
                            'Note: This is a test over PBMCs datasets collected from three different donors. In the '
                            'close future, it will be replaced with Case vs. Control datasets.',
                            style={'color': 'white', 'margin': '30px'}),
                    ], style={'border-style': 'solid', 'border-color': 'red'}),
                ],
                width=12,
                style={"width": "100%"}
            ),
        ],
        style={"margin": 50}
    ),

    dbc.Row([
        dbc.Col([
            html.H1("DEG without Cell Specification", style={'color': 'white'}),
            html.H2("Go to page:", style={'color': 'white'}),
            dcc.Input(id="all_cell_goto-page-input", type="number", style={'font-size': 'large'}),
            html.Br(),
            html.Br(),
            all_DEG_grid,
        ],
            width=12,
            style={"margin": 50})
    ]),

    dbc.Row([
        dbc.Col([
            html.H1("DEG with Cell Specification", style={'color': 'white'}),
            html.H2("Go to page:", style={'color': 'white'}),
            dcc.Input(id="specific_cell_goto-page-input", type="number", style={'font-size': 'large'}),
            html.Br(),
            html.Br(),
            specific_cell_DEG_grid,
        ],
            width=12,
            style={"margin": 50})
    ]),

#    dbc.Row(
#        [
#            dbc.Col(
#                [
#                    html.Img(src='assets/figures/X_mdeSCSA_cellType_annotation.png', style={'width': '44%'}),
#                    html.Img(src='assets/figures/tracksplotrank_genes_group_per_donor.png', style={'width': '45%'}),
#                ]),
#        ],
#        style={"margin": 50, 'text-align': 'center'}
#    ),
#
#    dbc.Row(
#        [
#            dbc.Col(
#                [
#                    html.Img(src='assets/figures/X_mdeDonor_cells.png', style={"width": "30%"}),
#                    html.Img(src='assets/figures/CellType_Donor_ratio.png', style={"width": "30%"}),
#                ])
#        ],
#        style={"margin": 50, 'text-align': 'center'}
#    ),

])


@callback(
    Output("all_cell_DEG_grid", "paginationGoTo"),
    Input("all_cell_goto-page-input", "value"),
)
def update_page_size(goto):
    if goto is None:
        return no_update
    # grid page numbers start at zero
    return goto - 1


@callback(
    Output("specific_cell_DEG_grid", "paginationGoTo"),
    Input("specific_cell_goto-page-input", "value"),
)
def update_page_size(goto):
    if goto is None:
        return no_update
    # grid page numbers start at zero
    return goto - 1


if __name__ == "__main__":
    app.run_server(debug=True)

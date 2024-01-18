from dash import Dash, html, dcc, Input, Output, no_update, callback
import dash_ag_grid as dag
import pandas as pd

df = pd.read_csv("DEG_all.csv")

app = Dash(__name__, serve_locally=True)

grid = dag.AgGrid(
    className="ag-theme-alpine-dark",
    id="quickstart-grid",
    rowData=df.to_dict("records"),
    columnDefs=[{"field": i, 'text-align': 'center'} for i in df.columns],
    defaultColDef={"resizable": True, "sortable": True, "filter": True, "Width": 125, },
    columnSize="sizeToFit",
    # dashGridOptions={"pagination": True, 'paginationPageSize': 12, "domLayout": "autoHeight"},
    dashGridOptions={"pagination": True, 'AutoPaginationPageSize': True},
    style={"height": 600}
)

app.layout = html.Div(
    [
        html.H1("BioNet Project", style={'color': 'white'}),
        html.H2("Results for Differential Gene Analysis", style={'color': 'white'}),
        html.P('This is a test over PBMCs datasets collected from three different donors. In the close future, '
               'it will be replaced with Case vs. Control datasets.', style={'color': 'white'}),
        html.H2("Go to page:", style={'color': 'white'}),
        dcc.Input(id="goto-page-input", type="number", style={'font-size': 'large'}),
        html.Br(),
        html.Br(),
        grid,
    ], style={"margin": 50},
)


@callback(
    Output("quickstart-grid", "paginationGoTo"),
    Input("goto-page-input", "value"),
)
def update_page_size(goto):
    if goto is None:
        return no_update
    # grid page numbers start at zero
    return goto - 1


if __name__ == "__main__":
    app.run_server(debug=True)

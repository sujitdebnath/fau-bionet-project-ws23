import pandas as pd
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo

app_ui = ui.page_fluid(
    ui.input_text("txt", "Enter the text to display below:", "delete me"),
    ui.row(
        ui.output_text("text"),
        ui.input_file("file1", "Choose CSV File", accept=[".csv"], multiple=False),
    ),
    ui.output_table("summary"),
)


def server(input: Inputs, output: Outputs, session: Session):
    @output
    @render.text
    def text():
        return input.txt()


    @output
    @render.table
    def summary():
        file = input.file1()
        d = pd.read_csv(file[0]["datapath"])
        print(d)
        # print(file[0]["datapath"])

app = App(app_ui, server)
app.run()
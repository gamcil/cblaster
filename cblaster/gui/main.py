"""A basic GUI for cblaster."""

import sys

import PySimpleGUI as sg

from cblaster import __version__
from cblaster import main
from cblaster.gui import search, makedb, citation


sg.theme("Lightgrey1")


def Column(layout, scrollable=False):
    return sg.Column(
        layout,
        scrollable=scrollable,
        size=(540, 520),
        vertical_scroll_only=True
    )


def run_cblaster(values):
    """Handles conversion of PySimpleGUI values to cblaster parameters.

    - Know which workflow tab we're on (search or makedb)
    - search
        - Know which search mode tab we're on
        - if remote, use entrez query, database, RID
        - if local, jdb and database

    Args:
        values (dict): Dictionary of values from PySimpleGUI.
    """
    if values["cblaster_tabs"] == "Search":
        args = dict(
            query_file=values["query_file"],
            query_ids=values["query_ids"],
            session_file=values["session_file"],
            mode=values["search_mode"],
            gap=int(values["gap"]),
            unique=int(values["unique"]),
            min_hits=int(values["min_hits"]),
            require=values["require"],
            min_identity=float(values["min_identity"]),
            min_coverage=float(values["min_coverage"]),
            max_evalue=float(values["max_evalue"]),
            recompute=values["recompute"],
        )

        if values["search_mode"] == "remote":
            args.update(
                database=values["database"],
                entrez_query=values["entrez_query"],
                rid=values["rid"]
            )
        else:
            args.update(
                database=values["dmnd_database"],
                json_db=values["json_db"]
            )

        if values["summary_gen"]:
            summary = sys.stdout

            if values["summary_text"]:
                summary = values["summary_text"]

            args.update(
                output=summary,
                output_human=values["summary_hr"],
                output_headers=values["summary_he"]
            )

        if values["binary_gen"]:
            args.update(
                binary=values["binary_text"],
                binary_human=values["binary_hr"],
                binary_headers=values["binary_he"]
            )

        if values["figure_gen"]:
            figure = values["figure_text"] if values["figure_text"] else True
            args.update(
                figure=figure,
                figure_dpi=values["figure_spin"],
                use_plotly=values["figure_plotly"]
            )

        main.cblaster(**args)

    elif values["cblaster_tabs"] == "Makedb":
        main.makedb(
            genbanks=values["makedb_genbanks"].split(";"),
            filename=values["makedb_filename"],
            indent=values["json_indent"]
        )
    else:
        raise ValueError("Expected 'Search' or 'Makedb'")


def cblaster_gui():
    layout = [
        [sg.Text("cblaster", font="Arial 18 bold", pad=(0, 0))],
        [sg.Text(f"v{__version__}", font="Arial 10", pad=(0, 0))],
        [sg.Text("Cameron Gilchrist, 2020", font="Arial 10", pad=(0, 0))],
        [sg.TabGroup([
            [sg.Tab("Search", [[Column(search.layout, scrollable=True)]])],
            [sg.Tab("Makedb", [[Column(makedb.layout)]])],
            [sg.Tab("Citation", [[Column(citation.layout)]])],
        ], enable_events=True, key="cblaster_tabs"
        )],
        [sg.Button("Start", key="start_button", button_color=["white", "green"]),
         sg.Button("Exit", key="exit_button", button_color=["white", "red"])],
    ]

    window = sg.Window(
        "cblaster",
        layout,
        size=(600, 700),
        element_padding=(5, 5),
        element_justification="center",
        finalize=True
    )

    while True:
        event, values = window.read()

        if event in (None, "exit_button"):
            break

        # Disable binary & summary table, figure options if not enabled
        for key in ("browse", "text", "hr", "he"):
            window[f"binary_{key}"].update(disabled=not values["binary_gen"])

        for key in ("browse", "text", "hr", "he"):
            window[f"summary_{key}"].update(disabled=not values["summary_gen"])

        for key in ("browse", "text", "spin", "plotly"):
            window[f"figure_{key}"].update(disabled=not values["figure_gen"])

        # Disable start button when on citation tab
        window["start_button"].update(
            disabled=values["cblaster_tabs"] not in ("Search", "Makedb")
        )

        if event:
            if event == "start_button":
                run_cblaster(values)

    window.close()


if __name__ == "__main__":
    cblaster_gui()

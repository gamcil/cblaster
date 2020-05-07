import PySimpleGUI as sg

from cblaster.gui.parts import TextLabel, Frame


sg.theme("Lightgrey1")

local_tab = sg.Tab("Local", [
    [TextLabel("JSON database"),
     sg.In(default_text="e.g. cblaster.json", size=(34, 1), key="json_db"),
     sg.FileBrowse(key="json_db")],
    [TextLabel("DIAMOND database"),
     sg.InputText(
         default_text="e.g. cblaster.dmnd",
         size=(34, 1),
         key="dmnd_database"
     ),
     sg.FileBrowse(key="dmnd_database")],
], key="local")

remote_tab = sg.Tab("Remote", [
    [TextLabel("Database"),
     sg.InputText(default_text="e.g. nr", key="database")],
    [TextLabel("Entrez query"),
     sg.InputText(
         default_text='e.g. "Aspergillus"[organism]',
         key="entrez_query"
     )],
    [TextLabel("Request Identifier (RID)"),
     sg.InputText(key="rid")],
], key="remote")

search_tabgroup = sg.TabGroup([[remote_tab], [local_tab]], key="search_mode")

search_frame = Frame(
    "Search",
    key="searching_frame",
    layout=[
        [sg.Text(
            "Specify the search mode and databases to be used in the cblaster run."
            " In remote mode, the database value should correspond to a BLAST"
            " database hosted by the NCBI. In local mode, the database arguments"
            " should refer to files generated using cblaster makedb.",
            size=(71, 4))],
        [search_tabgroup]
    ]
)

input_frame = Frame(
    "Input",
    key="input_frame",
    layout=[
        [sg.Text(
            "Specify the protein sequences that you want to search. These can"
            " be provided by either using a FASTA file or entering the NCBI"
            " accessions of sequences. Alternatively, a session file generated"
            " in a previous cblaster run can be loaded so that you do not have"
            " to repeat a search.",
            size=(71, 4),
        )],
        [TextLabel("File"),
         sg.InputText(size=(34, 1), key="query_file"),
         sg.FileBrowse(key="query_file")],
        [TextLabel("Sequence IDs"), sg.InputText(key="query_ids")],
        [TextLabel("Session file"), sg.InputText(key="session_file")],
    ],
)

clustering_frame = Frame(
    "Clustering",
    key="clustering_frame",
    layout=[
        [sg.Text(
             "Specify the conditions used when identifying clusters of hits on"
             " genomic scaffolds.",
             size=(71, 1))],
        [TextLabel("Max. intergenic gap (bp)"),
         sg.InputText(default_text="20000", key="gap")],
        [TextLabel("Min. unique query hits"),
         sg.InputText(default_text="3", key="unique")],
        [TextLabel("Min. hits in clusters"),
         sg.InputText(default_text="3", key="min_hits")],
        [TextLabel("Required sequences"),
         sg.InputText(key="require")],
    ],
)

filtering_frame = Frame(
    "Filtering",
    key="filtering_frame",
    layout=[
        [sg.Text("Specify hit quality thresholds used when filtering BLAST results")],
        [TextLabel("Max. e-value"),
         sg.InputText(default_text="0.01", key="max_evalue")],
        [TextLabel("Min. identity (%)"),
         sg.InputText(default_text="50", key="min_identity")],
        [TextLabel("Min. query coverage (%)"),
         sg.InputText(default_text="70", key="min_coverage")],
        [TextLabel("Recompute"),
         sg.InputText(key="recompute")],
    ],
)

summary_frame = Frame(
    "Summary table",
    key="summary_frame",
    layout=[
        [sg.Text(
            "This is the standard cblaster results table that is shown at the"
            " end of each run. To save this table to a file, pick a file path"
            " using the option below. If no path is provided, the table will"
            " be printed in the terminal.",
            size=(71, 3)
        )],
        [TextLabel("Generate summary table"),
         sg.Checkbox("", key="summary_gen", default=True, enable_events=True)],
        [TextLabel("Table formatting"),
         sg.Checkbox("Human-readable", default=True, key="summary_hr"),
         sg.Checkbox("Show headers", default=True, key="summary_he")],
        [TextLabel("Output file"),
         sg.In(key="summary_text", size=(34, 1)),
         sg.FileSaveAs(key="summary_browse")],
    ],
)

binary_frame = Frame(
    "Binary table",
    key="binary_frame",
    layout=[
        [sg.Text(
            "The binary table will give you an overview of the absence/presence of"
            " query genes in the hit clusters identified in the search. To generate"
            " this table, please provide a file name below.",
            size=(71, 3)
        )],
        [TextLabel("Generate binary table"),
         sg.Checkbox("", default=False, enable_events=True, key="binary_gen")],
        [TextLabel("Table formatting"),
         sg.Checkbox("Human-readable", key="binary_hr", disabled=True),
         sg.Checkbox("Show headers", key="binary_he", disabled=True)],
        [TextLabel("Output file"),
         sg.InputText(key="binary_text", disabled=True, size=(34, 1)),
         sg.FileSaveAs(key="binary_browse", disabled=True)],
    ],
)

figure_frame = Frame(
    "Figure",
    key="figure_frame",
    layout=[
        [sg.Text(
            "Generate a cblaster plot using matplotlib. This provides a visual"
            " representation of the binary table. If no file path is specified,"
            " the figure will be opened in the interactive viewer from the"
            " matplotlib library. Note that when the results"
            " set gets large, matplotlib tends to choke, in which case"
            " you can try using the plot.ly library instead, which will"
            " open an interactive plot in your web browser.",
            size=(71, 5)
        )],
        [TextLabel("Generate cblaster plot"),
         sg.Checkbox("", key="figure_gen", default=False, enable_events=True)],
        [TextLabel("Use plot.ly"), sg.Checkbox("", key="figure_plotly", disabled=True)],
        [TextLabel("Figure DPI"),
         sg.Spin(
            list(range(600)),
             initial_value=300,
             key="figure_spin",
             disabled=True
        )],
        [TextLabel("Output file"),
         sg.InputText(key="figure_text", disabled=True, size=(34, 1)),
         sg.FileSaveAs(key="figure_browse", disabled=True)],
    ],
)

layout = [
    [input_frame],
    [search_frame],
    [filtering_frame],
    [clustering_frame],
    [summary_frame],
    [binary_frame],
    [figure_frame],
]

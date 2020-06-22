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
        [TextLabel("Session file"),
         sg.InputText(size=(34, 1), key="session_file"),
         sg.FileBrowse(key="session_file")],
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
         sg.InputText(default_text="30", key="min_identity")],
        [TextLabel("Min. query coverage (%)"),
         sg.InputText(default_text="50", key="min_coverage")],
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
        [TextLabel("Delimiter"),
         sg.InputText(key="summary_delimiter", size=(34, 1))],
        [sg.Text(
            "Character used to delimit values in the summary table. If no delimiter,"
            " is specified, the table will be generated in human-readable format.",
            size=(71, 2)
        )],
        [TextLabel("Hide headers"), sg.Checkbox("", key="summary_hide_headers")],
        [sg.Text(
            "Hide all headers in the summary table. This includes organism and scaffold"
            " headers, as well as headers in the hit table.",
            size=(71, 2)
        )],
        [TextLabel("Decimal places"),
         sg.Spin(
            list(range(6)),
             initial_value=2,
             key="summary_decimals",
             disabled=True
         )],
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
        [TextLabel("Delimiter"),
         sg.InputText(key="binary_delimiter", disabled=True, size=(34,1))],
        [sg.Text(
            "Character used to delimit values in the binary table. If no delimiter,"
            " is specified, the table will be generated in human-readable format.",
            size=(71, 2)
        )],
        [TextLabel("Hide headers"),
         sg.Checkbox("", key="binary_hide_headers", disabled=True)],
        [sg.Text("Hide all headers in the binary table.")],
        [TextLabel("Key function"),
         sg.Drop(key="binary_key", disabled=True, values=("len", "sum", "max"))],
        [sg.Text(
            "This specifies a function used to compute the values given in the binary"
            " table. By default, this is 'len', which will calculate the 'length' of"
            " the list of hits for a given query sequence (i.e. cell counts). 'sum'"
            " or 'max' can be used to give the sum or max of hit attributes specified"
            " using the option below (e.g. cumulative identity).",
            size=(71, 4)
        )],
        [TextLabel("Hit attribute"),
         sg.Drop(
             key="binary_attr",
             disabled=True,
             values=("identity", "coverage", "bitscore", "evalue")
         )],
        [sg.Text(
            "This specifies the type of score value of a hit to use when computing"
            " cell values in the binary table. By default, percentage identity will"
            " be used.",
            size=(71, 2)
        )],
        [TextLabel("Decimal places"),
         sg.Spin(
            list(range(6)),
             initial_value=2,
             key="binary_decimals",
             disabled=True
         )],
        [TextLabel("Output file"),
         sg.InputText(key="binary_text", disabled=True, size=(34, 1)),
         sg.FileSaveAs(key="binary_browse", disabled=True)],
    ],
)

figure_frame = Frame(
    "Figure",
    key="figure_frame",
    layout=[
        [TextLabel("Generate cblaster plot"),
         sg.Checkbox("", key="figure_gen", default=False, enable_events=True)],
        [sg.Text(
            "This generates a visual representation of the binary table as a"
            " cluster heatmap. If no file path is specified,"
            " the figure will be served to an IP address and automatically opened,"
            " at which point the figure can be manipulated and saved as SVG."
            " If a file path is specified, a static HTML file will be generated at"
            " that path.",
            size=(71, 4)
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

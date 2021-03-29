import PySimpleGUI as sg

from cblaster.gui.parts import TextLabel, Frame, TEXT_WIDTH


sg.theme("Lightgrey1")

local_tab = sg.Tab("Local", [
    [TextLabel("DIAMOND database"),
     sg.InputText(
         default_text="e.g. cblaster.dmnd",
         size=(34, 1),
         key="dmnd_database"
     ),
     sg.FileBrowse(key="dmnd_database")],
    [TextLabel("Number of CPUs"),
     sg.InputText(key="cpus")],
    [sg.Text("If no value is supplied all available CPU's will be used.")]
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

    [TextLabel("Maximum hits"),
     sg.InputText(default_text=5000, key="max_hits")],
    [sg.Text("Maximum total hits to save from a remote BLAST search. Setting"
             " this value too low may result in missed hits/clusters.",
             size=(TEXT_WIDTH, 2))]
], key="remote")

hmm_tab = sg.Tab("Hmm", [
    [TextLabel("FASTA Database"),
     sg.InputText(
         default_text="e.g cblaster.fasta",
         key="fa database"
     ),
     sg.FileBrowse(key="fa database")],
    [TextLabel("Pfam database"),
     sg.InputText(
         default_text='directory with \'Pfam-A.hmm.dat.gz\' in it',
         key="pfam database"
     ),
     sg.FolderBrowse(key="pfam database")],
], key="hmm")

combi_local_tab = sg.Tab("Hmm and Local", [
    [TextLabel("FASTA Database"),
     sg.InputText(
         default_text="e.g cblaster.fasta",
         key="fa database cl",
         size=(34, 1)
     ),
     sg.FileBrowse(key="fa database cl")],
    [TextLabel("Pfam database"),
     sg.InputText(
         default_text='directory with \'Pfam-A.hmm.dat.gz\' in it',
         key="pfam database cl",
         size=(34, 1)
     ),
     sg.FolderBrowse(key="pfam database cl")],
    [TextLabel("DIAMOND database"),
     sg.InputText(
         default_text="e.g. cblaster.dmnd",
         size=(34, 1),
         key="dmnd_database cl"
     ),
     sg.FileBrowse(key="dmnd_database cl")],
    [TextLabel("Number of CPUs"),
     sg.InputText(default_text="1", key="cpus cl")],
], key="combi_local")

combi_remote_tab = sg.Tab("Hmm and remote", [
    [TextLabel("FASTA Database"),
     sg.InputText(
         default_text="e.g cblaster.fasta",
         key="fa database cr",
         size=(34, 1)
     ),
     sg.FileBrowse(key="fa database cr")],
    [TextLabel("Pfam database"),
     sg.InputText(
         default_text='e.g. directory with \'Pfam-A.hmm.dat.gz\' in it',
         key="pfam database cr",
         size=(34, 1)
     ),
     sg.FolderBrowse(key="pfam database cr")],
    [TextLabel("Database"),
     sg.InputText(default_text="e.g. nr", key="database cr")],
    [TextLabel("Entrez query"),
     sg.InputText(
         default_text='e.g. "Aspergillus"[organism]',
         key="entrez_query cr"
     )],
    [TextLabel("Request Identifier (RID)"),
     sg.InputText(key="rid cr")],
], key="combi_remote")

search_tabgroup = sg.TabGroup([[remote_tab], [local_tab], [hmm_tab], [combi_local_tab], [combi_remote_tab]],
                              key="search_mode")

search_frame = Frame(
    "Search",
    key="searching_frame",
    layout=[
        [sg.Text(
            "Specify the search mode and databases to be used in the cblaster "
            "run. In remote mode, the database value should correspond to a BLAST "
            "database hosted by the NCBI. In local mode, the database arguments "
            "should refer to files generated using cblaster makedb. When using any "
            "of the HMM modes, a local copy of the Pfam database will be stored at "
            "the indicated location or extracted from there. The FASTA database "
            "should refer to the FASTA file generated using cblaster makedb.",
            size=(TEXT_WIDTH, 6))],
        [search_tabgroup]
    ]
)

input_frame = Frame(
    "Input",
    key="input_frame",
    layout=[
        [sg.Text(
            "Specify the protein sequences that you want to search. These can be "
            "provided by either using a FASTA file or entering the NCBI accessions "
            "of sequences. When running any of the HMMs, at least one HMM profile "
            "has to be defined. Alternatively, a session file generated in a previous "
            "cblaster run can be loaded so that you do not have to repeat a search",
            size=(TEXT_WIDTH, 6),
        )],
        [TextLabel("File"),
         sg.InputText(size=(34, 1), key="query_file"),
         sg.FileBrowse(key="query_file")],
        [TextLabel("Sequence IDs"), sg.InputText(key="query_ids", default_text="e.g. space separated values")],
        [TextLabel("HMM profiles"), sg.InputText(key="query_profiles", default_text="e.g. space separated values")],
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
             size=(TEXT_WIDTH, 1))],
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
         sg.Checkbox("", key="recompute_gen", enable_events=True),
         sg.In(key="recompute_text", size=(28, 1), disabled=True, enable_events=True),
         sg.FileSaveAs(key="recompute_browse", disabled=True)],
        [sg.Text(
            "Recompute a previous search session using new thresholds. The "
            "filtered session will be written to the file specified by this argument. "
            "If this argument is specified with no value, the session will be "
            "filtered but not saved (e.g. for plotting purposes).",
            size=(TEXT_WIDTH, 3)
        )],
    ],
)

summary_frame = Frame(
    "Summary table",
    key="summary_frame",
    layout=[
        [sg.Text(
            "This is the standard cblaster results table that is shown at the "
            "end of each run. To save this table to a file, pick a file path using "
            "the option below. If no path is provided, the table will be printed"
            " in the terminal",
            size=(TEXT_WIDTH, 3)
        )],
        [TextLabel("Generate summary table"),
         sg.Checkbox("", key="summary_gen", default=True, enable_events=True)],
        [TextLabel("Delimiter"),
         sg.InputText(key="summary_delimiter", size=(34, 1))],
        [sg.Text(
            "Character used to delimit values in the summary table. If no delimiter,"
            " is specified, the table will be generated in human-readable format.",
            size=(TEXT_WIDTH, 2)
        )],
        [TextLabel("Hide headers"), sg.Checkbox("", key="summary_hide_headers")],
        [sg.Text(
            "Hide all headers in the summary table. This includes organism and scaffold"
            " headers, as well as headers in the hit table.",
            size=(TEXT_WIDTH, 2)
        )],
        [TextLabel("Decimal places"),
         sg.Spin(
            list(range(6)),
             initial_value=2,
             key="summary_decimals",
             disabled=True
         )],

        [TextLabel("Sort clusters"),  sg.Checkbox("", key="sort_clusters")],
        [sg.Text(
            "Sorts the clusters of the final output on score. This means that clusters"
            "of the same organism are not necessarily close together in the output.",
            size=(TEXT_WIDTH, 2)
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
            size=(TEXT_WIDTH, 3)
        )],
        [TextLabel("Generate binary table"),
         sg.Checkbox("", default=False, enable_events=True, key="binary_gen")],
        [TextLabel("Delimiter"),
         sg.InputText(key="binary_delimiter", disabled=True, size=(34,1))],
        [sg.Text(
            "Character used to delimit values in the binary table. If no delimiter,"
            " is specified, the table will be generated in human-readable format.",
            size=(TEXT_WIDTH, 2)
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
            size=(TEXT_WIDTH, 4)
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
            size=(TEXT_WIDTH, 2)
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
            " cluster heatmap. If this is not"
            " specified, the plot will be saved in a plot.html in the temporary folder,"
            " at which point the figure can be manipulated and saved as SVG."
            " If a file path is specified, a static HTML file will be generated at"
            " that path.",
            size=(TEXT_WIDTH, 4)
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

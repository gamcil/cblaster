import PySimpleGUI as sg

from cblaster.gui.parts import TextLabel, TEXT_WIDTH


sg.theme("Lightgrey1")


extract_clusters_frame = sg.Frame(
    "Extract Clusters",
    layout=[
        [sg.Text("This module allows you to extract clusters from a session.json file into a genbank file"
                 " for each separate cluster. The genbank files can be formatted to include qualifiers "
                 "to make them readable by bigscape.",
                 size=(TEXT_WIDTH, 5)
                 )],
        [TextLabel("Session file"),
         sg.InputText(size=(34, 1), key="extract_clusters_session"),
         sg.FileBrowse(key="extract_clusters_session")],
        [sg.Text(
            "A session file (.json) generated by a cblaster search.",
            size=(TEXT_WIDTH, 1)
        )],

        [TextLabel("Output directory"),
         sg.InputText(key="extract_clusters_output", size=(34, 1)),
         sg.FolderBrowse(key="extract_clusters_output")],
        [sg.Text(
            "Directory the extracted clusters will be saved in.",
            size=(TEXT_WIDTH, 1)
        )],

        [TextLabel("Prefix"),
         sg.InputText(key="prefix", size=(34, 1))
         ],
        [sg.Text(
            "Start of the name for each cluster file, the base name is cluster'clutser.number'",
            size=(TEXT_WIDTH, 1)
        )],

        [TextLabel("Output format"),
         sg.Drop(key="output format", default_value="genbank", values=("genbank", "bigscape"))],
        [sg.Text(
            "The format of the resulting files. The options are genbank and bigscape.",
            size=(TEXT_WIDTH, 1)
        )],

        [TextLabel("Clusters"), sg.InputText(key="clusters ec")],
        [sg.Text("Cluster numbers/ ranges provided by the summary file of the 'search' command. "
                 "For example to include clusters 1 to 4 use '1-4'. Multiple values can be"
                 " supplied separated by spaces.",
                 size=(TEXT_WIDTH, 3)
                 )],

        [TextLabel("Score threshold"), sg.InputText(key="score threshold ec")],
        [sg.Text("The minimum required score of a cluster in order to be extracted.",
                 size=(TEXT_WIDTH, 1)
                 )],

        [TextLabel("Organisms"), sg.InputText(key="organisms ec")],
        [sg.Text(
            "Organisms that extracted clusters must be from. These take the form"
            " of regular expression patterns and are therefore quite flexible."
            " You can provide more than one pattern."
            " For example, to extract sequences only from Aspergillus and Penicillium"
            " genomes, you might specify: 'Aspergillus.*' 'Penicillium.*'"
            " See the user guide for more examples. Multiple values can be supplied separated by"
            " spaces.",
            size=(TEXT_WIDTH, 5)
        )],

        [TextLabel("Scaffolds"), sg.InputText(key="scaffolds ec")],
        [sg.Text(
            "Scaffolds that extracted clusters must be on. These can be scaffold"
            " names or names AND coordinate ranges. For example, you could specify"
            " scaffold_1, which would retrieve ALL clusters on scaffold_1, or"
            " scaffold_1:10000-50000, which would retrieve only those from position"
            " 10000 to 50000. Multiple values can be supplied separated by spaces.",
            size=(TEXT_WIDTH, 5)
        )],

        [TextLabel("Maximum clusters"), sg.InputText(key="max clusters ec", default_text="50")],
        [sg.Text(
            "The maximum amount of clusters that will be extracted. Ordered on score.",
            size=(TEXT_WIDTH, 1)
        )],
    ],
    title_color="blue",
    font="Arial 10 bold",
    relief="flat",
)

layout = [[extract_clusters_frame]]

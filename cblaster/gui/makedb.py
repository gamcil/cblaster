import PySimpleGUI as sg

from cblaster.gui.parts import TextLabel


sg.theme("Lightgrey1")


makedb_frame = sg.Frame(
    "makedb",
    layout=[
        [sg.Text(
            "This workflow will allow you to create a local"
            " database for cblaster from a collection of genome"
            " files in GenBank or GFF3+FASTA formats. It will"
            " produce a formatted DIAMOND database (.dmnd)"
            " containing protein sequences in the supplied"
            " genomes, as well as a JSON file which stores"
            " genomic coordinates for these proteins. To use"
            " these files in a cblaster run, choose a local"
            " search, then supply the DIAMOND and JSON files"
            " using the Database and JSON Database fields,"
            " respectively.",
            size=(72, 6)
        )],
        [TextLabel("Select genome files"),
         sg.InputText(size=(34, 1), key="genbanks"),
         sg.FilesBrowse(key="makedb_genbanks")],
        [sg.Text(
            "Both the JSON and DIAMOND databases will take the name"
            " specified here, albeit with their own file type suffixes.",
            size=(73, 2)
        )],
        [TextLabel("Database name"),
         sg.InputText(key="makedb_filename", size=(34, 1)),
         sg.FileSaveAs(key="makedb_filename")],
        [TextLabel("JSON indent size"),
         sg.Spin(list(range(0, 21)), initial_value=0, key="json_indent")],
    ],
    title_color="blue",
    font="Arial 10 bold",
    relief="flat",
)

layout = [[makedb_frame]]

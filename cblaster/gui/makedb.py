import PySimpleGUI as sg

from cblaster.gui.parts import TextLabel, TEXT_WIDTH


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
            " genomes, as well as a sqlite3 file which stores"
            " genomic coordinates for these proteins. To use"
            " these files in a cblaster run, choose a local"
            " search, then supply the DIAMOND (.dmnd) file"
            " using the Database field.",
            size=(TEXT_WIDTH, 8)
        )],
        [TextLabel("Select genome files"),
         sg.InputText(size=(34, 1), key="genbanks"),
         sg.FilesBrowse(key="makedb_genbanks")],
        [sg.Text(
            "Select on or more genbank/gff3 files.",
            size=(TEXT_WIDTH, 1)
        )],

        [TextLabel("Database name"),
         sg.InputText(key="makedb_filename", size=(34, 1)),
         sg.FileSaveAs(key="makedb_filename")],
        [sg.Text(
            "Both the sqlite3 and DIAMOND databases will take the name"
            " specified here, albeit with their own file type suffixes.",
            size=(TEXT_WIDTH, 2)
        )],

        [TextLabel("Number of CPUs"),
         sg.InputText(key="cpus db")],
        [sg.Text("If no value is supplied all available CPU's will be used.")],

        [TextLabel("Batch size"),
         sg.InputText(key="batch size")],
        [sg.Text("Number of genome files to parse before saving them in the local"
                 " database. Useful when encountering memory issues with large/many"
                 " files. By default, all genome files will be parsed at once.",
                 size=(TEXT_WIDTH, 4))],

        [TextLabel("Force overwrite files"),
         sg.Checkbox("", default=False, enable_events=True, key="force")],
        [sg.Text("Overwrite pre-existing files, if any are present.")],
    ],
    title_color="blue",
    font="Arial 10 bold",
    relief="flat",
)

layout = [[makedb_frame]]

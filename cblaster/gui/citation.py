import PySimpleGUI as sg

from cblaster.gui.parts import TextLabel


sg.theme("Lightgrey1")


def get_citation(title, citation):
    return [
        [sg.Text(title, font="Arial 10 bold")],
        [sg.Multiline(citation, size=(71, 2))]
    ]


citation_frame = sg.Frame(
    "Citation",
    layout=[
        [sg.Text("If you found cblaster useful, please cite:")],
        *get_citation("cblaster", "Gilchrist et al. XX (2020)."),
        [sg.Text("cblaster makes use of the following tools:")],
        *get_citation(
            "DIAMOND",
            "Buchfink, B., Xie, C. & Huson, D. H."
            " Fast and sensitive protein alignment using DIAMOND."
            " Nat. Methods 12, 59–60 (2015).",
        ),
        *get_citation(
            "NCBI BLAST API",
            "Acland, A. et al."
            " Database resources of the National Center for Biotechnology Information."
            " Nucleic Acids Res. 42, 7–17 (2014)",
        )
    ],
    title_color="blue",
    font="Arial 10 bold",
    relief="flat",
)


layout = [[citation_frame]]

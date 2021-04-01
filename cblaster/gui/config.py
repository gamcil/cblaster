import PySimpleGUI as sg

from cblaster.gui.parts import TextLabel, TEXT_WIDTH


sg.theme("Lightgrey1")


config_frame = sg.Frame(
    "Config",
    layout=[
        [sg.Text(
            "This module will allow you to configure cblaster. "
            "At the moment, this is mostly used so that you can identify "
            "yourself to the NCBI before using their services. This is an "
            "anti-abuse measure to prevent people from submitting too many "
            "queries at a time (>3 per second, or >10 per second with an API key). "
            "This module must be run before you can do a remote cblaster search.",
            size=(TEXT_WIDTH, 6)
        )],
        [TextLabel("E-mail address"), sg.InputText(size=(34, 1), key="email")],
        [TextLabel("NCBI API Key"), sg.InputText(size=(34, 1), key="api_key")],
        [TextLabel("Maximum retries"), sg.InputText(size=(34, 1), key="max_tries")],
    ],
    title_color="blue",
    font="Arial 10 bold",
    relief="flat",
)

layout = [[config_frame]]

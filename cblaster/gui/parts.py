"""Common components shared between layouts."""


import PySimpleGUI as sg


TEXT_WIDTH = 65


def SectionLabel(text):
    return sg.Text(text, justification="l", font="Arial 12 bold")


def TextLabel(text):
    return sg.Text(
        f"{text}:",
        justification="l",
        size=(22, 1),
        font="Arial 10 bold"
    )


def Frame(text, key, layout):
    return sg.Frame(
        text,
        layout=layout,
        key=key,
        font="Arial 10 bold",
        title_color="blue",
        relief="flat",
    )

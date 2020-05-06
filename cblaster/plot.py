import numpy as np
import scipy
from matplotlib import pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree


def ply(session):
    """Plot a session using Plotly.

    Use this when data is large since matplotlib will choke.
    """

    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import plotly.figure_factory as ff

    names, scafs, counts, identities = session.form_matrices(html=True)

    counts = np.array(counts)
    identities = np.array(identities)

    fig = make_subplots(rows=1, cols=2)

    dendro = ff.create_dendrogram(identities, orientation="right")
    leaves = dendro["layout"]["yaxis"]["ticktext"]
    leaves = list(map(int, leaves))

    for item in dendro["data"]:
        item.pop("xaxis", None)
        item.pop("yaxis", None)
        fig.append_trace(item, 1, 1)

    identities = identities[leaves, :]
    names = [names[i] for i in leaves]
    scafs = [scafs[i] for i in leaves]

    heatmap = go.Heatmap(
        y=[int(v) for v in dendro["layout"]["yaxis"]["tickvals"]],
        z=identities,
        xaxis="x2",
        colorscale="Blues",
        colorbar={"title": "Identity (%)", "len": 0.2},
        ygap=1,
        xgap=1
    )

    fig.add_trace(heatmap, row=1, col=2)

    yvals = [int(v) for v in dendro["layout"]["yaxis"]["tickvals"]]

    fig.update_layout(
        height=15*len(yvals),
        showlegend=False,
        title_text="cblaster",
        autosize=True,
        margin=dict(l=50, r=50, b=100, t=100, pad=4),
        template="simple_white"
    )
    fig.update_yaxes(matches="y", range=[yvals[0] - 10, yvals[-1] + 10], automargin=True)
    fig.update_yaxes(
        tickvals=yvals,
        ticktext=names,
        row=1,
        col=1,
        showline=False,
        ticks="",
        side="right",
        automargin=True,
    )
    fig.update_xaxes(
        domain=[0, 0.2],
        row=1,
        col=1,
        ticks="",
        showline=False,
        showticklabels=False,
        automargin=True,
    )
    fig.update_xaxes(
        domain=[0.6, 1],
        tickvals=list(range(len(session.queries))),
        ticktext=session.queries,
        side="top",
        row=1,
        col=2,
        ticks="",
        showline=False,
        automargin=True,
    )
    fig.update_yaxes(
        tickvals=yvals,
        ticktext=scafs,
        showline=False,
        row=1,
        col=2,
        ticks="",
        automargin=True,
    )

    fig.show(config={"toImageButtonOptions": {"width": None, "height": None}})


def plot(session, figure=None, dpi=300, show_counts=False):
    """Plot a cblaster Session using `matplotlib`.

    `matplotlib` will produce nicer figures than `plot.ly`, but will choke when the
    session becomes large. When the `Session` object contains a large number of organisms,
    the `ply()` function should be used.

    Parameters
    ----------
    session: cblaster.models.Session
        A `Session` object containing `cblaster` results
    figure: str
        Path to write generated figure to. This should include a file extension, which
        matplotlib will use to automatically determine the file type to save in.
    dpi: int
        Resolution of generated figure in dots per inch (DPI)
    show_counts: bool
        Draw hit counts over the hit heatmap
    """

    rcParams["savefig.dpi"] = dpi

    names, scafs, counts, identities = session.form_matrices()

    counts = np.array(counts)
    identities = np.array(identities)

    # Rough initial figure size calculations
    # qname = longest query name
    qname = max(len(q) for q in session.queries) * 0.06
    width = 6 + 0.45 * len(session.queries) + qname
    height = max(2, 0.3 * len(names) + qname)

    fig, (dendro, matrix) = plt.subplots(
        1,
        2,
        figsize=(width, height),
        gridspec_kw={"width_ratios": [1, 0.3 * len(session.queries)]},
    )

    if len(names) > 1:
        # Plot dendrogram
        Y = linkage(identities, method="ward")
        Z = dendrogram(Y, orientation="left", ax=dendro, labels=names, leaf_font_size=9)

        # Hide borders
        for spine in dendro.spines.values():
            spine.set_visible(False)

        dendro.set_xticks([])

        # Plot matrix
        index = Z["leaves"]
        counts = counts[index, :]
        identities = identities[index, :]
        identities[identities == 0.0] = np.nan
        scafs = [scafs[i] for i in index]
        im = matrix.matshow(
            identities,
            cmap="Blues",
            aspect="auto",
            origin="lower",
            clim=(0, 100)
        )

        # Set yticklabels to scaffold locations
        matrix.set_yticks(range(len(scafs)))
        matrix.set_yticklabels(scafs, fontsize=9)
    else:
        # Only have one result
        name = f"{names[0]} {scafs[0]}"
        im = matrix.matshow(
            identities,
            cmap="Blues",
            aspect="auto",
            origin="lower",
            clim=(0, 100)
        )
        dendro.set_visible(False)
        matrix.set_yticks([0, 1])
        matrix.set_yticklabels([name], fontsize=9)

    # Annotate with counts
    if show_counts:
        colours = ["black", "white"]
        textkw = dict(fontsize=6, va="center", ha="right")

        for i in range(counts.shape[0]):
            for j in range(counts.shape[1]):
                # Choose black or white, depending if identity is greater than 80%
                # False -> 0, True -> 1, then use as index
                textkw["color"] = colours[int(im.norm(identities[i, j]) > 0.8)]
                im.axes.text(j, i, counts[i, j], **textkw)

    # Hide actual tick marks
    matrix.tick_params(axis="x", which="both", bottom=False)
    matrix.tick_params(axis="x", which="minor", top=False)
    matrix.tick_params(axis="y", which="both", left=False)

    # Set xticklabels to query protein IDs
    matrix.set_xticks(range(identities.shape[1]))
    matrix.set_xticklabels(session.queries, rotation=30, ha="left", fontsize=9)

    # Set minor ticks to be between 'cells', set white grid for separation
    matrix.set_xticks(np.arange(identities.shape[1] + 1) - 0.5, minor=True)
    matrix.set_yticks(np.arange(identities.shape[0] + 1) - 0.5, minor=True)
    matrix.grid(which="minor", color="w", linestyle="-", linewidth=2)

    # Set dividers on both axes so colorbar can be used without changing figure size
    div1 = make_axes_locatable(dendro)
    div2 = make_axes_locatable(matrix)
    cax1 = div1.append_axes("bottom", size=0.1, pad=0.1)
    cax2 = div2.append_axes("bottom", size=0.1, pad=0.1)
    cax1.axis("off")

    # Create colourbar
    cbar = plt.colorbar(im, cax=cax2, orientation="horizontal")
    cbar.ax.set_xlabel("Identity (%)", va="bottom", labelpad=12)

    plt.tight_layout()

    if not figure:
        fig.canvas.mpl_connect("resize_event", resize)
        plt.show()
    else:
        plt.savefig(figure)


def resize(event):
    plt.tight_layout()
    plt.gcf().canvas.draw()

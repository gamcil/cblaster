"""Plot cblaster results."""


import numpy as np
import scipy
from matplotlib import pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import dendrogram, linkage


def plot(session, figure=None, dpi=300, show_counts=False):
    """Plot a cblaster Session."""

    rcParams["savefig.dpi"] = dpi

    names, scafs, counts, identities = session.form_matrices()

    counts = np.array(counts)
    identities = np.array(identities)

    # Rough initial figure size calculations
    # qname = longest query name
    qname = max(len(q) for q in session.queries) * 0.06
    width = 6 + 0.45 * len(session.queries) + qname
    height = max(2, 0.2 * len(names) + qname)

    fig, (dendro, matrix) = plt.subplots(
        1,
        2,
        figsize=(width, height),
        gridspec_kw={"width_ratios": [1, 0.3 * len(session.queries)]},
    )

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
        identities, cmap="Blues", aspect="auto", origin="lower", clim=(0, 100)
    )

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

    # Set yticklabels to scaffold locations
    matrix.set_yticks(range(len(scafs)))
    matrix.set_yticklabels(scafs, fontsize=9)

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

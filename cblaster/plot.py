"""Plot cblaster results."""


from matplotlib import pyplot as plt, rcParams

import numpy as np

import scipy
from scipy.cluster.hierarchy import linkage, dendrogram


def plot(session, figure=None, dpi=300, show_counts=False):
    """Plot a cblaster Session."""

    rcParams["savefig.dpi"] = dpi

    names, scafs, counts, identities = session.form_matrices()

    fig, (dendro, matrix) = plt.subplots(1, 2, figsize=(10, 5))

    counts = np.array(counts)
    identities = np.array(identities)

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
    scafs = [scafs[i] for i in index]
    im = matrix.matshow(identities, cmap="Blues", aspect="auto", origin="lower")

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

    # Create colourbar
    cbar = fig.colorbar(im, ax=matrix, shrink=0.5)
    cbar.ax.set_ylabel("Identity (%)", rotation=-90, va="bottom")

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

    # Run resize() every time the viewer is resized
    plt.tight_layout()

    if not figure:
        fig.canvas.mpl_connect("resize_event", resize)
        plt.show()
    else:
        plt.savefig(figure)


def resize(event):
    plt.tight_layout()
    plt.gcf().canvas.draw()

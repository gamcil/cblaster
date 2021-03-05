Plotting extracted clusters using ``plot_clusters``
===================================================

By default, the visualisation offered by ``cblaster`` shows only a heatmap of query hits
per result cluster. While this is very useful for quickly identifying patterns in large
datasets, we generally still want to see how these clusters compare in a more
biologically relevant way.

The ``plot_clusters`` module allows you to do precisely this. Given a session and some
filters to choose specific clusters (exactly like in the ``extract_clusters`` module),
this module will automatically extract the clusters, then generate an interactive
visualisation showing each cluster to-scale using ``clinker`` (doi: 10.1093/bioinformatics/btab007,
https://github.com/gamcil/clinker).

Example usage
-------------

Minimum working example:

::

        $ cblaster plot_clusters session.json

Plot clusters 1-10 and cluster 25 (these numbers can be found in the summary file of the 'search' command):

::

        $ cblaster plot_clusters session.json -c 1-10 25 -o plot.html

Plot only from specific organisms (regular expressions):

::

        $ cblaster plot_clusters session.json -or "Aspergillus.*" "Penicillium.*" -o plot.html

Plot only clusters from a specific range on scaffold_123 and all clusters on scaffold_234 (note: assumes unique scaffold names):

::

        $ cblaster plot_clusters session.json -sc scaffold_123:1-80000 scaffold_234 -o plot.html

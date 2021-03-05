Extracting GenBank files from session files using ``extract_clusters``
======================================================================

A common next step after a ``cblaster`` search is to retrieve the identified gene
clusters so we can perform additional analysis. ``cblaster`` provides the
``extract_clusters`` module precisely for this purpose, allowing you to generate GenBank
files of specific gene clusters directly from a session file.
This works for sessions from both remote and local searches: for remote searches,
clusters are downloaded directly from the NCBI, and in local searches, from the SQL
database generated using the ``makedb`` module.

Example usage
-------------

Extract all clusters from a session (can take a long time for remote searches with many
results):

::

        $ cblaster extract_clusters session.json -o example_directory

Extract clusters 1-10 and cluster 25 (these numbers can be found in the summary file of the 'search' command):

::

        $ cblaster extract_clusters session.json -c 1-10 25 -o example_directory

Extract clusters only from specific organisms (regular expressions):

::

        $ cblaster extract_clusters session.json -or "Aspergillus.*" "Penicillium.*" -o example_directory

Extract clusters only from a specific range on scaffold_123 and all clusters on scaffold_234 (note: expects unique scaffold names):

::

        $ cblaster extract_clusters session.json -sc scaffold_123:1-80000 scaffold_234 -o example_directory

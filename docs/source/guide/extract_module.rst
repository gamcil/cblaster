Retrieving hit sequences with the ``extract`` module
====================================================

After a search has been performed, it can be useful to retrieve sequences matching a certain query for further analyses (e.g. sequence comparisons for phylogenies).
This is easily accomplished using ``cblaster``'s ``extract`` module.

This module takes a ``cblaster`` session file as input, then extracts sequences matching any filters you have specified.
If no filters are specified, ALL hit sequences will be extracted.
However, that's probably not too useful, so instead we could extract all hit sequences matching a query sequence:

::

        $ cblaster extract session.json -q "Query1"

By default, only sequence names are extracted.
This is because ``cblaster`` stores no actual sequence data for hit sequences during it's normal search workflow, only their coordinates.
However, sequences can automatically be retrieved from the NCBI by specifying the ``-d`` or ``--download`` argument.
``cblaster`` will then write them, in FASTA format, to either the command line or a file.
For example, we can do the same command as above, but retrieve the sequences and write them to ``output.fasta`` like so:

::

        $ cblaster extract session.json -q "Query1" -d -o output.fasta

Note that the ``-o`` or ``--output`` argument has been used here; this will write any results from the ``extract`` module to the specified file.

You can also provide multiple names of query sequences:

::

        $ cblaster extract session.json -q Query1 Query2 Query3 ...

Note, however, that all extracted sequences will be written to the same file.

The ``extract`` module can also filter based on the organism or scaffold that each hit sequence is on.
The organism filter uses regular expression patterns based on organism names.
Multiple patterns can be provided, and are additive (i.e. any organism matching any of the patterns will be saved).
For example, you could filter a search session of all fungal organisms on NCBI for only those sequences from *Aspergillus* or *Penicillium* species like so:

::

        $ cblaster extract fungi.json -or "Aspergillus.*" "Penicillium.*"

Note that patterns should be enclosed in quotation marks in order to be read in correctly.

The scaffold filter is less flexible, capable of matching exact scaffolds or scaffold ranges.
For example, to extract hit sequences on a scaffold, ``scaffold_1``, from position 10000 to 23000:

::

        $ cblaster extract session.json -sc "scaffold_1:10000-23000"

Like the organism filter, multiple scaffolds and/or scaffold ranges can be provided and they are additive.

By default, source information is added to each sequence name, for example:

::

        sequence [organism=Source organism] [scaffold=scaffold_1:123-456]

This can be turned off using the ``-no`` or ``--names_only`` argument.

Finally, the `extract` module can also generate delimited table files, for easy importing into spreadsheet programs.
For example, to generate a comma-delimited table (CSV file), simple provide the ``-de`` or ``--delimiter`` argument:

::

        $ cblaster extract session.json ... -de ","

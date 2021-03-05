.. _makedb_module:

Creating local sequence databases with the ``makedb`` module
============================================================

The ``makedb`` module is used to generate the databases used in local ``cblaster`` searches from genome files.
``makedb`` only takes two arguments: the genome files being used to build the databases, and some name to use when saving them.
For example, generating a database from a set of genomes is as simple as:

::

        $ cblaster makedb one.gbk two.gbk three.gbk four.gbk myDb

This will read in each GenBank file, then generate the files:

================ ==================================================================
``myDb.sqlite3`` Local database used for looking up genomic context of hits
``myDb.dmnd``    DIAMOND sequence search database
``myDb.fasta``   All protein sequences parsed from genomes; used for HMMER searches
================ ==================================================================

``cblaster`` can also build databases from GFF3 files; however, currently the FASTA sequence must be embedded within the GFF3 (i.e. under a ``##FASTA`` directive).
Typically it is easiest to have all your genome files within a folder and use a wildcard to avoid having to type every file name, like so:

::

        $ cblaster makedb genomes/*.gbk myDb

The shell will expand this automatically into a command that is functionally equivalent to the previous one. 
However, on Windows, we have run into some issues with this behaviour. For windows, instead use the command:

::

        $ cblaster makedb (ls *.gbk | \% FullName) myDb

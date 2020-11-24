Running a ``cblaster`` search using the ``search`` module
=========================================================

cblaster can also perform fully local searches, forgoing the need for any interaction
with NCBI whatsoever.
To do this, cblaster builds a JSON 'database' file from a list of GenBank files,
and generates a DIAMOND database of all protein sequences in this database.
Then, cblaster searches can be run against the created DIAMOND database, and genomic
context obtained from the JSON file.
cblaster provides the module, ``makedb``, which handles the creation of these databases
in one easy command:

::

  $ cblaster makedb folder/*.gbk mydatabase

This command will read in all files ending in ``.gbk`` inside ``folder``, then create
the DIAMOND and cblaster databases to ``mydatabase.dmnd`` and ``mydatabase.json``,
respectively.

Your new database can then be searched like so:

::

  $ cblaster search -mode local -qf query.fasta -db mydatabase.dmnd -jdb mydatabase.json

All other search parameters used in your remote searches can also be used here.

Obtaining genomes for a local database
--------------------------------------
As a working example, we'll download all RefSeq *Aspergillus* genomes from the NCBI and
create a local cblaster database.
The easiest way to do this is via the NCBI's batch assembly download option.

1. Search the NCBI Assembly database for *Aspergillus* genomes

.. image:: /_static/search.png
  :width: 500
  :alt: Search for Aspergillus assemblies

2. Click 'Download Assemblies', select 'Genomic GenBank format (.gbff)' and click 'Download'

.. image:: /_static/download.png
  :width: 400
  :alt: Download 'Genomic GenBank format (.gbff)' files

3. Extract all GenBank files

::

  $ tar -xf genome_assemblies_genome_gb.tar
  $ mv ncbi-genomes-2020-06-23/*.gz -t myfolder/
  $ pigz -d myfolder/*.gz

NCBI generates a tar archive containing the genomes. I first untar the file, then copy
all compressed GenBank files (ending with ``.gz``) to a folder of my choosing and
decompress them using pigz.

4. Create the DIAMOND and cblaster databases using ``makedb``

::

  $ cblaster makedb myfolder/*.gbff mydatabase
  [14:45:50] INFO - Starting cblaster makedb
  [14:45:50] INFO - Parsing 40 files...
  [14:45:50] INFO - 1. myfolder/GCF_000002655.1_ASM265v1_genomic.gbff
  [14:45:53] INFO - 2. myfolder/GCF_000002715.2_ASM271v1_genomic.gbff
  [14:45:55] INFO - 3. myfolder/GCF_000002855.3_ASM285v2_genomic.gbff
  [14:45:58] INFO - 4. myfolder/GCF_000006275.2_JCVI-afl1-v2.0_genomic.gbff
  ...
  [14:58:02] INFO - Writing FASTA file with database sequences: mydatabase.faa
  [14:58:02] INFO - Building DIAMOND database: mydatabase.dmnd
  [14:58:12] INFO - Building JSON database: mydatabase.json
  [14:58:59] INFO - Done.

5. Run `cblaster` against the newly created databases

::

  $ cblaster search -m local -qf query.fa -db mydatabase.dmnd -jdb mydatabase.dmnd

Alternatively, you could use a tool like `ncbi-genome-download`_ to retrieve the sequences
from the command line.

.. _`ncbi-genome-download`: https://github.com/kblin/ncbi-genome-download

Note that there is absolutely no requirement that the genomes come from NCBI; if you
have a local collection of GenBank files, you can easily make the cblaster databases
with those too.

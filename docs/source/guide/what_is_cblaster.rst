
What is ``cblaster``
====================

``cblaster`` is a tool for finding clusters of co-located sequences using BLAST searches.

By providing ``cblaster`` with amino acid sequences of interest, it can search a sequence database to find instances of genes encoding related proteins clustered on a specific scaffold. 
The software was developed as a spiritual successor to software such as `MultiGeneBlast <http://multigeneblast.sourceforge.net/>`_ and, to aid in the discovery of natural product biosynthetic gene clusters in bacteria and fungi.

It can, of course, be used to search for any group of clustered genes.

``cblaster`` offers a number of useful functions:


	* ``cblaster search`` allows users to query a FASTA file containing protein sequences against the NCBI database or a local database to identify instances of collocation.  A number of graphical and textual outputs are available for quick analysis or use in downstream applications.
	* ``cblaster makedb`` enables the user to rapidly create local databases from GenBank files for future searches.
	* ``cblaster gne`` can be used to analyse the effect of varying the intergenic distances on the search output. This can be used to evaluate the robustness of outputs and gain insights into the distribution of gene cluster sizes.
	* ``cblaster extract`` allows the user to quickly produce FASTA files containing the protein sequences of homologues of interest for downstream analysis.

Detailed information on how to install ``cblaster`` and use the above modules can be found in this guide.


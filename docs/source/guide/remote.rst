.. _remote:

Remotely search NCBI BLAST databases
====================================

At a minimum, a search could look like one of the following:

::

  $ cblaster search --query_file query.fasta 
  $ cblaster search --query_ids QBE85649.1 QBE85648.1 QBE85647.1 QBE85646.1 ...
 
This will launch a remote search against the non-redundant (nr) protein database,
retrieve and parse the results, then report any blocks of hits.
If you wanted to search against a different NCBI database, you can provide one using the
``--database`` argument, like so:

::

  $ cblaster search -qf query.fasta -db refseq_protein

Note that the short forms of ``--query_file`` (``-qf``) and ``--database`` (``-db``)
have been used here. Most commands in cblaster have a short form, which can be found
in the help menus (``cblaster search -h``).

By default, hits are only reported if they are above 30% identity and 50% query
coverage, and have an e-value below 0.01.
If we wanted to be stricter, we could change those values with their corresponding
arguments, for example:

::

  $ cblaster search -qf query.fasta --min_identity 70 --min_coverage 90 --evalue 0.001

You can also pass in NCBI search queries using `-eq / --entrez_query` to pre-filter
the target database, which can result in vastly reduced run-times and more
specific results. For example, to only search against *Aspergillus* sequences:

::

  $ cblaster search -qf query.fasta --entrez_query "Aspergillus"[ORGN]

Look here_ for a full description of Entrez search terms.

.. _here: https://www.ncbi.nlm.nih.gov/books/NBK49540/

Another useful feature of cblaster is the ability to retrieve results of a previous
remote BLAST search. Every search started on NCBI is automatically assigned a unique
request identifier (RID) which remains valid for up to 24 hours after the completion of
that search. cblaster is able to take an RID using the ``--rid`` argument, like so:

::

  $ cblaster search -qf query.fasta --rid WHS0UGYJ015

This will lookup the RID, retrieve the results if it is valid, and resume the rest of
the cblaster pipeline from there. 

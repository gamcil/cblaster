Miscellaneous functions
=======================

Accessing help dialogues
------------------------

Every module in the ``cblaster`` command line interface has a useful help dialogue which details the arguments you may specify.
To do this, simply type \texttt{-h/--help} as the sole argument after any command.
For example:

::

	  $ cblaster -h
	  $ cblaster search -h
	  $ cblaster gne -h
	  $ cblaster makedb -h
	  $ cblaster extract -h

Getting the current ``cblaster`` version.
To check the version of cblaster simply enter the command:

::

	$ cblaster --version

JSON indent level
-----------------

``cblaster`` uses JSON files at several stages in its pipelines (i.e. search sessions, context databases).
By default, these files are saved with no indentation level for the purpose of lowering file size. This means that all data is stored on a single line and is hard to read for anything but computers (e.g. human beings!).
In most cases, this is perfectly fine; however, sometimes you may wish to investigate something manually within a search session and would require some level of indenting to make the file readable.
This is achieved by using the ``-i`` or ``--indent`` argument.
The ``-i`` argument belongs to the first level of the command line interface, and therefore must be used directly after ``cblaster`` in the command line string.
So, this will **not** work:

::

	  $ cblaster search -qf query.fasta -s session.json -i 2

But **this will**:
	
::
	  $ cblaster -i 2 search -qf query.fasta -s session.json

This will set the indentation level of a given session to 2, meaning that for every new line in the file, 2 spaces will be drawn per indentation level.
For example, a session file with no indent (truncated) looks like this:

::
	{"queries": ["BuaB", "BuaC", "BuaD", "QBE85644.1", "BuaE", "BuaF", "BuaG",
	"QBE85648.1", "BuaA"], "params": {"mode": "remote", "database": "nr",
	"min_identity": 30, "min_coverage": 50, "max_evalue": 0.01, "query_file":
	"bua.faa", "rid": "RAV3P2F3014"}, "organisms": [{"name": "Aspergillus sp.
	CLMG-2019a", "strain": "FRR 5400", "scaffolds": ...

The session with indent 2 will look like this:

::
	{
	  "queries": [
	    "BuaB",
	    "BuaC",
	    "BuaD",
	    "QBE85644.1",
	    "BuaE",
	    "BuaF",
	    "BuaG",
	    "QBE85648.1",
	    "BuaA"
	  ],
	  "params": {
	    "mode": "remote",
	    "database": "nr",
	    "min_identity": 30,
	    "min_coverage": 50,
	    "max_evalue": 0.01,
	    "query_file": "bua.faa",
	    "rid": "RAV3P2F3014"
	  },
	  "organisms": [
	    {
	    	"name": "Aspergillus sp. CLMG-2019a",
	    	"strain": "FRR 5400",
    		"scaffolds": [{ ... }]
    		...
	
	Much more readable!
	Note though that, particularly in sessions with lots of results, this comes with a significant increase in file size.
	


Installation and Quick Start
============================

Installing Python for Windows
-----------------------------

In order to use ``cblaster`` you will first need to install Python on your computer.
Go to python.org/downloads/ and download the latest version of Python (3.8.5 at the time of writing).
This will initiate the download of the Python installer (python-x.x.x.exe).


Locate the installer in your downloads folder, run the program and follow the wizard.
**Make sure you tick the box `Add Python x.x to PATH'.**
This ensures that ``cblaster`` is available for you to use directly from the terminal.
It is not selected by default and you will have to do this step manually if you do not check the box here.


Once the installer has finished, you can try to run Python in PowerShell to verify that it has been installed correctly.
To open PowerShell in Windows, open a folder, Shift+Right click and select the option `Open PowerShell window here...'.
With PowerShell open, type \texttt{python} and press enter.
If python has installed correctly, the interactive shell should be launched, and the version number should be displayed in the console like so:

::
	
	Python 3.8.5 (default, Jul 20 2020, 17:41:41)
	Type "help", "copyright", "credits" or "license" for more information.
	>>> 
	
The Python package installer tool, ``pip`` is installed alongside Python.
This is necessary to install ``cblaster``, so verify that it is installed by typing ``pip`` in PowerShell as above.

Installing DIAMOND
------------------

``cblaster`` uses ``DIAMOND`` to perform local searches.
This can be freely obtained from \href{http://www.diamondsearch.org/index.php}{diamondsearch.org}.

Installing, uninstalling and updating cblaster
----------------------------------------------

To install ``cblaster``, simply input the command:

::

	pip install cblaster
	
This will install ``cblaster`` as well as all of its dependencies.

Should you decide to uninstall cblaster, this can also be done using pip:

::

	pip uninstall cblaster

Note: If a new version of cblaster is available, you can simply uninstall and reinstall the module as above to access the newer version.

Running your first search
-------------------------
Once ``cblaster`` has been installed, running a search is as simple as providing a collection of query amino acid sequences in a FASTA file, like so:

::

	  cblaster search -qf query.fasta

Visualisations can be generated using the ``-p`` or ``--plot`` argument:

::

	  cblaster search -qf query.fasta -p plot.html

Note: if no file name is provided, the plot will be dynamically served using Python's built in HTTP server.
The plot will be exactly the same, but it will not generate a static HTML file that can be shared around.

Search sessions can be saved for later re-use using the ``-s`` or ``--session`` argument:

::
	  
	  cblaster search -qf query.fasta -p plot.html -s session.json
	
Note: a session is saved as a JavaScript Object Notation (JSON) format file.
This is essentially just a dump of all the code objects, as well as search parameters, used during a ``cblaster`` search.
If you provide a pre-existing session file, ``cblaster`` will attempt to load it **instead** of performing a new search.

That is all you need to know about the basic usage of ``cblaster``.
However, there are many more ways to tweak and run the program to suit your needs which are further explored in the following sections.


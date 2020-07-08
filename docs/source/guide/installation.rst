.. _installation:


Installation
============

This section of the documentation covers the installation of cblaster.


Python version
--------------

`cblaster` is written using Python 3, and should work with any version above 3.3.

Dependencies
------------

These packages are automatically installed when installing cblaster:

- requests_ is used to send/retrieve data to/from the NCBI's public APIs
- numpy_ and scipy_ are used when generating data for cblaster visualisations
- PySimpleGUI_ is used as the graphical toolkit for the GUI
- genome2json_ parses GenBank/GFF genome files and converts them to the
  more workable JSON format, and is used when generating the genome 'database' for local
  cblaster searches

Other dependencies
------------------

- diamond_ is the search tool used in local cblaster searches

Installation
------------

1. (Optional) Create a new virtual environment

.. code-block:: python

        python3 -m virtualenv venv
        source venv/bin/activate

This will create (and activate) a sandboxed environment where you can install
Python packages separately to those available on your system. This isn't necessarily
required, but is recommended.

2. Install cblaster

The easiest way to obtain cblaster is to install it directly from PyPI using `pip`:

.. code-block:: sh

        pip install cblaster

This will install cblaster, as well as all of its required dependencies.
If you are not using a virtual environment, provide the ``--user`` flag to the above
statement. Alternatively, you could clone the cblaster repository from GitHub and
install it like so:

.. code-block:: sh

        git clone https://www.github.com/gamcil/cblaster
        cd cblaster
        pip install .

This will download the latest version of cblaster and install it from the downloaded
folder, rather than from PyPI.

cblaster should now be available directly on your terminal:

::

        $ cblaster -h
        usage: cblaster [-h] [--version] [-d] [-i INDENT] {gui,makedb,search,gne} ...

        cblaster is a tool for finding clusters of homologous proteins. Type -h/--help after either subcommand for full
        description of available arguments.

        positional arguments:
          {gui,makedb,search,gne}
            gui                 Launch cblaster GUI
            makedb              Generate JSON/diamond databases from GenBank files
            search              Start a local/remote cblaster search
            gne                 Perform gene neighbourhood estimation

        optional arguments:
          -h, --help            show this help message and exit
          --version             show program's version number and exit
          -d, --debug           Print debugging information
          -i INDENT, --indent INDENT
                                Total spaces to use as indent in JSON file (def. None)

        Cameron Gilchrist, 2020


.. _requests: https://requests.readthedocs.io/en/master/
.. _numpy: https://numpy.org/
.. _scipy: https://scipy.org/
.. _PySimpleGUI: https://pysimplegui.readthedocs.io/en/latest/
.. _genome2json: https://github.com/gamcil/genome2json
.. _diamond: https://github.com/bbuchfink/diamond

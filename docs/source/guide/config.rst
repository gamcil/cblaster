.. _config_module:

Pre-search configuration using the ``config`` module
====================================================

The NCBI requires that you provide some identification before using their
services in order to prevent abuse. This can be an e-mail address, or more recently,
an API key (https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).

You can use the ``config`` module to set these parameters for ``cblaster`` searches (you'll only have to do this once!).
This module will save a file, ``config.ini``, wherever your operating system stores configuration
files (for example, in Linux it will be saved in ~/.local/config/cblaster).
When you run remote searches in ``cblaster``, it will first check to see if it can find
this file, and then if an e-mail address or API key is saved; if they are not found,
``cblaster`` will throw an error.

To set an e-mail address:

::

        $ cblaster config --email "foo@bar.com"


...or an API key:

::

        $ cblaster config --api_key <your API key>


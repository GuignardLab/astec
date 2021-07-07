
Command line interfaces
=======================
   
Data organization
-----------------

It is assumed that there will be one directory per experiment. This
directory contains the acquired data, but will also contain the result
data as depicted below.
	  
.. code-block:: none

   $ /path/to/experiment/
   ├── RAWDATA/
   │  └── ...
   ├── FUSE
   │  └── ...
   ├── SEG
   │  └── ...
   └── POST
      └── ...




Principle
---------


.. _cli-common-options:

Command line interfaces common options
--------------------------------------

\label{sec:cli:common}


``-h, --help``
   prints a help message
   
``-p file, --parameters file``
   set the parameter file to be parsed
   
``-e path, --embryo-rep path``
   set the
   ``path`` to the directory where the
   ``RAWDATA/`` directory is located.
   Can also be given in the parameter file by the variable ``PATH_EMBRYO``.

``-k, --keep-temporary-files``
   allows to keep the temporary files. Not to be routinely used.

``-f, --force``
   forces execution, even if (temporary) result files
   are already existing

``-v, --verbose``
   increases verboseness (both at console and in the log file)

``-nv, --no-verbose``
   no verboseness

``-d, --debug``
   increases debug information (in the log file)

``-nd, --no-debug``
   no debug information

``-pp, --print-param``
   print parameters in console and exit. A parameter file has to be provided (``-p`` option). Allows to check the parameters that will be used before any processing; it is also a means to have access to the whole parameter list. 

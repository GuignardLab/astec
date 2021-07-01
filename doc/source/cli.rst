-----------------------------------
User guide: command line interfaces
-----------------------------------

Data organization
=================

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



Command Line Interface
======================

.. toctree::
   :maxdepth: 3
   :caption: Contents:

Astec
-----
.. automodule:: astec.bin.astec_astec
   :members:

Embryo properties
-----------------
.. automodule:: astec.bin.astec_embryoproperties
   :members:

Fuse
----
.. automodule:: astec.bin.astec_fuse
   :members:

Intra-registration
------------------
.. automodule:: astec.bin.astec_intraregistration
   :members:

Manual correction
-----------------
.. automodule:: astec.bin.astec_manualcorrection
   :members:

MARS
----
.. automodule:: astec.bin.astec_mars
   :members:

Naming
------
.. automodule:: astec.bin.astec_naming
   :members:

Post correction
---------------
.. automodule:: astec.bin.astec_postcorrection

------------
Installation
------------
   
Only works for Linux or MacOs systems.

User installation (without git)
===================================

Requires conda.

1. Create a conda environment (here the environment is named astec-test)

.. code-block:: bash

   conda create -n astec-test -c morpheme -c conda-forge astec
	
2. Activate the built conda environment
   
.. code-block:: bash

   conda activate astec-test

User installation (with git)
================================

Requires conda and git.

1. Astec code can be found at gitlab.inria.fr/astec/astec. It can be downloaded with

.. code-block:: bash
		    
   git clone https://gitlab.inria.fr/astec/astec.git
	
It creates an astec directory.

2. Create a conda environment named astec

.. code-block:: bash
		    
   cd astec
   conda env create -f pkg/env/astec.yaml	

3. Activate the built conda environment

.. code-block:: bash
		    
   conda activate astec

Developer installation (with git)
======================================

Requires conda and git.

1. Astec code can be found at gitlab.inria.fr/astec/astec. It can be downloaded with

.. code-block:: bash

   git clone https://gitlab.inria.fr/astec/astec.git

It creates an astec directory.

2. Create a conda environment named astec

.. code-block:: bash

   cd astec
   conda env create -f pkg/env/astec-dev.yaml
   
3. Activate the built conda environment

.. code-block:: bash

   conda activate astec-dev
   
4. Install astec package for use

.. code-block:: bash   

    python -m pip install -e .
    
The -e option install the package in "editable" mode, this is want you want if you aim at contributing to the astec project. This last command has to be repeated (within the conda environment every time the astec code has been modified).

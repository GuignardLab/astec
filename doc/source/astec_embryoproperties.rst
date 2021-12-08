
.. _cli-embryoproperties:

``astec_embryoproperties``
==========================

``astec_embryoproperties`` can be used either to extract cell properties as well as cell lineage from a co-registered image sequence or to handle a property file (pkl or xml).



``astec_embryoproperties`` additional options
---------------------------------------------

The following options are available:

``-i files ...``
   input files (pkl or xml) to be read

``-o files ...``
   output files (pkl or xml) to be read

``-c files ...``
   files (pkl or xml) to be compared to those given by ``-i``
   
``-feature features ...``
   features to be extracted from the input files, that are to be
   written in the output files. Features have to be chosen in
   'lineage',  'h_min', 'volume', 'surface', 'sigma', 
   'label_in_time', 'barycenter', 'fate', 'fate2',
   'fate3', 'fate4', 'all-cells', 'principal-value',
   'name', 'contact', 'history', 'principal-vector',
   'name-score', 'cell-compactness'
   
``-property features ...``
   same as ``-feature``
   
``--diagnosis``
   performs some test on the read properties
   
``--diagnosis-minimal-volume DIAGNOSIS_MINIMAL_VOLUME``
   displays all cells with volume smaller than ``DIAGNOSIS_MINIMAL_VOLUME``

``--diagnosis-items DIAGNOSIS_ITEMS``
   minimal number of items to be displayed

``-write-selection, --write-selection``
   convert xml selections into morphonet files

``-fate, --compute-fate``
   delete previous fates ('fate', 'fate2', 'fate3' and 'fate4') and recompute 'fate4'

``--print-content``
   print the keys of the input file(s) (read as python dictionary)

``--print-keys``
   same as ``--print-content``

``--print-types``
   print types of read features (for debug purpose)

   

Extracting properties from a co-registered image sequence
---------------------------------------------------------

When a parameter file is passed after the ``-p`` option, ``astec_embryoproperties`` will compute image sequence properties.
Computing cell related informations as well as the lineage tree
requires that the (post-corrected) segmentation images have already
been co-registered (with ``astec_intraregistration`` see section
:ref:`cli-intraregistration`).  
``astec_embryoproperties`` will parse the
``INTRAREG/INTRAREG_<EXP_INTRAREG>/`` directory, and will compute the
properties from the images in the ``POST/POST_<EXP_POST>/``
sub-directory, if existing, else of from the ``SEG/SEG_<EXP_SEG>/``
sub-directory. 



Embryo properties output data
-----------------------------

The results are stored in the ``POST/POST_<EXP_POST>/`` or
``SEG/SEG_<EXP_SEG>/`` sub-directory under the
``INTRAREG/INTRAREG_<EXP_INTRAREG>`` where
``<EXP_INTRAREG>`` is the value of the variable
``EXP_INTRAREG`` (its default value is ``'RELEASE'``).  
The resulting properties will be stored in the same directory than the images they are issued. It will be stored as a pickle python file, and also as a XML file. Both files contain exactly the same information.

According that the ``POST/POST_<EXP_POST>/`` sub-directory exists (that post-corrected segmentation images have been co-registered), 3 files will be created, named after ``<EN>``


.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── INTRAREG/
   │  └── INTRAREG_<EXP_INTRAREG>/
   │     ├── ...
   │     ├── POST/
   │     │   └── POST_<EXP_POST>/
   │     │       ├── <EN>_intrareg_post_lineage.pkl
   │     │       ├── <EN>_intrareg_post_lineage.txt
   │     │       └── <EN>_intrareg_post_lineage.xml
   │     ...
   ...

The computed information are

``all_cells``
   All the cell identifiers. Each cell (in a segmentation image) has a
   given label (ranging from 2 and above, 1 being used for the
   background) in each image. To uniquely identify a cell in the
   sequence, it has been given an unique identifier computed by $i *
   1000 + c$, $i$ and $c$ denoting respectively the image index
   (ranging in [``<begin>``, ``<end>``]) and the cell label.
   
``cell_barycenter``
   Cell center of mass (in voxel coordinates)
   
``cell_contact_surface``
   For each cell, give for each neighboring cell the contact
   surface. The sum of these contact surfaces is the cell surface.
   
``cell_principal_vectors``
   The cell principal vectors are issued from the diagonalization of
   the cell covariance matrix (in voxel unit).
   
``cell_principal_values``
   The cell principal value are issued from the diagonalization of the
   cell covariance matrix (in voxel unit).
   
``cell_volume``
   Cell volume (in voxel unit)

``cell_compactness``
   The cell compactness is defined by :math:`\mathcal{C}
   =\frac{\sqrt[3]{\mathcal{V}}}{\sqrt[2]{\mathcal{S}}}` where
   :math:`\mathcal{V}` is the volume of the cell and
   :math:`\mathcal{S}` is its surface. 

``cell_surface``
   Cell surface (in pixel unit). For this computation, is mandatory
   that the co-registered images are isotropic (the same voxel size
   along the 3 dimensions X, Y, and Z). 

``cell_lineage``

The text file ``<EN>_intrareg_post_lineage.txt`` contains diagnosis information about the sequence. It lists

* the cell with the smallest sizes as well as the ones with the
  largest sizes 
* the cell with a weird lineage: cells without a mother cell, or cells
  without daughter cells or having more than 2 daughter cells 
* cells having a small intersection with its mother cell with respect
  to either the mother cell volume or the cell volume.  

Note that a property file may contain some other information that can be computed either by ``astec_embryoproperties`` 
(e.g. with the ``--compute-fate`` option) or by other means.




Handling property files
-----------------------

``astec_embryoproperties`` can also help managing property files.

* Converting from ``xml`` to ``pkl`` and  the other way around.
  
  .. code-block:: bash

     $ astec_embryoproperties -i file.pkl -o file.xml

  convert the pickle file ``file.pkl`` into the ``xml`` file  ``file.xml``

* Converting the lineage information from either an ``xml``
  or an ``pkl`` file to a ``tlp`` [#]_ file for lineage visualization
  
  .. code-block:: bash

     $ astec_embryoproperties -i file.pkl -o file.tlp

  convert the pickle file ``file.pkl`` into the ``tlp`` file  ``file.tlp``

* Merging files.

  .. code-block:: bash

     $ astec_embryoproperties -i file1.pkl file2.xml ... filen.pkl -o merge.xml merge.pkl

  will merge the files  ``file1.pkl``,  ``file2.xml`` , ...,
  ``filen.pkl`` (note that they can be either xml or pkl) and write
  the result both in ``xml`` and ``pkl`` formats.
  
* Extracting properties.

  .. code-block:: bash
		  
     $ astec_embryoproperties -i file.pkl -feature volume surface -o file.xml

  will extract the cell volume and surface information from the
  pickle file ``file.pkl`` and write them into the xml file
  ``file.xml``. 

* Comparing property files may help to view changes and/or correction between two property files

  .. code-block:: bash

     $ astec_embryoproperties -i file.pkl -c file_to_be_compared_to.pkl

  compare the two files ``file.pkl`` and ``file_to_be_compared_to.pkl``. 
  The comparison is made on all common properties (according it has been implemented). 
  The ``-feature`` option allows to select the features to be compared.

* Assessing a property file

  .. code-block:: bash

     $ astec_embryoproperties -i file.pkl --diagnosis

  will run some test/diagnosis on some properties (only a few features are tested). It may help 
  at detecting errors either in the segmented images or in the property file.
  The ``-feature`` option allows to select the features to be tested.


.. _cli-embryoproperties-diagnosis:

Diagnosis on property file
--------------------------

Apart reporting diagnosis in the console and in the log file, 
``selection`` properties (in the ``morphonet`` sense) may be added in the output property file 
(if specified thanks to the ``-o`` option) that can be also written as ``selection`` files 
(if the option ``-write-selection`` is used).

Contact surfaces
^^^^^^^^^^^^^^^^

This diagnosis checks whether 
there are branches with large contact surface distance between consecutive cells
(large means larger than ``maximal_contact_distance``, see section :ref:`cli-parameters-diagnosis`).
A ``morphonet`` selection (of ``float`` type) is created whose values are the distance value
for both consecutive points that participate in the distance calculation.

Lineage
^^^^^^^
This diagnosis checks whether the lineage is well-formed. 

A ``morphonet`` selection (of ``selection`` type) is created whose values are

* ``10`` for the first cell of lineage trees starting after the first time point,
* ``20`` for cells with multiple mother cells,
* ``30`` for the last cell of branches ending before the last time point,
* ``40`` for dividing cells with more than 2 daughter cells,
* ``50`` for the first cell of short non-terminal branches 
  (short means less than the value of ``minimal_length``, see section :ref:`cli-parameters-diagnosis`).

Name
^^^^
This diagnosis checks whether the cell names obey to Conklin's syntax :cite:p:`conklin:1905aa`.

A ``morphonet`` selection (of ``selection`` type) is created whose values are

* ``10`` for cells that are not named but have a non-dividing predecessor that is named
* ``20`` for cells that are named but have an unamed predecessor
* ``30`` for cells that are named but differently than their non-dividing predecessor 
* ``40`` for named cells that come after a dividing cell, with a name that does not follow
  Conklin's syntax
* ``50`` for cells that come after a dividing cell, with a name that follows
  Conklin's syntax, but with a sibling without name
* ``60`` for cells that come after a dividing cell, with a name that follows
  Conklin's syntax, but with a sibling that has the same name
* ``70`` for cells that come after a dividing cell, with a name that follows
  Conklin's syntax, but with a sibling name that does not follow
  Conklin's syntax
* ``80`` for other errors


Volume 
^^^^^^

This diagnosis checks whether 

* there are small cell volume 
  (small means less than ``minimal_volume``, see section :ref:`cli-parameters-diagnosis`).

* there are branches with a large volume variation
  (large means larger than ``maximal_volume_variation``, see section :ref:`cli-parameters-diagnosis`).
  Branch volume variation is computed by
  :math:`100 * \frac{\max_{t} v(c_t) - \min_{t} v(c_t)}{\mathrm{med}_t v(c_t)}` where :math:`v(c_t)` is the volume
  of the cell :math:`c_t` and :math:`\mathrm{med}` is the median value.

* there are branches with large volume derivatives
  (large means larger than ``maximal_volume_derivative``, see section :ref:`cli-parameters-diagnosis`).
  The volume derivative along a branch is calculated as 
  :math:`100 * \frac{v(c_{t+1}) - v(c_{t})}{v(c_{t})}` where :math:`t` denotes the successive acquisition time points.
  A ``morphonet`` selection (of ``float`` type) is created whose values are the absolute value of derivative
  for both points (:math:`v(c_{t})` and :math:`v(c_{t+1})`) that participate in the derivative calculation.





.. [#] Tulip is a Data Visualization Software, see `tulip.labri.fr <http://tulip.labri.fr/>`_

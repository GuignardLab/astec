.. _cli-manual-correction:

``astec_manualcorrection``
===========================



Manual correction overview
--------------------------

``astec_manualcorrection`` aims at correcting segmentation errors, either for a single image, as the one produced by 
``astec_mars`` (see :ref:`cli-mars`), or a series of images, as the ones produced by ``astec_astec`` (see :ref:`cli-astec-astec`).

The seeded watershed is likely to produce segmentation errors, even with a careful choice of parameters. 
It is advised to set the parameters to favour over-segmentations insted of under-segmentations since the 
former are much more easier to correct, which is the purpose of ``astec_manualcorrection``. 

When dealing with a single image produced by ``astec_mars``,
the segmentation error correction could also be done at the seed correction step (see :ref:`sec-cli-mars-seed-correction`) 
where seeds can be fused, or new seeds can be added (by specifying point coordinate):
see the ``seed_edition_file`` variable for details in section :ref:`cli-parameters-seed-edition`.

Here, a under-segmentated cell can be splitted by changing the parameters controling the search for seeds
by :math:`h`-minima computation for the indicated cells, ie
``watershed_seed_hmin_min_value``, ``watershed_seed_hmin_max_value``, and ``watershed_seed_hmin_delta_value``.
The :math:`h`-minima computation start with a :math:`h` of ``watershed_seed_hmin_max_value`` until two seeds are found 
inside the cell of interest. Then, a watershed is performed to get the final splitting of the cell.

See section :ref:`cli-parameters-astec-manualcorrection` for a detailed presentation of the mapping file,
to specify cells to be fused or cells to be splitted.

For time series, newly formed divisions are propagated along the time series and embryo properties are updated (only cell lineage 
and cell volume are updated, other properties, that may become wrong, are deleted). 
The division propagation is done as in step 1 of ``astec_astec`` (see :ref:`cli-astec-step-1`)

1. Segmentation of previous time point is deformed onto the current time point.

2. Cells are eroded to get two seeds inside the cell to be divided.

3. Segmentation (watershed-based) from the deformed seeds.

This propagation only changes the progeny of the cell(s) to be splitted. Other cells remain unchanged. 

Since pre-processed images (see :ref:`cli-astec-pre-processing`) are used either for seed computation (as in :ref:`cli-astec-step-2`) or for watershed segmentation, it is adviced to keep the pro-processed images by setting ``keep_reconstruction`` to ``True``.



Manual correction / output data
-------------------------------

Single segmentation file correction 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The historical used of ``astec_manualcorrection`` was to correct the segmentation image issued from
``astec_mars`` to provide a good initialization (a good first segmentation image)
to ``astec_astec``. Hence, the file issued from ``astec_mars``, namely ``<EN>_mars_t<begin>.inr``,
was transformed into a file with a different name, namely``<EN>_seg_t<begin>.inr``.

The result can then be stored in the same sub-directory
``SEG/SEG_<EXP_SEG>`` under the
``/path/to/experiment/`` directory where ``<EXP_SEG>`` is the value of the variable ``EXP_SEG`` (its
default value is '``RELEASE``').


.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── SEG/
   │  └── SEG_<EXP_SEG>/
   │     ├── <EN>_mars_t<begin>.inr
   │     ├── <EN>_seg_t<begin>.inr
   │     ├── LOGS/
   │     └── RECONSTRUCTION/
   ...

However, it is also possible to change the name of the result directory, as detailed below.  
   

Time series segmentation files correction 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A time series of segmentation images can also be corrected. Here the segmentation files
from the sub-directory ``SEG/SEG_<EXP_SEG_FROM>`` are corrected into the ``SEG/SEG_<EXP_SEG_TO>``.
The property file is updated Only the cell lineage and volumes are updated, other properties are deleted 
and should be recomputed with ``astec_embryoproperties`` once the corrections are done.

.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── SEG/
   │  ├── SEG_<EXP_SEG_FROM>/
   │  │  ├── <EN>_seg_lineage.xml
   │  │  ├── <EN>_seg_t<begin>.inr
   │  │  ├── ...
   │  │  ├── <EN>_seg_t<end>.inr
   │  │  ├── LOGS/
   │  │  └── RECONSTRUCTION/
   │  └─ SEG_<EXP_SEG_TO>/
   │     ├── <EN>_seg_lineage.xml
   │     ├── <EN>_seg_t<begin>.inr
   │     ├── ...
   │     ├── <EN>_seg_t<end>.inr
   │     ├── LOGS/
   │     └── RECONSTRUCTION/
   ...


Segmentation correction parameters
----------------------------------

``astec_manualcorrection`` parses a correction file whose name is given by the variable ``manualcorrection_file``.
The syntax of this file is described in section :ref:`cli-parameters-astec-manualcorrection`.
See also the 
`tutorial section <https://astec.gitlabpages.inria.fr/astec-tutorial/astec_tutorial.html#correction-of-the-first-time-point-segmentation>`_
for an other example.



``astec-postcorrection``
========================



Post-correction overview
------------------------

The Astec segmentation procedure yields a series of segmented
images :math:`\{S^{\star}_t\}_t`, where each segmented
image :math:`S^{\star}_{t+1}` takes advantage of the knowledge of the
previously segmented image :math:`S^{\star}_t` to decide,
at the cell level, whether a cell division may
occur. However, there are still segmentation errors, that
can be detected from the study of the overall lineage
(see :cite:p:`guignard:tel-01278725` (section 2.3.3.7, page 74)
and :cite:p:`guignard:hal-02903409` (supp. mat.)).

As suggested by its name, the post-correction will try to a posteriori correct the segmentation resulting from the ``astec_astec`` stage (see section :ref:`cli-astec-astec`). 

The post-correction is made of the following steps.

1. Lineage pruning: it goes through the end branches (a branch does
   not have any cell division; an end branch finishes either at the
   end of the sequence or the cell vanishes between two time points)
   of the lineage tree. Some lineage end branches are deleted (see
   section :ref:`cli-postcorrection-pruning` for details), meaning
   that the corresponding cells are fused with other cells of the
   embryo. 

2. Division postponing: some divisions are postponed.



Post-correction / input data
----------------------------

Input data are the result of the ``astec_astec`` stage (see section :ref:`cli-astec-astec`) and will be searched in the directory ``SEG/SEG_<EXP_SEG>/`` (see section :ref:`cli-astec-output-data`).

.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── SEG/
   │  └── SEG_<EXP_SEG>/
   │     ├── <EN>_seg_lineage.xml.
   │     ├── <EN>_seg_t<begin>.inr
   │     ├── <EN>_seg_t<:math:`\ldots`>.inr.
   │     ├── <EN>_seg_t<end>.mha.
   │     ├── LOGS/
   │     └── ...
   ...



Post-correction / output data
-----------------------------

The results are stored in sub-directories ``POST/POST_<EXP_POST>`` under the
``/path/to/experiment/`` directory where where ``<EXP_POST>`` is the
value of the variable ``EXP_POST`` (its default value is ``'RELEASE'``). 

.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── POST/
   │  └── POST_<EXP_POST>/
   │     ├── <EN>_post_lineage.xml.
   │     ├── <EN>_post_t<begin>.inr
   │     ├── <EN>_post_t<:math:`\ldots`>.inr.
   │     ├── <EN>_post_t<end>.mha.
   │     ├── LOGS/
   │     └── ...
   ...


The image format to be used (here ``mha``) is given by the variable
``result_image_suffix``, while the lineage format to be used (here
``xml``) is given by the variable ``result_lineage_suffix``. 




.. _cli-postcorrection-pruning:

Step 1: lineage pruning
-----------------------

Bifurcations of the lineage tree correspond to cell division, while branches (between two bifurcations or between a bifurcation and a leaf) corresponds to the lifespan of a cell.
The purpose of this step is to detect suspicious end branches (terminating by a leaf) that may correspond to an over-segmentation error. 

An end branch is candidate for deletion if

* either it terminates before the last time point (it corresponds then to a cell without daughter cell in the next time point), 
* or the volume of its last cell is too small (threshold given by the variable ``postcorrection_volume_minimal_value``).


An end branch candidate for deletion is deleted if

* either it is too short (threshold given by the variable
  ``postcorrection_lifespan_minimal_value``), 
* or (if the variable ``postcorrection_test_early_division`` is set to
  ``True``) either its sister branch (which may not be an end branch)
  or its mother branch is too short, meaning that there are two
  divisions too close, (thresholds still given by the variable
  ``postcorrection_lifespan_minimal_value``), 
* or if the Pearson correlation coefficient between the volumes of the
  candidate end branch and its sister branch is less than
  -``postcorrection_correlation_threshold``, meaning that the volumes 
  are anti-correlated (typically the volumes of the candidate end branch
  are decreasing while the ones of the sister branch are increasing,
  indicating a fake division detection).  



Step 2: division postponing
---------------------------

      
* ``postcorrection_volume_minimal_value``
  branch ending with leaf cell below this value are candidate for deletion. Expressed in voxel unit.
* ``postcorrection_lifespan_minimal_value``
* ``postcorrection_test_early_division``
* ``postcorrection_test_volume_correlation``
* ``postcorrection_correlation_threshold``
* ``postcorrection_lineage_diagnosis``
  performs a kind of diagnosis on the lineage before and after the post-correction.








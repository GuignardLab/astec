``astec_manual_correction``
===========================



Manual correction overview
--------------------------

The seeded watershed is likely to produce segmentation errors, even with a careful choice of parameters. It is advised to set the parameters to favour over-segmentations insted of under-segmentations since the former are much more easier to correct, which is the purpose of ``astec_manualcorrection``. 
Note that the segmentation error correction could also be done at the
seed correction step of the ``astec_mars`` stage, see 
:doc:`astec_mars`.


``astec_manual_correction`` additional options
----------------------------------------------

The following options are available:

``-i input_image``
   set the ``input_image`` file to be corrected. Allows to skip the
   automated naming of files.  

``-o output_image``
   set the resulting ``ouput_image`` file to be saved. Allows to skip
   the automated naming of files.
   
``-m mapping_file``
   set the ``mapping_file`` to be used for the correction.

``-nsc smallest_cells``
   set the number of the smallest cells to be displayed after
   correction. The smallest cells are the most likely to be issued
   from an over-segmentation. 

``-nlc largest_cells``
   set the number of the largest cells to be displayed after
   correction. The largest cells are the most likely to be issued from
   an under-segmentation.   



Manual correction / output data
-------------------------------

The results are stored in sub-directories
``SEG/SEG\_<EXP\_SEG>`` under the
``/path/to/experiment/`` directory where ``<EXP\_SEG>`` is the value of the variable ``EXP\_SEG`` (its
default value is '``RELEASE``').
``<EN>\_seg\_t<begin>.inr`` is the correction of the segmentation image ``<EN>\_mars\_t<begin>.inr``.

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
   

   
Segmentation correction parameters
----------------------------------

``astec_manual_correction`` parses a correction file whose name is given by the variable ``mancor_mapping_file``. The syntax of this file is very simple. Lines beginning with ``#`` are ignored (and can be used to insert comments in the files). Non-empty lines should contain two numbers separated by a space, and ``astec_manual_correction`` will replace the first number by the second in the segmentation file.

E.g. a cell ``c`` is recognized to be over-segmented, and then is
represented by two labels, says 9 and 10. Thus the line

.. code-block:: none

   10 9		

will replace all 10's by 9's in the segmentation image,  thus ``c`` will only be represented by 9's after correction. See also the 
`tutorial section <https://astec.gitlabpages.inria.fr/astec-tutorial/astec_tutorial.html#correction-of-the-first-time-point-segmentation>`_
for an other example.



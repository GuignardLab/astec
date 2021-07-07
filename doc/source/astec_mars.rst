
.. _cli-mars:

``astec_mars``
==============



Mars method overview
--------------------

The name ``mars`` comes from :cite:p:`fernandez:hal-00521491` where ``MARS`` is the acronym of *multiangle image acquisition, 3D reconstruction and cell segmentation*.

This method aims at producing a segmentation of a membrane cell image (e.g.  a fused image) into a segmentation image. This segmentation image is a integer-valued image where each integer labeled an unique cell in the image. By convention, '1' is the background label, while cells have labels greater than 2. It is made of the following steps:


1. Pre-processing of the input image to produce the input seed image for seed computation.
This is described in section \ref{sec:cli:input:image:preprocessing}. The parameters that governed the pre-processing are described in section \ref{sec:cli:parameters:preprocessing} and prefixed by ``seed_``.

2. Seed extraction through the computation of the `h`-minima of the input seed image

3. Eventually seed correction

4. Pre-processing of the input image to produce the input membrane image for the seeded watershed.
This is described in section \ref{sec:cli:input:image:preprocessing}. The parameters that governed the pre-processing are described in section \ref{sec:cli:parameters:preprocessing} and prefixed by ``membrane_``.

5. A seeded watershed.



Mars / output data
------------------

The results are stored in sub-directories
``SEG/SEG_<EXP_SEG>`` under the
``/path/to/experiment/`` directory where where ``<EXP_SEG>`` is the
value of the variable ``EXP_SEG`` (its default value is ``'RELEASE'``). 

.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── SEG/
   │  └── SEG_<EXP_SEG>/
   │     ├── <EN>_mars_t<begin>.inr
   │     ├── LOGS/
   │     └── RECONSTRUCTION/
   ...


.. _cli-mars-input-data:

Steps 1 and 4: input image pre-processing
-----------------------------------------

The input image (typically the fused image representing the cell
membranes/walls) can be pre-processed before use in the seeded
watershed.The pre-processing can be different for the seed input image
(the one that will be used to extract the seeds) and the membrane
input image (the one that will be used as the height image for the
seeded watershed). Details about the pre-processing can be found in
section \ref{sec:cli:input:image:preprocessing}. 

Default settings are

.. code-block:: python
		
  intensity_transformation = 'Identity'
  intensity_enhancement = None

meaning that the original fused image is used for both
inputs. Different pre-processing can be done. E.g.

.. code-block:: python
		
  seed_intensity_transformation = 'Identity'
  membrane_intensity_transformation = 'normalization_to_u8'
  intensity_enhancement = None

comes to use the original image for the seed extraction, but its
normalization into 8 bits (1 byte) as the height image for the seeded
watershed. 

If the input image is transformed before segmented, the transformed
images can be saved in the directory
``SEG/SEG_<EXP_SEG>/RECONSTRUCTION/`` if the value of the variable
``keep_reconstruction`` is set to ``True``. 



Step 2: seed extraction
-----------------------

The seed extraction is made of the following steps:

1. Gaussian smoothing of the input image, the gaussian standard deviation being given by the variable ``seed_sigma``.

2. Extraction of the `h`-minima of the previous image, `h`  being given by the variable ``seed_hmin``.

3. Hysteresis thresholding (and labeling)  of the `h`-minima image, with
   a high threshold equal to ``seed_high_threshold`` (default is `h`)
   and and a low threshold equal to 1. It then only selects the
   `h`-minima that have an actual depth of `h`. 



.. _sec-cli-mars-seed-correction:

Step 3: seed correction
-----------------------

Several rounds of correction of the computed seeds can be done. At each round, different seeds can be assigned the same label (and this will fuse the further reconstructed cells) or new seeds (each new seed is a single voxel) can be added. See the \option{seed_edition_files} variable for details.

When correcting seeds, it is advised to launch ``astec_mars``  with the ``-k`` option. Indeed, temporary files, as the seed image, are kept in a temporary directory located in the ``SEG/SEG_'EXP_SEG'/`` directory and then re-used, and not recomputed at each ``astec_mars`` use.




Step 5 : seeded watershed
-------------------------

Given the seeds, the watershed is performed on the smoothed input
membrane image (gaussian standard deviation being given by the
variable ``membrane_sigma``).



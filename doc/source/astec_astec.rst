.. role:: python(code)
   :language: python



.. _cli-astec-astec:

``astec_astec``
===============

The name ``astec`` comes from the Phd work of L. Guignard :cite:p:`guignard:tel-01278725` where ``ASTEC`` is the acronym of *adaptive segmentation and tracking of embryonic cells*.

This method aims at producing a segmentation of each membrane cell image  (e.g. a fused image) temporal sequence. This is a method of segmentation by propagation: it used the segmentation of the previous timepoint (say :math:`t-1`) to constraint the segmentation at the aimed timepoint (say :math:`t`).



Astec method overview
---------------------

Astec principle is to guide the segmentation of the :math:`t` timepoint image :math:`I_{t}` with the segmentation :math:`S^\star_{t-1}` of the :math:`t-1` timepoint image :math:`I_{t-1}`.

1.  A first segmentation of :math:`I_{t}`, :math:`\tilde{S}_{t}`, is computed by a seeded watershed, where the seeds are the eroded cells of :math:`S^\star_{t-1}`, projected onto :math:`I_{t}`. By construction, no cell division can occur.
2. :math:`h`-minima are computed over a range of :math:`h` values. Studying the numbers of :math:`h`-minima located in each cell of :math:`\tilde{S}_{t}` gives an indication whether there might be a cell division or not. From this study a seed image :math:`Seeds_{t}` is computed, and then a new segmentation image :math:`\hat{S}_{t}`.
3. Some potential errors are detected by checking whether there is a significant volume decrease from a cell of :math:`S^\star_{t-1}` and its corresponding cells in :math:`\hat{S}_{t}`. For such cells, seeds may be recomputed, as well as the :math:`\hat{S}_{t}` segmentation.
4. It may occur, in previous step, that some cell from :math:`S^\star_{t-1}` correspond to 3 cells in :math:`\hat{S}_{t}`. This step aims at correcting this.
5. Some other potential errors are detected by checking whether there is a significant volume decrease from a cell of :math:`S^\star_t` and its corresponding cells in :math:`\hat{S}_{t}` due to a background invasion. For such cells, morphosnakes :cite:p:`marquez-neil:pami:2014` are computed to try to recover cell loss.
6. The morphosnake correction operated in previous step may invade the background too much. This step aims at correcting this.






.. _cli-astec-output-data:

Astec / output data
-------------------

The results are stored in sub-directories
``SEG/SEG_<EXP_SEG>`` under the
``/path/to/experiment/`` directory where where ``<EXP_SEG>`` is the value of the variable ``EXP_SEG`` (its
default value is '``RELEASE``'). 

.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── SEG/
   │  └── SEG_<EXP_SEG>/
   │     ├── <EN>_mars_t<begin>.inr
   │     ├── <EN>_seg_lineage.xml
   │     ├── <EN>_seg_t<begin>.inr
   │     ├── <EN>_seg_t<:math:`\ldots`>.inr
   │     ├── <EN>_seg_t<end>.inr
   │     ├── LOGS/
   │     └── RECONSTRUCTION/
   ...



Input image pre-processing
--------------------------

The input image (typically the fused image representing the cell membranes/walls) can be pre-processed before use in the astec stage (as for ``astec_mars``, see section :ref:`cli-mars-input-data`)
The pre-processing can be different for the 

* the seed input image (the one that will be used to compute the :math:`h`-minima),
* the membrane input image (the one that will be used as the height image for the seeded watersheds), and
* the morphosnake input image (the one that will be used to define the morphosnake energy).


Pre-processing parameters, described in section \ref{sec:cli:parameters:preprocessing}, and prefixed respectively by ``seed_``, ``membrane_`` and ``morphosnake_`` allow to tune these pre-processing.
Hence, the lines
\begin{verbatim}
seed_intensity_transformation = 'Identity'
membrane_intensity_transformation = 'normalization_to_u8'
morphosnake_intensity_transformation = 'Identity'
intensity_enhancement = None
\end{verbatim}
come to choose the original image for both the seed extraction and the morphosnake stage, but its normalization on 8 bits for the seeded watershed (this corresponds to the choice of the historical version of astec).


Step 1: :math:`\tilde{S}_{t}`
-----------------------------

A first segmentation of :math:`I_{t}`, :math:`\tilde{S}_{t}`, is computed by a seeded watershed, where the seeds are built from the eroded cells of :math:`S^\star_{t-1}`.

* ``previous_seg_method`` = ``'erode_then_deform'``
  
   The cells of :math:`S^\star_{t-1}` are first eroded, yielding the image :math:`S^e_{t-1}`, then this image is mapped onto :math:`I_{t}` frame thanks to the transformation :math:`\mathcal{T}_{t-1 \leftarrow t}`, resulting in the eroded seed image :math:`S^e_{t-1 \leftarrow t} = S^e_{t-1} \circ \mathcal{T}_{t-1 \leftarrow t}`. This is the historical astec behavior.

* ``previous_seg_method`` = ``'deform_then_erode'`` 

   :math:`S^\star_{t-1}` is first mapped onto :math:`I_{t}` frame thanks to the transformation :math:`\mathcal{T}_{t-1 \leftarrow t}`, resulting in the image :math:`S^\star_{t-1 \leftarrow t} = S^\star_{t-1} \circ \mathcal{T}_{t-1 \leftarrow t}`. Cells of :math:`S^\star_{t-1 \leftarrow t}` are then eroded to get :math:`S^e_{t-1 \leftarrow t}`


This seed image, :math:`S^e_{t-1 \leftarrow t}`, plus the membrane input image are used as input for a seeded watershed, and yield :math:`\tilde{S}_{t}`. 
By construction, no cell division can occur in :math:`\tilde{S}_{t}` with respect to :math:`S^\star_{t-1}`.

If the variable ``propagation_strategy`` is set to 
``'seeds_from_previous_segmentation'``, 
the segmentation propagation stops and :math:`\tilde{S}_{t}` is the final result. 



Step 2: :math:`\hat{S}_{t}`
---------------------------

The :math:`h`-minima are computed in the seed input image for a range of :math:`h \in [h_{min}, h_{max}]`, with a step of :math:`\delta h`.

:math:`h_{min}`, :math:`h_{max}` and :math:`\delta h` are set respectively by the variables
``watershed_seed_hmin_min_value``,
``watershed_seed_hmin_max_value``, and
``watershed_seed_hmin_delta_value``.

For a given cell of :math:`\tilde{S}_{t}`, if there is no cell
division, and if the :math:`h`-minima are well detected, ther should
be only one :math:`h`-minima included in the cell for all values of :math:`h`.
However, if a cell division occurs, there should be mostly
two :math:`h`-minima included in the cell.  Then, the study of the
number of :math:`h`-minima strictly included allows to decide whether
a cell division has occur (see :cite:p:`guignard:tel-01278725`,
:cite:p:`guignard:hal-02903409` for details). 

This step results in the image :math:`\hat{S}_{t}`.

If the variable ``propagation_strategy`` is set to 
``'seeds_selection_without_correction'``, the segmentation propagation
stops and :math:`\hat{S}_{t}` is the final result.  



Steps 3 and 4: volume checking
------------------------------

Some potential errors are detected by checking whether there is a
large volume decrease from a cell of :math:`S^\star_{t-1}` and its
corresponding cells in :math:`\hat{S}_{t}`. For such cells, seeds are
recomputed, as well as the :math:`\hat{S}_{t}` segmentation. 

It may occur, in this step, that some cell from :math:`S^\star_{t-1}`
correspond, after correction, to 3 cells in :math:`\hat{S}_{t}`. A
second step aims at correcting this. 



Steps 5 and 6: morphosnake correction
-------------------------------------

This step is performed if ``morphosnake_correction`` is set to :python:`True`.

Some other potential errors are detected by checking whether there is
a significant volume decrease from a cell of :math:`S^\star_t` and its
corresponding cells in :math:`\hat{S}_{t}` due to a background
invasion. For such cells, morphosnakes :cite:p:`marquez-neil:pami:2014`
are computed to try to recover cell loss. 





.. _cli-intraregistration:

``astec_intraregistration``
==============================



Intra-registration procedure overview
-------------------------------------

The sequence intra-registration procedure can be done either after the
fusion step, or after the (post-)segmentation step. It aims at

* compensating for the eventual motion of the imaged sample with respect to the microscope
* resampling the fusion and/or the segmentation images into a common frame/geometry, so they can better be compared, and
* building 2D+t images made of 2D sections from either the  fusion
  and/or the segmentation images, so that the quality of the fusion
  and/of the tracking step can be visually assessed. 

The intra-registration procedure is made of the following steps:

   1. Co-registration of pairs of successive fused images (section
   :ref:`cli-intraregistration-coregistration`). This yields the
   transformations :math:`T_{t+1 \leftarrow t}`. Fused images
   are located in ``<EMBRYO>/FUSE/FUSE_<EXP_FUSE>``: the parameter
   ``EXP_FUSE``is either set in the parameter file or is set at
   ``RELEASE``. This step may be long. 

   2. Composition of transformations issued from the co-registration
   step. This step computes the transformations :math:`T_{ref \leftarrow
   t}`. towards a reference image ``ref`` given by the parameter
   ``intra_registration_reference_index``.

   3. Computation of the *template* image (section
   :ref:`cli-intraregistration-template`). This *template*
   image dimension are computed so that the useful information of all
   resampled images fits into it. Useful information can be issued from
   either the fused sequence, the segmentation sequence or the
   post-segmentation sequence. It is indicated by the
   ``intra_registration_template_type`` which value can be either
   ``'FUSION'``, ``'SEGMENTATION'``, or
   ``'POST-SEGMENTATION'``. This step may be long.

   4. Resampling of either the fused or the segmentation images
   (section :ref:`cli-intraregistration-resampling`). Note that
   changing the parameters for this step will not require to re-compute
   the first steps.

   5. Extraction of 2D+t images from the resampled sequences (section
   :ref:`cli-intraregistration-movies`). Note that changing the
   parameters for this step (i.e. requiring extra movies) will not
   require to re-compute the first steps, with an eventual exception
   for the resampling step.

   6. Computation of a maximum image from the resampled images (section
   :ref:`cli-intraregistration-maximum`). Computing the maximum over  
   the resampled fusion images may be useful to define a common cropping
   area for the sequence.
   Note that changing the parameters for this step will not require 
   to re-compute the first steps. 


  
``astec_intraregistration`` additional options
----------------------------------------------

The following options are available:

``-t file``
   set the resampling transformation file for the reference image (see
   section :ref:`cli-intraregistration-template`)
   
``-a string``
   set the resampling transformation angles for the reference image
   (see section :ref:`cli-intraregistration-template`) 



Intra-registration / input data
-------------------------------

Intra-registration / multichannel acquisition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The co-registration transformations are computed from one series of
fused images, issued from the
``/path/to/experiment/FUSE/FUSE_<EXP_FUSE>`` directory. 

In case of multi-channel acquisition, all fused image directories
listed in the ``EXP_FUSE`` variable will be transformed by the
transformations computed on the *first* fused image directory of the
list. 

As detailed in section :ref:`cli-intraregistration-coregistration`,
specifying

.. code-block:: python

   EXP_FUSE = ['MEMBRANES', 'NUCLEI']

in the parameter file implies that co-registrations will be computed on
the fused images from ``FUSE/FUSE_MEMBRANES/``, but both fused image
series will be transformed. 

The same stands for segmentation and post-segmentation series:
multiple directories can be specified in either ``EXP_SEG`` or
``EXP_POST``. 



Intra-registration / output data
--------------------------------

The results are stored in sub-directories
``INTRAREG/INTRAREG_<EXP_INTRAREG>`` under the
``/path/to/experiment/`` directory where ``<EXP_INTRAREG>`` is the 
value of the variable ``EXP_INTRAREG`` (its default value is ``'RELEASE'``). 

.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── INTRAREG/
   │  └── INTRAREG_<EXP_INTRAREG>/
   │     ├── CO-TRSFS/
   │     ├── [FUSE/]
   │     ├── LOGS/
   │     ├── [MAXIMUM/]
   │     ├── [MOVIES/]
   │     ├── [POST/]
   │     ├── [SEG/]
   │     └── TRSFS_t<begin>-<end>/
   ...

Output data are of two kinds: image series (fused images, segmentation
images, post-corrected segmentation images) can be resampled in the
same common geometry (also known as the *template*), see
section :ref:`cli-intraregistration-resampling`, and 3D (ie 2D+t)
images of the evolution (with respect to time) of one section (XY, XZ,
or YZ) of the images of the series can be built, see
section :ref:`cli-intraregistration-movies`.



.. _cli-intraregistration-coregistration:

Step 1: co-registration
-----------------------

Default registration parameters for the co-registration are set by:

.. code-block:: python

   # intra_registration_compute_registration = True
   # intra_registration_transformation_type = 'rigid'
   # intra_registration_transformation_estimation_type = 'wlts'
   # intra_registration_lts_fraction = 0.55
   # intra_registration_pyramid_highest_level = 6
   # intra_registration_pyramid_lowest_level = 3
   # intra_registration_normalization = True

Computed transformations are stored in
``INTRAREG/INTRAREG_<EXP_INTRAREG>/CO-TRSFS``.It may be advised to set
the pyramid lowest level value to some higher value to speed up the
co-registrations (recall that all pairs of successive images will be
co-registered, i.e.

.. code-block:: python

   intra_registration_pyramid_lowest_level = 4

Co-registration are computed using the fused images of
``/path/to/experiment/FUSE/FUSE_<EXP_FUSE>``. If
``EXP_FUSE`` is a list of strings (ie indicates a list a
directories) rather than a single string, the fused image from the
first directory are used for the co-registration computation.

Typically, if there are several fused series (eg, in case of
multi-channel acquisition) as in


.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── FUSE/
   │  ├── FUSE_MEMBRANES/
   │  │   └── ...
   │  └── FUSE_NUCLEI/
   │      └── ...
   ...

Specifying

.. code-block:: python
		
   EXP_FUSE = ['MEMBRANES', 'NUCLEI']

in the parameter file implies that co-registrations will be done on
the fused images from ``FUSE/FUSE_MEMBRANES/``.



.. _cli-intraregistration-template:

Step 3: template building
-------------------------

.. code-block:: python
		
   # intra_registration_reference_index = None
   # intra_registration_reference_resampling_transformation_file = None
   # intra_registration_reference_resampling_transformation_angles = None
   #
   # intra_registration_template_type = "FUSION"
   # intra_registration_template_threshold = None
   # intra_registration_margin = None
   #
   # intra_registration_resolution = 0.6
   #
   # intra_registration_rebuild_template = False


* The ``intra_registration_reference_index`` allows to choose the reference image (the one which remains still, i.e. up to a translation), by default it is the first image image of the series (associated to ``begin``). 
  However, it may happen that this image has to be reoriented to fit the user's expectation. The resampling transformation\footnote{The resampling transformation is the one that goes from the destination image towards the input image.}, that re-orient the reference image, can then be given and will be applied to the whole series.

   * ``intra_registration_reference_resampling_transformation_file`` can be given a resampling transformation file name.
   * ``intra_registration_reference_resampling_transformation_angles`` can be given a string describing the successive rotations (with respect to the frame axis) to be applied. E.g. the string ``"X 30 Y 50"`` defines a resampling transformation equal to :math:`R_X(30) \circ R_Y(50)` where :math:`R_X(30)` is a rotation of 30 degrees around the X axis and :math:`R_Y(50)` is a rotation of 50 degrees around the Y axis.


* Depending on ``intra_registration_template_type`` (``'FUSION'``,
  ``'SEGMENTATION'`` or ``'POST-SEGMENTATION'``), the two latter
  assume obviously that the segmentation has been done), the
  *template* image can be built either after the fusion or the
  segmentation images. If no threshold is given by
  ``intra_registration_template_threshold``, the built template will
  be large enough to include all the transformed fields of view (in this
  case, the template is the same whatever
  ``intra_registration_template_type`` is).

  If ``intra_registration_template_type='FUSION'`` (respectively
  ``'SEGMENTATION'`` and ``'POST-SEGMENTATION'``),  the template
  is built from the images of the first directory indicated by
  ``EXP_FUSE`` (respectively
  ``EXP_SEG`` and ``EXP_POST``) in case of
  ``EXP_FUSE`` contains a list of strings.

  If a threshold is given, the built template will be large enough to
  include all the transformed points above the threshold. E.g., the
  background is labeled with either '1' or '0' in segmentation images,
  then a threshold of '2' ensures that all the embryo cells will not be
  cut by the resampling stage.  In this case, adding an additional
  margin (with ``intra_registration_margin``) to the template could be a good idea for visualization
  purpose. 

* Specifying  using a different resolution for the drift-compensated
  series than the ``target_resolution`` (the resolution of the fused
  images) allows to decrease the resampled images volume. This can be
  achieved by setting ``intra_registration_resolution`` to the desired  value  (default is 0.6).

* Last, co-registrations may have been computed during a first
  computation, fused images being used to compute the template. However,
  if a subsequent segmentation has been conducted, a smaller template
  is likely to be computed (with the segmentation images to build the
  template), without recomputing the co-registration. This is the
  purpose of the variable
  ``intra_registration_rebuild_template``.
  If set to ``True``, it forces to recompute the template as well
  as the transformations from the co-registrations (that are not
  re-computed). Obviously, resampling as well as 2D+t movies are also
  re-generated.


As an example, building a *template* image after the segmentation images can be done with

.. code-block:: python
		
   # intra_registration_reference_index = None
  intra_registration_template_type = "SEGMENTATION"
  intra_registration_template_threshold = 2
  # intra_registration_resolution = 0.6
  intra_registration_margin = 10

Computed transformations from the *template* image as well as the
*template* image itself are stored in
``INTRAREG/INTRAREG<EXP_INTRAREG>/TRSFS_t<F>-<L>/`` where ``<F>`` and
``L`` are the first and the last index of the series (specified by
``begin`` and ``end`` from the parameter file). 



.. _cli-intraregistration-resampling:

Step 4: resampling fusion/segmentation images
---------------------------------------------

The resampling of the fused and/or segmentation images are done
depending on the value of the following variables (here commented). Resampling is done
either if the following parameters are set to ``True`` or if movies
are requested to be computed (section :ref:`cli-intraregistration-movies`).

.. code-block:: python
		
   # intra_registration_resample_fusion_images = True
   # intra_registration_resample_segmentation_images = False
   # intra_registration_resample_post_segmentation_images = False

This default behavior implies that the fusion images will be resampled
while the segmentation and the post-corrected segmentation images are not.


Resampled images will be stored in the
``INTRAREG/INTRAREG_<EXP_INTRAREG/>`` directory, with the same
hierarchy than under ``/path/to/experiment``. E.g. 

.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── FUSE/
   │  ├── FUSE_1/
   │  │   └── ...
   │  └── FUSE_2/
   │      └── ...
   ...

Specifying

.. code-block:: python
		
   EXP_FUSE = ['1', '2']

in the parameter file causes the resampling of both fused image series
(``FUSE/FUSE_1/`` and ``FUSE/FUSE_2/``)


.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── FUSE/
   │  ├── FUSE_1/
   │  │   └── ...
   │  └── FUSE_2/
   │      └── ...
   ├── INTRAREG/
   │  └── INTRAREG_<EXP_INTRAREG>/
   │     ├── CO-TRSFS/
   │     ├── FUSE/
   │     │   ├── FUSE_1/
   │     │   │   └── ...
   │     │   └── FUSE_2/
   │     │       └── ...
   │     ...
   ...

The same behavior stands for ``EXP_SEG`` and  ``EXP_POST``.



.. _cli-intraregistration-movies:

Step 5: 2D+t movies
-------------------

For either visual assessment or illustration purposes, 2D+t (i.e. 3D) images can be built from 2D sections extracted from the resampled temporal series. This is controlled by the following parameters:

.. code-block:: python
		
   # intra_registration_movie_fusion_images = True
   # intra_registration_movie_segmentation_images = False
   # intra_registration_movie_post_segmentation_images = False
   
   # intra_registration_xy_movie_fusion_images = [];
   # intra_registration_xz_movie_fusion_images = [];
   # intra_registration_yz_movie_fusion_images = [];

   # intra_registration_xy_movie_segmentation_images = [];
   # intra_registration_xz_movie_segmentation_images = [];
   # intra_registration_yz_movie_segmentation_images = [];
   
   # intra_registration_xy_movie_post_segmentation_images = [];
   # intra_registration_xz_movie_post_segmentation_images = [];
   # intra_registration_yz_movie_post_segmentation_images = [];

If ``intra_registration_movie_fusion_images`` is set to ``True``, a
movie is made with the  XY-section located at the middle of each
resampled fusion image (recall that, after resampling, all images have
the same geometry). Additional XY-movies can be done by specifying the
wanted Z values in ``intra_registration_xy_movie_fusion_images``. E.g.

.. code-block:: python

   intra_registration_xy_movie_fusion_images = [100, 200];

will build two movies with XY-sections located respectively at Z values of 100 and 200. The same stands for the other orientation and for the resampled segmentation images.

Movies will be stored in the
``INTRAREG/INTRAREG_<EXP_INTRAREG>/MOVIES/`` directory, with the same
hierarchy than under ``/path/to/experiment``. E.g., 

.. code-block:: python

   EXP_FUSE = ['1', '2']

in the parameter file results in

.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── FUSE/
   │  ├── FUSE_1/
   │  │   └── ...
   │  └── FUSE_2/
   │      └── ...
   ├── INTRAREG/
   │  └── INTRAREG_<EXP_INTRAREG>/
   │     ├── CO-TRSFS/
   │     ├── FUSE/
   │     │   ├── FUSE_1/
   │     │   │   └── ...
   │     │   └── FUSE_2/
   │     │       └── ...
   │     ├── MOVIES/
   │     │   └── FUSE/
   │     │       ├── FUSE_1/
   │     │       │   └── ...
   │     │       └── FUSE_2/
   │     │           └── ...
   │     ...
   ...

The same behavior stands for ``EXP_SEG`` and  ``EXP_POST``.



.. _cli-intraregistration-maximum:

Step 6: 3D maximum over the 3D+t sequence
-----------------------------------------

To set a cropping area valid for the whole resampled sequence, a maximum image can be built from the resampled temporal series. This is controlled by the following parameters:

.. code-block:: python
		
   # intra_registration_maximum_fusion_images = False
   # intra_registration_maximum_segmentation_images = False
   # intra_registration_maximum_post_segmentation_images = False

If ``intra_registration_maximum_fusion_images`` is set to ``True``, a maximum image is computed over the sequence of resampled fusion images (recall that, after resampling, all images have the same geometry). The value of a voxel in this maximum image is the maximum value (over time) of this voxel in the sequence.

The maximum image will be stored in the
``INTRAREG/INTRAREG_<EXP_INTRAREG>/MAXIMUM/`` directory, with the same
hierarchy than under ``/path/to/experiment``. E.g., 

.. code-block:: python
		
   EXP_FUSE = ['1', '2']

in the parameter file results in

.. code-block:: none

   /path/to/experiment/
   ├── ...
   ├── FUSE/
   │  ├── FUSE_1/
   │  │   └── ...
   │  └── FUSE_2/
   │      └── ...
   ├── INTRAREG/
   │  └── INTRAREG_<EXP_INTRAREG>/
   │     ├── CO-TRSFS/
   │     ├── FUSE/
   │     │   ├── FUSE_1/
   │     │   │   └── ...
   │     │   └── FUSE_2/
   │     │       └── ...
   │     ├── MAXIMUM/
   │     │   └── FUSE/
   │     │       ├── FUSE_1/
   │     │       │   └── ...
   │     │       └── FUSE_2/
   │     │           └── ...
   │     ...
   ...

The same behavior stands for ``EXP_SEG`` and  ``EXP_POST``.

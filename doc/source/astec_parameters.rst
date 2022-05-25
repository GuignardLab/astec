Parameters
==========

The different command line interfaces, or CLIs, (``astec_fuse``, ``astec_mars``, etc.) requires a parameter file (which is nothing but a ``python`` file) that contains both information on the experiment (path to the experiment directory, on the sub-directory names -- see section :ref:`cli-parameters-data`) as well as specific parameters for the CLIs.



Prefixed parameters
-------------------

Some of the parameter sets are said to be *prefixed*, such as the two sets of pre-processing parameters for the ``astec_mars`` CLI (see section :ref:`cli-parameters-mars`). Indeed, the pre-processing can be set differently for the seed input image and the membrane input image (eg see section :ref:`cli-mars`).

Prefixing parameters allows to either set *all* the parameters with the same name together or set them *independently*.

As exemplified in section :ref:`cli-mars-input-data`, 
the parameter file lines (where the variables are not prefixed)

.. code-block:: python
		
   intensity_transformation = 'normalization_to_u8'
   intensity_enhancement = None

will set the corresponding pre-processing parameters for both the seed and the membrane image pre-processing. However, using prefixes, as in the lines

.. code-block:: python
		
   seed_intensity_transformation = 'Identity'
   membrane_intensity_transformation = 'normalization_to_u8'
   intensity_enhancement = None

allows to set them independently.

This mechanism is designed to simplify the parameter file, but may have undesired consequences. Indeed, using the basic variable names of the registration parameters (see section :ref:`cli-parameters-registration`)  for the ``astec_astec`` CLI will change all registration parameters included in the pre-processing parameters.

To check whether the parameters have been set correctly, one can either use
the ``--print-param`` CLI option (see section :ref:`cli-common-options`) beforehand, and/or to a posteriori check the used parameters in the ``log`` file. 
  
  

Common parameters
-----------------

``begin``
  first time point to be processed (``astec_fuse``, 
  ``astec_astec`` or ``astec_postcorrection``) 
  or single time point to be processed
  (``astec_mars`` or ``astec_manualcorrection``).

``end``
  last time point to be processed (``astec_fuse``, 
  ``astec_astec`` or ``astec_postcorrection``).

``delta``
  interval between two time points to be processed. Set to 1 by default.
  Fragile.
  
``raw_delay``
  Delay to be added to the time points to build the file names. 
  Fragile.
  
``time_digits_for_filename``
  number of digits used to build the file names.
  
``time_digits_for_cell_id``
  number of digits used to define unique cellule id. in lineage file.
  The unique id of cell :math:`c` at time :math:`t` is
  :math:`t \times 10^d + c` where :math:`d` is set by
  ``time_digits_for_cell_id``.  

``default_image_suffix``:
  used for both the result and the temporary data.

    * ``'inr'``: Inrimage format, kept for historical reasons.
    * ``'mha'``: MetaImage format, readable by Fiji.
    * ``'tif'``: not advised, since the tiff format does not allow
      to keep the voxel size along the z direction (aka spacing),
      at least in a standardized way. 
    * ``'nii'``: Nifti format, compatible with Omero.

  Gzipped image files (with the additional extension ``'.gz'`` are also readable.
  
``result_image_suffix``:
  used for both the result data.

``result_lineage_suffix``:

    6. ``'pkl'``: Pickle file
    * ``'xml'``: Xml file


.. _cli-parameters-data:

Data organisation parameters
----------------------------

``DIR_LEFTCAM_STACKONE`` 
  see section :ref:`cli-fuse-input-data`, 
  see figures :numref:`fig-data-rawdata-1`, 
  :numref:`fig-data-rawdata-2`, 
  and :numref:`fig-data-rawdata-3`.

``DIR_LEFTCAM_STACKONE_CHANNEL_2`` 
  see section :ref:`cli-fuse-input-data`

``DIR_LEFTCAM_STACKONE_CHANNEL_3`` 
  see section :ref:`cli-fuse-input-data`

``DIR_LEFTCAM_STACKZERO`` 
  see section :ref:`cli-fuse-input-data`,
  see figures :numref:`fig-data-rawdata-1`, :numref:`fig-data-rawdata-2`, 
  and :numref:`fig-data-rawdata-3`.

``DIR_LEFTCAM_STACKZERO_CHANNEL_2`` 
  see section :ref:`cli-fuse-input-data`

``DIR_LEFTCAM_STACKZERO_CHANNEL_3`` 
  see section :ref:`cli-fuse-input-data`

``DIR_RAWDATA`` 
  see section :ref:`cli-fuse-input-data`,
  see figures :numref:`fig-data-rawdata-1`, :numref:`fig-data-rawdata-2`, 
  and :numref:`fig-data-rawdata-3`.

``DIR_RAWDATA_CHANNEL_2`` 
  see section :ref:`cli-fuse-input-data`

``DIR_RAWDATA_CHANNEL_3`` 
  see section :ref:`cli-fuse-input-data`

``DIR_RIGHTCAM_STACKONE`` 
  see section :ref:`cli-fuse-input-data`,
  see figures :numref:`fig-data-rawdata-1`, :numref:`fig-data-rawdata-2`, 
  and :numref:`fig-data-rawdata-3`.

``DIR_RIGHTCAM_STACKONE_CHANNEL_2`` 
  see section :ref:`cli-fuse-input-data`

``DIR_RIGHTCAM_STACKONE_CHANNEL_3`` 
  see section :ref:`cli-fuse-input-data`

``DIR_RIGHTCAM_STACKZERO`` 
  see section :ref:`cli-fuse-input-data`,
  see figures :numref:`fig-data-rawdata-1`, :numref:`fig-data-rawdata-2`, 
  and :numref:`fig-data-rawdata-3`.

``DIR_RIGHTCAM_STACKZERO_CHANNEL_2`` 
  see section :ref:`cli-fuse-input-data`

``DIR_RIGHTCAM_STACKZERO_CHANNEL_3`` 
  see section :ref:`cli-fuse-input-data`

``EN``: 
  the so-called *embryo* name. 
  All files will be named after this name.
  E.g. see section :ref:`cli-fuse-output-data`
  and figure :numref:`fig-data-fusion`.

``EXP_FUSE``:
  String (``str`` type) or list (``list`` type) of
  strings. It indicates what are the fused images directories, 
  of the form ``<PATH_EMBRYO>/FUSE/FUSE_<EXP_FUSE>``.

  .. code-block:: python

    EXP_FUSE = 'exp1'
    EXP_FUSE = ['exp1', 'exp2']

  are then both valid. 
  Default value of ``EXP_FUSE`` is ``'RELEASE'``.
  See section :ref:`cli-fuse-output-data`,
  see figure :numref:`fig-data-fusion`.

``EXP_FUSE_CHANNEL_2`` 
  see section :ref:`cli-fuse-output-data`

``EXP_FUSE_CHANNEL_3`` 
  see section :ref:`cli-fuse-output-data`

``PATH_EMBRYO``: 
  path to the *experiment*.
  If not present, the current directory is used.
  See section :ref:`cli-fuse-input-data`,
  see figures :numref:`fig-data-rawdata-1`, :numref:`fig-data-rawdata-2`, 
  :numref:`fig-data-rawdata-3`, and :numref:`fig-data-fusion`

``acquisition_leftcam_image_prefix``  
  see section :ref:`cli-fuse-input-data`,
  see figures :numref:`fig-data-rawdata-1`, :numref:`fig-data-rawdata-2`, 
  and :numref:`fig-data-rawdata-3`.

``acquisition_rightcam_image_prefix``  
  see section :ref:`cli-fuse-input-data`,
  see figures :numref:`fig-data-rawdata-1`, :numref:`fig-data-rawdata-2`, 
  and :numref:`fig-data-rawdata-3`.



.. code-block:: none
  :caption: Typical organisation of mono-channel data.
  :name: fig-data-rawdata-1

   ``<PATH_EMBRYO>``/
   ├── ``<DIR_RAWDATA>``/
   │  ├── ``<DIR_LEFTCAM_STACKZERO>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  ├── ``<acquisition_leftcam_image_prefix>001.zip``
   │  │  └── ...
   │  ├── ``<DIR_RIGHTCAM_STACKZERO>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  ├── ``<acquisition_leftcam_image_prefix>001.zip``
   │  │  └── ...
   │  ├── ``<DIR_LEFTCAM_STACKONE>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  ├── ``<acquisition_leftcam_image_prefix>001.zip``
   │  │  └── ...
   │  └── ``<DIR_RIGHTCAM_STACKONE>``/
   │     ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │     ├── ``<acquisition_leftcam_image_prefix>001.zip``
   │     └── ...
   ...

.. code-block:: none
  :caption: Typical organisation of multi-channel data.
  :name: fig-data-rawdata-2

   ``<PATH_EMBRYO>``/
   ├── ``<DIR_RAWDATA>``/
   │  ├── ``<DIR_LEFTCAM_STACKZERO>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_RIGHTCAM_STACKZERO>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_LEFTCAM_STACKONE>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_RIGHTCAM_STACKONE>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_LEFTCAM_STACKZERO_CHANNEL_2>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_RIGHTCAM_STACKZERO_CHANNEL_2>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_LEFTCAM_STACKONE_CHANNEL_2>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  └── ``<DIR_RIGHTCAM_STACKONE_CHANNEL_2>``/
   │     ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │     └── ...
   ...

.. code-block:: none
  :caption: Alternative organisation of multi-channel data.mono-channel data.
  :name: fig-data-rawdata-3

   ``<PATH_EMBRYO>``/
   ├── ``<DIR_RAWDATA>``/
   │  ├── ``<DIR_LEFTCAM_STACKZERO>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_RIGHTCAM_STACKZERO>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_LEFTCAM_STACKONE>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  └── ``<DIR_RIGHTCAM_STACKONE>``/
   │     ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │     └── ...
   ├── ``<DIR_RAWDATA_CHANNEL_2>``/
   │  ├── ``<DIR_LEFTCAM_STACKZERO>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_RIGHTCAM_STACKZERO>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  ├── ``<DIR_LEFTCAM_STACKONE>``/
   │  │  ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │  │  └── ...
   │  └── ``<DIR_RIGHTCAM_STACKONE>``/
   │     ├── ``<acquisition_leftcam_image_prefix>000.zip``
   │     └── ...
   ...

.. code-block:: none
  :caption: Typical organisation of fused images.
  :name: fig-data-fusion

   ``<PATH_EMBRYO>``/
   ├── ``<DIR_RAWDATA>``/
   │  └── ...
   ├── ``<FUSE>``/
   │  └── ``FUSE_<EXP_FUSE>``/
   │     ├── ``<EN>_fuse_t000.<result_image_suffix>``
   │     ├── ``<EN>_fuse_t001.<result_image_suffix>``
   │     └── ...
   ...
   




.. _cli-parameters-ace:

Ace parameters
--------------

Ace stand for *Automated Cell Extractor*. [G[L]]ACE methods aim at detecting and enhancing membranes in a 3D images (see also section :ref:`cli-input-image-preprocessing-membrane`). 

1. Hessian-based detection of 2-D manifolds, computation of a center-membrane image.  
2. Thresholding of the center-membrane image to get a binary image.
3. Reconstruction of a membrane images from the binary image through tensor voting.
   

``sigma_membrane``
  this is the gaussian sigma that is used to compute image derivatives
  (in real units), for the Hessian-based detection of 2-D manifolds.
  
``hard_thresholding``
  ``True`` or ``False``.
  If set to ``True``, a hard threshold (set by variable 
  ``hard_threshold``) is used instead of an automated threshold.
  
``hard_threshold``

``manual``
  ``True`` or ``False``.
  By default, this parameter is set to False. If failure, 
  (meaning that thresholds are very bad, meaning that the binarized 
  image is very bad), set this parameter to True and relaunch the
  computation on the test image. If the method fails again, "play" 
  with the value of ``manual_sigma`` ... and good luck.
  
``manual_sigma``
  Axial histograms fitting initialization parameter for the computation 
  of membrane image binarization axial thresholds (this parameter is 
  used if ``manual`` is set to ``True``). 
  One may need to test different values of 
  ``manual_sigma``. We suggest to test values between 5 and 25 
  in case of initial failure. Good luck.

``sensitivity``
  Membrane binarization parameter.
  Use larger values (smaller than or equal to 1.0) to increase 
  the quantity of binarized membranes to be used for tensor voting.

``sigma_TV``
  Parameter which defines the voting scale for membrane structures 
  propagation by tensor voting method (real coordinates). This parameter
  should be set between :math:`3 \mu m` (little cells) and 
  :math:`4.5 \mu m` (big gaps in the binarized membrane image).
  
``sigma_LF``:
  Additional smoothing parameter for reconstructed image 
  (in real coordinates).
  It seems that the default value = :math:`0.9 \mu m` is 
  ok for standard use.
  
``sample``:
  Set the fraction (in [0, 1]) of the binarized membranes further
  used for tensor voting.  
  It allows tensor voting computation speed optimisation (do not 
  touch if not bewared): the more sample, the higher the cost.

``sample_random_seed``
  Drawing a sample from the binarized membranes (see parameter 
  ``sample``) is a stochastic process. Setting this parameter 
  to some ``int`` value allows to make this stochastic process 
  reproducible.
  
``bounding_box_dilation``

``default_image_suffix``



.. _cli-parameters-morphosnake:

Morphosnake parameters
----------------------

``dilation_iterations``
  dilation of the cell bounding box for computation purpose.

``iterations``
  maximal number of morphosnake iterations.

``delta_voxel``: 
  error on voxel count to define a stopping criteria.

``energy``
  * ``'gradient'``: uses the same formula as in 
    :cite:p:`marquez-neil:pami:2014`, as in the historical 
    astec version. But seems to be a poor choice.
  * ``'image'``: uses directly the image as the energy image.
    
``smoothing``:
  internal parameter for the morphosnake.

``balloon``:
  internal parameter for the morphosnake.

``processors``: number of processors used for the 
  morphosnake correction.

``mimic_historical_astec``:
  ``True`` or ``False``. 
  If set to ``True``, same implementation than the historical 
  astec version. Kept for comparison purpose.



.. _cli-parameters-preprocessing:
  
Preprocessing parameters
------------------------

The input image may be pre-processed before being used as

* either the *membrane* image (ie the height image) for watershed segmentation,
* or the *seed* image (ie the image with which the regional minima are
  computed), 
* or the *morphosnake* image (ie the image with which the morphosnake energy is
  computed).
  
For more details, see section :ref:`cli-input-image-preprocessing`.

* Ace parameters:
  see section :ref:`cli-parameters-ace`.


* ``intensity_prenormalization``: possible values are
  * ``'identity'``
  * ``'normalization_to_u8'``
  * ``'normalization_to_u16'``
  Performs a global robust normalization of the input image, prior to other
  pre-processing. The
  intensity value corresponding to the min percentile is set
  to 0, while the intensity value corresponding to the max
  percentile is set either to 255 (u8) or 4000 (u16). In-between
  values are linearly interpolated.
  Should be left to 'identity' for integer-encoded images.
  It has been introduced for real-encoded images.

  It is governed by the variables:
  * ``prenormalization_max_percentile``
  * ``prenormalization_min_percentile``

* ``intensity_transformation``:
  set the (histogram based) intensity transformation of the original image
  (see section :ref:`cli-input-image-preprocessing-histogram`)

  * ``None``: no intensity transformation of the original image is used
    to pre-process the input image. 
  * ``'identity'``: the input image is used without any transformation.
  * ``'normalization_to_u8'``: the input image (usually encoded on 16 bits)
    is normalized onto 8 bits. The values corresponding to percentiles
    given by the  variables ``normalization_min_percentile`` and
    ``normalization_max_percentile`` are mapped respectively on 0 and 255. 
  * ``'cell_normalization_to_u8'``: same principle than
    ``'normalization_to_u8'`` but values mapped on 0 and 255 are computed
    on a cell basis (cells are the ones of :math:`S^{\star}_{t-1} \circ
    \mathcal{T}_{t-1 \leftarrow t}` -- see cite:p:`guignard:tel-01278725` for
    notations --, ie the segmentation obtained for the previous time point
    :math:`t-1` and deformed onto the frame at the current time point
    :math:`t`). This can be used only with ``astec_astec`` (section
    :ref:`cli-astec-astec`).
	  
    This feature has been added for tests, but has not demonstrated yet
    any benefit.
    
* ``intensity_enhancement``:
  set the membrane enhancement transformation of the original image
  (see section :ref:`cli-input-image-preprocessing-membrane`)

  * ``None``: no membrane enhancement of the original image is used to
    pre-process the input image.
  * ``'GACE'``: stands for *Global Automated Cell Extractor*. It tries to
    reconstructed a membrane image through a membrane detector, an automated
    thresholding and a tensor voting step. The automated thresholding is
    computed once for the whole image.
  * ``'GLACE'``: stands for *Grouped Local Automated Cell Extractor*. It differs
    from one step from ``GACE``: the threshold of extrema image is not computed
    globally (as in ``GACE``), but one threshold is computed per cell of
    :math:`S^{\star}_{t-1} \circ \mathcal{T}_{t-1 \leftarrow t}`, from the
    extrema values of the cell bounding box.
    This can be used only with ``astec_astec`` (section
    :ref:`cli-astec-astec`).

* ``outer_contour_enhancement``:
  This feature has been added for tests, but has not demonstrated yet
  any benefit.

* ``reconstruction_images_combination``:
  * ``'addition'``
  * ``'maximum'``

* ``cell_normalization_min_method``:
  set the cell area where is computed the percentile value that 
  will give the :math:`0` value in the normalized image
  * ``'cell'``
  * ``'cellborder'``
  * ``'cellinterior'``

* ``cell_normalization_max_method``:
  set the cell area where is computed the percentile value that 
  will give the :math:`255` value in the normalized image
  * ``'cell'``
  * ``'cellborder'``
  * ``'cellinterior'``

* ``normalization_min_percentile``
* ``normalization_max_percentile``
* ``cell_normalization_sigma``:
  the ``'cell_normalization_to_u8'`` method computes a couple
  :math:`(I_{min}, I_{max})` for each cell of
  :math:`S^{\star}_{t-1} \circ \mathcal{T}_{t-1 \leftarrow t}`,
  yielding discontinuities in the :math:`I_{min}` and :math:`I_{max}`
  from cell to cell. To normalize the whole image, images of :math:`I_{min}`
  and :math:`I_{max}` are built and then smoothed with a gaussian kernel
  (sigma given by the variable ``cell_normalization_sigma``.

* ``intensity_transformation``

  * ``'identity'``: intensities of the input image are left unchanged
  * ``'normalization_to_u8'``: intensities of the input image are
    globally normalized into [0,255] 
  * ``'normalization_to_u16'``: intensities of the input image are
    globally normalized into [0,2047] 
  * ``'cell_normalization_to_u8'``: intensities of the input image are
    normalized into [0,255] on a cell-based manner, those cells 
    being the cells at :math:`t-1` deformed on :math:`t`. Works only in the 
    propagation segmentation stage. Fragile

* Registration parameters 
  (see section :ref:`cli-parameters-registration`) prefixed 
  by ``linear_registration_``
* Registration parameters 
  (see section :ref:`cli-parameters-registration`) prefixed 
  by ``nonlinear_registration_``
* ``keep_reconstruction``:
  ``True`` or ``False``. If set to ``True``, 
  pre-processed images are kept in a ``RECONSTRUCTION/`` directory.



.. _cli-parameters-registration:

Registration parameters
-----------------------

* ``compute_registration``
* ``pyramid_highest_level``:
  highest level of the pyramid image for registration.
  Registration is done hierarchically with a pyramid of images. At 
  each pyramid level, image dimensions are divided by 2.
  Setting this variable to 6 means that registration starts with images 
  whose dimensions are 1/64th of the original image.
* ``pyramid_lowest_level``:
  lowest level of the pyramid image for registration. Setting it
  to 0 means that the lowest level is with the image itself.
  Setting it to 1 or even 2 allows to gain computational time.
* ``gaussian_pyramid``
* ``transformation_type``
* ``elastic_sigma``
* ``transformation_estimation_type``
* ``lts_fraction``
* ``fluid_sigma``
* ``normalization``


.. _cli-parameters-seed-edition:

Seed edition parameters
-----------------------

* ``seed_edition_dir``:
* ``seed_edition_file``:
  if run with ``'-k'``, temporary files, including the computed 
  seeds are kept into a temporary directory, and can be corrected in
  several rounds
  
  .. code-block:: python

    seed_edition_file = [['seeds_to_be_fused_001.txt', 'seeds_to_be_created_001.txt'], \
                     ['seeds_to_be_fused_002.txt', 'seeds_to_be_created_002.txt'],
                     ...
                     ['seeds_to_be_fused_00X.txt', 'seeds_to_be_created_00X.txt']]

  Each line of a ``seeds_to_be_fused_00x.txt`` file contains 
  the labels to be fused, e.g. "10 4 2 24". A same label can be found 
  in several lines, meaning that all the labels of these lines will be 
  fused. Each line of ``seeds_to_be_created_00x.txt`` contains 
  the coordinates of a seed to be added.


  
.. _cli-parameters-watershed:

Watershed parameters
--------------------

* ``seed_sigma``:
  gaussian sigma for smoothing of initial image for seed extraction
  (real coordinates).
* ``seed_hmin``:
  :math:`h` value for the extraction of the :math:`h`-minima,
* ``seed_high_threshold``:
  regional minima thresholding. 
* ``membrane_sigma``:
  gaussian sigma for smoothing of reconstructed image for image 
  regularization prior to segmentation
  (real coordinates).


.. _cli-parameters-fuse:

``astec_fuse`` parameters
-------------------------

* ``acquisition_orientation``:
  image orientation (``'right'`` or ``'left'``)
  gives the rotation (with respect to the Y axis) of the left camera 
  frame of stack #0 to be aligned with the the left camera 
  frame of stack #1.
  
  * ``'right'``: +90 degrees
  * ``'left'``: -90 degrees
    
  See section :ref:`cli-fuse-important-parameters`.
* ``acquisition_mirrors``:
  mirroring of the right camera image along the X-axis.
  Right camera images may have to be mirrored along the X-axis 
  to be aligned with the left camera images.
  
  * ``True``: +90 degrees
  * ``False``: -90 degrees
    
  Since it should depend on the apparatus,
  this parameter should not change for all acquisitions 
  performed by the same microscope.
  See section :ref:`cli-fuse-important-parameters`.
* ``acquisition_resolution``:
  acquisition voxel size
  e.g.
  
  .. code-block:: python
		  
    raw_resolution = (.21, .21, 1.)

  see section :ref:`cli-fuse-important-parameters`.
* ``acquisition_stack0_leftcamera_z_stacking``:
  see ``acquisition_leftcamera_z_stacking``.
* ``acquisition_stack1_leftcamera_z_stacking``:
  see ``acquisition_leftcamera_z_stacking``.
* ``acquisition_slit_line_correction``:
  ``True`` or ``False``.
  See section :ref:`cli-fuse-overview`.
* ``target_resolution``:
  isotropic voxel size of the fusion result (fused images).
  See section :ref:`cli-fuse-output-data`.
  
* ``fusion_strategy``:

  * ``'direct-fusion'``: each acquisition is linearly 
    co-registered with the first acquisition (stack #0, left camera). 
    Used registration parameters are the ones  prefixed by 
    ``fusion_preregistration_`` and 
    ``fusion_registration_``.
    Then weights and images are transformed thanks to the 
    computed transformations.
  * ``'hierarchical-fusion'``: from the couple 
    (left camera, right camera), each stack is reconstructed (with the
    registration parameters prefixed by 
    ``fusion_preregistration_`` and 
    ``fusion_registration_``), following the same scheme than 
    the direct fusion but with only 2 images. 
    Then stack #1 is (non-)linearly co-registered with stack #0 with the
    registration parameters prefixed by 
    ``fusion_stack_preregistration_`` and 
    ``fusion_stack_registration_``.
    Images and weights associated with stack\#1 are then (non-)linearly 
    transformed. Finally a weighted linear combination gives the result.

  See section :ref:`cli-fuse-image-registration`.

* ``acquisition_cropping``:
  ``True`` or ``False``. If set to ``True``, 
  the acquisitions stacks are cropped before fusion along the X and Y directions.
  See section :ref:`cli-fuse-raw-data-cropping`.
* ``acquisition_z_cropping``:
  ``True`` or ``False``. If set to ``True``, 
  the acquisitions stacks are cropped before fusion along the Z direction.

* ``acquisition_cropping_margin_x_0``:
  extra margin for the left side of the X direction.
* ``acquisition_cropping_margin_x_1``:
  extra margin for the right side of the X direction.
* ``acquisition_cropping_margin_y_0``:
  extra margin for the left side of the Y direction.
* ``acquisition_cropping_margin_y_1``:
  extra margin for the right side of the Y direction.
* ``acquisition_cropping_margin_z_0``:
  extra margin for the left side of the Z direction.
* ``acquisition_cropping_margin_z_1``:
  extra margin for the right side of the Z direction.
* ``acquisition_cropping_margin_x``: 
  allows to set both ``acquisition_cropping_margin_x_0`` and
  ``acquisition_cropping_margin_x_1``
* ``acquisition_cropping_margin_y``: 
  allows to set both ``acquisition_cropping_margin_y_0`` and
  ``acquisition_cropping_margin_y_1``
* ``acquisition_cropping_margin_z``: 
  allows to set both ``acquisition_cropping_margin_z_0`` and
  ``acquisition_cropping_margin_z_1``
* ``acquisition_cropping_margin``: 
  allows to set the six margin variables.

* ``raw_crop`` same as ``acquisition_cropping``

* Registration parameters
  (see section :ref:`cli-fuse-image-registration`) prefixed 
  by ``fusion_preregistration_``
* Registration parameters 
  (see section :ref:`cli-fuse-image-registration`) prefixed 
  by ``fusion_registration_``
* Registration parameters 
  (see section :ref:`cli-fuse-image-registration`) prefixed 
  by ``fusion_stack_preregistration_``
* Registration parameters 
  (see section :ref:`cli-fuse-image-registration`) prefixed 
  by ``fusion_stack_registration_``
  
* ``xzsection_extraction``:
  ``True`` or ``False``.
  Setting ``xzsection_extraction`` to ``True`` allows 
  to extract XZ-sections of the 4 co-registered stacks as well as 
  the weighting function images. It provides an efficient way to check
  whether the ``acquisition_leftcamera_z_stacking`` variable 
  was correcly set.
  See section :ref:`cli-fuse-nonlinear-registration`
  
* ``fusion_cropping``:
  ``True`` or ``False``. If set to ``True``, 
  the fusion result is cropped along X and Y directions.
  see section :ref:`cli-fuse-fused-data-cropping`
* ``fusion_z_cropping``:
  ``True`` or ``False``. If set to ``True``, 
  the fusion result is cropped along the Z direction.

* ``fusion_cropping_margin_x_0``
* ``fusion_cropping_margin_x_1``
* ``fusion_cropping_margin_y_0``
* ``fusion_cropping_margin_y_1``
* ``fusion_cropping_margin_z_0``
* ``fusion_cropping_margin_z_1``
* ``fusion_cropping_margin_x``:
  allows to set both ``fusion_cropping_margin_x_0``
  and ``fusion_cropping_margin_x_1``
* ``fusion_cropping_margin_y``:
  allows to set both ``fusion_cropping_margin_y_0``
  and ``fusion_cropping_margin_y_1``
* ``fusion_cropping_margin_z``:
  allows to set both ``fusion_cropping_margin_z_0``
  and ``fusion_cropping_margin_z_1``
* ``fusion_cropping_margin``:
  allows to set the six margin variables.
  
* ``acquisition_leftcamera_z_stacking``:
  allows to set both ``acquisition_stack0_leftcamera_z_stacking`` 
  and ``acquisition_stack1_leftcamera_z_stacking``.
  Gives the order of stacking of in the Z direction

  * ``'direct'``: from the high-contrasted images (small values of z)
    to the fuzzy/blurred ones (large values of z)
  * ``'inverse'``: the other way around.

  See section :ref:`cli-fuse-important-parameters`.

* ``fusion_weighting``: 
  set the weighting function for the weighted sum of the registered
  acquisition stacks (for all channels to be processed).

  * ``'uniform'``: uniform (or constant) weighting, it comes 
    to the average of the resampled co-registered stacks
  * ``'ramp'``: the weights are linearly increasing or 
    decreasing along the Z axis
  * ``'corner'``: the weights are constant in a corner portion 
    of the stack, defined by two diagonals in the XZ-section
  * ``'guignard'``: original historical weighting function, 
    described in Leo Guignard's Phd thesis :cite:p:`guignard:tel-01278725`, 
    that puts more weight to sections close to the camera and take
    also account the traversed material.

  See section :ref:`cli-fuse-nonlinear-registration`.
* ``fusion_weighting_channel_1``:
  set the weighting function for the weighted sum of the registered
  acquisition stacks for the first channel (in case of multi-channel
  acquisition).
* ``fusion_weighting_channel_2``:
  set the weighting function for the weighted sum of the registered
  acquisition stacks for the second channel (in case of multi-channel
  acquisition).
* ``fusion_weighting_channel_3``:
  set the weighting function for the weighted sum of the registered
  acquisition stacks for the third channel (in case of multi-channel
  acquisition).


The following parameters are kept for backward compatibility:

* ``fusion_crop`` 
  same as ``fusion_cropping``
* ``fusion_margin_x_0``
  same as ``fusion_cropping_margin_x_0``
* ``fusion_margin_x_1``
  same as ``fusion_cropping_margin_x_1``
* ``fusion_margin_y_0``
  same as ``fusion_cropping_margin_y_0``
* ``fusion_margin_y_1``
  same as ``fusion_cropping_margin_y_1``
* ``fusion_xzsection_extraction`` 
  same as ``xzsection_extraction``
* ``raw_crop`` 
  same as ``acquisition_cropping``
* ``raw_margin_x_0``
  same as ``acquisition_cropping_margin_x_0``
* ``raw_margin_x_1``
  same as ``acquisition_cropping_margin_x_1``
* ``raw_margin_y_0``
  same as ``acquisition_cropping_margin_y_0``
* ``raw_margin_y_1``
  same as ``acquisition_cropping_margin_y_1``
* ``raw_mirrors`` 
  same as ``acquisition_mirrors``
* ``raw_ori`` 
  same as ``acquisition_orientation``
* ``raw_resolution`` 
  same as ``acquisition_resolution``

* ``begin`` 
  see section :ref:`cli-fuse-important-parameters`
* ``delta``
* ``end`` 
  see section :ref:`cli-fuse-important-parameters`
* ``fusion_weighting``
* ``fusion_weighting_channel_1``
* ``fusion_weighting_channel_2``
* ``fusion_weighting_channel_3``
* ``raw_delay``



``astec_intraregistration`` parameters
--------------------------------------

These parameters are prefixed by ``intra_registration_``.

* Registration parameters 
  (see section :ref:`cli-fuse-image-registration`) 

* ``reference_index``:
  defines the still image after transformation compositions it will 
  only translated, except if ``reference_transformation_file`` 
  or ``reference_transformation_angles`` are set.
  See section :ref:`cli-intraregistration-template`.
* ``reference_transformation_file``:
  resampling transformation to be applied to the reference image 
  (and to the whole serie) after transformation compositions.
  See section :ref:`cli-intraregistration-template`.
* ``reference_transformation_angles``:
  list of rotations wrt the X, Y,or Z axis that defines the resampling
  transformation.

  .. code-block:: python

     reference_transformation_angles = 'X 30 Y 50'

  represents a rotation of 30 degree around the X axis followed by a 
  rotation of 50 degrees around the Y axis.
  
  Beware: rotation composition depends on the order, so 
  ``'X 30 Y 50'`` is not equivalent to ``'Y 50 X 30'``.

* ``template_type``
* ``template_threshold``
* ``margin``
* ``resolution``
* ``rebuild_template``:
  ``True`` or ``False``.
  If set to ``True``, force to recompute the template as well as the
  transformations from existing co-registrations (that are not
  re-computed). It is useful when a first intra-registration has been
  done with only the fusion images: a second intra-registration with
  the segmentation images as template can be done without recomputing 
  the co-registrations.
* ``sigma_segmentation_images``
* ``resample_fusion_images``
* ``resample_segmentation_images``
* ``resample_post_segmentation_images``
* ``movie_fusion_images``
* ``movie_segmentation_images``
* ``movie_post_segmentation_images``
* ``xy_movie_fusion_images``
* ``xz_movie_fusion_images``
* ``yz_movie_fusion_images``
* ``xy_movie_segmentation_images``
* ``xz_movie_segmentation_images``
* ``yz_movie_segmentation_images``
* ``xy_movie_post_segmentation_images``
* ``xz_movie_post_segmentation_images``
* ``yz_movie_post_segmentation_images``
* ``maximum_fusion_images``
* ``maximum_segmentation_images``
* ``maximum_post_segmentation_images``



.. _cli-parameters-mars:

``astec_mars`` parameters
-------------------------

These parameters are prefixed by ``mars_``.

* ``first_time_point``:
  first time point to be segmented by the mars method.
  Overrides the value of the ``begin`` variable.
* ``last_time_point``:
  last time point to be segmented by the mars method.
* Watershed parameters 
  (see section :ref:`cli-parameters-watershed`)
* Seed edition parameters
  (see section :ref:`cli-parameters-seed-edition`)
* Preprocessing parameters
  (see section :ref:`cli-parameters-preprocessing`)
  prefixed by ``seed_``
* Preprocessing parameters
  (see section :ref:`cli-parameters-preprocessing`)
  prefixed by ``membrane_``


.. _cli-parameters-manualcorrection:

``astec_manualcorrection`` parameters
-------------------------------------

* Diagnosis parameters 
  (see section :ref:`cli-parameters-diagnosis`)

* Astec parameters
  (see section :ref:`cli-parameters-astec`)

* ``first_time_point``:
  first time point to be corrected.
  Overrides the value of the ``begin`` variable.
* ``last_time_point``:
  lats time point to be corrected.
* ``input_image``:
  defines the input file names (to be used when correcting
  other files than the ``astec_mars`` output file.
* ``output_image``:
  defines the output file names (to be used when correcting
  other files than the ``astec_mars`` output file.
* ``manualcorrection_dir``:
  path to directory where to find the mapping file.
  * ``manualcorrection_file``:
  path to mapping file for manual correction of a segmentation (ie label)
  image. See above the syntax of this file.
  
  * 1 line per label association
  * background label has value 1
  * the character ``\#`` denotes commented lines 

  Example of ``mapping_file``:

  .. code-block:: none

     # a line beginning by '#' is ignored (comment)

     # lines with only numbers concern changes for the first time point of the time series
     # or the only time point when correcting the segmentation of the first time point
     # - one single number: label of the cell to be divided at the first time point
     # - several numbers: labels of the cells to be fused
     # Hence

     8
     # means that cell of label 8 have to be splitted

     9 2 7
     # means that cells of label 9, 7, and 2 have to be fused

     30 1 
     # means that cell of label 30 have to be fused with the background (of label 1)

     # lines beginning by 'timevalue:' concern changes for the given time point
     # - 'timevalue:' + one single number: label of the cell to be splitted
     # - 'timevalue:' + several numbers: labels of the cells to be fused
     # Note there is no space between the time point and ':'

     8: 7
     # means that cell of label 7 of time point 8 have to be splitted

     10: 14 12 6
     # means that cells of label 14, 12 and 6 of time point 10 have to be fused

     # lines beginning by 'timevalue-timevalue:' concern changes for the given time point range
     # - 'timevalue-timevalue:' + several numbers: labels of the cells to be fused

     10-12: 14 16
     # means that cells of label 14 and 16 of time point 10 have to be fused
     # their offspring will be fused until time point 12


.. _cli-parameters-astec:

``astec_astec`` parameters
--------------------------

These parameters are prefixed by ``astec_``.

* Watershed parameters 
  (see section :ref:`cli-parameters-watershed`)
* Preprocessing parameters
  (see section :ref:`cli-parameters-preprocessing`)
  prefixed by ``seed_``
* Preprocessing parameters
  (see section :ref:`cli-parameters-preprocessing`)
  prefixed by ``membrane_``
* Preprocessing parameters
  (see section :ref:`cli-parameters-preprocessing`)
  prefixed by ``morphosnake_``
* Morphosnake parameters
  (see section :ref:`cli-parameters-morphosnake`)
* ``propagation_strategy``:

  * ``'seeds_from_previous_segmentation'``
  * ``'seeds_selection_without_correction'``
 
* ``previous_seg_method``:
  how to build the seeds :math:`S^e_{t-1 \leftarrow t}` 
  for the computation of :math:`\tilde{S}_{t}`

  * ``'deform_then_erode'``: :math:`S^{\star}_{t-1}` is transformed
    towards :math:`I_t` frame through :math:`\mathcal{T}_{t-1 \leftarrow t}`,
    and then the cells and the background  are eroded.
  * ``'erode_then_deform'``: historical method. The cells 
    and the background of :math:`S^{\star}_{t-1}` are eroded, and
    then transformed
    towards :math:`I_t` frame through :math:`\mathcal{T}_{t-1 \leftarrow t}`.

* ``previous_seg_erosion_cell_iterations``:
  set the cell erosion size for :math:`S^e_{t-1 \leftarrow t}` computation. 
* ``previous_seg_erosion_background_iterations``:
  set the background erosion size for :math:`S^e_{t-1 \leftarrow t}`
  computation. 
* ``previous_seg_erosion_cell_min_size``:
  size threshold. Cells whose size is below this threshold will 
  be discarded seeds in :math:`S^e_{t-1 \leftarrow t}` 

* ``watershed_seed_hmin_min_value``:
  set the :math:`h_{min}` value of the :math:`[h_{min}, h_{max}]` interval.
* ``watershed_seed_hmin_max_value``:
  set the :math:`h_{max}` value of the :math:`[h_{min}, h_{max}]` interval.
* ``watershed_seed_hmin_delta_value``
  set the :math:`\delta h` to go from one :math:`h` to the next in the 
  :math:`[h_{min}, h_{max}]` interval.
  
* ``background_seed_from_hmin``:
  ``True`` or ``False``. 
  Build the background seed at time point :math:`t` by cell propagation.  
* ``background_seed_from_previous``:
  ``True`` or ``False``.
  Build the background seed at time point :math:`t` by using the background 
  seed from :math:`S^e_{t-1 \leftarrow t}`. 
  Fragile. 
  
* ``seed_selection_tau``:
  Set the :math:`\tau` value for division decision (seed selection step).
 
* ``minimum_volume_unseeded_cell``:
  Volume threshold for cells without found seeds in the seed 
  selection step. Cells with volume (in :math:`\tilde{S}_t`) whose size is below
  this threshold and for which no seed was found are discarded.

* ``volume_ratio_tolerance``:
  Ratio threshold to decide whether there is a volume decrease 
  (due to the background) for morphosnake correction.
* ``volume_ratio_threshold``:
  Ratio threshold to decide whether there is a large volume decrease
  for segmentation consistency checking.
* ``volume_minimal_value``:
  Size threshold for seed correction step. For a given cell at time 
  point :math:`t-1`, if the corresponding cell(s) at time point :math:`t` has(ve) 
  volume below this threshold, they are discarded (and the cell at time 
  point :math:`t-1` has no lineage.
  
* ``morphosnake_correction``:
  ``True`` or ``False``. 
* ``outer_correction_radius_opening``





``astec_postcorrection`` parameters
-----------------------------------

These parameters are prefixed by ``postcorrection_``.

* ``volume_minimal_value``
  branch ending with leaf cell below this value are candidate
  for deletion. Expressed in voxel unit.
* ``lifespan_minimal_value``
* ``test_early_division``
* ``test_volume_correlation``
* ``correlation_threshold``
* ``test_postponing_division``
* ``postponing_correlation_threshold``
* ``postponing_minimal_length``
* ``postponing_window_length``
* ``lineage_diagnosis``
  performs a kind of diagnosis on the lineage before and after
  the post-correction.



.. _cli-parameters-contact-surface:

Contact surface parameters
--------------------------

* ``cell_contact_distance``
  Defines the contact surface similarity. Contact surface vectors are normalized before
  comparison (by the l1-norm, so percentages of the total surface are compared). Possible values are:

  * ``l1_distance``: sum of absolute value of coordinate difference (or Manhattan distance).
  * ``l2_distance``: euclidean distance.

  This measure is normalized into [0, 1]: 0 means perfect equality, 1 means total dissimilarity.



.. _cli-parameters-diagnosis:

Diagnosis parameters
--------------------

These parameters are prefixed by ``diagnosis_``.

* Contact surface parameters 
  (see section :ref:`cli-parameters-contact-surface`)

* ``minimal_volume``: for diagnosis on cell volume.
  Threshold on cell volume. Snapshot cells that have a volume below this threshold are displayed.

* ``maximal_volume_variation``: for diagnosis on cell volume.
  Threshold on volume variation along branches. Branches that have a volume variation
  above this threshold are displayed.
  The volume variation along a branch is calculated as 
  :math:`100 * \frac{\max_{t} v(c_t) - \min_{t} v(c_t)}{\mathrm{med}_t v(c_t)}` where :math:`v(c_t)` is the volume
  of the cell :math:`c_t` and :math:`\mathrm{med}` is the median value.

* ``maximal_volume_derivative``: for diagnosis on cell volume.
  Threshold on volume derivative along branches. 
  Time points along branches that have a volume derivative
  above this threshold are displayed.
  The volume derivative along a branch is calculated as 
  :math:`100 * \frac{v(c_{t+1}) - v(c_{t})}{v(c_{t})}` where :math:`t` denotes the successive acquisition time points.

* ``items``: if strictly positif, number minimal of items (ie cells) to be displayed in diagnosis.

* ``minimal_length``: for diagnosis on lineage.
  Threshold on branch length. Branches that have a length below this threshold are displayed.

* ``maximal_contact_distance``: for diagnosis on cell contact surface. Threshold on cell contact surface distance 
  along branches. Time points along branches that have a cell contact surface distance 
  above this threshold are displayed (recall that the distance is in [0, 1]).



.. _cli-parameters-contact-atlas:

``astec_contact_atlas`` parameters
----------------------------------

These parameters are prefixed by ``atlas_``.

* Diagnosis parameters 
  (see section :ref:`cli-parameters-diagnosis`)


* ``atlasFiles``: list of atlas files. An atlas file is a property file that contains lineage,
  names, and contact surfaces for an embryo.

* ``referenceAtlas``: reference atlas. Use for time alignment of atlases. If not provide, the first atlas of
  ``atlasFiles`` is used as reference. Warning, the reference atlas has to be in ``atlasFiles`` list also.

* ``outputDir``: output directory where to write atlas-individualized output files,
  ie morphonet selection files or figure files.

* ``write_selection``: write morphonet selection file on disk.

* ``add_symmetric_neighborhood``: if ``True``, add the symmetric neighborhood as additional exemplar. It means
  that left and right embryo hemisphere are considered together.

* ``differentiate_other_half``: if ``True``, differentiate the cells of the symmetric half-embryo.
  If 'False', consider all the cells of the symmetric half-embryo as a single cell.
  This option has been introduced for test purpose. Please do not consider changing 
  its default value (``True``)

* ``use_common_neighborhood``: the same cell has different neighbors from an atlas to the other.
  If 'True' build and keep an unique common neighborhood (set of neighbors) for all atlases by
  keeping the closest ancestor for neighboring cells. Eg, if a division has occurred in some
  embryos and not in others, daughter cells will be somehow fused so that all neighborhoods only
  exhibit the parent cell. 
  Please do not consider changing its default value (``True``).

* ``name_delay_from_division``: Delay from the division to extract the neighborhooods.
  0 means right after the division.

* ``confidence_delay_from_division``: Delay from the division to compute a name confidence with
  respect to all the atlases (the ones from ``atlasFiles``).
  0 means right after the division.
  The naming confidence is computed from the ``N`` closest atlases: for a given division
  division-to-division distances are computed between the named atlas (to be assessed) and each atlas of 
  ``atlasFiles`` and the ``N`` with the smallest distances are retained.
  The confidence measure is the average of the differences between the division-to-division distance for 
  the wrong choice (names are inverted) and the division-to-division distance for the right choice.

  If there are less atlases for the given division than ``confidence_atlases_nmin``, the confidence is
  not computed. Else ``N`` is the maximum of ``confidence_atlases_nmin`` and the percentage of atlases 
  ``confidence_atlases_percentage``.

* ``confidence_atlases_nmin``: minimum number of atlases to compute the naming
  confidence.

* ``confidence_atlases_percentage``: percentage of atlases (for a given division) to compute the naming
  confidence. The actual number of atlases used is the maximum value of thsi percentage and
  ``confidence_atlases_nmin``.

* ``delay_from_division``: set both ``name_delay_from_division`` and ``confidence_delay_from_division``.

* ``division_contact_similarity``: How to compare two division patterns (a division is considered here
  as the concatenation of the contact surface vectors of the 2 daughter cells). Choices are:

  * ``distance``: the distance type is given by ``cell_contact_distance`` 
    (see section :ref:`cli-parameters-contact-surface`).
    Distances are normalized between 0 (perfect match) and 1 (complete mismatch).
  * ``probability``: 1-(division probability) is used to keep the same meaning
    for the 0 and 1 extremal values. Probabilities are built with the distance
    ``cell_contact_distance``. This is kept for test purposes and should not be used.

* ``diagnosis_properties``: ``True`` or ``False``. Performs some diagnosis when reading an additional 
  property file into the atlases. Incrementing the verboseness ('-v' in the command line) may give more details.

* ``division_permutation_proposal``: If ``True``, will propose some daughters switches in the atlases. 
  For a given division, a global score is computed as the sum of all pairwise division similarity. 
  A switch is proposed for an atlas if it allows to decrease this global score.

* ``dendrogram_cluster_distance``: cluster distance used to build dendrograms.
  Dendrograms are used either for diagnosis purpose (if ``diagnosis_properties`` is set to ``True``)
  or to generate figures (if ``generate_figure`` is set to ``True``)
  See `scipy.cluster.hierarchy.linkage documentation <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html>`_. Choices are:

  * ``'single'``
  * ``'complete'``
  * ``'average'``
  * ``'weighted'``
  * ``'centroid'``
  * ``'median'``
  * ``'ward'``

* ``generate_figure``: if ``True``, generate python files in (prefixed by ``figures_``) that generate figures.

* ``figurefile_suffix``: suffix used to named the above python files as well as the generated figures.


.. _cli-parameters-contact-naming:

``astec_contact_naming`` parameters
-----------------------------------

These parameters are prefixed by ``naming_``.

* ``astec_contact_atlas`` parameters
  (see section :ref:`cli-parameters-contact-atlas`)

* ``inputFile``: input property file to be named. Must contain lineage and contact surfaces
  as well as some input names (one time point should be entirely named).

* ``outputFile``: output property file.

* ``selection_method``: decision method to name the daughters after a division. Distances are computed
  between the couple of daugthers to be named as well as the couple of switched daughters and 
  all the atlas divisions.

  * ``mean``: choose the couple of names that yield the minimal average distance over all the atlases.
  * ``minimum``: choose the couple of names that yield a minimal distance.
    It comes to name after the closest atlas (for this division).
  * ``sum``: same as ``mean``.

* ``testFile``: input property file to be tested (must include cell names), for instance for
   leave-one-out test.
   64-cells time point is searched, cell names at other time points are delated and the embryo 
   is entirely renamed from this given time point. Comparison between new names and actual ones are reported.
   If given, ``inputFile`` is ignored.


* ``test_diagnosis``: if ``True``, some diagnosis are conducted on the property file to be tested.




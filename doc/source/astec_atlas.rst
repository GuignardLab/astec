.. role:: python(code)
   :language: python

.. _cli-astec-atlas:

``astec_atlas``
===============

``astec_atlas`` can be used to assess the quality/coherency of a set of already 
named ascidian embryos (or atlases) as well to point out potential corrections.
It assumes that the set (or mathematically speaking the vector) of valued surface
contact of a named cell can be used as its signature.

``astec_atlas`` additional options
----------------------------------

The following options are available:

``-write-selection, --write-selection``
   write out morphonet selection files


``astec_atlas`` principle
-------------------------

Eg, the cell :math:`c^{R}_i` (the cell :math:`c_i` of atlas :math:`R`) is represented by the 
vector
:math:`c^{R}_i = \left( \begin{array}{c} 
s^{R}_{i,1} \\
\vdots \\
s^{R}_{i,j} \\
\vdots
\end{array}
\right)`.

To account for size differences (between embryos, or between time points within an embryo), normalized  surface
contact vector should be used for computation (parameter ``cell_normalization``, see section :ref:`cli-parameters-atlas`).
A distance between cells comes to a L1 distance 
between two vectors.

From the cell-to-cell distance, a division-to-division similarity can be built, a division being represented by 
the couple of daughter cells (extracted at the distance ``delay_from_division`` from the division).

To enrich division exemplars, the symmetric division (of the other hemi-embryo) can be symmetrized
(ie ``a7.0002_`` will be changed in ``a7.0002*``). This is governed by ``add_symmetric_neighborhood``.

Thus, to assess the quality of a set of atlases, a typical parameter file may be

.. code-block:: python

   atlasFiles = []
   atlasFiles += ['/path_to_reference_embryos/Astec-pm1.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm3.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm4.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm5.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm7.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm8.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm9.pkl']
   #
   # how to select cells
   #
   add_symmetric_neighborhood = True
   use_common_neighborhood = True
   delay_from_division = 3
   # 
   # how to compute distances
   #
   cell_normalization = 'global'
   #
   # how to extract/display information
   #
   atlas_diagnosis = True
   division_diagnosis = True
   division_permutation_proposal = True
   generate_figure = True
   figurefile_suffix = 'some_suffix'

* ``atlas_diagnosis`` and ``division_diagnosis``may be quite verbose. It may be adviced to set them to ``True`` when
  introducing a new atlas, but not when using a set of already curated atlases.
  Two kinds of diagnosis are conducted.

  * ``atlas_diagnosis``

    * on each single atlas/reference file, the ``name`` and the ``contact`` properties are assessed. 
      Such diagnosis can also be conducted with ``astec_embryoproperties``
      (see section :ref:`cli-embryoproperties`)
    * on the population of division neighborhoods:

  * ``division_diagnosis``

    * pairwise disagreements: for a given cell and every couple of reference embryos, 
      the distance of the two divisions (one per reference) is compared to the distance
      of one division compared to the other being switched. If the later is prefered, 
      it is denoted as a disagreement.
    * linkage/dendrogram analysis: it is checked whether adding the switched divisions
      to the set of divisions changes the largest value of cluster distance in
      a dendrogram. If yes, it also suggest that some divisions may be switched.
      Individualized morphonet selection files are written (if ``write_selection`` is set to ``True``)
      in the ``outputDir`` directory.


* ``division_permutation_proposal`` may propose to switch the daughter names of some divisions. It calculates
  whether a name switch result in a global score improvement, and, if so, proposes the switch.
  It is somehow computationally costly.
  Individualized ``morphonet`` selection files are written (if ``write_selection`` is set to ``True``)
  in the ``outputDir`` directory.

* ``generate_figure`` will generate python files (in the ``outputDir`` directory) 
  that, when executed, generates some figures.
  It is somehow computationally costly.


Section :ref:`cli-parameters-astec-atlas` provides a view on all the parameters.



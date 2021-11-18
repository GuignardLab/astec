.. role:: python(code)
   :language: python

.. _cli-astec-contact-atlas:

``astec_contact_atlas``
=======================

``astec_contact_atlas`` can be used to assess the quality/coherency of a set of already 
named embryos (or atlases) as well to point out potential corrections.
It assumes that the set (or mathematically speaking the vector) of valued surface
contact of a named cell can be used as its signature.

Eg, the cell :math:`c^{R}_i` (the cell :math:`c_i` of atlas :math:`R`) is represented by the 
vector
:math:`c^{R}_i = \left( \begin{array}{c} 
s^{R}_{i,1} \\
\vdots \\
s^{R}_{i,j} \\
\vdots
\end{array}
\right)`.

To account for embryo size difference, normalized  surface
contact vector are used for computation (contact surfaces are divided by the total cell surface, 
so values represented the fraction of the total cell surface). A distance between cells comes to a distance 
between two vectors (or the modulus of the difference of the two vectors).
This distance can be set by ``cell_contact_distance`` 
(see section :ref:`cli-parameters-contact-surface`).

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
   delay_from_division = 0
   # 
   # how to compute distances
   #
   cell_contact_distance = 'l1-distance'
   division_contact_similarity = 'distance'
   #
   # how to extract/display information
   #
   diagnosis_properties = True
   daughter_switch_proposal = True
   generate_figure = True
   figurefile_suffix = 'some_suffix'

* ``diagnosis_properties`` may be quite verbose. It may be adviced to set it to ``True`` when
  introducing a new atlas, but not when using a set of already curated atlases.

* ``daughter_switch_proposal`` may propose to switch the names of some divisions. It calculates
  whether a name switch result in a global score improvement, and, if so, proposes the switch.
  It is somehow computationally costly.

* ``generate_figure`` will generate python files that, when executed, generates some figures.
  It is somehow computationally costly.


Section :ref:`cli-parameters-contact-atlas` provides a view on all the parameters.



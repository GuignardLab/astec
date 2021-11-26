.. role:: python(code)
   :language: python

.. _cli-astec-contact-naming:

``astec_contact_naming``
========================

``astec_contact_naming`` can be used either to name an embryo (for which the lineage and the contact surfaces are known), or to evaluate the naming on an already named embryo.

Section :ref:`cli-parameters-contact-naming` provides a view on all the parameters.

``astec_contact_naming`` additional options
-------------------------------------------

The following options are available:

``-write-selection, --write-selection``
   write out morphonet selection files

Naming an embryo
----------------

To name an embryo, the minimal parameter file has to contain:

* the input property file, containing the lineage, the contact surfaces, and some names 
  (typically a 64-cells time point has to be named),

* the output property file,

* the list of atlases (already named embryos).

.. code-block:: python

   inputFile = 'property_file_partially_named.xml'
   outputFile = 'property_file_named_after_the_atlases.xml'

   atlasFiles = []
   atlasFiles += ['/path_to_reference_embryos/Astec-pm1.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm3.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm4.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm5.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm7.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm8.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm9.pkl']

Section :ref:`cli-parameters-contact-naming` provides a view on all the parameters. Most of the parameters concern 
the atlas building and are shared with ``astec_contact_atlas`` (see section :ref:`cli-astec-contact-atlas` and 
section :ref:`cli-parameters-contact-atlas`).


Assessing the naming procedure on an already named embryo
---------------------------------------------------------

If the input file is given to the ``testFile`` variable (see below), it is assumed to be already named.

.. code-block:: python

   testFile = 'property_file_already_named.xml'
   outputFile = 'property_file_named_after_the_atlases.xml'

   atlasFiles = []
   atlasFiles += ['/path_to_reference_embryos/Astec-pm1.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm3.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm4.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm5.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm7.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm8.pkl']
   atlasFiles += ['/path_to_reference_embryos/Astec-pm9.pkl']

Then a 64-cells time point is searched and the embryo is renamed from
this given time point. Comparison between new names and actual ones
are reported.
A selection  (in the ``morphonet`` sense), named ``selection_test_errors`` is added to the output file where

* ``100`` indicates the future naming errors (cells that will divide, whose daughter names will be switched) 

* ``200`` indicates the sister cells whose name has been switched

* ``255`` indicates the other errors

Calling ``astec_contact_naming`` with the ``-extract-selection`` option allows to write selection files (in the ``morphonet`` sense) together with the ouput file.


.. role:: python(code)
   :language: python

.. _cli-astec-contact-naming:

``astec_contact_naming``
========================

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



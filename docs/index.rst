
.. image:: https://www.earlham.ac.uk/sites/default/files/images/pages/Brand/EI-double-large-trans-highquality.png
  :width: 200
  :alt: Alternative text

|

############################################
EMotEP: Earlham Motif Extrapolation pipeline
############################################

EMotEP is a Snakemake based pipeline developed to lift-over regulatory motif
annotations drawn from ChIP-seq data. The pipeline is designed to use data from
well annotated model species, and lift-over to species which have annotated
genomes but no regulatory motif annotations.

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Documentation:

   description
   installation
   usage
   output

************
Dependencies
************

EMotEP uses the following softwares:

Snakemake: `Köster J, Rahmann S. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics. 2018 Oct 15;34(20):3600. <https://academic.oup.com/bioinformatics/article/28/19/2520/290322>`_

Bedtools: `Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. <https://academic.oup.com/bioinformatics/article/26/6/841/244688>`_

****************
Acknowledgements
****************

EMotEP would not be useable without

Gene transcription regulation Database:
`GTRD: a database on gene transcription regulation—2019 update <https://doi.org/10.1093/nar/gky1128>`_
`Database Website <https://gtrd.biouml.org/>`_

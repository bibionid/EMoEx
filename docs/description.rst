###########
Description
###########

.. toctree::
   :maxdepth: 2
   :caption: Desc:

EMotEP is a Snakemake based pipeline developed to lift-over regulatory motif
annotations drawn from ChIP-seq data. The pipeline is designed to use data from
well annotated model species, and lift-over to species which have annotated
genomes but no regulatory motif annotations.

The pipeline uses a mixture of custom python scripts and established softwares
to define putative regulatory regions.

Taking the target species genome and annotation as input, ``cisRegion.py`` creates
an annotation of a user defined region upstream of each gene in the annotation
provided. The distribution of these regions is plotted by ``cisRegion_displot.py``
and then the sequences extracted from the input genome using ``bedtools getfasta``

|

•••••••
Scripts
•••••••
The custom scripts used in the pipeline are described below:

|

.. autoprogram:: workflow.scripts.cisRegion:parser
   :prog: cisRegion.py

.. autoprogram:: workflow.scripts.cisRegion_displot:parser
   :prog: cisRegion_displot.py

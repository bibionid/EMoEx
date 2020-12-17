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

The pipeline uses a mixture of custom python and shell scripts along with
established softwares to define putative regulatory annotations in cis regions
of protein coding genes.

The pipeline will take target and reference genomes and uses the
``sequencePrep.sh`` and ``annotationPrep.sh`` scripts to parse this information
from Ensembl or NCBI formats into those required for later EMotEP scripts.

Taking the target species genome and annotation as input, ``cisRegion.py``
creates an annotation of a user-defined region upstream of each gene in the
annotation provided. The distribution of these regions is plotted by
``cisRegion_displot.py`` and then the sequences extracted from the input genome
using ``bedtools getfasta``.

For the well annotated models species, ChIP-seq data collected by the
`GTRD Database <gtrd.biouml.org>`_ is parsed from .bidBed to .bed by
``GTRD_parse_bigBed.py``. From this point, ``bedtools window`` is used to
retrieve genes within a 10kb window on up and downstream of each binding event.
``bedtools closest`` is then used to retrieve the closest genic annotation up
and downstream of each binding event. Binding events which fall within gene
space are removed at this point. ``GTRD_parse_closest.py`` is used to filter the
remaining annotaions to select a single gene per binding event. This gene is the
gene with the TSS closest to the binding event coordinates.
``bedtools getfasta`` is then used to extract the sequences under these binding
event annotations, and ``GTRD_parse_reformat.py`` is used to transform from
fasta format to a three column format required by the following scripts.

|

*******
Scripts
*******

The custom scripts used in the pipeline are described below:

.. autoprogram:: workflow.scripts.cisRegion:parser
    :prog: cisRegion.py

.. autoprogram:: workflow.scripts.cisRegion_displot:parser
    :prog: cisRegion_displot.py

.. autoprogram:: workflow.scripts.GTRD_parse_bigBed:parser
    :prog: GTRD_parse_bigBed.py

.. autoprogram:: workflow.scripts.GTRD_parse_closest:parser
    :prog: GTRD_parse_closest.py

.. autoprogram:: workflow.scripts.GTRD_parse_reformat:parser
    :prog: GTRD_parse_reformat.py

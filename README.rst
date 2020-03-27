####################################
MoCha Oncogenic MOI Annotator - MOMA
####################################
The MoCha Oncogenic MOI Annotator (MOMA) utility is a sequencing platform
agnostic tool used to annotate variants from a NGS sequencing assay and classify
the variants as Mutations of Interest (MOIs) or Variants of Unknown Significance
(VuS).  These classifications are based on data from annotating the variants
with `Annovar`_ and mapping them to `OncoKB`_. 

In addition, simple filtering is performed on the output data to remove variants
that are not clinically relevant.  We remove variants that are above a set
population frequency as determined from `GnomAD`_, `ExAC`_, and `1000G`_, as
well as non-coding variants and synonymous variants. The remaining calls, then
are mapped to OncoKB and the variant's Oncogenicity and Oncogenic Effect
annotation is added to the variant call if there is a match.

MOMA is designed to be compatible with VCF files from all NGS platforms either
natively or with the addition of a simple helper script to modfy the input to be
compatible with the annotation pipeline.  MOMA can also run starting from a
generic MAF file.  In fact, the first steps of the pipeline when starting from a
VCF file is to annotate the data and generate a MAF file that the tool will use
for downstream processing. As shown below, the spirit of the tool is to be able
to accomodate any kind of data, using helper scripts to staget he data in a way
that can be easily processed through the rest of the tool.

********
Workflow
********

.. image:: docs/moma_workflow.png
   :height: 825px
   :width: 637px

************
Installation
************
Installation details can be found in the documentation PDF file or online at the 
`MOMA Read The Docs Page`_.

.. _Annovar: https://doc-openbio.readthedocs.io/projects/annovar/en/latest/
.. _OncoKB: https://www.oncokb.org/
.. _GnomAD: https://gnomad.broadinstitute.org/
.. _ExAC: https://exac.broadinstitute.org
.. _1000G: https://www.internationalgenome.org/1000-genomes-browsers/
.. _MOMA Read the Docs Page: https:// mocha-oncogenic-moi-annotator.readthedocs.io

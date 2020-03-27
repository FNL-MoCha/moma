.. index:: 
    MOMA; Technical Details
    Technical Details

.. _technical-details:

#################
Technical Details
#################
This section describes a little more details on the inner workings of MOMA,
including the various filtering and annotation steps and the necessary resources
to generate the final report.  It is helpful to consult the flow diagram at the 
beginning of this document (:ref:`moma-workflow`) to get a better sense of how
the data flows and what is happening at each step.

Briefly, we can break the flow of data down into a few components:

    1. VCF Parsing and preparation for Annovar
    2. Annovar annotation
    3. MAF file generation
    4. MOMA MOI / VuS Annotation
    5. Reporting

As the flow diagram shows, there certainly are other steps and components that
occur during a run.  But, these 5 are the main components leading to a final
result. 

*********************
Pipeline Organization
*********************
The MOMA pipeline is constructed of a series of Perl, Python, and BASH scripts
that all work together to work the data through the flow described above.  These
scripts can be categorized as helper scripts used to modify the data in a way
the pipeline can use it, pipeline scripts that shunt the data through the
various paths in the workflow, and analytical scripts that generate data. Most
of these scripts are all stored within the ``scripts/`` directory in the package
root.  


********************
MOI Annotation Rules
********************
Mutations of interest (MOIs) are defined as variants that have an oncogencity
rating per OncoKB annotation. Per OncoKB's website, "OncoKB contains detailed 
information about specific alterations in 671 cancer genes."  

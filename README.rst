#####################
AMG-232 Documentation
#####################

The AMG-232 Reporter is an Ion Torrent Suite plugin developed for version 5.2 of 
the system. The utility will read in a Ion Torrent Variant Caller (TVC) file, and 
generate a report of TP53 variants that are non-synonymous, pathogenic variants
not found in the human population.  Abscence of these kinds of alterations will
make a patient suitable for the AMG-232 Study, and is the rationale for this 
tool.

Requirements
************
The following packages and tools are required to run this plugin:
    - Python 3
    - vcftools
    - vcfExtractor.pl
    - JSON::XS
    - Annovar

As indicated above, Annovar is employed by this plugin to do variant annotation,
and as such some databases need to be loaded.  The following are the databases 
used by this plugin:

Annovar Databases:
    - hg19_clinvar_20170905
    - hg19_cosmic85 (will have to be custom built)
    - hg19_popfreq_all_20150413
    - hg19_dbnsfp35a
    - hg19_ensGene, hg19_ensGeneMrna, and hg19_knownGene (from standard Annovar
      build).
    - hg19_refGene, hg19_refGeneMrna, and hg19_refGeneVersion (from standard 
      Annovar build).

Installation
************
Installation follows the standard Ion Torrent Plugin installation method.  One
can just navigate to the Plugins page in their Ion Torrent Browser, click the 
*Install or Upgrade Plugin* button, and upload the pluging ZIP file.  Additionally
one can just copy this whole git repository to their ``/results/plugins/`` 
directory and then *Rescan Plugins for Changes* from the Ion Torrent Browser 
plugins page.

Once installation is complete, there is a small setup that must be performed 
from the command line.  Since the required Annovar databases are 39G in size,
they can not be included in this repository.  Instead, you must download the 
databases following the guidelines in the `Annovar Documentation 
<http://annovar.openbioinformatics.org/en/latest/user-guide/download/>`_  Once 
the files are downloaded, they should be decompressed and stored in a directory
called ``annovar_db`` in the plugin (full path: 
``AMG232_Reporter/resource/annovar_db/``.


Running the Utility
*******************
This plugin can be run from the Torrent Suite GUI just as any other plugin.  At 
this time, there are no configurable options.  Note that the **Torrent Variant 
Caller** must have been run on any specimens for which you want to generate a 
AMG-232 Reporter report.  If there are no VCF files from TVC (i.e. TVC has not
yet been run), the plugin will complete with an error.  Consult the error logs
for information.

Plugin Output
*************
The plugin will generate a block report indicating whether or not there were
TP53 variants found in the specimen.  Clicking on the *Barcode Name* of the
specimen you are interested in will bring you to a full webpage that outlines
the details of any TP53 variants present, as well as providing links to a CSV
report of the data and a ZIP file of intermediate VCF and Annovar files that 
can be used for troubleshooting and further data analysis.  

Troubleshooting
***************
Running the managing the plugin is fairly straightforward and simple, with little
intervention needed.  If, for some reason, the plugin completes with an error, 
always click the dropdown array next to the plugin status and click on the 
**View Plugin Log** link.  The reason for the failure should be outlined in this
logfile.  

..note:
    The plugin was intended to run on a DNA specimen that has been registered in
    a standard run plan.  Any deviation from that will likely generate an error.


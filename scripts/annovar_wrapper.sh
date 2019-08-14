#!/bin/bash
# Wrapper script to launch the Annovar pipeline.
VERSION='1.5.20190814'

PLUGIN_DIR=$(dirname $(readlink -f $0) | sed 's/\/scripts//')
ANNOVAR_ROOT="${PLUGIN_DIR}/lib/annovar/"
ANNOVAR_DB="${PLUGIN_DIR}/resource/annovar_db/"

function usage() {
    scriptname=$(basename $0)
    echo
    echo "$scriptname - v$VERSION"
    echo "Wrapper script to help run Annovar on VCF file."
    echo 
    echo "USAGE: $scriptname <VCF> [outpath]"
    exit
}

vcf=$1
if [[ -z $vcf ]]; then
    echo "ERROR: you must input a VCF file!"
    exit 1
elif [[ $vcf == '-h' || $vcf == '--help' ]]; then
    usage
elif [[ ! -e $vcf ]]; then 
    echo "ERROR: VCF file '$vcf' can not be found!"
    exit 1
fi

outpath=$2
if [[ -z $outpath ]]; then
    outpath=$vcf
fi

# Annovar cmd
# NOTE: This takes much, much longer if you use the threading option.  Just run
# as single thread.
$ANNOVAR_ROOT/table_annovar.pl \
    -buildver hg19 \
    -polish \
    -remove \
    -nastring . \
    -outfile $outpath \
    -protocol trunc_refGene,cosmic89_noEnst,avsnp142,dbnsfp35a,clinvar_20190305,popfreq_all_20150413 \
    -operation g,f,f,f,f,f \
    -argument '-hgvs -splicing_threshold 5,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' \
    -vcfinput $vcf \
    $ANNOVAR_DB \

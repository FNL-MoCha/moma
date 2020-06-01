#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 06/2019 - D Sims
################################################################################
"""
MoCha Oncogenic Mutation Annotator Tool (MOMA)
This utility will take a VCF file from MoCha pipelines including the Oncomine
Cancer Assay (OCA), TSO500 ctDNA assay, and our current whole exome sequencing
(WES) assay, and generate a report that contains all coding variants that are
non-synonymous and not likely SNPs and adds OncoKB annotations to determine
which of them are oncogenic and not.  In addition to the SNVs and Indels, we can
also add copy number variation (CNV) and fusion data to the output in one master
report. 
"""

import os
import re
import sys
import csv 
import shutil
import argparse
import subprocess

from pprint import pformat # noqa
from textwrap import indent # noqa
from operator import itemgetter
from collections import defaultdict

from lib import logger
from lib import utils

# Globals
global verbose, debug, quiet
verbose = False
debug = False
quiet = False

# Global Paths
output_root = os.getcwd()
outdir_path = os.path.abspath(output_root)
package_root = os.path.dirname(__file__)
scripts_dir = os.path.join(package_root, 'scripts')
lib = os.path.join(package_root, 'lib')
resources = os.path.join(package_root, 'resource')

vfile = os.path.join(package_root, '_version.py')
with open(vfile) as fh:
    verline = fh.readline().rstrip('\n').split('=')[1]
version = verline.lstrip(' ').replace('\'', '')

# Want to sort the results and need a natural sort algorithm. Module natsort is
# not available by default on some systems, so do a local import.
sys.path.insert(0, lib)
from natsort import natsorted

# Default resource files.
oncokb_cnv_file = os.path.join(resources, 'moma_cnv_lookup.tsv')
oncokb_fusion_lookup = os.path.join(resources, 'moma_fusion_genes.tsv')
blacklist_file = os.path.join(resources, 'blacklisted_vars.txt')
reference = os.path.join(resources, 'hg19.fasta.gz')

# Default Thresholds
pfreq = '0.01'   # Max pop frequency < 1%
camp = '4'       # 95% confidence interval >= 4 copies for amplification
closs = '1'      # 5% confidence interval <= 1 copy for deletion.
ct_amp = '1.15'  # Fold Change default for TSO500
ct_loss = '0.85' # Fold change default for TSO500
freads = '250'   # Minimum fusion reads > 250 for call.
fvaf = '0.01'    # Fusion VAF for TSO500


def get_args():
    global verbose, quiet

    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'vcf', 
        metavar="<VCF File>",
        help='VCF file on which to run the analysis. VCF files must be derived '
            'from the Ion Torrent TVC plugin.'
    )

    # Mandatory args
    required_args = parser.add_argument_group('Required Arguments')
    required_args.add_argument(
        '-s', '--source',
        required=True,
        choices=['oca', 'wes', 'tso500', 'genmaf'],
        help='Type of data being loaded. Will determine what steps and script '
            'are necessary to process these data.'
    )

    # Arguments used for filtering data
    filter_args = parser.add_argument_group('Filtering Arguments')
    filter_args.add_argument(
        '-g', '--genes', 
        metavar="<gene>", 
        default="all",
        help='Gene or comma separated list of genes to report. Use the string '
        '"all" to remove the gene filter and report data for all genes. DEFAULT: '
            '%(default)s.'
    )

    filter_args.add_argument(
        '--cu',
        metavar='FLOAT <amp threshold>',
        help='Threshold for calling amplifications. Will be in ploidy space '
            'for Oncomine data, and fold change for TSO500. Only valid with the '
            '`-CNV` option selected. Default: %(default)s.'
    )
    filter_args.add_argument(
        '--cl', 
        metavar='FLOAT <loss threshold>',
        help='Threshold for calling deletions. Will be in ploidy space for '
            'Oncomine data, and fold change for TSO500. Only valid with the '
            '`-CNV` option selected. Default: %(default)s.'
    )
    filter_args.add_argument(
        '--reads',
        metavar='INT <reads>',
        default=freads, 
        help='Threshold for reporting fusions. This filter arg is used for '
             'Oncomine results only, and is nly valid with the `-Fusion` '
             'option selected. Default: %(default)s.'
    )
    filter_args.add_argument(
        '--fvaf',
        metavar='FLOAT <fusion allele frequency<',
        default=fvaf,
        help='Variant Allele Frequency (VAF) threshold for fusion reads. This '
            'filter is only applicable to the TSO500 assay, which does not '
            'report read counts like Oncomine. Default: %(default)s.'
    )
    filter_args.add_argument(
        '-p', '--popfreq',
        metavar='<pop_frequency>',
        type=float,
        default=pfreq,
        help='Population frequency threshold above which variants will be '
        'filtered out. Default: %(default)s.'
    )

    parser.add_argument(
        '-C', '--CNV',
        action='store_true',
        help='Add CNV data to the output. CNV reporting is off by default.'
    )
    parser.add_argument(
        '-F', '--Fusion',
        action='store_true',
        help='Add Fusion data to the output. Fusion reporting is off by default.'
    )
    parser.add_argument(
        '-n', '--name', 
        metavar="<sample_name>", 
        help='Sample name to use for output.'
    )
    parser.add_argument(
        '-o', '--outdir', 
        metavar='<output_directory>',
        help='Directory to which the output data should be written. DEFAULT: '
            '<sample_name>_out/'
    )
    parser.add_argument(
        '-k', '--keep',
        action='store_true',
        help='Keep all intermediate files and do not delete. Useful for '
            'troubleshooting.'
    )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress output to the command line.  Will still write to log '
             'file.'
    )
    parser.add_argument(
        '-V', '--Verbose',
        action='store_true',
        help='Output more log messages and general stderr messages.'
    )
    parser.add_argument(
        '-r', '--rave',
        metavar='<protocol;treatment_id;specimen_id>',
        help='Generate a CSV file that can be uploaded into Rave Web Reporting '
            'to support a clinical trial. The protocol number should be "10231" '
            'in general, the treatment ID will be the site identifier (e.g. '
            '"IA004-0004"), and the specimen ID will be the same as what the '
            'lab used for sequencing (e.g. "10231-6D4410K1-1"). Since we are '
            'using semicolons to delimit, the whole string must be enclosed in '
            'quotes.'
    )
    parser.add_argument(
        '-v', '--version', 
        action='version',
        version='%(prog)s - v' + version
    )
    '''
    NOTE:
    Undocumented arg to help with testing and debugging. This will skip the
    Annovar step so that you can focus on just processing the data. Instead
    of supplying a VCF and running Annovar, you merely need to submit the
    Annovar file as though it were a VCF instead.
    '''
    parser.add_argument(
        '--noanno',
        action='store_true',
        help=argparse.SUPPRESS
    )
    args = parser.parse_args()

    if args.Verbose:
        verbose = True

    if args.quiet:
        sys.stderr.write("Running in quiet mode.")
        quiet = True

    # If we want CNV data, need to get a different threshold depending on the
    # assay type.
    if args.source == 'oca':
        args.cu = camp if not args.cu else args.cu
        args.cl = closs if not args.cl else args.cl

    elif args.source == 'tso500':
        args.cu = ct_amp if not args.cu else args.cu
        args.cl = ct_loss if not args.cl else args.cl

    if args.rave:
        args.rave = __validate_rave_ids(args.rave)
    return args

def get_name_from_vcf(vcf):
    name = ''
    with open(vcf) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                elems = line.split()
                try:
                    name = elems[9]
                except IndexError:
                    # We don't have a name field in this VCF for some reason, so
                    # just use the VCF filename.
                    name = vcf.rstrip('.vcf')
    return name

def simplify_vcf(vcf, outdir, source):
    """
    Use the `simplify_vcf.pl` script to remove reference and NOCALLs from the 
    input VCF. Return a simplified VCF containing only 1 variant per line, and 
    with only the critical VAF and coverage info.  Return the resultant simple
    VCF filename for downstream processing.
    """
    new_name = '{}_simple.vcf'.format(os.path.join(outdir, vcf.rstrip('.vcf')))
    if source == 'ion':
        cmd = [os.path.join(scripts_dir, 'simplify_vcf.pl'), '-f', new_name, vcf]
        status = run(cmd, 'simplify the Ion VCF', silent=not verbose)
    elif source == 'mutect':
        cmd = [os.path.join(scripts_dir, 'mutect2_parser.pl'), '-o', new_name,
                vcf]
        status = run(cmd, 'simplify the MuTect2 VCF', ret_data=True, 
        #  status = run(cmd, 'simplify the MuTect2 VCF',
            silent=not verbose)
        for line in status[2:]:
            log.write_log(None, line)
    if status and source != 'mutect':
        log.write_log("error", "Could not run `simplify_vcf.pl`. Exiting.")
        sys.exit(1)
    else:
        with open(new_name) as fh:
            num_vars = len([line for line in fh if line.startswith('chr')])
        return num_vars, new_name

def run_annovar(simple_vcf, sample_name):
    """
    Run Annovar on the simplified VCF to generate an annotate dataset that can
    then be filtered by gene. Return the resultant Annovar .txt file for 
    downstream processing.
    """
    file_path = f'{outdir_path}/{sample_name}'
    cmd = [os.path.join(scripts_dir, 'annovar_wrapper.sh'), simple_vcf,
            file_path]
    status = run(cmd, 'annotate VCF with Annovar', silent=not verbose)

    if status:
        sys.exit(1)
    else:
        # Rename the files to be shorter and cleaner
        annovar_txt_out = os.path.abspath('%s.hg19_multianno.txt' % file_path)
        annovar_vcf_out = os.path.abspath('%s.hg19_multianno.vcf' % file_path)
        for f in (annovar_vcf_out, annovar_txt_out):
            new_name = f.replace('hg19_multianno', 'annovar')
            os.rename(f, new_name)
        return new_name

def generate_report(annovar_data, genes, outfile):
    """
    Process the OncoKB annotated Annovar data to filter out data by gene, 
    population frequency, and any other filter. Use `wanted_fields` to order the
    output.
    """
    with open(outfile, 'w') as outfh:
        csv_writer = csv.writer(outfh, lineterminator="\n", delimiter=",")

        wanted_fields = ('Chromosome', 'Start_Position', 'Reference_Allele',
            'Tumor_Seq_Allele2', 'i_TumorVAF', 't_ref_count', 't_alt_count', 
            'Existing_variation', 'Hugo_Symbol', 'RefSeq', 'HGVSc', 
            'HGVSp_Short', 'Exon_Number', 'Variant_Classification',
            'SIFT', 'PolyPhen', 'Clinvar_Significance', 'Clinvar_Review_Status', 
            'MOI_Type', 'Oncogenicity', 'Effect')

        csv_writer.writerow(('Chr', 'Pos', 'Ref', 'Alt', 'VAF', 'Ref_Reads',
            'Alt_Reads', 'Variant_Id', 'Gene', 'Transcript', 'CDS',
            'AA', 'Location', 'Function', 'Sift', 'Polyphen', 
            'Clinvar_Significance', 'Clinvar_Review_Status', 'MOI_Type', 
            'Oncogenicity', 'Effect'))

        for var in natsorted(annovar_data.keys(), 
                key=lambda k: (k.split(':')[1], k.split(':')[2])):
            if genes and annovar_data[var]['Hugo_Symbol'] not in genes:
                continue
            try:
                csv_writer.writerow([annovar_data[var][x] for x in wanted_fields])
            except:
                print("ERROR with entry:")
                utils.pp(var)
                raise

def run(cmd, task, ret_data=False, silent=True):
    """
    Generic subprocess runner. Can return results if `ret_data` set to True
    """
    data = []
    with subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT,
            bufsize=1,
            universal_newlines=False
        ) as proc:

        if not silent:
            log.write_log('info', 'Messages from %s:' % task)

        for line in proc.stdout:
            if ret_data:
                data.append(line.decode('utf-8'))
            else:
                if not silent:
                    log.write_log('', line.decode("utf-8"))

    if proc.returncode != 0:
        log.write_log('error', 'An error has occurred while trying to %s' %
                task)
        #  logger.write_log('unformatted', '\n'.join(proc.stdout))
        raise subprocess.CalledProcessError(proc.returncode, ' '.join(proc.args))
    return data

def __load_blacklist():
    with open(blacklist_file) as fh:
        return [line.split()[0] for line in fh]

def __get_avg_pop(annovar_data, population):
    """
    Read in a set of population frequency data, and return the average.
    """
    freq_data = [annovar_data[x] for x in annovar_data if population in x]
    # For some reason, the Annovar population frequency data can show up as '0.'
    # instead of just '0', and the dot screws stuff up.  
    freq_data = list(map(lambda x: x.split('.')[0] if x == '0.' else x,
        freq_data))

    freqs = [float(x) for x in freq_data if x.replace('.','',1).isdigit()]
    return (sum(freqs) / len(freqs)) if freqs else '.'

def __get_vaf(varid, vafstr, source):
    if source == 'oca':
        return vafstr.split(';')[0].split('=')[1]
    elif source == 'wes': 
        return '{:.2f}'.format(float(vafstr.split(':')[2]) * 100.00)
    elif source == 'tso500':
        ref_reads, alt_reads = vafstr['AD'].split(',')
        return '{:.4f}'.format(float(alt_reads) / (float(ref_reads) +
            float(alt_reads)) * 100.00)

def run_var_sig(vcf, source, outfile):
    """
    Run the `calculate_deamination_score.pl` script to get some mutational
    signature data, along with Ts/Tv ratio and Deamination Score (C>T at other
    sites compared to C>T at CpG sites).
    """
    platform_map = {
        'oca'        : 'ion',
        'tso500'     : 'illumina',
        'wes'        : 'illumina',
        'wes_mutect' : 'illumina'
    }
    cmd = [os.path.join(scripts_dir, 'calc_tstv_deam.py'), '-s',
        platform_map[source], '-r', reference, '-o', outfile, vcf]
    status = run(cmd, 'Calculating Ts/Tv ratio, deamination score, and '
        'generating SBS-6 data.', ret_data=True, silent=not verbose)
    log.write_log('info', "Results from SBS script:")

    # Capture the output we want from the script, avoiding error messages.
    wanted_data = [ x.rstrip('\n') for x in status if x.startswith('\t') or 
            x.startswith('F')]

    for i in wanted_data:
        log.write_log(None, i)

    results = {}
    for line in wanted_data:
        if line.startswith('\t '):
            key, val = (line.replace(' ','').split(':'))
            results[key.replace('\t', '')] = val

    tstv = results.get('Ts/Tv', None)
    deam = results.get('DeamScore', None)

    if float(deam) > 1.7:
        log.write_log('warn', f'Deamination score is high ({deam}). Sample may '
                'be damaged.')
    if float(tstv) > 4:
        log.write_log('warn', f'Ts/Tv ratio is high ({tstv}). Sample may '
            'be damaged.')

def oncokb_annotate(annovar_file, source, popfreq):
    """
    Create a temp file that our TSO500, OncoKB driven annotator can use and
    annotate using those rules.  Those annotations will be added to the dataset
    for filtering and reporting downstream.
    """
    # Make a truncated MAF file of the annovar data that can be read by OncoKB
    # annotator.
    log.write_log('info', 'Converting Annovar data to a truncated MAF file.')
    trunc_maf = annovar_file.replace('.annovar.txt', '.maf')
    annovar2maf_cmd = [os.path.join(scripts_dir, 'annovar2maf.pl'), '-o', 
            trunc_maf, annovar_file]
    run(annovar2maf_cmd, 'Annovar2MAF', silent=not verbose)

    # Use TSO500 to annotate the truncated MAF file.
    log.write_log('info', 'Running MOMA')
    annotated_maf = trunc_maf.replace('.maf', '.oncokb.maf')
    oncokb_annot_cmd = (
        os.path.join(scripts_dir, 'moma.pl'), 
        '-a', 'oncokb',
        '-p', 'exac:%f' % popfreq, 
        '--no-trim', 
        '-V', 
        '-o', annotated_maf, 
        '-l', '/dev/null',
        trunc_maf
    )

    summary_report = run(oncokb_annot_cmd, 'MoCha OncoKB MOI Annotator', 
            ret_data=True, silent=not verbose)

    # Always want the summary report from MOMA to make it to at least the log
    # file, even if we have asked for non-verbose run.
    log.write_log('info', 'MOMA Summary Report:')
    for line in summary_report:
        log.write_log('', line)

    # Run final filters and prepare data for report.
    log.write_log('info', 'Running final filters and preparing data for report.')
    return (__filter_oncokb_data(annotated_maf, source))

def __filter_oncokb_data(annotated_maf, source):
    """
    Filter out unwanted calls and prepare data for final report.
    """
    final_data = defaultdict(dict)
    blacklist_vars = __load_blacklist()

    # Load the OncoKB Annotated data first...
    with open(annotated_maf) as amaf_fh:
        header = amaf_fh.readline().rstrip('\n').split('\t')

        for line in amaf_fh:
            #  data = line.rstrip('\n').split('\t')
            var_data = dict(zip(header, line.rstrip('\n').split('\t')))
            varid = ':'.join(
                itemgetter(*['Hugo_Symbol', 'Chromosome','Start_Position', 
                    'Reference_Allele', 'Tumor_Seq_Allele2'])(var_data)
            )

            # Check to see if variant is in blacklist
            if varid in blacklist_vars:
                log.write_log('note', 'Removing variant {}:{} because it is '
                    'blacklisted.'.format(varid, var_data['HGVSp_Short']))
                continue

            # Remove variants with low VAF if not running ctDNA panel.
            if source != 'tso500' and float(var_data['i_TumorVAF']) < 0.01:
                log.write_log('note', 'Removing variant {} because VAF '
                    'is <1% and this is not ctDNA assay data.'.format(varid))
                continue

            # Run the Routbort Rule; Filter out any variants with alt reads
            # < 25.
            if source != 'tso500' and float(var_data['t_alt_count']) < 25:
                log.write_log('note', 'Removing variant {} because Alt '
                    'reads are less than 25 (alt reads: {}).'.format(
                        varid, var_data['t_alt_count']))
                continue

            # Change the output format for some fields.
            var_data['i_TumorVAF'] = '{:.2f}'.format(
                    float(var_data['i_TumorVAF']) * 100.00)
            var_data['Exon_Number'] = 'exon{}'.format(
                    var_data['Exon_Number'].split('/')[0])


            # Add the HGVS ids to each. Remap the keys for the function, just to
            # keep the function more generic for potential use in other bits.
            func_keys = ('gene', 'chr', 'pos', 'ref', 'alt', 'cds', 'aa',
                    'transcript')
            current_keys = ('Hugo_Symbol', 'Chromosome', 'Start_Position',
                'Reference_Allele', 'Tumor_Seq_Allele2', 'HGVSc', 'HGVSp',
                'RefSeq')
            map_data = dict(zip(func_keys, itemgetter(*current_keys)(var_data)))

            var_data.update(dict(zip(('HGVSg_full', 'HGVSc_full', 'HGVSp_full'),
                utils.gen_hgvs(map_data))))
            final_data[varid] = var_data

    #  utils.pp(final_data)
    #  utils.__exit__(432,'')
    return dict(final_data)

def gen_cnv_report(cnv_data_file, cu, cl, genelist, outfile, source):
    """
    If required, run an Oncomine CNV report to get that data from the VCF file
    as well.
    """
    cnv_data = []
    parsed_cnv_results = []
    if source == 'oca':
        cmd = (os.path.join(scripts_dir, 'get_cnvs.pl'), '--cu', cu, '--cl', cl,
                cnv_data_file)
        parsed_cnv_results = run(cmd, "Generating Oncomine CNV results data",
            True, silent=not verbose)
    elif source == 'tso500':
        cmd = (os.path.join(scripts_dir, 'tso500_cnvs.pl'), cnv_data_file)
        parsed_cnv_results = run(cmd, "Generating TSO500 CNV results data",
                True, silent=not verbose)
    elif source == 'wes':
        logger.write_log('error', 'CNV data acquisition is not yet implemented '
                'for WES panels yet. Skipping this step.')
        return None

    header = parsed_cnv_results.pop(0).rstrip('\n').split(',')

    # If we have results, then generate a report.  Else, bail out.
    if parsed_cnv_results:
        for r in parsed_cnv_results:
            data = dict(zip(header, r.rstrip('\n').split(',')))

            # Fix measurements to accomodate multiple assays.
            if source == 'oca':
                cn            = float(data['Raw_CN'])
                lower_reading = float(data['CI_05'])
                upper_reading = float(data['CI_95'])
            elif source == 'tso500':
                cn = lower_reading = upper_reading = float(data['FoldChange'])

            #  cnv_str = '{},{}'.format(data['Gene'], data['Raw_CN'])
            cnv_str = '{},{}'.format(data['Gene'], cn) 

            if lower_reading >= float(cu):
                data.update({'var_type' : 'Amplification'})
            elif upper_reading <= float(cl):
                data.update({'var_type' : 'Deletion'})

            # Remove any genes not in our filter list if there is one.
            if genelist and data['Gene'] not in genelist:
                log.write_log('debug', 'Filtering out {} because it is not '
                        'in the requested list.'.format(cnv_str))
                continue

            # Filter events outside of our thresholds.
            if not __cnv_filter(lower_reading, upper_reading, cu, cl):
                log.write_log('debug', 'Filtering out {} because it is not '
                    'an copy number event within the threshold (copy loss={}; '
                    'copy_gain={}).'.format(cnv_str, cl, cu))
                continue

            cnv_data.append(data)

    if not cnv_data:
        log.write_log('info', 'No copy number amplifications or deletions '
                'found in this specimen.')
    else:
        log.write_log('info', 'Found {} copy number events in the '
            'specimen.'.format(len(cnv_data)))
        log.write_log('info', 'Annotating with OncoKB lookup.')
        __oncokb_cnv_and_fusion_annotate(cnv_data, oncokb_cnv_file, 'cnv')

    # Set up column headers for different assays.
    assay_elems = {
        'oca' : ('Chr', 'Start', 'End', 'Gene', 'CN', 'CI_05', 'CI_95', 'MAPD',
            'Oncogenicity', 'Effect'),
        'tso500' : ('Chr', 'Start', 'End', 'Gene', 'FoldChange', 'Oncogenicity',
            'Effect')
    }

    with open(outfile, 'w') as outfh:
        csv_writer = csv.writer(outfh, lineterminator="\n", delimiter=",")
        wanted = assay_elems[source]

        csv_writer.writerow(wanted)
        if cnv_data:
            for elem in cnv_data:
                data = [elem[x] for x in wanted]
                csv_writer.writerow(data)
        else:
            outfh.write("No CNVs found.\n")

    return cnv_data

def __cnv_filter(lower_reading, upper_reading, cu, cl):
    if lower_reading > float(cu):
        return True
    elif upper_reading < float(cl):
        return True
    return False

def __oncokb_cnv_and_fusion_annotate(data, lookup_file, datatype):
    okb_lookup = {}
    oncokb_fields = ('Alteration', 'Oncogenicity', 'Effect')

    with open(lookup_file) as fh:
        _ = fh.readline()
        for line in fh:
            elems = line.rstrip('\n').split('\t')
            okb_lookup[elems[0]] = elems[1:]

    for d in data:
        gene = d['Gene'] if datatype == 'cnv' else d['Driver_Gene']
        if gene in okb_lookup and d['var_type'] == okb_lookup[gene][0]:
            okb_annot = dict(zip(oncokb_fields, okb_lookup[gene]))
        else:
            okb_annot = dict(zip(oncokb_fields, ['.', '.', '.']))
        d.update(okb_annot)

def gen_fusion_report(vcf, reads, genelist, filename):
    """
    Generate a report of Oncomine Fusion data. Start with a 100 reads filter,
    but then let this script handle the final filter in order to record those
    filtered out in the logs.
    """
    
    '''

    if source == 'oca':
        cmd = (os.path.join(scripts_dir, 'get_cnvs.pl'), '--cu', cu, '--cl', cl,
                cnv_data_file)
        parsed_cnv_results = run(cmd, "Generating Oncomine CNV results data",
            True, silent=not verbose)
    elif source == 'tso500':
        cmd = (os.path.join(scripts_dir, 'tso500_cnvs.pl'), cnv_data_file)
        parsed_cnv_results = run(cmd, "Generating TSO500 CNV results data",
                True, silent=not verbose)
    elif source == 'wes':
        logger.write_log('error', 'CNV data acquisition is not yet implemented '
                'for WES panels yet. Skipping this step.')
        return None

    header = parsed_cnv_results.pop(0).rstrip('\n').split(',')
    '''


    cmd = [
        os.path.join(scripts_dir, 'get_fusions.pl'), 
        '--reads', '100',
        vcf
    ]

    if genelist:
        cmd[3:3] = ['--gene', ','.join(genelist)]

    fusion_data = []
    fusion_results = run(cmd, "Generating Oncomine Fusion results data", True,
            silent=not verbose)
    header = fusion_results.pop(0).rstrip('\n').split(',')

    if fusion_results:
        for f in fusion_results:
            data = dict(zip(header, f.rstrip('\n').split(',')))
            fusion_str = '{}.{}.{}'.format(data['Fusion'], data['Junction'],
                data['ID'])

            if int(data['Read_Count']) < int(reads):
                log.write_log('debug', 'Filtering out {} because number of '
                    'reads is not high enough ({}).'.format(fusion_str,
                    data['Read_Count']))
                continue
            data['var_type'] = 'Fusions'
            fusion_data.append(data)

    if not fusion_data:
        log.write_log('info', 'No Fusion events found in this specimen.')
    else:
        log.write_log('info', 'Found {} Fusion events in the '
            'specimen.'.format(len(fusion_data)))
        log.write_log('info', 'Annotating the captured fusions with OncoKB '
            'lookup.')
        __oncokb_cnv_and_fusion_annotate(fusion_data, oncokb_fusion_lookup,
        'fusion')

    with open(filename, 'w') as outfh:
        csv_writer = csv.writer(outfh, lineterminator="\n", delimiter=",")
        wanted = ('Fusion', 'Junction', 'ID', 'Driver_Gene', 'Partner_Gene', 
            'Read_Count', 'Oncogenicity', 'Effect')
        csv_writer.writerow(wanted)
        if fusion_data:
            for elem in fusion_data:
                data = [elem[x] for x in wanted]
                csv_writer.writerow(data)
        else:
            outfh.write("No Fusions found.\n")

    return fusion_data

def combine_reports(sample_name, mut_report, cnv_report, fusion_report, path):
    final_report_name = '{}_moma_report_{}.csv'.format(sample_name,
            utils.today())

    log.write_log('info', 'Collating all variant data into a single report')

    report_path = os.path.join(path, final_report_name)
    with open(report_path, "w") as outfh:
        outfh.write(f':::  SNV / Indel Report for {sample_name}  :::\n')
        var_data = __read_report(mut_report)
        if len(var_data) == 1:
            outfh.write(f'{var_data[0]}\nNo SNVs or Indels found.\n')
        else:
            outfh.write('\n'.join(var_data))

        if cnv_report is not None:
            outfh.write(f'\n\n:::  CNV Data Report for {sample_name}  :::\n')
            var_data = __read_report(cnv_report)
            outfh.write('\n'.join(var_data))

        if fusion_report is not None:
            outfh.write(f'\n\n::: Fusion Data Report for {sample_name}  :::\n')
            var_data = __read_report(fusion_report)
            outfh.write('\n'.join(var_data))
        outfh.write('\n')

    return report_path

def run_moma2rave(moma_report, path, rave_args):
    """
    Generate a CSV file that can be uploaded into Rave.
    """
    pnum, tpid, spid = rave_args.split(';')
    
    msg = ('\tProtocol Number:      {}\n\tTreatment Patient ID: {}\n\tSpecimen '
        'ID:          {}\n'.format(pnum, tpid, spid))
    log.write_log('unformatted', msg)

    outfile = os.path.join(path, moma_report.replace('moma_report', 'rave_report'))
    cmd = [os.path.join(scripts_dir, 'moma2rave.py'), '-p', pnum, '-t', 
            tpid, '-s', spid, '-o', outfile, moma_report]
    status = run(cmd, 'Generate a Rave compatible CSV file', silent=not verbose)
    if status:
        sys.exit(1)

def __validate_rave_ids(ids):
    ''' Validate the IDs passed to the ``rave`` arg. '''
    elems = ids.split(';')

    if elems[0] != '10231':
        sys.stderr.write('WARN: Protocol number must always be 10231 for '
            'now (may be more protocols later). Changing to match.\n')
        elems[0] = '10231'

    trid_regex = r'[A-Z]{2}\d{3}-\d{4}'
    if not re.fullmatch(trid_regex, elems[1]):
        sys.stderr.write(f'ERROR: {elems[1]} is not a valid treatment patient '
            'ID string.\n')
        sys.exit(1)

    spid_regex = r'^' + elems[0] + '-\w+-\d$'
    if not re.fullmatch(spid_regex, elems[2]):
        sys.stderr.write(f'ERROR: {elems[2]} is not a valid specimen ID.\n')
        sys.exit(1)
    return ';'.join(elems)

def __read_report(report):
    with open(report) as fh:
        return [line.rstrip('\n') for line in fh]

def __verify_vcf(vcf, source, noanno=False):
    '''
    Make sure that the VCF loaded matches the source string so that we can make
    sure that we're parsing it correctly.

    oca: ##source="tvc 5.12-11 (35c9114f) - Torrent Variant Caller"
    tso500: ##source=Pisces 5.2.11.62
    wes: ##source=Platypus_Version_0.7.9.5
    wes: ##joinAdjacentSNPs # Rajesh reommends this approach...but no always
    reliable?
    wes: ##SentieonCommandLine # Seems we're always using this pipeline.

    '''
    id_strings = {
        'oca'    : '##source="tvc',
        'tso500' : '##source=Pisces',
        'wes'    : '##SentieonCommandLine'
    }

    # Let anything pass as a correct file if we use the `--noanno` option
    if noanno: return True 

    with open(vcf) as fh:
        header = [line for line in fh if line.startswith('#')]

    if any(id_strings[source] in x for x in (header)):
        # We can have two different WES piplines, tumor only (haplotype caller)
        # or matched tumor / normal pair (MuTect2 somatic caller).  Need to know
        # which is which so that we can handle the variant columns
        # appropriately. Check to see if the tumor defined in the  VCF header 
        # to know if we have MuTect2 or Haplotype Caller run; we'll only have a
        # normal if a paired normal was run and not the "panel of normals".
        if source == 'wes':
            if any(x.startswith('##tumor_sample') for x  in header):
                source = 'wes_mutect'
        return source

    # If we didn't get a match, throw an error.
    log.write_log('error', f'Source type "{source}" does not appear to be '
        'correct for this VCF. Check the source value.')
    sys.exit(1)

def __jsonify_moma_data(sample_name, snv_data, cnv_data, fusion_data, outfile):
    # Create a JSON blob of data that we can use for databasing and whatnot.

    parsed_data = {sample_name : {'snv_indels' : snv_data, 'cnvs' : cnv_data,
        'fusions' : fusion_data}}
    #  utils.print_json(parsed_data)
    utils.make_json(outfile=outfile, data=parsed_data)


def main(vcf, data_source, sample_name, genes, popfreq, get_cnvs, cu, cl, 
        get_fusions, fusion_threshold, outdir, keep_intermediate_files, noanno, 
        rave_report):

    global outdir_path, debug, verbose, log

    if sample_name is None:
        sample_name = get_name_from_vcf(vcf)

    # Set up the logger.
    if verbose:
        debug = True
    level = 'debug' if debug else 'error'
    logfile = 'moma_reporter_{}_{}.log'.format(sample_name, utils.today())

    # By default logger is set to quiet. If we want both logfile and stdout
    # output, then change this arg.
    log = logger.Logger(loglevel=level, colored_output=True, dest=logfile,
            quiet=quiet)
            
    # Get some usage stats
    cmd = [os.path.join(scripts_dir, 'usage.sh')]
    run(cmd, 'Collect usage information')

    if not os.path.exists(vcf):
        log.write_log('error', f'VCF file {vcf} can not be processed! '
                'No such file or directory!')
        sys.exit(1)

    # Write out welcome and header
    pipeline_version = ''
    with open(os.path.join(package_root, '_version.py')) as fh:
        pipeline_version = fh.readline().rstrip("'\n'").split("'")[1]
    welcome_str = ('{0}\n:::::  Running the MoCha Oncogenic Mutation Annotator '
        '(MOMA) pipeline - v{1}  :::::\n{0}\n'.format('='*97, 
            pipeline_version))
    log.write_log('header', welcome_str)

    # Verify the VCF source is correct for the VCF that we've loaded.
    data_source = __verify_vcf(vcf, data_source, noanno)

    if debug or noanno:
        log.write_log('unformatted', '\n\u001b[33m{decor}  TEST VERSION!  {decor}'
            '\u001b[0m\n\n'.format(decor='='*25))

    log.write_log('info', f'Processing sample: {sample_name}')
        
    if noanno:
        log.write_log('unformatted', 
            '\u001b[33m\t\t===>  Running without Annovar annotation  <===\u001b[0m\n\n')
        annovar_file = vcf
        outdir_path = os.getcwd()

        log.write_log('debug', 
                f'Selected destination dir for data is: {outdir_path}')
    else:
        if outdir is None:
            outdir_path = os.path.join(outdir_path, '%s_out' % sample_name)
        else:
            outdir_path = os.path.join(outdir_path, outdir)
        log.write_log('debug', 
                f'Selected destination dir for data is: {outdir_path}')

        if not os.path.exists(outdir_path):
            os.mkdir(os.path.abspath(outdir_path), 0o755)

        # Simplify the VCF
        if data_source == 'oca':
            log.write_log('info', 'Simplifying the VCF file.')
            num_vars, simple_vcf = simplify_vcf(vcf, outdir_path, 'ion')
            log.write_log('info', 'Done simplifying VCF file. There are %s '
                'variants in this specimen.' % str(num_vars))
        elif data_source == 'wes_mutect':
            log.write_log('info', 'MuTect2 Tumor / Normal run detected. Need '
                'to remove normal data from the VCF.')
            num_vars, simple_vcf = simplify_vcf(vcf, outdir_path, 'mutect')
            log.write_log('info', 'Done simplifying VCF file. There are %s '
                'somatic variants in this specimen.' % str(num_vars))
        else:
            log.write_log('info', 'VCF already ready for processing with '
                'Annovar.')
            simple_vcf = vcf

        # Run Mutational Signature QC to determine Ts/Tv and deamination score.
        log.write_log('info', 'Generating Ts/Tv, Deamination Score, and '
            'SBS-6 information for sample.')
        outfile = f'{outdir_path}/{sample_name}_sbs_metrics.csv'
        sbs_input_vcf = simple_vcf if data_source != 'oca' else vcf
        #  run_var_sig(vcf, data_source, outfile)
        run_var_sig(sbs_input_vcf, data_source, outfile)

        # Annotate the vcf with ANNOVAR.
        log.write_log('info', 'Annotating the simplified VCF file with '
             'Annovar.')
        annovar_file = run_annovar(simple_vcf, sample_name)

    # Implement the MOMA here.
    log.write_log('info', 'Running OncoKB Filter rules on data to look for '
            'oncogenic variants.')
    oncokb_annotated_data = oncokb_annotate(annovar_file, data_source, popfreq) 

    # Generate a filtered CSV file of results for the report.
    mutation_report = os.path.join(outdir_path, 
        f'{sample_name}_mocha_snv_indel_report.csv')
    log.write_log('info', f'Writing report data to {mutation_report}.')
    generate_report(oncokb_annotated_data, genes, mutation_report)

    # Generate a CNV report if we're asking for one.
    cnv_report = None
    oncokb_cnv_data = []
    if get_cnvs:
        log.write_log('info', 
            f"Generating a CNV report for sample '{sample_name}'.\n")
        cnv_report = os.path.join(outdir_path, 
            f'{sample_name}_mocha_cnv_report.csv')
        log.write_log('info', f'Writing report data to {cnv_report}.')

        # We need different files depending on the assay.
        cnv_data_file = ''
        if data_source == 'oca':
            cnv_data_file = vcf
        elif data_source == 'tso500':
            cnv_data_file = '{}.cnv.fc.txt'.format(vcf.replace('.vcf', ''))

        if not os.path.exists(cnv_data_file):
            log.write_log('error', f'CNV file "{cnv_data_file}" can not be found.'
                'Make sure the file is in the same location as the VCF.')
            sys.exit(1)

        oncokb_cnv_data = gen_cnv_report(cnv_data_file, cu, cl, genelist, 
                cnv_report, data_source)
        log.write_log('info', 'Done with CNV report.')

    # Generate a Fusions report if we've asked for one.
    fusion_report = None
    oncokb_fusion_data = []
    if get_fusions:

        units = '%' if data_source == 'tso500' else ' reads'

        log.write_log('info', 'Generating a Fusions report for sample {}. Using '
            'threshold of {}{}.'.format(sample_name, fusion_threshold, units))
        fusion_report= os.path.join(outdir_path,
            f'{sample_name}_mocha_fusion_report.csv')
        log.write_log('info', f'Writing report data to {fusion_report}.')

        if data_source == 'oca':
            fusion_data_file = vcf
        elif data_source == 'tso500':
            fusion_data_file = '{}.fusion.txt'.format(vcf.replace('.vcf', ''))

        oncokb_fusion_data = gen_fusion_report(fusion_data_file, 
            fusion_threshold, genelist, fusion_report)
        log.write_log('info', "Done with Fusions report.")

    # Combine the three reports into one master report.
    moma_report = combine_reports(sample_name, mutation_report, cnv_report, 
            fusion_report, outdir_path)

    # If we need a Rave report, generate one now.
    if rave_report:
        log.write_log('info', 'Generating a CSV file that can be uploaded '
                'into Theradex Rave.')
        run_moma2rave(moma_report, outdir_path, rave_report)

    # Generate a JSON blob of data for databasing and other purposes
    json_out = os.path.join(outdir_path, f'{sample_name}_moma_data.json')
    __jsonify_moma_data(
        sample_name, 
        snv_data=oncokb_annotated_data,
        cnv_data=oncokb_cnv_data, 
        fusion_data=oncokb_fusion_data, 
        outfile=json_out
    )

    # Clean up and move the logfile to data dir.
    contents = [os.path.join(outdir_path, f) for f in os.listdir(outdir_path)]
    if not keep_intermediate_files:
        for f in contents:
            if any(f.endswith(x) for x in ('truncmaf', 'avinput', 'oncokb.maf', 
                'simple.vcf', 'snv_indel_report.csv', 'cnv_report.csv', 
                'fusion_report.csv')):
                    log.write_log('debug', f'Removing {f}.')
                    os.remove(f)

    # Move our logfile into the output dir now that we're done.
    shutil.move(os.path.abspath(logfile), os.path.join(outdir_path, logfile))

    log.write_log('info', 'MoCha OncoKB Annotator Pipeline completed '
         'successfully! Data can be found in %s.' % outdir_path)

    if debug:
        sys.stderr.write('\u001b[33m{decor}\u001b[0m\n'.format(decor='='*67))
        sys.exit()

if __name__ == '__main__':
    args = get_args()

    """
    Check the environment to make sure we have everything we need before we get
    started.
    """
    try:
        subprocess.check_call([os.path.join(scripts_dir, 'env_check.py')])
    except subprocess.CalledProcessError:
        sys.exit(1)

    if not os.path.exists(os.path.join(lib, 'annovar')):
        sys.stderr.write("ERROR: Annovar has not been added to this package. "
            "You must first install Annovar and add this binary to the 'lib' "
            "directory.\n")
        sys.exit(1)
    elif not os.path.exists(os.path.join(resources, 'annovar_db')):
        sys.stderr.write("ERROR: The required Annovar database files are not "
            "found. You must get a copy of the appropriate database files and "
            "put them in the `resources` dir.\n")
        sys.exit(1)

    if args.genes == 'all':
        genelist = []
    else:
        genelist = args.genes.split(',')

    # Centralize this into one metric rather than trying to sort through 
    # different ones depending on the platform.
    fusion_threshold = 0
    if args.source == 'oca':
        fusion_threshold = int(args.reads)
    elif args.source == 'tso500':
        fusion_threshold = float(args.fvaf)

    main(args.vcf, args.source, args.name, genelist, args.popfreq, args.CNV, 
            args.cu, args.cl, args.Fusion, fusion_threshold, args.outdir,
            args.keep, args.noanno, args.rave)

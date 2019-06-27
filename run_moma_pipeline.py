#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline wrapper script for the MoCha Oncomine Clinical Annotation Tool (MOMA)
plugin. 
"""

import os
import sys
import csv 
import argparse
import subprocess

from pprint import pprint as pp # noqa
from pprint import pformat
from textwrap import indent
from operator import itemgetter
from collections import defaultdict
from natsort import natsorted

from lib import logger
from lib import utils # noqa

version = '1.5.20190625-dev0'

# Globals
output_root = os.getcwd()
package_root = os.path.dirname(__file__)
scripts_dir = os.path.join(package_root, 'scripts')
resources = os.path.join(package_root, 'resource')
lib = os.path.join(package_root, 'lib')

debug = True

logger = logger.Logger(loglevel='debug', colored_output=True,
        dest='moma_reporter_{}.log'.format(utils.today()))

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'vcf', 
        metavar="<VCF File>",
        help='VCF file on which to run the analysis. VCF files must be derived '
            'from the Ion Torrent TVC plugin.'
    )
    parser.add_argument(
        '-g', '--genes', 
        metavar="<gene>", 
        default="TP53",
        help='Gene or comma separated list of genes to report. Use the string '
        '"all" to remove the gene filter and report data for all genes. DEFAULT: '
            '%(default)s.'
    )
    parser.add_argument(
        '-n', '--name', 
        metavar="<sample_name>", 
        help='Sample name to use for output.'
    )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress output to the command line.  Will still write to log '
             'file.'
    )
    parser.add_argument(
        '-o', '--outdir', 
        metavar='<output_directory>',
        help='Directory to which the output data should be written. DEFAULT: '
            '<sample_name>_out/'
    )
    parser.add_argument(
        '-v', '--version', 
        action='version',
        version='%(prog)s - v' + version
    )
    args = parser.parse_args()

    if args.quiet is False:
        logger.quiet = False

    if args.genes == 'all':
        args.genes = None

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

def simplify_vcf(vcf, outdir):
    """
    Use the `simplify_vcf.pl` script to remove reference and NOCALLs from the 
    input VCF. Return a simplified VCF containing only 1 variant per line, and 
    with only the critical VAF and coverage info.  Return the resultant simple
    VCF filename for downstream processing.
    """
    new_name = '{}_simple.vcf'.format(os.path.join(outdir, vcf.rstrip('.vcf')))
    cmd = [os.path.join(scripts_dir, 'simplify_vcf.pl'), '-f', new_name, vcf]
    status = run(cmd, 'simplify the Ion VCF')
    if status:
        logger.write_log("error", "Could not run `simplify_vcf.pl`. Exiting.")
        sys.exit(1)
    else:
        with open(new_name) as fh:
            num_vars = len([line for line in fh if line.startswith('chr')])
        return num_vars, new_name

def run_annovar(simple_vcf):
    """
    Run Annovar on the simplified VCF to generate an annotate dataset that can
    then be filtered by gene. Return the resultant Annovar .txt file for 
    downstream processing.
    """
    cmd = [os.path.join(scripts_dir, 'annovar_wrapper.sh'), simple_vcf]
    status = run(cmd, 'annotate VCF with Annovar')
    if status:
        sys.exit(1)
    else:
        # Rename the files to be shorter and cleaner
        annovar_txt_out = os.path.abspath('%s.hg19_multianno.txt' % simple_vcf)
        annovar_vcf_out = os.path.abspath('%s.hg19_multianno.vcf' % simple_vcf)
        for f in (annovar_vcf_out, annovar_txt_out):
            new_name = f.replace('vcf.hg19_multianno', 'annovar')
            os.rename(f, new_name)
        return new_name

def generate_report(annovar_data, genes, outfile, outdir):
    """
    Process the OncoKB annotated Annovar data to filter out data by gene, 
    population frequency, and any other filter.
    """
    outfh = open(os.path.join(outdir, outfile), 'w')
    csv_writer = csv.writer(outfh, lineterminator="\n", delimiter=",")
    wanted_fields = ('Chromosome', 'Start_Position', 'Reference_Allele',
            'Tumor_Seq_Allele2', 'vaf', 'Hugo_Symbol', 'Transcript_ID', 'HGVSc',
            'HGVSp_Short', 'Exon_Number', 'Variant_Classification', 'sift', 
            'polyphen', 'CLNSIG', 'CLNREVSTAT', 'MOI_Type', 
            'Oncogenicity', 'Effect')

    csv_writer.writerow(('Chr', 'Pos', 'Ref', 'Alt', 'VAF', 'Gene',
        'Transcript', 'CDS', 'AA', 'Location', 'Function', 'Sift', 'Polyphen',
        'Clinvar_Significance', 'Clinvar_Review_Status', 'MOI_Type', 
        'Oncogenicity', 'Effect'))

    for var in natsorted(annovar_data.keys(), 
            key=lambda k: (k.split(':')[0], k.split(':')[1])):
        if genes is None:
            csv_writer.writerow([annovar_data[var][x] for x in wanted_fields])
        else:
            if annovar_data[var]['Hugo_Symbol'] in genes:
                csv_writer.writerow([annovar_data[var][x] for x in wanted_fields])

    logger.write_log('info', 'MoCha OncoKB Annotator Pipline completed '
         'successfully! Data can be found in %s.' % outdir)

def run(cmd, task):
    """
    Generic subprocess runner.
    """
    # TODO: We can either capture all output with either only a stdout->PIPE
    # call or stdout -> PIPE plus stderr -> STDOUT (to mix the two).  Doesn't
    # seem to be a way to handle the two together while still getting a real
    # time buffer flush.  This is really not set up too ideally!
    with subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT,
            bufsize=1,
            universal_newlines=False
        ) as proc:
        for line in proc.stdout:
            logger.write_log('', line, raw=True)
    '''
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, bufsize=1)
    msg, err = proc.communicate()
    '''
    if proc.returncode != 0:
        logger.write_log('error', 'An error has occurred while trying to %s' %
                task)
        logger.write_log('', subprocess.CalledProcessError(proc.returncode,
            ' '.join(proc.args)))
        raise subprocess.CalledProcessError(proc.returncode, ' '.join(proc.args))
    return 0

def __marry_oncokb_data(annovar_file, annotated_maf):
    """
    Add the oncokb annotations to the annovar data so that we can use this in
    the report generation step later on.  
    """
    final_data = defaultdict(dict)

    with open(annotated_maf) as amaf_fh:
        header = amaf_fh.readline().rstrip('\n').split('\t')
        for line in amaf_fh:
            data = line.rstrip('\n').split('\t')
            varid = ':'.join(itemgetter(*[1,2,4,5])(data))
            final_data[varid] = dict(zip(header, data))

    with open(annovar_file) as annovar_fh:
        header = annovar_fh.readline().rstrip('\n').split('\t')
        for line in annovar_fh:
            data = line.rstrip('\n').split('\t')
            varid = ':'.join(itemgetter(*[0,1,3,4])(data))
            if varid in final_data:
                new = {
                    'sift'       : __translate(data[13]),
                    'polyphen'   : __translate(data[16]),
                    'vaf'        : __get_vaf(data[115]),
                    'CLNREVSTAT' : data[84],
                    'CLNSIG'     : data[85],
                }
                final_data[varid].update(new)
    return dict(final_data)

def __translate(string):
    sp_conversion = {
        'T' : 'Tolerated',
        'D' : 'Damaging',
        'B' : 'Benign',
        'P' : 'Probably Damaging',
        '.' : '-'
    }
    return sp_conversion.get(string, '???')

def __get_vaf(vafstr):
    return vafstr.split(';')[0].split('=')[1]

def oncokb_annotate(annovar_file):
    """
    Create a temp file that our TSO500, OncoKB driven annotator can use and
    annotate using those rules.  Those annotations will be added to the dataset
    for filtering and reporting downstream.
    """
    # Make a truncated MAF file of the annovar data that can be read by OncoKB
    # annotator.
    logger.write_log('info', 'Converting Annovar data to a truncated MAF file.')
    trunc_maf = annovar_file.replace('.txt', '.truncmaf')
    annovar2maf_cmd = [os.path.join(scripts_dir, 'annovar2maf.pl'), '-o', 
            trunc_maf, annovar_file]
    run(annovar2maf_cmd, 'Annovar2MAF')

    # Use TSO500 to annotate the truncated MAF file.
    logger.write_log('info', 'Running MOMA')
    annotated_maf = trunc_maf.replace('.annovar.truncmaf', '.oncokb.tsv')
    oncokb_annot_cmd = (
        os.path.join(scripts_dir, 'moma.pl'), 
        '-a', 'oncokb',
        '-p', 'exac:0.05', 
        '--no-trim', 
        '-V', 
        '-o', annotated_maf, 
        '-l', '/dev/null',
        trunc_maf
    )
    run(oncokb_annot_cmd, 'TSO500 MOI Annotator')

    # Marry the new OncoKB annotations with the original Annovar Data.
    logger.write_log('info', 'Marrying OncoKB data to Annovar data.')
    return(__marry_oncokb_data(annovar_file, annotated_maf))

def main(vcf, sample_name, genes, outdir, quiet):
    # Create an output directory based on the sample_name
    if sample_name is None:
        sample_name = get_name_from_vcf(vcf)
    logger.write_log('info', f'Sample Name is: {sample_name}')
        
    outdir_path = os.path.abspath(output_root)
    if outdir is None:
        outdir_path = os.path.join(outdir_path, '%s_out' % sample_name)
    else:
        outdir_path = os.path.join(outdir_path, outdir)
    logger.write_log('debug', 
            f'Selected destination dir for data is: {outdir_path}')

    if not os.path.exists(outdir_path):
        os.mkdir(os.path.abspath(outdir_path), 0o755)

    # Simplify the VCF
    logger.write_log('info', 'Simplifying the VCF file.')
    num_vars, simple_vcf = simplify_vcf(vcf, outdir_path)
    logger.write_log('info', 'Done simplifying VCF file. There are %s variants'
            ' in this specimen.' % str(num_vars))

    # Annotate the vcf with ANNOVAR.
    # TODO: Comment back in.
    logger.write_log('info', 'Annotating the simplified VCF file with Annovar.')
    annovar_file = run_annovar(simple_vcf)

    # TODO: delete this.
    #  annovar_file = os.path.join(outdir_path,
    #  'MSN31633_v2_MSN31633_RNA_v2_simple.annovar.txt')

    # Implement the MOMA here.
    logger.write_log('info', 'Running OncoKB Filter rules on data to look for '
            'oncogenic variants.')
    oncokb_annotated_data = oncokb_annotate(annovar_file)

    # Generate a filtered CSV file of results for the report.
    filename = f'{sample_name}_mocha_moi_report.csv' 
    logger.write_log('info', f'Writing report data to {filename}.')
    generate_report(oncokb_annotated_data, genes, filename, outdir_path)

    utils.__exit__(315)
    # Clean up and move the logfile to data dir.
    contents = [os.path.join(outdir_path, f) for f in os.listdir(outdir_path)]
    for f in contents:
        if any(f.endswith(x) for x in ('truncmaf', 'avinput')):
            logger.write_log('debug', f'Removing {f}.')
            os.remove(f)


if __name__ == '__main__':
    args = get_args()

    pipeline_version = ''
    with open(os.path.join(package_root, '_version.py')) as fh:
        pipeline_version = fh.readline().rstrip("'\n'").split("'")[1]
    welcome_str = ('{0}\n:::::  Running the MoCha OncoKB Annotation (MOMA) '
         'pipeline - Version {1}  :::::\n{0}\n'.format('='*88, pipeline_version))
    logger.write_log('header', welcome_str)

    if debug:
        args_dict = pformat(vars(args))
        logger.write_log('debug', "Args passed to the script:\n{}".format(
            indent(args_dict, '    ')))

    main(args.vcf, args.name, args.genes, args.outdir, args.quiet)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline wrapper script for the MoCha Oncomine Clinical Annotation Tool (MOMA)
plugin. 
"""

import os
import sys
import csv 
import shutil
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

version = '1.8.20190708-dev1'

# Globals
output_root = os.getcwd()
package_root = os.path.dirname(__file__)
scripts_dir = os.path.join(package_root, 'scripts')
resources = os.path.join(package_root, 'resource')
lib = os.path.join(package_root, 'lib')
oncokb_cnv_file = os.path.join(package_root, 'resource', 'oncokb_cnv_lookup.tsv')
oncokb_fusion_lookup = os.path.join(package_root, 'resource',
        'oncokb_fusion_genes.tsv')

debug = True

logfile = 'moma_reporter_{}.log'.format(utils.today())
logger = logger.Logger(loglevel='debug', colored_output=True, dest=logfile)

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
        default="all",
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
        '-C', '--CNV',
        action='store_true',
        help='Add CNV data to the output. CNV reporting is off by default.'
    )
    parser.add_argument(
        '--cu',
        metavar='INT <5% CI>',
        default='4',
        help='Threshold for calling amplifications. Default: %(default)s.'
    )
    parser.add_argument(
        '--cl', 
        metavar='INT <95% CI>',
        default='1',
        help='Threshold for calling deletions. Default: %(default)s.'
    )
    parser.add_argument(
        '-F', '--Fusion',
        action='store_true',
        help='Add Fusion data to the output. Fusion reporting is off by default.'
    )
    parser.add_argument(
        '--reads',
        metavar='INT <reads>',
        default='250',
        help='Threshold for reporting fusions. Default: %(default)s.'
    )
    parser.add_argument(
        '-k', '--keep',
        action='store_true',
        help='Keep all intermediate files and do not delete. Useful for '
            'troubleshooting.'
    )
    parser.add_argument(
        '-v', '--version', 
        action='version',
        version='%(prog)s - v' + version
    )
    args = parser.parse_args()

    if args.quiet is False:
        logger.quiet = False

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

def generate_report(annovar_data, genes, outfile):
    """
    Process the OncoKB annotated Annovar data to filter out data by gene, 
    population frequency, and any other filter.
    """
    with open(outfile, 'w') as outfh:
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
            if genes and annovar_data[var]['Hugo_Symbol'] not in genes:
                continue
            csv_writer.writerow([annovar_data[var][x] for x in wanted_fields])

def run(cmd, task, ret_data=False):
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
        logger.write_log('info', 'Messages from %s:' % task)
        for line in proc.stdout:
            if ret_data:
                data.append(line.decode('utf-8'))
            else:
                logger.write_log('', line.decode("utf-8"))

    if proc.returncode != 0:
        logger.write_log('error', 'An error has occurred while trying to %s' %
                task)
        #  logger.write_log('', subprocess.CalledProcessError(proc.returncode))
        #  logger.write_log('', subprocess.CalledProcessError(proc.returncode,
            #  ' '.join(proc.args)))
        raise subprocess.CalledProcessError(proc.returncode, ' '.join(proc.args))

    return data

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
    run(oncokb_annot_cmd, 'MoCha OncoKB MOI Annotator')

    # Marry the new OncoKB annotations with the original Annovar Data.
    logger.write_log('info', 'Marrying OncoKB data to Annovar data.')
    return(__marry_oncokb_data(annovar_file, annotated_maf))

def gen_cnv_report(vcf, cu, cl, genelist, outfile):
    """
    If required, run an Oncomine CNV report to get that data from the VCF file
    as well.
    """
    cnv_data = []
    cmd = (
        os.path.join(scripts_dir, 'get_cnvs.pl'), 
        '--cu', cu,
        '--cl', cl,
        vcf
    )
    cnv_results = run(cmd, "Generating Oncomine CNV results data", True)
    header = cnv_results.pop(0).rstrip('\n').split(',')

    # If we have results, then generate a report.  Else, bail out.
    if cnv_results:
        for r in cnv_results:
            data = dict(zip(header, r.rstrip('\n').split(',')))
            cnv_str = '{},{}'.format(data['Gene'], data['Raw_CN'])
            if float(data['CI_05']) >= float(cu):
                data.update({'var_type' : 'Amplification'})
            elif float(data['CI_95']) <= float(cl):
                data.update({'var_type' : 'Deletion'})

            if genelist and data['Gene'] not in genelist:
                logger.write_log('debug', 'Filtering out {} because it is not '
                        'in the requested list.'.format(cnv_str))
                continue
            if not __cnv_filter(data, cu, cl):
                logger.write_log('debug', 'Filtering out {} because it is not '
                    'an copy number event within the threshold (5% CI={}; '
                    '95% CI={}).'.format(cnv_str, data['CI_05'], data['CI_95']))
                continue
            cnv_data.append(data)

    if not cnv_data:
        logger.write_log('info', 'No copy number amplifications or deletions '
                'found in this specimen.')
    else:
        logger.write_log('info', 'Found {} copy number events in the '
            'specimen.'.format(len(cnv_data)))
        logger.write_log('info', 'Annotating with OncoKB lookup.')
        __oncokb_cnv_and_fusion_annotate(cnv_data, oncokb_cnv_file, 'cnv')

    with open(outfile, 'w') as outfh:
        csv_writer = csv.writer(outfh, lineterminator="\n", delimiter=",")
        wanted = ('Chr', 'Start', 'End', 'Gene', 'CN', 'CI_05', 'CI_95', 'MAPD',
                'Oncogenicity', 'Effect')
        csv_writer.writerow(wanted)
        if cnv_data:
            for elem in cnv_data:
                data = [elem[x] for x in wanted]
                csv_writer.writerow(data)
        else:
            outfh.write("No CNVs found.\n")

def __cnv_filter(data, cu, cl):
    if float(data['CI_05']) > int(cu):
        return True
    elif float(data['CI_95']) < int(cl):
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

    #  okb_annot = {} 
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
    cmd = [
        os.path.join(scripts_dir, 'get_fusions.pl'), 
        '--reads', '100',
        vcf
    ]

    if genelist:
        cmd[3:3] = ['--gene', ','.join(genelist)]

    fusion_data = []
    fusion_results = run(cmd, "Generating Oncomine Fusion results data", True)
    header = fusion_results.pop(0).rstrip('\n').split(',')

    if fusion_results:
        for f in fusion_results:
            data = dict(zip(header, f.rstrip('\n').split(',')))
            fusion_str = '{}.{}.{}'.format(data['Fusion'], data['Junction'],
                data['ID'])

            if int(data['Read_Count']) < int(reads):
                logger.write_log('debug', 'Filtering out {} because number of '
                    'reads is not high enough ({}).'.format(fusion_str,
                    data['Read_Count']))
                continue
            data['var_type'] = 'Fusions'
            fusion_data.append(data)

    if not fusion_data:
        logger.write_log('info', 'No Fusion events found in this specimen.')
    else:
        logger.write_log('info', 'Found {} Fusion events in the '
            'specimen.'.format(len(fusion_data)))
        logger.write_log('info', 'Annotating the captured fusions with OncoKB '
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

def combine_reports(sample_name, mut_report, cnv_report, fusion_report, path):
    final_report_name = '{}_moi_report_{}.csv'.format(sample_name,
            utils.today())

    logger.write_log('info', 'Collating all variant data into a single report')

    with open(os.path.join(path, final_report_name), "w") as outfh:
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
            #  if (var_data) == 1:
                #  outfh.write(f'{var_data[0]}\nNo CNVs found.\n')
            #  else:
                #  outfh.write('\n'.join(var_data))

        if fusion_report is not None:
            outfh.write(f'\n\n::: Fusion Data Report for {sample_name}  :::\n')
            var_data = __read_report(fusion_report)
            outfh.write('\n'.join(var_data))
            #  if len(var_data) == 1:
                #  outfh.write(f'{var_data[0]}\nNo Fuions found.\n')
            #  else:
                #  outfh.write('\n'.join(var_data))

def __read_report(report):
    with open(report) as fh:
        return [line.rstrip('\n') for line in fh]

def main(vcf, sample_name, genes, get_cnvs, cu, cl, get_fusions, reads,
        outdir, quiet, keep_intermediate_files):

    pipeline_version = ''
    with open(os.path.join(package_root, '_version.py')) as fh:
        pipeline_version = fh.readline().rstrip("'\n'").split("'")[1]
    welcome_str = ('{0}\n:::::  Running the MoCha OncoKB Annotation (MOMA) '
         'pipeline - Version {1}  :::::\n{0}\n'.format('='*88, pipeline_version))
    logger.write_log('header', welcome_str)

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
    logger.write_log('info', 'Annotating the simplified VCF file with Annovar.')
    annovar_file = run_annovar(simple_vcf)

    # Implement the MOMA here.
    logger.write_log('info', 'Running OncoKB Filter rules on data to look for '
            'oncogenic variants.')
    oncokb_annotated_data = oncokb_annotate(annovar_file)

    # Generate a filtered CSV file of results for the report.
    mutation_report = os.path.join(outdir_path, 
        f'{sample_name}_mocha_snv_indel_report.csv')
    logger.write_log('info', f'Writing report data to {mutation_report}.')
    generate_report(oncokb_annotated_data, genes, mutation_report)

    # Generate a CNV report if we're asking for one.
    cnv_report = None
    if get_cnvs:
        logger.write_log('info', 'Generating a CNV report for sample {}. Using '
                'thresholds 5% CI = {}, 95% CI = {}'.format(sample_name, cu,
                    cl))
        cnv_report = os.path.join(outdir_path, 
            f'{sample_name}_mocha_cnv_report.csv')
        logger.write_log('info', f'Writing report data to {cnv_report}.')
        gen_cnv_report(vcf, cu, cl, genelist, cnv_report)
        logger.write_log('info', 'Done with CNV report.')

    # Generate a Fusions report if we've asked for one.
    fusion_report = None
    if get_fusions:
        logger.write_log('info', 'Generating a Fusions report for sample {}. '
                'Using threshold of {} reads.'.format(sample_name, reads))
        fusion_report= os.path.join(outdir_path,
            f'{sample_name}_mocha_fusion_report.csv')
        logger.write_log('info', f'Writing report data to {cnv_report}.')
        gen_fusion_report(vcf, reads, genelist, fusion_report)
        logger.write_log('info', "Done with Fusions report.")

    # Combine the three reports into one master report.
    combine_reports(sample_name, mutation_report, cnv_report, fusion_report, 
            outdir_path)

    # Clean up and move the logfile to data dir.
    contents = [os.path.join(outdir_path, f) for f in os.listdir(outdir_path)]
    if not keep_intermediate_files:
        for f in contents:
            if any(f.endswith(x) for x in ('truncmaf', 'avinput')):
                logger.write_log('debug', f'Removing {f}.')
                os.remove(f)
    # Move our logfile into the output dir now that we're done.
    shutil.move(os.path.abspath(logfile), os.path.join(outdir_path, logfile))

    logger.write_log('info', 'MoCha OncoKB Annotator Pipline completed '
         'successfully! Data can be found in %s.' % outdir_path)

if __name__ == '__main__':
    args = get_args()

    args_dict = pformat(vars(args))
    logger.write_log('debug', "Args passed to the script:\n{}".format(
        indent(args_dict, '    ')))

    if args.genes == 'all':
        genelist = []
    else:
        genelist = args.genes.split(',')

    main(args.vcf, args.name, genelist, args.CNV, args.cu, args.cl, 
            args.Fusion, args.reads, args.outdir, args.quiet, args.keep)

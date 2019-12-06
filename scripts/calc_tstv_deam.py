#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Re-write of the `calculate_deamination_score.pl` script to run faster since we
# have a direct Cython implementation of samtools (`pysam`).  This speeds the
# generation of these data up at least 3x!  The best solution will be to re-code
# in C directly, but this is a great step along the way.
# TODO:
#  - Can we speed this up with parallelization for the context generation?
#  - 
# 11/26/2019 - D Sims
################################################################################
"""
Read in an Ion Torrent or Illumina VCF and output single base substitution
metrics indicating the counts for the 6 general base substitutions (i.e. C>A,
C>T, C>G, T>A, T>C, T>G), along with the transition : transversion ratio (Ts/Tv)
and a deamination score (defined as the number of C>T changes at non-CpG sites
divided by the number of C>T changes at CpG sites).
"""

import sys
import os
import re
import argparse
#import shutil
import pysam
import subprocess

from pprint import pprint as pp # noqa


# Need to import the natsort module, which is not commonly installed.  So, point
# to a local copy of the module.
lib_path = os.path.join(os.path.dirname(__file__), '../', 'lib')
sys.path.insert(0, lib_path)
import natsort


version = '1.0.120319'
global quiet

reference=os.path.join(os.path.abspath(os.path.dirname(__file__)), '..', 
        'resource', 'hg19.fasta.gz')

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'vcfs',
        nargs='+',
        metavar='<VCFs>',
        help='VCF(s) file(s) to process.'
    )
    parser.add_argument(
        '-s', '--source',
        required=True,
        metavar='<VCF source>',
        choices=['ion', 'illumina'],
        help='Source of the data.  Allowed values are either "ion" or '
            '"illumina"'
    )
    parser.add_argument(
        '-n', '--numprocs',
        type = int,
        metavar='<num_procs>',
        default = 0,
        help="Number of parallel processes to run if parallel processing. If "
            "you wish to avoid parallel processing, enter 0. Default: "
            "%(default)s"
    )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Do not print extra messages to stdout. Only print the final '
            'result table.'
    )
    parser.add_argument(
        '-r', '--reference',
        metavar='<reference>',
        default=reference,
        help='Reference genome to use for getting sequence contexts. Default: '
            '%(default)s'
    )
    parser.add_argument(
        '-o', '--outfile',
        metavar='<outfile>',
        help='Write data to file instead of just printing to stdout.'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s - v' + version
    )
    args = parser.parse_args()
    return args

def parse_vcf(vcf, source):
    if not quiet:
        print(f"Processing {vcf}...")
    if source == 'ion':
        return run_ion_parser(vcf)
    else:
        return run_illumina_parser(vcf)

def run_illumina_parser(vcf):
    vcf_data = []

    # Rajesh says that the correct sample name will always be in the VCF file's
    # name, and that we can rely on that to map the VCF when a paired tumor /
    # normal is sequenced.
    sample_name = os.path.basename(vcf).split('.')[0]

    tumor_index = 0
    with open(vcf) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                elems = line.rstrip('\n').split('\t')
                for i, e in enumerate(elems):
                    if elems[i] == sample_name:
                        tumor_index = i
                        break
            elif line.startswith('chr'):
                if 'AD' not in line: continue
                elems = line.strip('\n').split('\t')
                if any(',' in x for x  in (elems[3], elems[4])): continue
                if elems[6] in ('PASS', '.'):
                    if './.' in elems[tumor_index]: continue
                    vaf, ref_cof, alt_cov = __get_vaf(elems, tumor_index)
                    if float(vaf) < 1: continue
                    vcf_data.append(["{}:{}".format(elems[0], elems[1]),
                        "{}>{}".format(elems[3], elems[4])])

    return sample_name, vcf_data

def run_ion_parser(vcf):
    sample_name = __get_sample_name(vcf)
    vcf_data = []

    cmd = [
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
            "simplify_vcf.pl"), 
        "-o", "csv", vcf
    ]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    data = out.decode('ascii').split('\n')

    for d in data:
        elems = d.split(',')
        if elems and elems[0].startswith('chr'):
            if float(elems[3]) < 1:
                continue
            else:
                vcf_data.append([elems[0], '{}>{}'.format(elems[1], elems[2])])
    return sample_name, vcf_data


def proc_vcfs_parallel(vcfs, source):
    sys.stderr.write('not yet implemented!')
    sys.exit(1)

def get_var_metrics(vcf_data, sample_name):
    global quiet 

    tot_vars = 0
    indels = 0

    # Get a copy of all contexts to fill for this sample. 
    sbs_96 = contexts.copy()
    
    for data in vcf_data:
        coord, change = data
        tot_vars += 1

        # Count and get rid of indels. 
        if len(change) > 3:
            indels += 1
            continue

        chrom, pos = coord.split(':')
        query = '{}:{}-{}'.format(
            'chr%s' % chrom.lstrip('chr'), 
            int(pos) - 1,
            int(pos) + 1)

        if any(change.startswith(x) for x in ('C', 'T')):
            context = get_context(query)
        else:
            context = __rev_comp(get_context(query))
            comp = __rev_comp(change)
            change = '>'.join(reversed(comp))
        first, second, third = list(context)
        c_string = '{}[{}]{}'.format(first, change, third)

        if c_string not in sbs_96:
            sys.stderr.write(f"Can't find {c_string} in the 96 patterns hash.\n")
        else:
            sbs_96[c_string] += 1

    if not quiet:
        print(f'\tFound {tot_vars} variants in {sample_name}:')
        print("\t  SNVs:          {}".format(tot_vars - indels))
        print(f"\t  Indels:        {indels}")
        print("\tComputing SBS-6 and variant motifs...", end='')

    tstv, var_motifs, transitions, transversions = get_motifs(sbs_96)
    deam_score = '{:.2f}'.format(int(var_motifs['C>T_at_other']) /
            int(var_motifs['C>T_at_CpG']))

    if not quiet:
        print("Done!")
        print(f'\t  Transitions:   {transitions}')
        print(f'\t  Transversions: {transversions}')
        print(f'\t  Ts/Tv:         {tstv}')
        print(f'\t  Deam Score:    {deam_score}')

    return tot_vars, indels, tstv, deam_score, var_motifs, sbs_96

def proc_vcfs(vcfs, source):
    base_data = {}
    keys = ['tot_vars', 'num_indels', 'titv', 'deam_score', 'sbs_6', 'sbs_96']
    
    for vcf in vcfs:
        sample_name, data = parse_vcf(vcf, source)
        base_data[sample_name] = dict(zip(keys, get_var_metrics(data, sample_name)))

    return base_data

def get_motifs(var_contexts):
    transitions = 0
    transversions = 0
    motifs = {
        'C>A'          : 0,
        'C>T'          : 0,
        'C>G'          : 0,
        'T>A'          : 0,
        'T>C'          : 0,
        'T>G'          : 0,
        'C>T_at_CpG'   : 0,
        'C>T_at_other' : 0
    }

    regex = re.compile(r'([ACTG])\[([ACTG]>[ACTG])\]([ACTG])')
    for context in var_contexts:
        (first_base, change, last_base) = re.search(regex, context).group(1,2,3)
        motifs[change] += var_contexts[context]
        if change == 'C>T':
            transitions += var_contexts[context]
            if context.endswith('G'):
                motifs['C>T_at_CpG'] += var_contexts[context]
            else:
                motifs['C>T_at_other'] += var_contexts[context]
        elif change == 'T>C':
            transitions += var_contexts[context]
        else:
            transversions += var_contexts[context]

    tstv = '{:.3f}'.format(transitions / transversions)
    return tstv, motifs, transitions, transversions


def get_context(query):
    return pysam.faidx(reference, query).split('\n')[1]

def __rev_comp(string):
    comp = {'A' : 'T', 'C' : 'G', 'G': 'C', 'T' : 'A'}
    new_str = ''
    bases = [char for char in string if char in ('A', 'C', 'T', 'G')]
    for b in bases:
        new_str += comp[b]
    return ''.join(reversed(new_str))

def __get_sample_name(vcf):
    with open(vcf) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                return line.rstrip('\n').split('\t')[-1]

def __get_vaf(vcf_elems, info_index):
    if info_index is None:
        info_index = 9
    vcf_format = dict(zip(vcf_elems[8].split(':'), vcf_elems[info_index].split(':')))

    ref_reads, alt_reads = vcf_format['AD'].split(',')
    try:
        vaf = '{:.2f}'.format(
            (float(alt_reads) / (float(alt_reads) + float(ref_reads))) * 100)
    except:
        sys.stderr.write("WARNING: Can not get VAF from alt/ref reads for this "
            "entry:\n")
        pp(format)
        sys.stderr.write("Trying to get allele frequency from format field "
            "entry.\n")
        for k,v in vcf_format.items():
            if 'AF' in k:
                vaf = v
                break
        if vaf is None:
            return 0
    return vaf, ref_reads, alt_reads

def __make_context_hash():
    contexts = {}
    bases = ('A', 'C', 'T', 'G')
    mut_classes = ('[C>A]', '[C>G]', '[C>T]', '[T>A]', '[T>C]', '[T>G]')
    for first_base in bases:
        for c in mut_classes:
            for third_base in bases:
                contexts[f'{first_base}{c}{third_base}'] = 0
    return contexts

def __get_width(elems):
    colwidth = 0
    for e in elems:
        colwidth = len(e) if len(e) > colwidth else colwidth
    return colwidth + 2

def print_results(data, outfile):
    if outfile:
        sys.stderr.write(f"Writing output to {outfile}.\n")
        outfh = open(outfile, 'w')
    else:
        outfh = sys.stdout

    sbs_order = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G', 'C>T_at_CpG',
            'C>T_at_other']
    field_widths = {
        'sample'       : __get_width(data.keys()),
        'tot_vars'     : 9,
        'indels'       : 7,
        'tstv'         : 6,
        'C>A'          : 5,
        'C>G'          : 5,
        'C>T'          : 5,
        'T>A'          : 5,
        'T>C'          : 5,
        'T>G'          : 5,
        'C>T_at_CpG'   : 11,
        'C>T_at_other' : 12,
        'deam_score'   : 20,
    }
    header = ['Sample', 'Tot_Vars', 'Indels', 'Ts/Tv'] + sbs_order + ['C>T_Ratio']

    fline = ('{:{sp}<{sample}} {:{sp}<{tot_vars}} {:{sp}<{indels}} {:{sp}<{tstv}} '
        '{:{sp}<{C>A}} {:{sp}<{C>G}} {:{sp}<{C>T}} {:{sp}<{T>A}} {:{sp}<{T>C}} '
        '{:{sp}<{T>G}} {:{sp}<{C>T_at_CpG}} {:{sp}<{C>T_at_other}} '
        '{:{sp}<{deam_score}}\n')

    outfh.write(fline.format(*header, sp=' ', **field_widths))

    for sample in natsort.natsorted(data):
        sbs_data = [data[sample]['sbs_6'][x] for x in sbs_order]
        outdata = [
            sample,
            data[sample]['tot_vars'], 
            data[sample]['num_indels'],
            data[sample]['titv'],
            *sbs_data,
            data[sample]['deam_score']]
        outfh.write(fline.format(*outdata, sp=' ', **field_widths))

def main(vcfs, source, numprocs, quiet, reference, outfile):
    base_data = {}

    if numprocs > 0:
        base_data = proc_vcfs_parallel(vcfs, source)
    else:
        base_data = proc_vcfs(vcfs, source)

    print_results(base_data, outfile)

if __name__ == '__main__':
    contexts = __make_context_hash()

    args = get_args()
    global quiet
    quiet = args.quiet

    main(args.vcfs, args.source, args.numprocs, quiet, args.reference, 
            args.outfile)

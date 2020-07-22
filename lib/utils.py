# -*- coding: utf-8 -*-
# Misc global functions for the package
import os
import sys
import json
import gzip
import shutil
import inspect
import datetime

from termcolor import cprint
from pprint import pprint


def __exit__(line=None, msg=None, color=None):
    """
    Exit, outputting a colored, formatted, line indicating where we stopped in
    the script with an optional message telling us why we stopped.  Useful for
    debugging and writing new code.
    """
    filename, lineno, function = inspect.stack()[1][1:4]
    outcolor = 'on_green' if color is None else 'on_%s' % color

    line = lineno if line is None else line
    output = ('Script "{}" stopped in `{}()` at line: {} with message: '
            '"{}".'.format(os.path.basename(filename), function, line, msg))
    sys.stderr.write('\n')
    cprint(output, "white", outcolor, attrs=['bold'], file=sys.stderr)
    sys.exit()

def today():
        return datetime.datetime.today().strftime('%Y%m%d')

def pp(data):
    pprint(data, stream=sys.stderr)

def gen_hgvs(var_elems):
    """
    Generate a valid HGVS genome, coding change, and amino acid change term for
    a variant. Requires a dict of variant elems that must include the terms:

        gene : Offical HUGO gene symbol
        chr  : Chromosome; can include the string 'chr' or not.
        pos  : Start position where the variant occurs.
        ref  : Reference base(s). Entry should be '-' for insertions.
        alt  : Alternative base(s). Entry should be '-' for deletions.
        cds  : Coding sequence change (e.g. c.123G>A)
        aa   : Amino acid change. Can be long or short version (e.g. p.V600E 
               or p.Val600Glu).
        transcript : Transcript ID used for the call.

    Return data will be tuple of HGVSg, HGVSc, and HGVSp changes.
    """
    # Based on GRCh37.p13
    chr_to_refseq = {
        'chr1' : 'NC_000001.10', 'chr2' : 'NC_000002.11', 'chr3' : 'NC_000003.11',
        'chr4' : 'NC_000004.11', 'chr5' : 'NC_000005.9', 'chr6' : 'NC_000006.11',
        'chr7' : 'NC_000007.13', 'chr8' : 'NC_000008.10', 'chr9' : 'NC_000009.11',
        'chr10' : 'NC_000010.10', 'chr11' : 'NC_000011.9', 'chr12' : 'NC_000012.11',
        'chr13' : 'NC_000013.10', 'chr14' : 'NC_000014.8', 'chr15' : 'NC_000015.9',
        'chr16' : 'NC_000016.9', 'chr17' : 'NC_000017.10', 'chr18' : 'NC_000018.9',
        'chr19' : 'NC_000019.9', 'chr20' : 'NC_000020.10', 'chr21' : 'NC_000021.8',
        'chr22' : 'NC_000022.10', 'chrX' : 'NC_000023.10', 'chrY' : 'NC_000024.9',
    }

    # Make sure we have the "chr" string.
    var_elems['chr'] = 'chr{}'.format(var_elems['chr'].lstrip('chr'))

    hgvs_g = '{}:g.{}{}>{}'.format(
        chr_to_refseq[var_elems['chr']],
        var_elems['pos'],
        var_elems['ref'],
        var_elems['alt']
    )
    hgvs_c = '{}({}):{}'.format(
        var_elems['transcript'],
        var_elems['gene'],
        var_elems['cds']
    )
    hgvs_p = '{}({}):{}'.format(
        var_elems['transcript'],
        var_elems['gene'],
        var_elems.get('aa', 'p.?')
    )
    return hgvs_g, hgvs_c, hgvs_p

def make_json(*, outfile, data, sort=True):
    with open(outfile, 'w') as fh:
        json.dump(data, fh, sort_keys=sort, indent=4)

def print_json(data):
    print(json.dumps(data, sort_keys=True, indent=4))

def read_json(jfile):
    with open(jfile) as fh:
        return json.load(fh)

def gunzip_vcf(vcf):
    """
    Decompress gzipped and bgzipped (the gzip lib seems to be OK with either)
    VCF files so that we can use them in MOMA.
    """
    gunzipped_file = vcf.replace('.gz', '')
    with gzip.open(vcf, 'rb') as gz_fh:
        with open(gunzipped_file, 'wb') as fh:
            shutil.copyfileobj(gz_fh, fh)
    return gunzipped_file


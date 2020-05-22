#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Simple setup script for pulling in the rest of the data and whatnot needed to
# get the package running. This won't globally work for every scenario.  But, it
# will at least make all of the iterations I need to do a bit easier!
#
# 5/22/2020 - D Sims
################################################################################
"""
Simple setup script for getting MOMA up and running on a new system.  
"""
import os
import sys
import shutil
import argparse

from pprint import pprint as pp # noqa

version = '1.0.052220'

quiet = True

package_root = os.path.dirname(__file__)
lib_dir = os.path.join(package_root, 'lib')
resource_dir = os.path.join(package_root, 'resource')

def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--annovar',
        metavar='<path_to_annovar>',
        required=True,
        help='Location of Annovar program directory to copy into package. '
            'Should contain the 6 perl scripts that make up Annovar at least. '
            'The database files will be imported elsewhere.'
    )
    parser.add_argument(
        '--annovar_db',
        metavar='<path_to_annovar_db_files>',
        required=True,
        help='Location of Annovar database files to copy into package.'
    )
    parser.add_argument(
        '--hg19_ref',
        metavar='<path_to_hg19.fasta>',
        required=True,
        help="Location of the GRCh37 (hg19) reference FASTA file. It is "
            "recommended to use a bgzipped version of this file (or can use "
            "gzip if you manually index it), as the file size is smaller, but "
            "either will work. Also no need to include the index file (unless "
            "using gzip instead of bgzip), as this will be created the first "
            "time it's used."
    )
    parser.add_argument(
        '-q', '--quiet',
        action = 'store_true',
        help='Suppress output messages.'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version="%(prog)s - v" + version
    )
    args = parser.parse_args()
    
    global quiet
    quiet = args.quiet

    return args

def __format_message_block(msg, pos=None):
    if pos == 'first':
        sys.stdout.write('  {}\n    {}\n'.format('-' * 75, msg))
        sys.stdout.flush()
    elif pos == 'last':
        sys.stdout.write('    {}\n  {}\n'.format(msg, '-' * 75))
        sys.stdout.flush()
    else:
        sys.stdout.write('    {}\n'.format(msg))
        sys.stdout.flush()

def set_up_annovar(annovar_source):
    __format_message_block('Setting up local copy of Annovar...', 'first')

    annovar_dest = os.path.join(package_root, 'lib', 'annovar')
    annovar_manifest = ['variants_reduction.pl', 'table_annovar.pl', 
        'annotate_variation.pl', 'coding_change.pl', 'retrieve_seq_from_fasta.pl',
        'convert2annovar.pl']

    try:
        os.mkdir(annovar_dest)
        
    except FileExistsError:
        # Already have an Annovar dir for some reason; make room for new one.
        if quiet is False:
            sys.stdout.write('Already have an Annovar dir for some reason. Making '
                'room for a new one.\n')
        shutil.rmtree(annovar_dest)
        os.mkdir(annovar_dest)

    annovar_files = [os.path.join(os.path.abspath(annovar_source), f) for f in
            os.listdir(annovar_source)]

    if list(map(os.path.basename, annovar_files)) != annovar_manifest:
        sys.stderr.write('ERROR: Expected to find the following Annovar '
            'files, but some or all are missing from import dir:\n')
        pp(annovar_files)
        sys.exit(1)

    for f in annovar_files:
        dest = os.path.join(annovar_dest, os.path.basename(f))
        if quiet is False:
            sys.stdout.write(f'Copying {f} to {dest}...\n')
        shutil.copyfile(f, dest)

    __format_message_block('Done!', 'last')

def set_up_annovar_db(db_source):
    # TODO: Create progress bar hook for this since it takes so long.
    __format_message_block('Setting up local copy of MOMA Annovar Database. '
        'This will take a while...', 'first')

    db_dest_dir = os.path.join(package_root, 'resource', 'annovar_db')

    db_manifest = [
     'hg19_dbnsfp35a.txt',
     'hg19_avsnp142.txt.idx',
     'hg19_gnomad_exome.txt',
     'hg19_popfreq_all_20150413.txt.idx',
     'hg19_trunc_refGeneMrna.fa',
     'hg19_cosmic89_noEnst.txt',
     'hg19_gnomad_exome.txt.idx',
     'hg19_avsnp142.txt',
     'hg19_cytoBand.txt',
     'hg19_knownGene.txt',
     'hg19_clinvar_20190305.txt',
     'hg19_clinvar_20190305.txt.idx',
     'hg19_dbnsfp35a.txt.idx',
     'hg19_popfreq_all_20150413.txt',
     'hg19_trunc_refGene.txt']

    try:
        os.mkdir(db_dest_dir)
        
    except FileExistsError:
        # Already have an Annovar database dir for some reason; make room for new one.
        if quiet is False:
            sys.stdout.write('Already have an Annovar database dir for some '
                'reason. Making room for a new one.\n')
        shutil.rmtree(db_dest_dir)
        os.mkdir(db_dest_dir) 

    db_files = [os.path.join(os.path.abspath(db_source), f) for f in
            os.listdir(db_source)]

    if list(map(os.path.basename, db_files)) != db_manifest:
        sys.stderr.write('ERROR: Expected to find the following database '
            'files, but some or all are missing from import dir:\n')
        pp(db_files)
        sys.exit(1)

    for f in db_files:
        dest = os.path.join(db_dest_dir, os.path.basename(f))
        if quiet is False:
            sys.stdout.write(f'Copying {f} to {dest}...\n')
            sys.stdout.flush()
        shutil.copyfile(f, dest)

    __format_message_block('Done!', 'last')

def set_up_hg19(hg19_ref):
    __format_message_block('Copying hg19 reference to resource dir.', 'first')

    if not hg19_ref.endswith('gz'):
        __format_message_block("Using uncompressed hg19 reference. In future, "
            "try using a bgzip compressed file; it's smaller!")

    dest_dir = os.path.join(package_root, 'resource')

    try:
        shutil.copyfile(os.path.abspath(hg19_ref), 
            os.path.join(dest_dir, os.path.basename(hg19_ref)))
    except FileNotFoundError:
        sys.stderr.write("ERROR: hg19 ref submitted can not be found.\n")
        sys.exit()

    __format_message_block('Done!', 'last')

def main(annovar, annovar_db, hg19_ref):
    sys.stdout.write("<|{0}  Setting up MOMA  {0}|>\n".format('='*30))

    set_up_annovar(annovar)
    set_up_annovar_db(annovar_db)
    set_up_hg19(hg19_ref)

    sys.stdout.write(

    """
    Done!

    All resource files have been installed and are
    ready for use. See user guide for additional program modules and
    libraries required. It is also recommended to run some of the VCF
    files in the 'test' directory to ensure your installation is
    working correctly.
    """
    )
    sys.stdout.write('\n<|{}|>\n'.format('='*75))

if __name__ == '__main__':
    args = get_args()
    main(args.annovar, args.annovar_db, args.hg19_ref)


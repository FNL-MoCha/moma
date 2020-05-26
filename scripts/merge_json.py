#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quickie script to merge the JSON objects from MOMA into one if we want just one 
blob.
"""
import sys
import os
import argparse

lib_path = os.path.join(os.path.dirname(__file__), '..', 'lib')
sys.path.insert(0, lib_path)
import utils


version = '1.0.052620'

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'json_files', 
        nargs='+',
        metavar='<JSON files>',
        help='JSON files to merge'
    )
    parser.add_argument(
        '-o', '--outfile',
        metavar="<output file>",
        default='out.json',
        help="File to which the data should be written. The default is "
            "'%(default)s'"
    )
    parser.add_argument(
        '-v', '--version', 
        action='version',
        version=f"%(prog)s - v{version}"
    )
    return parser.parse_args()

def main(json_files, outfile):
    jdata = {}
    for f in json_files:
        jdata.update(utils.read_json(f))
    sys.stdout.write(f'Writing data to {outfile}.\n')
    utils.make_json(data=jdata, outfile=outfile)

if __name__ == '__main__':
    args = get_args()
    main(args.json_files, args.outfile)

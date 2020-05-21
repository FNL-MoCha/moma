#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that the environment has the proper packages and libraries to run - for 
people that don't check the manual for installation instructions!
"""
import sys
import subprocess
import importlib.util
import shutil

perl_mods = ['Data::Dump', 'Text::CSV', 'Sort::Versions', 'Log::Log4perl']
python_mods = ['pysam'] # Also natsort, but I decided to package it with this 
packages = ['vcftools', 'samtools', 'bedtools']

version = '1.0.052120'

def check_perl_mods(mods):
    for mod in mods:
        cmd = ['perl', f'-M{mod}', '-e', "\'1\'"]
        try:
            subprocess.check_call(cmd, stdout=subprocess.DEVNULL, 
                stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            __throw_error('Perl module', mod)

def check_python_mods(mods):
    for mod in mods:
        spec = importlib.util.find_spec(mod)
        if spec is None:
            __throw_error('Python3 module', mod)

def check_packages(pkgs):
    for p in pkgs:
        if shutil.which(p) is None:
            __throw_error("program", p)

def __throw_error(item, elem):
    sys.stderr.write(f'\nERROR: Can not import {item} "{elem}". '
        'Please install this item before running MOMA.\n')
    sys.exit(1)

def main():
    check_perl_mods(perl_mods)
    check_python_mods(python_mods)
    check_packages(packages)

    #  sys.stdout.write('All modules and packages are installed and available. '
        #  'Continuing to run MOMA.\n')

if __name__ == '__main__':
    main()

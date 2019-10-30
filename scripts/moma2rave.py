#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read in a MOMA Report and output a CSV file that can be read into Theradex Rave
for the NCORP-10231 trial.
"""
import sys
import re
import csv
import argparse

from pprint import pprint as pp # noqa

version = '1.0.103019'

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'input_file',
        metavar="<input_file>",
        help='MOMA MOI report file to parse.'
    )
    parser.add_argument(
        '-p', '--protocol_number',
        metavar='<protocol_number>',
        default='10231',
        type=str,
        help='Protocol number required by Rave. Default: %(default)s.'
    )
    parser.add_argument(
        '-t', '--treatment_id',
        metavar='<treatemnt_patient_ID>',
        required = True,
        type=str,
        help='Treatment Patient ID that is required for registration into RAVE.'
    )
    parser.add_argument(
        '-s', '--specimen_id',
        metavar='<specimen_id>',
        required=True,
        type=str,
        help='Specimen ID required by RAVE.'
    )
    parser.add_argument(
        '-o', '--outfile',
        metavar="<output_file>",
        help="File to which the data should be written.  By default data just "
            "written to stdout."
    )
    parser.add_argument(
        '-v', '--version',
        action = 'version',
        version = '%(prog)s - v' + version
    )
    args = parser.parse_args()

    # Validate that the input protocol number, treatment patient id, and
    # specimen id strings appear to be valid based on pattern match.  
    if args.protocol_number != '10231':
        sys.stderr.write("ERROR: Protocol number should be `10231` only.\n")
        sys.exit(1)
    #  print('')

    #  print('treatment patient id: {} (type: {})'.format(args.treatment_id,
        #  type(args.treatment_id)))
    id_type = 'treatment_patient_id'
    trid_regex = r'^[A-Z]{2}\d{3}-\d{4}'
    validate_id(trid_regex, args.treatment_id, id_type)
    #  print('')

    #  print('specimen id:          {} (type: {})'.format(args.specimen_id,
        #  type(args.specimen_id)))
    id_type = 'specimen_id'
    spid_regex = r'^' + args.protocol_number + '-[0-9A-Z]{8}-\d$'
    validate_id(spid_regex, args.specimen_id, id_type)

    return args

def validate_id(regex, string, id_type):
    """
    Make sure that the ID string that's been passed in follows the correct
    format and is likely the correct string (to also check for data entry
    errors).
    """
    #  print(f"regex: {regex}")
    if not re.fullmatch(regex, string):
        sys.stderr.write(f"ERROR: `{string}` is not a valid {id_type} string.\n")
        sys.exit(1)

def read_input(moi_file):
    """
    Read data into chunks, one for each variant type, and store in a dict for
    printing later on.
    """
    with open(moi_file) as fh:
        lines = fh.read()
    # split data into chunks based on a blank line in between sections. The
    # first section will be SNV data, the second will be CNV data, and the last
    # will be Fusion data.
    snv_data, cnv_data, fusion_data = re.split(r'(?:\r?\n){2,}', lines)

    data = {}
    data['SNV']    = parse_data(snv_data)
    data['CNV']    = parse_data(cnv_data)
    data['Fusion'] = parse_data(fusion_data)

    # If we got `None` for any of these categories, then we didn't have a
    # variant for that type. Remove from the dict.
    for var_type in ('SNV', 'CNV', 'Fusion'):
        if data[var_type] is None:
            del data[var_type]
    return data

def parse_data(data):
    """
    Parse the variant lines into dicts that can be used later for printing.
    """
    parsed_data = []
    lines = data.split('\n')
    # Return None if there is no data for a particular category.
    if len(lines) < 2 or lines[2].startswith('No '):
        return None
    for l in lines[2:]:
        # Have to dump a trailing newline
        if l.split():
            parsed_data.append(line2dict(lines[1], l))
    return parsed_data

def line2dict(header_line, line):
    """
    Read in the header line and the data for each variant type, and return a
    dict with the data keyed by the Theradex header elems suitable for printing
    later.
    """
    theradex_key_map = {
        'Chr' : 'Chromosome', 'Pos' : 'Position', 'Ref' : 'Reference Allele', 
        'Alt' : 'Alternative Allele', 'VAF' : 'Variant Allele Frequency', 
        'Ref_Reads' : 'Reference Coverage', 'Alt_Reads' : 'Alternative Coverage',
        'COSMIC_Id' : 'Variant ID', 'CDS' : 'Coding Sequence', 'AA' : 'Amino Acids',
        'CN' : 'Copies', 'Read_Count' : 'Read Counts'
    }
    header = next(csv.reader(header_line.splitlines(), delimiter=','))
    theradex_header = [theradex_key_map.get(x, x) for x in header]
    return dict(zip(theradex_header, 
        next(csv.reader(line.splitlines(), delimiter=','))))

def print_data(data, pnum, tid, sid, outfile):
    if outfile is None:
        # Write to STDOUT
        outfh = sys.stdout
    else:
        outfh = open(outfile, 'w')

    csv_writer = csv.writer(outfh, delimiter=',')
    header = ['Protocol Number', 'Treatment patient ID', 'Specimen ID',
        'Mutation Type', 'Gene', 'Chromosome', 'Position', 'Reference Allele', 
        'Alternative Allele', 'Transcript', 'Coding Sequence', 'Amino Acids', 
        'Variant ID', 'Variant Allele Frequency', 'Coverage', 
        'Reference Coverage', 'Alternative Coverage', 'Copies', 'Read Counts',
        'Function']
    csv_writer.writerow(header)

    out_data = []

    for var_type in ('SNV', 'CNV', 'Fusion'):
        if var_type in data:
            for var in data[var_type]:
                var_data = {x : var.get(x, '') for x in header}
                var_data.update({'Mutation Type' : var_type})
                
                if var_type == 'SNV':
                    # Need to get total coverage info and function into the data
                    coverage = str(int(var_data['Reference Coverage']) +
                            int(var_data['Alternative Coverage']))
                    var_data.update({'Coverage' : coverage})
                elif var_type == 'CNV':
                    # If CNV, have to add effect
                    var_data.update({'Function' : var['Effect']})

                elif var_type == 'Fusion':
                    # If Fusion, need to get the Fusion ID instead of gene, and
                    # need to add effect.
                    var_data.update({
                        'Gene' : var['Driver_Gene'], 
                        'Variant ID' : '%s.%s' % (var['Fusion'], var['Junction']),
                        'Function' : var['Effect']
                    })
                var_data.update({
                    'Protocol Number' : pnum,
                    'Treatment patient ID' : tid,
                    'Specimen ID' : sid,
                })
                out_data.append(var_data)
        else:
            # No data for this variant type. Pad out the section.
            var_data = {x :  '' for x in header}
            var_data.update({
                'Mutation Type' : var_type,
                'Protocol Number' : pnum,
                'Treatment patient ID' : tid,
                'Specimen ID' : sid,
            })
            out_data.append(var_data)

    for var in out_data:
        csv_writer.writerow([var.get(x, '') for x in header])

def main(protocol_number, treatment_id, specimen_id, input_file, outfile):
    data = read_input(input_file)
    print_data(data, protocol_number, treatment_id, specimen_id, outfile)

if __name__ == '__main__':
    args = get_args()
    main(args.protocol_number, args.treatment_id, args.specimen_id,
            args.input_file, args.outfile)

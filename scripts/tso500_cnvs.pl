#!/usr/bin/env perl
# -*- coding: utf-8 -*-
# Simple, quick script to pull out CNV data from a TSO500 assay, and get the
# data ready for integration into the rest of the MOMA dataset.
#
# 2/19/2020 - D Sims
################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long;
use File::Basename;
use Data::Dump;

my $version = "v1.0.021920";
my $scriptdir = dirname($0);

# Gene reference file.
my $gene_ref_file = "$scriptdir/../resource/gene_reference.csv";
die "ERROR: Can not locate gene reference file '$gene_ref_file'! ",
        "Check that this file is in your package.\n" unless -e $gene_ref_file;
my $gene_data = __read_ref($gene_ref_file);

my $scriptname = basename($0);
my $description = <<"EOT";
Read in the TSO CNV datafile and generate a dataset that can be incorporated 
into the MOMA data.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <TSO500 *.cnv.fc.txt file>
    Program Options:
        -v, --version      Version information
        -h, --help         Print this help information
EOT

my $help;
my $ver_info;

GetOptions( 
    "version|v"     => \$ver_info,
    "help|h"        => \$help,
) or die $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version if $ver_info;
 
# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "ERROR: No CNV data file passed to script!\n";
    print "$usage\n";
    exit 1;
}

################------ END Arg Parsing and Script Setup ------#################
my $input_file = shift;
my @cnv_data;

open(my $fh, "<", $input_file);
chomp(my $header_line = readline($fh));
my @header = split(/\t/, $header_line);

while(<$fh>) {
    my %tmp;
    chomp(my @elems = split(/\t/));
    @tmp{@header} = @elems;

    my $gene = $tmp{'Gene'};
    # if the gene is in the file, then it's an OncoKB curated gene; add the rest
    # of the info.  Otherwise, not an OncoKB annotated gene, and we should
    # remove from this analysis.
    if (exists $gene_data->{$gene}) {
        $tmp{'Chr'}   = $gene_data->{$gene}{'Chromosome'};
        $tmp{'Start'} = $gene_data->{$gene}{'Start'};
        $tmp{'End'}   = $gene_data->{$gene}{'End'};
    } else {
        next;
    }
    push(@cnv_data, \%tmp);
}
close $fh;

my @wanted = qw(Sample Chr Start End Gene FoldChange Status);
print join(',', @wanted), "\n";

for my $cnv (@cnv_data) {
    print join(',', @{$cnv}{@wanted}), "\n";
}

sub __read_ref {
    my $ref_file = shift;
    my %data;

    open(my $fh, "<", $ref_file);
    my $vline = readline($fh);
    my $hline = readline($fh);
    chomp(my @header = split(/,/, $hline));

    while(<$fh>) {
        my %tmp;
        chomp(my @elems = split(/,/));
        @tmp{@header} = @elems;
        $data{$tmp{'Hugo_Symbol'}} = \%tmp;
    }
    return \%data;
}

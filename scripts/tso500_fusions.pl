#!/usr/bin/env perl
# -*- coding: utf-8 -*-
# Pull out fusions from TSO500 data for MOMA integration.
#
# 2/20/20 - D Sims
################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long;
use File::Basename;
use Data::Dump;
#use Sort::Versions;

# TODO: Do we need this?
use constant TRUE  => 1;
use constant FALSE => 0;
use constant DEBUG => 0;

my $scriptname = basename($0);
my $version = "1.0.022020";

my $vaf = 0.01; 
my $outfile;
my $help;
my $ver_info;

my $description = <<"EOT";
Read in the TSO500 Fusions file and generate a dataset that can be incorporated
into a MOMA report.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <TSO500 *.fusions.txt file>

Options
    -V, --VAF         Only report fusions above this threshold (DEFAULT: 
                      $vaf).
    -o, --output      Write output to file.
    -v, --version     Display version information.
    -h, --help        Display this help text.
EOT

GetOptions( 
    "VAF|V=f"       => \$vaf,
    "outfilel|o=s"  => \$outfile,
    "help|h"        => \$help,
    "version|v"     => \$ver_info,
);

sub version { print "$scriptname - version: $version\n\n" }
sub help { version(); print "$description\n$usage\n\n" }

help() and exit if $help;
version() and exit if $ver_info;

# Check we have some files to process
if ( @ARGV < 1 ) {
    print "ERROR: You did not load the TSO500 fusions file.\n";
    print $usage;
    exit 1;
}

# Set up output
my $out_fh;
if ( $outfile ) {
    open( $out_fh, ">", $outfile );
} else {
    $out_fh = \*STDOUT;
}

########======================  END ARG Parsing  ======================########
my $input_file = shift;
# my @results;

my $return_data = proc_fusions_file(\$input_file);

dd $return_data;
exit;

=cut
my @header = qw (Fusion Junction ID Read_Count Driver_Gene Partner_Gene);
print {$out_fh} join(',', @header), "\n"; 

# Generate and print out the final results table(s)
for my $fusion ( sort{ versioncmp($a, $b) } keys %$return_data ) {
    # Some basic filters.
    next if ($gene 
        and ! grep { $return_data->{$fusion}->{'DRIVER'} eq $_ } @genes_list);
    next if (grep { $fusion eq $_ } qw(Non-targeted Novel) and $novel_filter);
    next if $return_data->{$fusion}->{'COUNT'} < $reads;

    my ($name, $junct, $id ) = split(/\|/, $fusion);

    print {$out_fh} join(',', $name, $junct, $id, 
        @{$return_data->{$fusion}}{qw(COUNT DRIVER PARTNER)}), "\n";
}
=cut


sub proc_fusions_file {
    my $input_file = shift;

    my @drivers = qw(ABL1 ALK BCR BRAF CD74 EGFR ETV1 ETV4 ETV6 EWSR1 FGFR2 FGFR3
        NAB2 NTRK1 NTRK2 NUTM1 PAX3 PAX8 PPARG RET ROS1 TFE3 TMPRSS2); 

    my $sample_name;
    my %results;

    open( my $fh, "<", $$input_file );
    chomp(my $hline = readline($fh));
    my @header = split(/\t/, $hline);

    while (<$fh>) {
        my %tmp_data;

        chomp(my @elems = split(/\t/));
        @tmp_data{@header} = @elems;

        # Create a junction field.
        $tmp_data{'junction'} = sprintf('%s:%s::%s:%s', 
            @tmp_data{qw(Chr1 Pos1 Chr2 Pos2)});
        dd \%tmp_data;

=cut
        next if /^#/;
        my @data = split;

        # Get rid of FAIL and NOCALL calls to be more compatible with MATCHBox 
        # output.
        next if ($data[6] eq 'FAIL' || $data[6] eq 'NOCALL');

        if ( $data[7] =~ /SVTYPE=Fusion/ ) {
            my ( $name, $elem ) = $data[2] =~ /(.*?)_([12])$/;
            my ($count) = map { /READ_COUNT=(\d+)/ } @data;
            next if ( $count == 0 );
            my ($pair, $junct, $id) = split(/\./, $name);
            $id //= '-';

            # if ($id eq 'Non-Targeted' || $id eq 'Novel') {
                # next unless $novel;
            # }
            #
            my ($gene1, $gene2) = split(/-/, $pair);
            my $fid = join('|', $pair, $junct, $id);

            if ( $pair eq 'MET-MET' || $pair eq 'EGFR-EGFR' ) {
                $results{$fid}->{'DRIVER'} = $results{$fid}->{'PARTNER'} = $gene1;
            }
            elsif (grep {$_ eq $gene1} @drivers) {
                $results{$fid}->{'DRIVER'} = $gene1;
                $results{$fid}->{'PARTNER'} = $gene2;
            }
            elsif (grep {$_ eq $gene2} @drivers) {
                $results{$fid}->{'DRIVER'} = $gene2;
                $results{$fid}->{'PARTNER'} = $gene1;
            }
            else {
                $results{$fid}->{'DRIVER'} = 'UNKNOWN';
                $results{$fid}->{'PARTNER'} = "$gene1,$gene2";
            }
            $results{$fid}->{'COUNT'} = $count;
        }
=cut
    }
    return \%results;
}

#!/usr/bin/perl
# Light (very!) version of vcfExtractor just for the AMG-232 plugin. This will 
# parse out the NOCALLs and reference calls from the VCF, and explode the 
# multiple calls per line into 1 call per line.
################################################################################
use strict;
use warnings;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use List::Util qw( sum min max );
use Sort::Versions;
use Data::Dump;
use File::Basename;
use Time::Piece;

my $scriptname = basename($0);
my $version = "v1.0.091218";

my $description = <<"EOT";
Utility to read an Ion Torrent VCF file and remove reference calls, and NOCALLs, 
as well as, expand the data so that we have just 1 entry per line.  This will 
help create a VCF that can then more readily be passed into VEP or Annovar.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <input_vcf_file>
    -f, --filename  Name to use for the resultant output file. 
                    Default: <file>_flattened.vcf
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $debug_pos;  # undocumented.
my $out_filename;

GetOptions( 
    "filename|f=s"  => \$out_filename,
    "debug_pos=s"   => \$debug_pos, #undocumented.
    "version|v"     => \$ver_info,
    "help|h"        => \$help )
or die "\n$usage";

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version_info {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version_info if $ver_info;

# Check for vcftools; we can't run without it...for now.
if ( ! qx(which vcftools) ) {
    die "ERROR: Required package 'vcftools' is not installed on this system. ",
        "Install vcftools ('vcftools.sourceforge.net') and try again.\n";
}

# Implementing a debug position method to help with development.  This is an
# undocumented method that will allow for one to input a position and output
# the parsed hash of data and rest of method on only that position alone.
if ($debug_pos) {
    print '-'x25 . '  DEBUG  ' . '-'x25, "\n";
    print "Outputting data for position $debug_pos only.\n";
    print '-'x59, "\n";
}

# Make sure VCF has been passed to script.
if ( scalar( @ARGV ) < 1 ) {
    print "ERROR: No VCF file passed to script!\n\n";
    print $usage;
    exit 1;
}
my $input_vcf = shift;

# Set up output type
my $outfile;
($out_filename)
    ? ($outfile = $out_filename) 
    : (($outfile = $input_vcf) =~ s/\.vcf/_flattened.vcf/);

open(my $out_fh, ">", $outfile);

#########--------------- END Arg Parsing and validation ---------------#########

# Check VCF file and options to make sure they're valid
open ( my $vcf_fh, "<", $input_vcf);
my @header = grep { /^#/ } <$vcf_fh>;
die "ERROR: '$input_vcf' does not appear to be a valid VCF file or does not ",
    "have a header.\n" unless @header;
close $vcf_fh;

# Get the data from VCF Tools
my @wanted_fields = qw(%CHROM:%POS %REF %ALT %FILTER %INFO/FR %INFO/OID 
    %INFO/OPOS %INFO/OREF %INFO/OALT %INFO/OMAPALT --- --- [%GTR %AF %FRO %RO 
    %FAO %AO %DP]);

my $vcf_format = join('\t', @wanted_fields);
my @extracted_data = qx( vcf-query $input_vcf -f "$vcf_format\n" );

# Read in the VCF file data and create a hash
my $vcf_data = parse_data(\@extracted_data);

# Write the data out to a new VCF file.
make_vcf(\@header, $vcf_data, $out_fh);

sub parse_data {
    # Extract the VCF information and create a hash of the data.  
    my $data = shift;
    my %parsed_data;

    for ( @$data ) {
        chomp;
        my ( $pos, $ref, $alt, $filter, $reason, $oid, $opos, $oref, $oalt, 
            $omapalt, $func, $lod, $gtr, $af, $fro, $ro, $fao, $ao, 
            $dp ) = split( /\t/ );

        # Limit processing to just one position and output more metrics so
        # that we can figure out what's going on.
        if ($debug_pos) {
            next unless $pos eq $debug_pos;
            print_debug_output([split(/\t/)]);
        }

        # Clean up filter reason string
        $reason =~ s/^\.,//;

        # Filter out vars we don't want to print out later anyway.
        next if $reason eq "NODATA";
        
        # Sometimes not getting correct 'NOCALL' flag; manually set if need to.
        $filter = "NOCALL" if (($gtr =~ m|\./\.|) or ($reason eq 'REJECTION'));
        # Get rid of reference and NOCALLS.
        next if ($filter eq "NOCALL" or $gtr eq '0/0');
        # next if ( $nocall && $filter eq "NOCALL" );
        # next if ( $noref && $gtr eq '0/0' );

        # Create some arrays to hold the variant data in case we have MNP calls
        my @alt_array     = split( /,/, $alt );
        my @oid_array     = split( /,/, $oid );
        my @opos_array    = split( /,/, $opos );
        my @oref_array    = split( /,/, $oref );
        my @oalt_array    = split( /,/, $oalt );
        my @omapalt_array = split( /,/, $omapalt );
        my @fao_array     = split( /,/, $fao );
        my @ao_array      = split( /,/, $ao );

        for my $alt_index ( 0..$#alt_array ) {

            # NOTE: New vars for flagging.
            my $caller = 'tvc'; # set default as TVC and change later.

            my $alt_var = $alt_array[$alt_index];

            # Get the normalizedRef, normalizedAlt, and normalizedPos values
            # from the REF and ALT fields so that we can map the FUNC block.
            my @coords = split(/:/, $pos);
            my %norm_data = normalize_variant(\$ref, \$alt_var, $coords[1]);

            my @array_pos = grep { $omapalt_array[$_] eq $alt_var } 0..$#omapalt_array;
            for my $index ( @array_pos ) {
                my @var_data; # Temp holder array for collected variant data.
                (my $parsed_pos = $pos) =~ s/(chr\d+:).*/$1$norm_data{'normalizedPos'}/; 
                
                my $var_id = join( ":", $parsed_pos, $oref_array[$index], 
                    $oalt_array[$index] );
                my $cosid = $oid_array[$index];

                # Stupid bug with de novo and hotspot merge that can create two 
                # duplicate entries for the same variant but one with and one 
                # without a HS (also different VAF, coverage,etc). Try this to 
                # capture only HS entry.
                if ($parsed_data{$var_id}) {
                    # if data already there, check if it's from LIA and replace
                    if ($parsed_data{$var_id}->[-1] eq 'lia') {
                        delete $parsed_data{$var_id};
                    }
                    # Pedmatch has duplicate hotspots, both of which are 
                    # populating the output. If we already have a variant 
                    # entry and it's the PM version replace with COSMIC 
                    # version.
                    elsif ($parsed_data{$var_id}->[8] =~ /^PM/) {
                        delete $parsed_data{$var_id};
                    }

                    # If TVC Duplicate Hotspot bug, 
                    delete $parsed_data{$var_id} if ( $cosid eq '.');
                } 

                # Start the var string.
                push( @var_data,
                    $parsed_pos,
                    $norm_data{'normalizedRef'},
                    $norm_data{'normalizedAlt'},
                    $filter, 
                    $reason
                );

                # Check to see if call is result of long indel assembler and 
                # handle appropriately. 
                my ($vaf, $tot_coverage);
                if ( $fao_array[$alt_index] eq '.' ) {
                    # Result is from LIA
                    $caller = 'lia';
                    $tot_coverage = $ao_array[$alt_index] + $ro;
                    $vaf = vaf_calc( \$filter, \$dp, \$ro, 
                        \$ao_array[$alt_index] );
                    
                    push(@var_data, $vaf, $tot_coverage, $ro, 
                        $ao_array[$alt_index], $cosid);
                } else {
                    my @cleaned_fao_array = grep { $_ ne '.' } @fao_array;
                    $tot_coverage = sum( @cleaned_fao_array ) + $fro;
                    $vaf = vaf_calc( \$filter, \$tot_coverage, \$fro, 
                        \$fao_array[$alt_index] );
                    
                    push(@var_data, $vaf, $tot_coverage, $fro, 
                        $fao_array[$alt_index], $cosid );
                }
                
                # Filter out reference calls if we have turned on the noref 
                # filter. Have to leave the NOCALL calls if we have left those 
                # in, and have to deal with sub 1% VAFs for cfDNA assay.
                my $calc_vaf = $var_data[5];
                if ( $calc_vaf ne '.' ) {
                    next if $calc_vaf == 0;
                }
                # Load data into hash.
                $parsed_data{$var_id} = \@var_data;
            }
        }
    }
    
    # dd \%parsed_data;
    # exit;
    return \%parsed_data;
}

sub normalize_variant {
    my ($ref,$alt,$pos) = @_;
    my ($norm_ref, $norm_alt);

    my ($rev_ref, $rev_alt, $position_delta) = rev_and_trim($ref, $alt);
    ($norm_ref, $norm_alt, $position_delta) = rev_and_trim(\$rev_ref, 
        \$rev_alt);

    my $adj_position = $position_delta + $pos;
    return ( 'normalizedRef' => $norm_ref, 'normalizedAlt' => $norm_alt, 
        'normalizedPos' => $adj_position );
}

sub rev_and_trim {
    my ($ref, $alt) = @_;
    my $position_delta = 0;

    my @rev_ref = split(//, reverse($$ref));
    my @rev_alt = split(//, reverse($$alt));

    while (@rev_ref > 1 && @rev_alt > 1 && $rev_ref[0] eq $rev_alt[0]) {
        shift @rev_ref;
        shift @rev_alt;
        $position_delta++;
    }
    return (join('',@rev_ref), join('', @rev_alt), $position_delta);
}

sub vaf_calc {
    # Determine the VAF
    my ($filter, $tcov, $rcov, $acov) = @_;
    my $vaf;

    if ( $$filter eq "NOCALL" ) { 
        $vaf = '.';
    }
    elsif( $$filter eq "NODATA" || $$tcov == 0) {
        $vaf = 0;
    }
    else {
        $vaf = sprintf( "%.2f", 100*($$acov / $$tcov) );
    }
    return $vaf;
}

sub make_vcf {
    # Output data in VCF format instead of a pretty table.
    my ($header, $data, $out_fh) = @_;
    my $final_header = __make_header($header);

    my @vcf_lines;
    for (sort { versioncmp( $a, $b ) } keys %$data) {
        push(@vcf_lines, __make_vcf_line($$data{$_}));
    }

    print {$out_fh} $_ for @$final_header;
    print {$out_fh} $_ for @vcf_lines;
}

sub print_debug_output {
    # DEBUG: Can add position to filter and output a hash of data to help.
    my $data = shift;
    my @fields = qw(pos ref alt filter reason oid opos oref oalt omapalt func 
        lod gtr af fro ro fao ao dp);
    my %foo;

    @foo{@fields} = map{chomp; $_} @$data;

    print '='x25, "  DEBUG  ", "="x25, "\n";
    dd \%foo;
    print '='x59, "\n";
}

sub __make_header {
    my $header = shift;
    my @captured;
    no warnings;
    
    my @wanted = qw(##fileformat ##fileDate ##source ##reference ##sampleGender
        ##sampleDisease ##contig ##INFO=<ID=AF ##INFO=<ID=AO ##INFO=<ID=DP 
        ##INFO=<ID=RO ##FORMAT=<ID=GT ##FORMAT=<ID=DP );

    for my $line (@$header) {
        push(@captured, $line) if grep { $line =~ /$_/ } @wanted;
    }
    my $new_field = qq(##FORMAT=<ID=AD,Number=G,Type=Integer,Description="Allelic Depths of REF and ALT(s) in the order listed">\n);
    push(@captured, $new_field);
    push(@captured, $header->[-1]);

    # Update the date to today's
    my $today = localtime->strftime('%Y-%m-%d');
    ($captured[1] =  $captured[1]) =~ s/\d{8}/$today/;

    # Change AF field to VAF to help with MAF conflicts.
    s/ID=AF/ID=VAF/ for @captured;
    return \@captured;
}

sub __make_vcf_line {
    my $data = shift;
    my ($chr, $pos) = split(/:/, $data->[0]);
    my $vaf = $data->[5];
    my $gt;
    if ($vaf eq '.') {
        $gt = './.';
    }
    elsif ($vaf > 50) {
        $gt = '1/1';
    }
    else {
        $gt = '0/1';
    }

    my $info = sprintf('VAF=%s;DP=%s;RO=%s;AO=%s', @$data[5..8]);
    my $sample_data = join(':', $gt, "$data->[7],$data->[8]", $data->[6]);

    my $vcf_line = join("\t", $chr, $pos, $data->[9], $data->[1], $data->[2],
        '.', $data->[3], $info, 'GT:AD:DP', $sample_data) . "\n"; 
    return $vcf_line;
}

#!/usr/bin/env perl
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;
use Sort::Versions;

use constant TRUE  => 1;
use constant FALSE => 0;
use constant DEBUG => 0;

my $scriptname = basename($0);
my $version = "1.0.070319-dev";

my $novel_filter = FALSE;
my $gene;
my $reads = 100;
my $outfile;
my $help;
my $ver_info;

my $description = <<"EOT";
<<  Description  >>
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <vcf_file(s)>
    -n, --novel       Exclude 'Non-Targeted' fusions in the output (DEFAULT: 
                      Off).
    -r, --reads       Only report fusions above this threshold (DEFAULT: 
                      $reads).
    -g, --gene        Only output data for a specific driver gene or genes
                      separated by a comma. 

    -o, --output      Write output to file.
    -v, --version     Display version information.
    -h, --help        Display this help text.
EOT

GetOptions( 
    "novel|n"       => \$novel_filter,
    "reads|r=i"     => \$reads,
    "gene|g=s"      => \$gene,
    "outfilel|o=s"  => \$outfile,
    "help|h"        => \$help,
    "version|v"     => \$ver_info,
);

sub help { 
    printf "%s - %s\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
    exit;
}

sub version_info {
    printf "%s - %s\n", $scriptname, $version;
    exit;
}

help() if $help;
version_info() if $ver_info;

# Check we have some files to process
if ( @ARGV < 1 ) {
    print "ERROR: You must load at least one VCF file\n";
    print $usage;
    exit 1;
}

# Set up output
my $out_fh;
if ( $outfile ) {
    open( $out_fh, ">", $outfile ) 
        || die "Can't open the output file '$outfile' for writing: $!";
} else {
    $out_fh = \*STDOUT;
}

my @genes_list = map{uc($_)} split(/,/, $gene) if $gene;

#######======================  END ARG Parsing  ======================#######
my $vcf = shift;
# my @results;

my $return_data = proc_vcf(\$vcf);

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

sub proc_vcf {
    my $vcf = shift;
    # Version 1,2, and 3 drivers. Not all exist in the current version, but 
    # keep all for backward compatibility.
    my @drivers = qw(ABL1 AKT2 AKT3 ALK AR AXL BRAF BRCA1 BRCA2 CDKN2A EGFR 
        ERBB2 ERBB4 ERG ESR1 ETV1 ETV1a ETV1b ETV4 ETV4a ETV5 ETV5a ETV5d FGFR1
        FGFR2 FGFR3 FGR FLT3 JAK2 KRAS MDM4 MET MYB MYBL1 NF1 NOTCH1 NOTCH4 NRG1
        NTRK1 NTRK2 NTRK3 NUTM1 PDGFRA PDGFRB PIK3CA PPARG PRKACA PRKACB PTEN 
        RAD51B RAF1 RB1 RELA RET ROS1 RSPO2 RSPO3 TERT
    );

    my $sample_name;
    my %results;

    open( my $fh, "<", $$vcf );
    while (<$fh>) {
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
    }
    return \%results;
}

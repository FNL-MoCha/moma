#!/usr/bin/env perl
# Script to make a temporary simple VCF file for MOMA processing. It is
# necessary to modify the data in order to get a VAF, coverage, etc. in a
# format that MOMA can easily handle.
#
# 7/21/2020 - D Sims
################################################################################
use strict;
use warnings;
use autodie;

use Data::Dump;
use Getopt::Long;
use File::Basename;

my $version = '1.0.072120';
my $scriptname = basename($0);
my $description = <<"EOT";
Since the Strelka2 VCF files have both a germline and tumor sample in order to 
determine if a variant is somatic or germline, we need to parse out just the
tumor data so that MOMA can interpret the fields correctly. We also need to 
compute VAF for each entry as this is not easily obtained from the VCF entry.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF>

Options
    -o, --outfile    Write the results to a file rather than STDOUT.
    -v, --version    Print the version information and exit.
    -h, --help       Print this help text and exit.
EOT

my $help;
my $ver_info;
my $outfile;

GetOptions(
    'o|outfile=s'  => \$outfile,
    'h|help'       => \$help,
    'v|version'    => \$ver_info,
) or die $usage;


my $ver_string = "$scriptname - v$version\n";
if ($ver_info) {
    print $ver_string;
    exit;
} 
elsif ($help) {
    print "$ver_string\n";
    print "$usage\n";
    exit;
}

my $vcf = shift;
die "ERROR: You must input a VCF file file!\n" unless $vcf;

my ($sample_name, $extracted_data) = read_vcf(\$vcf);

my %filter_counts;
make_new_vcf($sample_name, $extracted_data, $outfile, \%filter_counts);

# Print the filter counts for record keeping purposes.
print "\nSample filtering metrics for $sample_name\n";
my $total = 0;
for (sort{ $a cmp $b } keys %filter_counts) {
    printf("\t%20s %7s\n", $_, $filter_counts{$_});
    $total += $filter_counts{$_};
}
printf("\t%s\n\t%20s %7s*\n", "-" x 28, "Total", $total);
print <<"EOT";

  *note: Number can be higher than total variants since multiple filter terms
         can be applied to one variant.
EOT

sub read_vcf {
    my $vcf = shift;
    my @vcf_lines;
    my $sample_name;

    open(my $fh, "<", $$vcf);
    while (<$fh>) {
        # Get the sample name, but dump the rest of the header.
        if (/^#/) {
            if (/^##cmdline/) {
                ($sample_name) = /--tumorBam (.*?)~WES.*/;
                $sample_name = basename($sample_name);
            }
            next;
        }
        
        # Don't keep mitochondrial variants.
        next unless /^chr[0-9XY]+\b/;
        chomp(my @elems = split(/\t/));
        
        my %format;
        @format{split(/:/, $elems[8])} = (split(/:/, $elems[10]));

        # Not really a format element, but needed for VAF calc.
        $format{'REF'} = $elems[3];
        $format{'ALT'} = $elems[4];

        my ($vaf, $alt_counts, $ref_counts);
        # if (length($elems[3]) > 1 || length($elems[4]) > 1) {
        if (length($elems[3]) != length($elems[4])) {
            # We have an indel. Previously tested if either ref or alt had more
            # than 1 bp, but with MNV, this is not feasible. Not sure how to
            # handle block subs, but let's keep it simple for now.
            ($vaf, $alt_counts, $ref_counts) = compute_vaf(\%format, 'indel');
        } else {
            ($vaf, $alt_counts, $ref_counts) = compute_vaf(\%format, 'snv');
        }
        next if $vaf eq 'None';
        # Very (!) crude genotype algorithm just to satisfy Annovar (fails
        # parsing without this, but not sure why).
        my $genotype;
        if ($vaf > 0.5) {
            $genotype = '1/1';
        }
        elsif ($vaf == 0) {
            $genotype = '0/0';
        } else {
            $genotype = '0/1';
        }
        my @var_data = (@elems[0..7], 'GT:AD:AF', 
            "$genotype:$ref_counts,$alt_counts:$vaf");
        push(@vcf_lines, \@var_data);
    }

    return $sample_name, \@vcf_lines;
}

sub compute_vaf {
    # To compute VAF, need to know whether we have an indel or snv, and then
    # need to pull the correct information out of the FORMAT field depending on
    # the variant type. See Strelka2 User Guide for "Somatic variant allele
    # Frequencies" for details.
    my ($format_field, $var_type) = @_;
    my ($vaf, $t1_ref, $t1_alt);

    local $SIG{__WARN__} = sub {
        my $msg = shift;
        print "Expected '$var_type' but got\n";
        dd $format_field;
        exit;
    };

    if ($var_type eq 'indel') {
        $t1_alt = (split(/,/, $format_field->{'TIR'}))[0];
        $t1_ref = (split(/,/, $format_field->{'TAR'}))[0];
    }
    elsif ($var_type eq 'snv') {
        my $first_ref = (split('', $format_field->{'REF'}))[0];
        my $first_alt = (split('', $format_field->{'ALT'}))[0];
        # $t1_alt = (split(/,/, $format_field->{"$format_field->{'ALT'}U"}))[0];
        # $t1_ref = (split(/,/, $format_field->{"$format_field->{'REF'}U"}))[0];
        $t1_alt = (split(/,/, $format_field->{"${first_ref}U"}))[0];
        $t1_ref = (split(/,/, $format_field->{"${first_alt}U"}))[0];
    }

    # We can sometimes get issues where there is a variant call in normal but
    # not tumor for some reason.  Catch and skip these entries by returning
    # "None" and then skipping in the calling loop.
    my $total_counts = $t1_ref + $t1_alt;
    return 'None' if $total_counts == 0;

    $vaf = sprintf("%0.4f", $t1_alt / $total_counts);
    return($vaf, $t1_alt, $t1_ref);
}

sub make_new_vcf {
    my ($sample_id, $data, $outfile, $filter_counts) = @_;

    my $outfh;
    if ($outfile) {
        print "Writing data to '$outfile'.\n";
        open($outfh, ">", $outfile);
    } else {
        $outfh = \*STDOUT;
    }

    # Print the new header to file.
    print {$outfh} $_ while <DATA>;

    no warnings; # needed to print the header with a '#' char in it. 
    my @header = qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
    use warnings;
    print {$outfh} join("\t", @header, $sample_id), "\n";

    # Now print out the rest of the variant data.
    for my $var (@$data) {
        if ($var->[6] =~ /PASS/) {
            print {$outfh} join("\t", @$var), "\n";
        } 
        my @filter_terms = split(/;/, $var->[6]);
        $filter_counts->{$_}++ for @filter_terms;
    }
}

# Simple VCF header just to (mostly) pass validation.
__DATA__
##fileformat=VCFv4.2
##source=MOMA
##reference=file:///home/dave/Dropbox/reference/ucsc_hg19/ucsc.hg19.fasta
##contig=<ID=chrM,length=16571>
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
##contig=<ID=chr3,length=198022430>
##contig=<ID=chr4,length=191154276>
##contig=<ID=chr5,length=180915260>
##contig=<ID=chr6,length=171115067>
##contig=<ID=chr7,length=159138663>
##contig=<ID=chr8,length=146364022>
##contig=<ID=chr9,length=141213431>
##contig=<ID=chr10,length=135534747>
##contig=<ID=chr11,length=135006516>
##contig=<ID=chr12,length=133851895>
##contig=<ID=chr13,length=115169878>
##contig=<ID=chr14,length=107349540>
##contig=<ID=chr15,length=102531392>
##contig=<ID=chr16,length=90354753>
##contig=<ID=chr17,length=81195210>
##contig=<ID=chr18,length=78077248>
##contig=<ID=chr19,length=59128983>
##contig=<ID=chr20,length=63025520>
##contig=<ID=chr21,length=48129895>
##contig=<ID=chr22,length=51304566>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chrY,length=59373566>
##contig=<ID=chr1_gl000191_random,length=106433>
##contig=<ID=chr1_gl000192_random,length=547496>
##contig=<ID=chr4_ctg9_hap1,length=590426>
##contig=<ID=chr4_gl000193_random,length=189789>
##contig=<ID=chr4_gl000194_random,length=191469>
##contig=<ID=chr6_apd_hap1,length=4622290>
##contig=<ID=chr6_cox_hap2,length=4795371>
##contig=<ID=chr6_dbb_hap3,length=4610396>
##contig=<ID=chr6_mann_hap4,length=4683263>
##contig=<ID=chr6_mcf_hap5,length=4833398>
##contig=<ID=chr6_qbl_hap6,length=4611984>
##contig=<ID=chr6_ssto_hap7,length=4928567>
##contig=<ID=chr7_gl000195_random,length=182896>
##contig=<ID=chr8_gl000196_random,length=38914>
##contig=<ID=chr8_gl000197_random,length=37175>
##contig=<ID=chr9_gl000198_random,length=90085>
##contig=<ID=chr9_gl000199_random,length=169874>
##contig=<ID=chr9_gl000200_random,length=187035>
##contig=<ID=chr9_gl000201_random,length=36148>
##contig=<ID=chr11_gl000202_random,length=40103>
##contig=<ID=chr17_ctg5_hap1,length=1680828>
##contig=<ID=chr17_gl000203_random,length=37498>
##contig=<ID=chr17_gl000204_random,length=81310>
##contig=<ID=chr17_gl000205_random,length=174588>
##contig=<ID=chr17_gl000206_random,length=41001>
##contig=<ID=chr18_gl000207_random,length=4262>
##contig=<ID=chr19_gl000208_random,length=92689>
##contig=<ID=chr19_gl000209_random,length=159169>
##contig=<ID=chr21_gl000210_random,length=27682>
##contig=<ID=chrUn_gl000211,length=166566>
##contig=<ID=chrUn_gl000212,length=186858>
##contig=<ID=chrUn_gl000213,length=164239>
##contig=<ID=chrUn_gl000214,length=137718>
##contig=<ID=chrUn_gl000215,length=172545>
##contig=<ID=chrUn_gl000216,length=172294>
##contig=<ID=chrUn_gl000217,length=172149>
##contig=<ID=chrUn_gl000218,length=161147>
##contig=<ID=chrUn_gl000219,length=179198>
##contig=<ID=chrUn_gl000220,length=161802>
##contig=<ID=chrUn_gl000221,length=155397>
##contig=<ID=chrUn_gl000222,length=186861>
##contig=<ID=chrUn_gl000223,length=180455>
##contig=<ID=chrUn_gl000224,length=179693>
##contig=<ID=chrUn_gl000225,length=211173>
##contig=<ID=chrUn_gl000226,length=15008>
##contig=<ID=chrUn_gl000227,length=128374>
##contig=<ID=chrUn_gl000228,length=129120>
##contig=<ID=chrUn_gl000229,length=19913>
##contig=<ID=chrUn_gl000230,length=43691>
##contig=<ID=chrUn_gl000231,length=27386>
##contig=<ID=chrUn_gl000232,length=40652>
##contig=<ID=chrUn_gl000233,length=45941>
##contig=<ID=chrUn_gl000234,length=40531>
##contig=<ID=chrUn_gl000235,length=34474>
##contig=<ID=chrUn_gl000236,length=41934>
##contig=<ID=chrUn_gl000237,length=45867>
##contig=<ID=chrUn_gl000238,length=39939>
##contig=<ID=chrUn_gl000239,length=33824>
##contig=<ID=chrUn_gl000240,length=41933>
##contig=<ID=chrUn_gl000241,length=42152>
##contig=<ID=chrUn_gl000242,length=43523>
##contig=<ID=chrUn_gl000243,length=43341>
##contig=<ID=chrUn_gl000244,length=39929>
##contig=<ID=chrUn_gl000245,length=36651>
##contig=<ID=chrUn_gl000246,length=38154>
##contig=<ID=chrUn_gl000247,length=36422>
##contig=<ID=chrUn_gl000248,length=39786>
##contig=<ID=chrUn_gl000249,length=38502>
##content=strelka somatic snv calls
##priorSomaticSnvRate=0.0001
##INFO=<ID=QSS,Number=1,Type=Integer,Description="Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal">
##INFO=<ID=TQSS,Number=1,Type=Integer,Description="Data tier used to compute QSS">
###INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description="Data tier used to compute QSS_NT">
##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref read-position in the tumor">
##INFO=<ID=SNVSB,Number=1,Type=Float,Description="Somatic SNV site strand bias">
##INFO=<ID=PNOISE,Number=1,Type=Float,Description="Fraction of panel containing non-reference noise at this site">
##INFO=<ID=PNOISE2,Number=1,Type=Float,Description="Fraction of panel containing more than one non-reference noise obs at this site">
##INFO=<ID=SomaticEVS,Number=1,Type=Float,Description="Somatic Empirical Variant Score (EVS) expressing the phred-scaled probability of the call being a false positive observation.">
##FILTER=<ID=LowEVS,Description="Somatic Empirical Variant Score (SomaticEVS) is below threshold">
##FILTER=<ID=LowDepth,Description="Tumor or normal sample read depth at this locus is below 2">
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed.">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Variant Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

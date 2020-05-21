#!/usr/bin/env perl
# Simple script to pull somatic calls from a MuTect2 run for MOMA processing.
#
# 5/4/2020 - D Sims
################################################################################
use strict;
use warnings;
use autodie;

use Data::Dump;
use Getopt::Long;
use File::Basename;

my $version = '1.0.050420';
my $scriptname = basename($0);
my $description = <<"EOT";
Since the MuTect2 VCF files have both a germline and tumor sample in order to 
determine if a variant is somatic or germline, we need to parse out just the
tumor data so that MOMA can interpret the fields correctly.  Use VCF Tools to 
generate a simplified VCF that can be read and processed the rest of the way.
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

unless (qx(which vcftools)) {
    die "ERROR: You must have VCF Tools installed and avaiable to your \$PATH!\n";
}

my $vcf = shift;
die "ERROR: You must input a VCF file file!\n" unless $vcf;

my $tumor_sample = get_tumor_id(\$vcf);
print "The tumor specimen ID is: $tumor_sample\n";

my @wanted_fields = qw(%CHROM %POS %ID %REF %ALT %QUAL %FILTER %INFO [%GTR %AD %AF]);
my $vcf_format = join('\t', @wanted_fields);
my @extracted_data = qx( vcf-query $vcf -c $tumor_sample -f "$vcf_format\n" );

make_new_vcf($tumor_sample, \@extracted_data, $outfile);

sub make_new_vcf {
    my ($sample_id, $data, $outfile) = @_;

    my $outfh;
    if ($outfile) {
        print "Writing data to '$outfile'.\n";
        open($outfh, ">", $outfile);
    } else {
        $outfh = \*STDOUT;
    }

    while (<DATA>) {
        print {$outfh} $_;
    }

    no warnings; # needed to print the header with a '#' char in it. 
    my @header = qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
    use warnings;
    print {$outfh} join("\t", @header, $sample_id), "\n";

    for my $var (@$data) {
        # Filter out non "PASS" variants as I'm told non-somatic calls will not
        # have a PASS string, which is the point of this algorithm. Not sure why
        # the calls aren't just filtered out of the VCF.
        if ($var =~ /PASS/) {
            chomp $var;
            my @data = split(/\t/, $var);
            my $format = "GT:AD:AF";
            my $fdata = join(":", splice(@data, -3));
            print {$outfh} join("\t", @data, $format, $fdata), "\n";
        }
    }
}

sub get_tumor_id {
    my $vcf = shift;
    my $tumor_id;

    open(my $fh, "<", $$vcf);
    while (<$fh>) {
        last unless $_ =~ /^#/;
        if ($_ =~ /^##tumor_sample=(.*?)$/) {
            return $1;
        }
    }

    # If we got here, we didn't have the appropriate information to know which 
    # is which; probably not a MuTect file.
    die "ERROR: no tumor / normal identification in VCF file. Is this data " .
        "really from a MuTect2 pipeline?\n";
}

# Simple VCF header just to pass validation.
__DATA__
##fileformat=VCFv4.2
##source=MOMA
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chrM,length=16571,assembly=hg19>
##contig=<ID=chr1,length=249250621,assembly=hg19>
##contig=<ID=chr2,length=243199373,assembly=hg19>
##contig=<ID=chr3,length=198022430,assembly=hg19>
##contig=<ID=chr4,length=191154276,assembly=hg19>
##contig=<ID=chr5,length=180915260,assembly=hg19>
##contig=<ID=chr6,length=171115067,assembly=hg19>
##contig=<ID=chr7,length=159138663,assembly=hg19>
##contig=<ID=chr8,length=146364022,assembly=hg19>
##contig=<ID=chr9,length=141213431,assembly=hg19>
##contig=<ID=chr10,length=135534747,assembly=hg19>
##contig=<ID=chr11,length=135006516,assembly=hg19>
##contig=<ID=chr12,length=133851895,assembly=hg19>
##contig=<ID=chr13,length=115169878,assembly=hg19>
##contig=<ID=chr14,length=107349540,assembly=hg19>
##contig=<ID=chr15,length=102531392,assembly=hg19>
##contig=<ID=chr16,length=90354753,assembly=hg19>
##contig=<ID=chr17,length=81195210,assembly=hg19>
##contig=<ID=chr18,length=78077248,assembly=hg19>
##contig=<ID=chr19,length=59128983,assembly=hg19>
##contig=<ID=chr20,length=63025520,assembly=hg19>
##contig=<ID=chr21,length=48129895,assembly=hg19>
##contig=<ID=chr22,length=51304566,assembly=hg19>
##contig=<ID=chrX,length=155270560,assembly=hg19>
##contig=<ID=chrY,length=59373566,assembly=hg19>
##contig=<ID=chr1_gl000191_random,length=106433,assembly=hg19>
##contig=<ID=chr1_gl000192_random,length=547496,assembly=hg19>
##contig=<ID=chr4_ctg9_hap1,length=590426,assembly=hg19>
##contig=<ID=chr4_gl000193_random,length=189789,assembly=hg19>
##contig=<ID=chr4_gl000194_random,length=191469,assembly=hg19>
##contig=<ID=chr6_apd_hap1,length=4622290,assembly=hg19>
##contig=<ID=chr6_cox_hap2,length=4795371,assembly=hg19>
##contig=<ID=chr6_dbb_hap3,length=4610396,assembly=hg19>
##contig=<ID=chr6_mann_hap4,length=4683263,assembly=hg19>
##contig=<ID=chr6_mcf_hap5,length=4833398,assembly=hg19>
##contig=<ID=chr6_qbl_hap6,length=4611984,assembly=hg19>
##contig=<ID=chr6_ssto_hap7,length=4928567,assembly=hg19>
##contig=<ID=chr7_gl000195_random,length=182896,assembly=hg19>
##contig=<ID=chr8_gl000196_random,length=38914,assembly=hg19>
##contig=<ID=chr8_gl000197_random,length=37175,assembly=hg19>
##contig=<ID=chr9_gl000198_random,length=90085,assembly=hg19>
##contig=<ID=chr9_gl000199_random,length=169874,assembly=hg19>
##contig=<ID=chr9_gl000200_random,length=187035,assembly=hg19>
##contig=<ID=chr9_gl000201_random,length=36148,assembly=hg19>
##contig=<ID=chr11_gl000202_random,length=40103,assembly=hg19>
##contig=<ID=chr17_ctg5_hap1,length=1680828,assembly=hg19>
##contig=<ID=chr17_gl000203_random,length=37498,assembly=hg19>
##contig=<ID=chr17_gl000204_random,length=81310,assembly=hg19>
##contig=<ID=chr17_gl000205_random,length=174588,assembly=hg19>
##contig=<ID=chr17_gl000206_random,length=41001,assembly=hg19>
##contig=<ID=chr18_gl000207_random,length=4262,assembly=hg19>
##contig=<ID=chr19_gl000208_random,length=92689,assembly=hg19>
##contig=<ID=chr19_gl000209_random,length=159169,assembly=hg19>
##contig=<ID=chr21_gl000210_random,length=27682,assembly=hg19>
##contig=<ID=chrUn_gl000211,length=166566,assembly=hg19>
##contig=<ID=chrUn_gl000212,length=186858,assembly=hg19>
##contig=<ID=chrUn_gl000213,length=164239,assembly=hg19>
##contig=<ID=chrUn_gl000214,length=137718,assembly=hg19>
##contig=<ID=chrUn_gl000215,length=172545,assembly=hg19>
##contig=<ID=chrUn_gl000216,length=172294,assembly=hg19>
##contig=<ID=chrUn_gl000217,length=172149,assembly=hg19>
##contig=<ID=chrUn_gl000218,length=161147,assembly=hg19>
##contig=<ID=chrUn_gl000219,length=179198,assembly=hg19>
##contig=<ID=chrUn_gl000220,length=161802,assembly=hg19>
##contig=<ID=chrUn_gl000221,length=155397,assembly=hg19>
##contig=<ID=chrUn_gl000222,length=186861,assembly=hg19>
##contig=<ID=chrUn_gl000223,length=180455,assembly=hg19>
##contig=<ID=chrUn_gl000224,length=179693,assembly=hg19>
##contig=<ID=chrUn_gl000225,length=211173,assembly=hg19>
##contig=<ID=chrUn_gl000226,length=15008,assembly=hg19>
##contig=<ID=chrUn_gl000227,length=128374,assembly=hg19>
##contig=<ID=chrUn_gl000228,length=129120,assembly=hg19>
##contig=<ID=chrUn_gl000229,length=19913,assembly=hg19>
##contig=<ID=chrUn_gl000230,length=43691,assembly=hg19>
##contig=<ID=chrUn_gl000231,length=27386,assembly=hg19>
##contig=<ID=chrUn_gl000232,length=40652,assembly=hg19>
##contig=<ID=chrUn_gl000233,length=45941,assembly=hg19>
##contig=<ID=chrUn_gl000234,length=40531,assembly=hg19>
##contig=<ID=chrUn_gl000235,length=34474,assembly=hg19>
##contig=<ID=chrUn_gl000236,length=41934,assembly=hg19>
##contig=<ID=chrUn_gl000237,length=45867,assembly=hg19>
##contig=<ID=chrUn_gl000238,length=39939,assembly=hg19>
##contig=<ID=chrUn_gl000239,length=33824,assembly=hg19>
##contig=<ID=chrUn_gl000240,length=41933,assembly=hg19>
##contig=<ID=chrUn_gl000241,length=42152,assembly=hg19>
##contig=<ID=chrUn_gl000242,length=43523,assembly=hg19>
##contig=<ID=chrUn_gl000243,length=43341,assembly=hg19>
##contig=<ID=chrUn_gl000244,length=39929,assembly=hg19>
##contig=<ID=chrUn_gl000245,length=36651,assembly=hg19>
##contig=<ID=chrUn_gl000246,length=38154,assembly=hg19>
##contig=<ID=chrUn_gl000247,length=36422,assembly=hg19>
##contig=<ID=chrUn_gl000248,length=39786,assembly=hg19>
##contig=<ID=chrUn_gl000249,length=38502,assembly=hg19>
##reference=file:///data/MoCha/patidarr/ref/ucsc.hg19.fasta
##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">
##FILTER=<ID=artifact_in_normal,Description="artifact_in_normal">
##FILTER=<ID=base_quality,Description="alt median base quality">
##FILTER=<ID=clustered_events,Description="Clustered events observed in the tumor">
##FILTER=<ID=contamination,Description="contamination">
##FILTER=<ID=duplicate_evidence,Description="evidence for alt allele is overrepresented by apparent duplicates">
##FILTER=<ID=fragment_length,Description="abs(ref - alt) median fragment length">
##FILTER=<ID=germline_risk,Description="Evidence indicates this site is germline, not somatic">
##FILTER=<ID=mapping_quality,Description="ref - alt median mapping quality">
##FILTER=<ID=multiallelic,Description="Site filtered because too many alt alleles pass tumor LOD">
##FILTER=<ID=panel_of_normals,Description="Blacklisted site in panel of normals">
##FILTER=<ID=read_position,Description="median distance of alt variants from end of reads">
##FILTER=<ID=str_contraction,Description="Site filtered due to contraction of short tandem repeat region">
##FILTER=<ID=strand_artifact,Description="Evidence for alt allele comes from one read direction only">
##FILTER=<ID=t_lod,Description="Tumor does not meet likelihood threshold">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad  mates are filtered)">
##INFO=<ID=ECNT,Number=1,Type=Integer,Description="Number of events in this haplotype">
##INFO=<ID=IN_PON,Number=0,Type=Flag,Description="site found in panel of normals">
##INFO=<ID=NLOD,Number=A,Type=Float,Description="Normal LOD score">
##INFO=<ID=N_ART_LOD,Number=A,Type=Float,Description="log odds of artifact in normal with same allele fraction as tumor">
##INFO=<ID=POP_AF,Number=A,Type=Float,Description="population allele frequencies of alt alleles">
##INFO=<ID=P_GERMLINE,Number=A,Type=Float,Description="Posterior probability for alt allele to be germ line variants">
##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
##INFO=<ID=TLOD,Number=A,Type=Float,Description="Tumor LOD score">

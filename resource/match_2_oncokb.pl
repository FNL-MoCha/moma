#!/usr/bin/env perl
# 7/30/2018 - D Sims
################################################################################
use strict;
use warnings;
use autodie;

use Getopt::Long qw(:config bundling auto_abbrev no_ignore_case);
use File::Basename;
use Data::Dump;

my $version = '1.1.061419';
my $scriptname = basename($0);
my $description = <<'EOT';
Starting with an NCI-MATCH BED file that has been converted to a tsv file using
Ion's `tvcutils prepare_hotspots` tool, and then further processed with Annovar,
generate a formatted dataset that can easily be merged into the OnocoKB dataset 
being use for variant annotation. 

One caveat to this process is that the Ion dataset won't have the same "Effect"
or "Oncogenicity" fields that we have in the OncoKB set, and so we will use the 
following rules:

    1. All variants will be defined as being "Likely Oncogenic" since they are 
       known MOIs and truly likely oncogenic.
    2. All variants in TSGs will be defined as "Loss-of-function", while all 
       variants in oncogenes will be defined as "Gain-of-function".  This is not
       100% correct, but will suffice for now.

Note that in addition to the required VCF file, a file containing a list of TSGs
that are part of the panel (one gene per line), as well as, the OncoKB file into 
which you are going to merge the data must be input.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] -T <TSG_file> -O <OncoKB_file> -I <Ion_HS.vcf>
    Required Files:
        -T, --TSG_file     File containing Tumor Suppressor Genes (TSGs) covered
                           by the panel.
        -O, --OncoKB_file  OncoKB allAnnotatedVariants file that we want to merge
                           with the Ion data.
        -I, --Ion_file     Ion Torrent Hotspots VCF that we want to convert.

    Other Options:
        -o, --outfile    Name of output file to which the data should be written
        (default: STDOUT).
        -v, --version    Print the version information and exit.
        -h, --help       Print this help text and exit.
EOT

my $help;
my $ver_info;
my $tsg_file;
my $oncokb_file;
my $hs_file;
my $outfile;

GetOptions( 
    'TSG_file|T=s'    => \$tsg_file,
    'OncoKB_file|O=s' => \$oncokb_file,
    'Ion_file|I=s'    => \$hs_file,
    'outfile|o=s'     => \$outfile,
    'version|v'       => \$ver_info,
    'help|h'          => \$help,
) or die "$usage\n";

sub help {
    print "$scriptname - $version\n\n$description\n$usage\n";
    exit;
}

sub version {
    print "$scriptname - $version\n";
    exit;
}

help if $help;
version if $ver_info;

die "ERROR: You must input a TSG file, an OncoKB file and an Ion VCF to ",
    "process!\n" unless ($tsg_file and $oncokb_file and $hs_file);

my $out_fh;
if ($outfile) {
    print "Writing output to $outfile.\n";
    open($out_fh, ">", $outfile);
} else {
    $out_fh = \*STDOUT;
}

################------ END Arg Parsing and Script Setup ------#################

# Get a TSG set to map.
open(my $tsg_fh, "<", $tsg_file);
my %tsgs = map{chomp; $_ => ''} <$tsg_fh>;
close $tsg_fh;

my %results;
my %oncokb_genes;
open(my $okb_fh, "<", $oncokb_file);
my $header = readline($okb_fh); # Dump the header.
while (<$okb_fh>) {
    chomp(my @fields = split(/\t/));

    # Have a few oddball genes that are not found in most datasets that I
    # have. Not sure their significance, so skip them for now, since not much in
    # OncoKB file at the moment.
    next if grep { $fields[0] eq $_ } qw(CIITA BCORL1 CXORF67);

    $results{"$fields[0]:$fields[2]"} = \@fields;
    $oncokb_genes{$fields[0]} = $fields[1]; #build gene / tscript lookup.

}
close $okb_fh;


# dd \%results;
# dd \%oncokb_genes;
# exit;

sub getent {
    my $candidates = shift;

    for my $var (@$candidates) {
        my ($gene, $tscript, $exon, $cds, $aa) = split(/:/, $var);

        # Check to make sure we have all fields: gene, transcript, exon, cds, aa.
        if (grep { ! defined $_ } qq($gene $tscript $exon $cds $aa)) {
            print "WARN: Fewer than expected 5 fields in entry. May be a problem\n";
            dd (split(/:/, $var));
            next;
        }

        $tscript =~ s/\.[0-9]+$//; # Allow for matching all but version string.
        if (exists $oncokb_genes{$gene}) {
            if ($oncokb_genes{$gene} =~ /$tscript/) {
                return join(':', ($gene, $oncokb_genes{$gene}, $exon, $cds, $aa));
            }
        } else {
            # Gene not in OncoKB, so let's use Oncomine
            open (my $stream, "-|", "get_transcript_id.pl $gene");
            my ($oca_script) = map { chomp; (split(/,/))[1] } <$stream>;
            if ($oca_script =~ /$tscript/) {
                return $var;
            }
        }
    }
    # If we made it here, then we couldn't find any match. Will have to resort
    # to the Oncomine Lookup.
    print "\nWARN: Can't find a suitable candidate for var:\n\t";
    dd $candidates;
    print "You may have to manually input this one.\n";
    return 0;
}

open(my $fh, "<", $hs_file);
while(<$fh>) {
    next unless /^chr/;
    chomp(my @fields = split(/\t/));

    # have a few intronic variants from BRCA that I think should be removed.
    next if $fields[5] eq 'intronic';

    # ANNOVAR will output all transcripts that it knows about, which means we'll
    # get multiple per gene in many cases.  For our case, we want to either use
    # the same transcript as OncoKB (highest priority) or, failing that, the
    # same as OCA. Filter out the entries that don't have a transcript that we
    # know about.

    my $selected_variant;
    if ($fields[9] eq '.') {
        # We have a variant that's a non-transcript variant (i.e. splicesite).
        my @var_data = map { "$fields[6]:$_:p.?" } split(';', $fields[7]);
        $selected_variant = getent(\@var_data);
    } else {
        $selected_variant = getent([split(',', $fields[9])]);
    } 
    my ($gene, $refseq, $exon, $hgvs_c, $hgvs_p) = split(/:/, $selected_variant);

    next unless $gene;

    # DEBUG: Check to see what you got if there is no hgvs_p...would be an edge
    # case.
    unless ($hgvs_p) {
        dd \@fields;
        print "results from `getent()`: $gene:$refseq:$exon:$hgvs_c\n";
        next;
    }

    my $effect;
    (exists $tsgs{$gene}) 
        ? ($effect = 'Loss-of-function')
        : ($effect = 'Gain-of-function');
    my $varid = "$gene:$hgvs_p";

    # push the Oncomine variant to the results hash unless it's already in there
    # from OncoKB
    $results{$varid} = [$gene, $refseq, $hgvs_p, 'Likely Oncogenic', 
        $effect, "$refseq($gene):$hgvs_p"] unless exists $results{$varid};
}

# Here lie the manually curated variants to output.
my %manually_curated = (
    'MYD88:p.S219C' => ['MYD88', 'NM_001172567.1', 'p.S219C', 'Likely Oncogenic',
        'Gain-of-function', 'NM_001172567.1(MYD88):p.S219C'],
);
@results{keys %manually_curated} = values %manually_curated;

# DEBUG: Check for dots (i.e. version strings) in each transcript to harmonize.
for (keys %results) {
    dd $results{$_} unless $results{$_}->[1] =~ /\./;
}

print {$out_fh} $header;
for my $var (sort { $results{$a}->[0] cmp $results{$b}->[0] } keys %results) {
    print {$out_fh} join("\t", @{$results{$var}}), "\n";
}

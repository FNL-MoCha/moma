#!/usr/bin/env perl
# Convert an Annovar TSV file to a truncated MAF file for use in other MAF 
# related tools and importation into cBioPortal.
#
# 6/17/2019 - D Sims
################################################################################
use strict;
use warnings;
use autodie;

use Getopt::Long qw(:config bundling auto_abbrev no_ignore_case);
use File::Basename;
use Term::ANSIColor;
use Data::Dump;


my $version = '1.1.042120';

use constant DEBUG => 0;
use constant TRUE  => 1;
use constant FALSE => 0;

my $scriptdir = dirname($0);
my $gene_reference_file = "$scriptdir/../resource/gene_reference.csv";

my $scriptname = basename($0);
my $description = <<"EOT";
Convert a standard Annovar file into a valid MAF file that can be used with 
other MAF analysis tools.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <Annovar File>

    Program Options
    -o, --outfile  Custom filename for the data. Default is <annovar_file>.maf
    -v, --version  Print the version information and exit.
    -h, --help     Print this help text and exit.
EOT

my $help;
my $ver_info;
my $outfile;
my $sample_name;

GetOptions (
    'sample|s=s'    => \$sample_name,
    "outfile|o=s"   => \$outfile,
    "help|h"        => \$help,
    "version|v"     => \$ver_info ) 
or die "\n$usage";

sub help {
    print "$scriptname - $version\n\n$description\n\n$usage\n";
    exit;
}

sub version_info {
    print "$scriptname - $version\n";
    exit;
}

help() if $help;
version_info() if $ver_info;

my $err   = colored("ERROR: ", "bold red on_black");
my $warn  = colored("WARN: ", "bold yellow on_black");
my $debug = colored("DEBUG: ", "bold cyan on_black");
my $info  = colored("INFO: ", "bold green on_black");

if (@ARGV < 1) {
    print "$err You must pass an Annovar file to this script!\n\n";
    print $usage;
    exit 1;
}
my $annovar_file = shift;
print "$warn is this an Annovar Output file?\n" unless $annovar_file =~ /txt$/;
if (! defined $sample_name) {
    ($sample_name = $annovar_file) =~ s/\.annovar\.txt//;
}

my $gene_reference_data = read_cantran($gene_reference_file);

my $outfh;
($outfile = $annovar_file) =~ s/\.annovar.txt/.maf/ unless $outfile;
open ($outfh, ">", $outfile);

#########--------------- END Arg Parsing and validation ---------------########
my $var_data = parse_annovar(\$annovar_file);
print_data($var_data, $outfh);


sub print_data {
    my ($data, $outfh) = @_;

    # Print the necessary MAF fields for the annotator
    my @header_order = qw(Hugo_Symbol Entrez_Gene_Id Chromosome Start_Position
        End_Position Strand Variant_Classification Variant_Type Reference_Allele
        Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode 
        Matched_Norm_Sample_Barcode Transcript_ID RefSeq HGVSc HGVSp HGVSp_Short
        Existing_variation Exon_Number Consequence t_depth t_ref_count
        t_alt_count i_TumorVAF n_depth n_ref_count n_alt_count SIFT PolyPhen
        Clinvar_Significance Clinvar_Review_Status ExAC_AF ExAC_AF_AFR 
        ExAC_AF_AMR ExAC_AF_EAS ExAC_AF_FIN ExAC_AF_NFE ExAC_AF_OTH ExAC_AF_SAS
        gnomAD_AF gnomAD_AFR_AF gnomAD_AMR_AF gnomAD_ASJ_AF gnomAD_EAS_AF 
        gnomAD_FIN_AF gnomAD_NFE_AF gnomAD_OTH_AF gnomAD_SAS_AF);

    print {$outfh} "#version 2.4\n";
    print {$outfh} join("\t", @header_order), "\n";

    for my $var (@$data) {
        print {$outfh} join("\t", @$var{@header_order}), "\n";
    }
}

sub parse_annovar {
    # Make an array of hashes to hold all the data.  Will be easier to work
    # with later. Get rid of synonymous and intronic variants.
    my $annovar_file = shift;
    my @data;
    my @not_in_okb_genes;

    open(my $fh, "<", $$annovar_file);
    my $header_line = readline($fh);
    chomp(my @tmp_header = split(/\t/, $header_line));

    my @header = map{ (split(/\./, $_))[0] } @tmp_header;

    # Add in header names for the rest of the "Otherinfo" columns
    @header = (@header, qw(otherinfo2 otherinfo3 vcf_chr vcf_pos vcf_varid 
        vcf_ref vcf_alt vcf_qual vcf_filter vcf_info vcf_format vcf_sample));

    my %var_data;
    while (my $line = <$fh>) {
        chomp(@var_data{@header} = split(/\t/, $line));

        # WES produces some interesting results.
        next if $var_data{'ExonicFunc'} =~ /unknown/i 
            && $var_data{'AAChange'} =~ /UNKNOWN/i;

        # DEBUG 
        # next unless $var_data{'Gene'} eq 'FANCI';
        # print "$debug Skipping some entries!\n";

        # NOTE: Just some basic obvious filters.
        my @unwanted_locations = qw(intronic ncrna upstream downstream utr);
        next if grep { $var_data{'Func'} =~ /$_/i } @unwanted_locations;
        next if $var_data{'ExonicFunc'} =~ /\bsynonymous/i;

        # NOTE: MYO18A / TIAF is a bit tough to work with since they overlap. As
        # of now, there is not much data in OncoKB, My Cancer Genome, ExAC, etc.
        # suggesting oncogenicity / actionability.  We'll skip these for now,
        # and revisit later if the data improve.
        if ($var_data{'Gene'} eq 'TIAF1;MYO18A') {
            print( "Skipping $var_data{'Gene'} since we don't ",
                "know about oncogenicity and mapping of these vars.\n");
            next;
        }

        # NOTE: Genes U2AF1 and U2AF1L5 completely overlap and we can get the
        # Annovar Gene1;Gene2 notation just like with MYO18A above.  Just keep
        # the U2AF1 entries since we have hotspots specifically for that gene.
        if ($var_data{'Gene'} eq 'U2AF1;U2AF1L5') {
            print("Got U2AF1;U2AF1L5 entry. Selecting U2AF1 from list.\n");
            $var_data{'Gene'} = 'U2AF1';
            my @tx_data;
            ($var_data{'AAChange'} ne '.') 
                ? (@tx_data = split(/,/, $var_data{'AAChange'}))
                : (@tx_data = split(/;/, $var_data{'GeneDetail'}));
            my @selected = grep { /U2AF1\b/ } @tx_data;
            $var_data{'AAChange'} = join(',', @selected);
        }

        # If we can't find the gene in the reference, then it is not in OncoKB
        # and we will not report on it. Store in an array to avoid reporing
        # duplicates.
        next if grep { $var_data{'Gene'} eq $_ } @not_in_okb_genes;

        if ( ! exists $gene_reference_data->{$var_data{'Gene'}} ) {
            push(@not_in_okb_genes, $var_data{'Gene'});
            next;
        }

        my $parsed_data = translate_annovar(\%var_data);
        push(@data, $parsed_data) if $parsed_data ne 0;
    }
    print "$info Total parsed and retained variants: " . scalar(@data) . "\n";
    return \@data;
}

sub translate_annovar {
    # Convert between MAF fields and Annovar fields. Expecting an array ref of
    # variant data from Annovar. If intronic and not a splicesite variant, then
    # dump the entry.  Otherwise, map the relevent fields to those in a MAF
    # file.
    my $var_data = shift;

    my %results;
    my %translator = (
        'Gene'       => 'Hugo_Symbol',
        'Chr'        => 'Chromosome',
        'Start'      => 'Start_Position',
        'End'        => 'End_Position',
        'Ref'        => 'Reference_Allele',
        'Alt'        => 'Tumor_Seq_Allele2',
        'ExonicFunc' => 'Consequence', #XXX: Need this? 
    );

    while (my ($k,$v) = each %translator) {
        $results{$v} = $var_data->{$k};
    }

    my $ret_data;
    if ($var_data->{'AAChange'} ne '.') {
        # We have a coding sequence variant.
        my @vars = split(/,/, $var_data->{'AAChange'});
        $ret_data = map_transcript(\@vars, 'coding', undef);
        $ret_data->{'Variant_Classification'} = map_consequence(
            $var_data->{'ExonicFunc'}, 'maf_term'
        );
        $ret_data->{'Consequence'} = map_consequence(
            $var_data->{'ExonicFunc'}, 'so_term'
        );
    } else {
        # We have a non-coding variant.
        my @vars = split(/;/, $var_data->{'GeneDetail'});
        if ($var_data->{'Func'} eq 'splicing') {
            # Have a splice variant.
            $ret_data = map_transcript(\@vars, 'splicesite', $var_data->{'Gene'});
        } else {
            $ret_data = map_transcript(\@vars, 'noncoding', $var_data->{'Gene'});
        }

        if ($ret_data->{'location'} eq '-') {
            $ret_data->{'location'} = $var_data->{'Func'};
        }
        $ret_data->{'Variant_Classification'} = map_consequence(
            $var_data->{'Func'}, 'maf_term'
        );
        $ret_data->{'Consequence'} = map_consequence(
            $var_data->{'Func'}, 'so_term'
        );
    }

    # Pare down COSMIC data to just COSMIC ID(s)
    my $cosstring;
    ($var_data->{'cosmic89_noEnst'} ne '.')
        ? ($cosstring = (split(/;/, $var_data->{'cosmic89_noEnst'}))[0])
        : ($cosstring = '.');
    my @cosids = map{ s/ID=//; split(/,/) } $cosstring;
    $var_data->{'cosids'} = join(',', @cosids);

    return FALSE if $ret_data->{'Transcript_ID'} eq '-';
    @results{keys %$ret_data} = values %$ret_data;

    # Add the ExAC Population data to the output, changing the header names to
    # make the MAF headers for mapping later.
    my @maf_exac_header     = qw(ExAC_AF ExAC_AF_AFR ExAC_AF_AMR ExAC_AF_EAS 
        ExAC_AF_FIN ExAC_AF_NFE ExAC_AF_OTH ExAC_AF_SAS);
    my @annovar_exac_header = qw(ExAC_ALL ExAC_AFR ExAC_AMR ExAC_EAS ExAC_FIN 
        ExAC_NFE ExAC_OTH ExAC_SAS);

    my @maf_gnomad_header= qw(gnomAD_AF gnomAD_AFR_AF gnomAD_AMR_AF gnomAD_ASJ_AF
        gnomAD_EAS_AF gnomAD_FIN_AF gnomAD_NFE_AF gnomAD_OTH_AF gnomAD_SAS_AF);
    my @annovar_gnomad_header = qw(gnomAD_exome_ALL gnomAD_exome_AFR gnomAD_exome_AMR
        gnomAD_exome_ASJ gnomAD_exome_EAS gnomAD_exome_FIN gnomAD_exome_NFE
        gnomAD_exome_OTH gnomAD_exome_SAS);

    my @pop_data = map { ($_ eq '.') ? '' : $_ } @{$var_data}{@annovar_exac_header};
    @results{@maf_exac_header} = @pop_data;

    @pop_data = map { ($_ eq '.') ? '' : $_ } @{$var_data}{@annovar_gnomad_header};
    @results{@maf_gnomad_header} = @pop_data;

    # Add in the other MAF required fields
    get_maf_fields($var_data, \%results);

    if (DEBUG) {
        print '-'x50, "\n";
        dd \%results;
        print '-'x50, "\n";
        # __exit__(__LINE__, "Finished parsing annovar data and translating.");
    }

    # dd \%results;
    # exit;
    return \%results;
}

sub get_maf_fields {
    # Add in the rest of the required MAF fields.
    my ($var_data, $collected_data) = @_;

    my $gene_data = $gene_reference_data->{$collected_data->{'Hugo_Symbol'}};

    $collected_data->{'Strand'}               = $gene_data->{'Strand'};
    $collected_data->{'Entrez_Gene_Id'}       = $gene_data->{'EntrezID'};
    $collected_data->{'Tumor_Seq_Allele1'}    = $var_data->{'Ref'};
    $collected_data->{'Tumor_Sample_Barcode'} = basename($sample_name);
    $collected_data->{'Matched_Norm_Sample_Barcode'} = 'NORMAL';

    $collected_data->{'Variant_Type'} = __get_var_type($var_data->{'Ref'},
        $var_data->{'Alt'});

    # Generate an "Existing Variation" field from COSMIC and dbSNP IDs.
    my $varid;
    my @ids = grep{ ! /\./ } ($var_data->{'avsnp142'}, $var_data->{'cosids'});
    $collected_data->{'Existing_variation'} = join(';', @ids);

    # Translate the SIFT and PolyPhen data.
    $collected_data->{'SIFT'}     = __translate_sift_polyphen($var_data, 'sift');
    $collected_data->{'PolyPhen'} = __translate_sift_polyphen($var_data, 'polyphen');

    $collected_data->{'Clinvar_Significance'} = $var_data->{'CLNSIG'};
    $collected_data->{'Clinvar_Review_Status'} = $var_data->{'CLNREVSTAT'};

    # Get tumor specific coverage information.
    my ($t_depth, $t_ref_count, $t_alt_count, $i_TumorVAF) = __get_cov_info($var_data);
    $collected_data->{'t_depth'}      = $t_depth     // '.';
    $collected_data->{'t_ref_count'}  = $t_ref_count // '.';
    $collected_data->{'t_alt_count'}  = $t_alt_count // '.';
    $collected_data->{'i_TumorVAF'}   = $i_TumorVAF  // '.';

    # There are some empty fields that are required for a valid MAF.
    my @empty = qw(n_depth n_ref_count n_alt_count);
    $collected_data->{$_} = '' for @empty;
}

sub map_transcript {
    # Map alteration to our primary transcript from the gene_reference file.
    # Also pick up the rest of the data we want from that resource.
    my ($vars, $type, $gene) = @_;

    my %fields = (
        'Hugo_Symbol'    => '-',
        'Transcript_ID'  => '-',
        'location'       => '-',
        'HGVSc'          => '-',
        'HGVSp'          => '',
        'HGVSp_Short'    => '-',
        'Exon_Number'    => '',
        'RefSeq'         => '-',
    );

    for my $candidate (@$vars) {
        my @elems = split(/:/, $candidate);
    
        local $SIG{__WARN__}  = sub {
            my $msg = shift;
            print "$warn Issue mapping transcript for the following entry:\n";
            print "\tGene: $gene; type: $type; candidate: $candidate\n";
            dd \@elems;
            print $msg;
            # print "Skipping this entry.\n";
            # print "Skipping for now....this is a temporary fix!\n";
            # exit 1;
        };
        
        if ($type eq 'noncoding' || $type eq 'splicesite') {
            my $tx_root = (split(/\./, $elems[0]))[0];
            if ($gene_reference_data->{$gene}{'RefSeq'} =~ /$tx_root/) {
                if ($type eq 'noncoding')  {
                    @fields{qw(Hugo_Symbol Transcript_ID RefSeq HGVSc)} = (
                        $gene, 
                        $gene_reference_data->{$gene}{'Isoform'}, 
                        $gene_reference_data->{$gene}{'RefSeq'}, 
                        $elems[1]
                    );
                } else {
                    @fields{qw(Hugo_Symbol Transcript_ID RefSeq location HGVSc)} = (
                        $gene, 
                        $gene_reference_data->{$gene}{'Isoform'}, 
                        $gene_reference_data->{$gene}{'RefSeq'}, 
                        @elems[1,2]
                    );
                }
                return \%fields;
            }
        } else {
            $gene = $elems[0];
            my $tx_root = (split(/\./, $elems[1]))[0];

            if ($gene_reference_data->{$gene}{'RefSeq'} =~ /$tx_root/) {
                $elems[1] = $gene_reference_data->{$gene}{'RefSeq'};
                my @field_order = qw(Hugo_Symbol RefSeq location HGVSc
                    HGVSp_Short);
                @fields{@field_order} = @elems;

                # Make some adjustments and add a little more data.
                $fields{'Transcript_ID'} = $gene_reference_data->{$gene}{'Isoform'};
                $fields{'HGVSp_Short'} =~ tr/X/*/;
                $fields{'HGVSp'} = convert_aa_change($fields{'HGVSp_Short'});
                if ($fields{'location'} =~ /^exon/) {
                    $fields{'Exon_Number'} = sprintf("%s/%s", 
                        $fields{'location'} =~ s/exon//r,
                        $gene_reference_data->{$gene}{'Num_Exons'});
                }
                return \%fields;
            }
        }
    }

    # If we got here, then we couldn't find a transcript that was in our DB.
    print '-'x50, "\n";
    print "$warn No suitable variant entries / transcripts found!\n";
    print "  Var candidate(s):\n\t  $gene => [\n";
    print "\t\t$_\n" for @$vars;
    print "\n\t  ]\n";
    print "  Primary Transcript: $gene_reference_data->{$gene}{'RefSeq'}\n"; 
    print "Skipping entry.\n";
    print '-'x50, "\n";
    return \%fields;
}

sub convert_aa_change {
    # Input the HGVSp_Short (the default from my script) and output the HGVSp
    # (long form) version for a MAF file.
    my $aa_change = shift;
    my @converted;
    my @chars = split('', $aa_change);
    push(@converted, (/[A-Z]/) ? __convert_aa($_) : $_) for @chars;
    return join('', @converted);
}

sub map_consequence {
    # Translate between annovar functional consequence, MAF consequence, and
    # sequence ontology (SO) term.
    my ($query, $outterm) = @_;
    my $map = {
         'frameshift deletion'               => [qw(Frame_Shift_Del frameshift_variant)],
         'frameshift insertion'              => [qw(Frame_Shift_Ins frameshift_variant)],
         'frameshift block substitution'     => [qw(Frame_Shift_Sub frameshift_variant)],
         'frameshift substitution'           => [qw(Frame_Shift_Sub frameshift_variant)],
         'nonframeshift deletion'            => [qw(In_Frame_Del inframe_deletion)],
         'nonframeshift insertion'           => [qw(In_Frame_Ins inframe_insertion)],
         'nonframeshift substitution'        => [qw(Missense_Mutation missense_variant)],
         'nonframeshift block substitution'  => [qw(In_Frame_Sub missense_variant)],
         'nonsynonymous SNV'                 => [qw(Missense_Mutation missense_variant)],
         'stopgain'                          => [qw(Nonsense_Mutation stop_gained)],
         'stoploss'                          => [qw(Nonstop_Mutation stop_lost)],
         'splicing'                          => [qw(Splice_Site splice_region_variant)],
         'UTR3'                              => [qw(3'UTR 3_prime_UTR_variant)],
         'UTR5'                              => [qw(5'UTR 5_prime_UTR_variant)],
         'intronic'                          => [qw(Intron intron_variant)],
         'ncRNA'                             => [qw(RNA non_coding_transcript_variant)],
    };

    if (! exists $map->{$query})  {
         print("$warn the Annovar Term '$query' does not map to any MAF ",
             "consequence!\n");
         return "UNK";
     } else {
        ($outterm eq 'maf_term') 
            ? return $map->{$query}[0]
            : return $map->{$query}[1];
    }
}

sub read_cantran {
    # Read in the refseq file and output a hash that we can use to filter data
    # to the primary transcript. 
    my $txfile = shift;
    die ("ERROR: Can not find the transcript file: '$txfile'!\n") unless -e $txfile;

    my %tx_db;
    open (my $fh, "<", $txfile);
    # return map{ chomp; split(/,/, $_, 2) } grep { ! /^#/ } <$fh>;
    my $tx_version = readline($fh);
    print "$info Gene and Transcript DB version: $tx_version";
    chomp(my @header = split(/,/, readline($fh)));
    while (<$fh>) {
        chomp(my @elems = split(/,/, $_));
        @{$tx_db{$elems[3]}}{@header} = @elems;
    }
    return \%tx_db;
}

sub __convert_aa {
    my $query = shift;

    my $single_letter = 'ACDEFGHIKLMNPQRSTVWY*';
    my @three_letter = qw(Ala Cys Asp Glu Phe Gly His Ile Lys Leu Met Asn Pro
        Gln Arg Ser Thr Val Trp Tyr Ter');

    my %lookup;
    (length($query) == 1)
        ? (@lookup{split('', $single_letter)} = @three_letter)
        : (@lookup{@three_letter} = split('', $single_letter));
    return $lookup{$query};
}

sub __get_var_type {
    my ($ref, $alt) = @_;

    if ($ref eq '-') {
        return 'INS';
    }
    elsif ($alt eq '-') {
        return 'DEL';
    }
    elsif (length($ref) == 1 and length($alt) == 1) {
        return 'SNP';
    }
    elsif (length($ref) == 2 and length($alt) == 2) {
        return 'DNP';
    }
    elsif (length($ref) == 3 and length($alt) == 3) {
        return 'TNP';
    }
    else {
        return "MNV";
    }
}

sub __translate_sift_polyphen {
    my ($data, $type) = @_;

    # Map the qbbreviations in the output to full terms. VEP adds some
    # confidence terms and scores to the output, but it's not clear how these
    # data are generated or if that level of detail is even useful.
    my $mapper = {
        'sift' => {
            'T' => 'tolerated',
            'D' => 'delterious',
        },
        'polyphen' => {
            'D' => 'Probably Damaging',
            'P' => 'Possibly Damaging',
            'B' => 'Benign',
        }
    };

    my $query;
    if ($type eq 'sift')  {
        $query = $data->{'SIFT_pred'};
    }
    elsif ($type eq 'polyphen') {
        $query = $data->{'Polyphen2_HDIV_pred'};
    }
    my $result = $mapper->{$type}{$query};
    $result //= '';
    return $result;
}

sub __get_cov_info {
    my $data = shift;

    my %var_info;
    @var_info{split(/:/, $data->{'vcf_format'})} = split(/:/, $data->{'vcf_sample'});

    my ($ref_count, $alt_count, $vaf);
    # For some reason, some WES entries do not have an AD field, so one can not
    # determine the ref and alt reads easily.  For now, try to back calculate
    # from AF and DP. May want to ultimately skip these entries?
    # if (! $var_info{'AD'}) { dd \%var_info; dd $data; exit; }
    # XXX
    if (! $var_info{'AD'}) { 
        print "$warn No 'AD' field for this entry. Attempting to calculate " . 
            "values from AF and DP.\n";
        if (! $var_info{'DP'} || ! $var_info{'AF'}) {
            print "$warn Can not get coverage info for this variant. Skipping.\n";
            return 0;
        }
        $vaf = $var_info{'AF'};
        $alt_count = $vaf * $var_info{'DP'};
        $ref_count = $var_info{'DP'} - $alt_count;
    } else {
        ($ref_count, $alt_count) = split(/,/, $var_info{'AD'});
        $vaf = sprintf("%.4f", $alt_count / ($alt_count + $ref_count));
    }
    return ($var_info{'DP'}, $ref_count, $alt_count, $vaf);
}

sub __exit__ {
    my ($lineno, $msg) = @_;
    $msg //= '';
    print colored("Stopped at $lineno with message: '$msg'\n", 
        'bold green on_black');
    exit;
}

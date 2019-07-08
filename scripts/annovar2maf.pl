#!/usr/bin/env perl
# Convert an Annovar TSV file to a MAF file for use in other MAF related tools.
# 6/17/2019 - D Sims
################################################################################
use strict;
use warnings;
use autodie;

use Getopt::Long qw(:config bundling auto_abbrev no_ignore_case);
use File::Basename;
use Term::ANSIColor;
use Data::Dump;


my $version = '0.7.070119';

use constant DEBUG => 0;
use constant TRUE  => 1;
use constant FALSE => 0;

my $scriptdir = dirname($0);
my $cantran = "$scriptdir/../resource/refseq.txt";

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

GetOptions (
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

my %transcript_db = read_cantran($cantran);

my $outfh;
($outfile = $annovar_file) =~ s/\.txt/.maf/ unless $outfile;
# print("$info Writing output data to $outfile\n");
open ($outfh, ">", $outfile);

#########--------------- END Arg Parsing and validation ---------------########

my $var_data = parse_annovar(\$annovar_file);
print_data($var_data, $outfh);

sub print_data {
    # my ($data, $header, $outfh) = @_;
    my ($data, $outfh) = @_;

    my $first_entry= $data->[0];

    # Print the necessary MAF fields for the annotator
    my @header_order = qw(Hugo_Symbol Chromosome Start_Position End_Position 
         Reference_Allele Tumor_Seq_Allele2 Transcript_ID HGVSc HGVSp_Short 
         Exon_Number Variant_Classification ExAC_AF ExAC_AF_AFR ExAC_AF_AMR 
         ExAC_AF_EAS ExAC_AF_FIN ExAC_AF_OTH ExAC_AF_SAS);


    print {$outfh} "#version 2.4\n";
    print {$outfh} join("\t", @header_order), "\n";

    for my $var (@$data) {
        local $SIG{__WARN__}  = sub {
            my $msg = shift;
            print "Issue with the following entry:\n";
            dd $var;
            print $msg;
            exit 1;
        };
        print {$outfh} join("\t", @$var{@header_order}), "\n";
    }
}

sub parse_annovar {
    # Make an array of hashes to hold all the data.  Will be easier to work
    # with later. Get rid of synonymous and intronic variants.
    my $annovar_file = shift;
    my @data;

    open(my $fh, "<", $$annovar_file);
    my $header_line = readline($fh);
    chomp(my @header = split(/\t/, $header_line));

    # Add in header names for the rest of the "Otherinfo" columns
    @header = (@header, qw(otherinfo2 otherinfo3 vcf_chr vcf_pos vcf_varid 
        vcf_ref vcf_alt vcf_qual vcf_filter vcf_info vcf_format vcf_sample));

    my %var_data;
    while (my $line = <$fh>) {
        chomp(@var_data{@header} = split(/\t/, $line));

        # XXX
        # next unless $var_data{'Gene.refGeneWithVer'} eq 'CDKN2A';
        # print "$debug Skipping some entries!\n";

        # NOTE: Just some basic obvious filters.
        next if $var_data{'Func.refGeneWithVer'} =~ /(intronic|ncrna)/i;
        next if $var_data{'ExonicFunc.refGeneWithVer'} =~ /\bsynonymous/i;

        # NOTE: MYO18A / TIAF is a bit tough to work with since they overlap. As
        # of now, there is not much data in OncoKB, My Cancer Genome, ExAC, etc.
        # suggesting oncogenicity / actionability.  We'll skip these for now,
        # and revisit later if the data improve.
        if ($var_data{'Gene.refGeneWithVer'} eq 'TIAF1;MYO18A') {
            print( "Skipping $var_data{'Gene.refGeneWithVer'} since we don't ",
                 "know about oncogenicity and mapping of these vars.\n");
             next;
         }

        my $parsed_data = translate_annovar(\%var_data);
        ($parsed_data eq 0)
            ? print "Skipping entry.\n"
            : push(@data, {%var_data, %$parsed_data});
    }
    print "Total parsed and retained variants: " . scalar(@data) . "\n";
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
        'Gene.refGene'  => 'Hugo_Symbol',
        'Chr'           => 'Chromosome',
        'Start'         => 'Start_Position',
        'End'           => 'End_Position',
        'Ref'           => 'Reference_Allele',
        'Alt'           => 'Tumor_Seq_Allele2',
    );

    while (my ($k,$v) = each %translator) {
        $results{$v} = $var_data->{$k};
    }

    my $ret_data;
    if ($var_data->{'AAChange.refGeneWithVer'} ne '.') {
        # We have a coding sequence variant.
        my @vars = split(/,/, $var_data->{'AAChange.refGeneWithVer'});
        $ret_data = map_transcript(\@vars, 'coding', undef);
        $ret_data->{'Variant_Classification'} = map_consequence(
            $var_data->{'ExonicFunc.refGeneWithVer'}
        );
    } else {
        # We have a non-coding variant.
        my @vars = split(/;/, $var_data->{'GeneDetail.refGeneWithVer'});
        if ($var_data->{'Func.refGeneWithVer'} eq 'splicing') {
            # Have a splice variant.
            $ret_data = map_transcript(\@vars, 'splicesite',
                $var_data->{'Gene.refGeneWithVer'});
        } else {
            $ret_data = map_transcript(\@vars, 'noncoding',
                $var_data->{'Gene.refGeneWithVer'});
        }
        if ($ret_data->{'Exon_Number'} eq '-') {
            $ret_data->{'Exon_Number'} = $var_data->{'Func.refGeneWithVer'};
        }
        $ret_data->{'Variant_Classification'} = map_consequence(
            $var_data->{'Func.refGeneWithVer'}
        );
    }
    return FALSE if $ret_data->{'Transcript_ID'} eq '-';
    @results{keys %$ret_data} = values %$ret_data;

    # Add the ExAC Population data to the output, changing the header names to
    # make the MAF headers for mapping later.
    my @maf_exac_header     = qw(ExAC_AF ExAC_AF_AFR ExAC_AF_AMR ExAC_AF_EAS 
        ExAC_AF_FIN ExAC_AF_NFE ExAC_AF_OTH ExAC_AF_SAS);
    my @annovar_exac_header = qw(ExAC_ALL ExAC_AFR ExAC_AMR ExAC_EAS ExAC_FIN 
        ExAC_NFE ExAC_OTH ExAC_SAS);
    my @pop_data = map { ($_ eq '.') ? '' : $_ } @{$var_data}{@annovar_exac_header};
    # @results{@maf_exac_header} = @{$var_data}{@annovar_exac_header};
    @results{@maf_exac_header} = @pop_data;


    if (DEBUG) {
        print '-'x50, "\n";
        dd \%results;
        print '-'x50, "\n";
        # __exit__(__LINE__, "Finished parsing annovar data and translating.");
    }
    return \%results;
}

sub map_transcript {
    my ($vars, $type, $gene) = @_;

    my %fields = (
        'Hugo_Symbol'    => '-',
        'Transcript_ID'  => '-',
        'Exon_Number'    => '-',
        'HGVSc'          => '-',
        'HGVSp_Short'    => '-',
    );

    for my $candidate (@$vars) {
        my @elems = split(/:/, $candidate);
        local $SIG{__WARN__}  = sub {
            my $msg = shift;
            print "Issue mapping transcript for the following entry:\n";
            print "Gene: $gene; type: $type; candidate: ";
            dd $candidate;
            dd \@elems;
            print $msg;
            exit 1;
        };
        
        if ($type eq 'noncoding' || $type eq 'splicesite') {
            my $tx_root = (split(/\./, $elems[0]))[0];
            if ($transcript_db{$gene} =~ /$tx_root/) {
                if ($type eq 'noncoding')  {
                    @fields{qw(Hugo_Symbol Transcript_ID HGVSc)} = (
                        $gene, $transcript_db{$gene}, $elems[1]
                    );
                } else {
                    @fields{qw(Hugo_Symbol Transcript_ID Exon_Number HGVSc)} = (
                        $gene, $transcript_db{$gene}, @elems[1,2]
                    );
                }
                return \%fields;
            }
        } else {
            $gene = $elems[0];
            my $tx_root = (split(/\./, $elems[1]))[0];
            if ($transcript_db{$gene} =~ /$tx_root/) {
                $elems[1] = $transcript_db{$gene};
                my @field_order = qw(Hugo_Symbol Transcript_ID Exon_Number HGVSc
                    HGVSp_Short);
                @fields{@field_order} = @elems;
                $fields{'HGVSp_Short'} =~ tr/X/*/;
                return \%fields;
            }
        }
    }

    # If we got here, then we couldn't find a transcript that was in our DB.
    print "WARNING: No suitable variant entries / transcripts found!\n";
    print "Input data:\n";
    print "Gene: $gene => ";
    dd $vars;

    return \%fields;
}

sub map_consequence {
    # Translate between annovar functional consequence and MAF consequences.
    my $term = shift;
    my %map = (
         'frameshift deletion'               => 'Frame_Shift_Del',
         'frameshift insertion'              => 'Frame_Shift_Ins',
         'framshift block substitution'      => "???",
         'nonframeshift deletion'            => 'In_Frame_Del',
         'nonframeshift insertion'           => 'In_Frame_Ins',
         'nonframeshift block substitution'  => '???',
         'nonsynonymous SNV'                 => 'Missense_Mutation',
         'stopgain'                          => 'Nonsense_Mutation',
         'stoploss'                          => 'Nonstop_Mutation',
         'splicing'                          => 'Splice_Site',
         'UTR3'                              => "3'UTR",
         'UTR5'                              => "5'UTR",
         'intronic'                          => 'Intron',
         'ncRNA'                             => 'RNA',
    );
    if (! exists $map{$term})  {
         print "WARNING: the Annovar Term does not map to any VEP consequence!\n";
         return "UNK";
     } else {
        return $map{$term};
    }
}

sub read_cantran {
    # Read in the refseq file and output a hash that we can use to filter data
    # to the primary transcript. 
     
    # NOTE: We have all genes, even rna and noncoding genes in this original
    # file, and we'll just keep them.  Also, there may be more than 1 transcript
    # per gene in the file (the canonical transcript changed for some genes in
    # some cases).  Keep all transcripts and filter on both later on.
    my $txfile = shift;
    die ("ERROR: Can not find the transcript file: '$txfile'!\n") unless -e $txfile;
    open (my $fh, "<", $txfile);
    return map{ chomp; split(/,/, $_, 2) } grep { ! /^#/ } <$fh>;
}

sub __exit__ {
    my ($lineno, $msg) = @_;
    print colored("Stopped at $lineno with message: '$msg'\n", 
        'bold green on_black');
    exit;
}

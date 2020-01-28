#!/usr/bin/env perl
# Read in a set of MOMA CSV reports, and ouput a single, collated file that
# contains all variant information by sample.  
#
# 1/27/2020 - D Sims
################################################################################
use strict;
use warnings;
use autodie;

use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use Text::CSV;
use Sort::Versions;
use Data::Dump;

my $version = '1.0.012820';
my $scriptname = basename($0);
my $description = <<"EOT";
Starting with a set of MOI reports from the MOMA tool, generate a single CSV
file of all samples' variants together than can be used for cohort reporting.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <MOMA Report(s)>

General Options
    -o, --outfile  Write results to a file rather than just STDOUT.
    -v, --version  Print the version information and exit.
    -h, --help     Print this help text and exit.

Filter Options
    -m, --mois_only  Only output data for MOIs; exclude VuS.
    -g, --genes      Only output results for a gene or set of genes. You can
                     enter multiple genes by separating each with a comma.
                     
    -l, --list       Instead of just one or a few genes to filter on, input a
                     file of all genes to be included in the output.
    -t, --type       Only output variants of these types.  Valid values are:
                     'snv', 'cnv', 'fusion', or 'all'.  You can choose multiple
                     types by comma separating each type.  'All' will allow all
                     types to be output.
EOT

my $help;
my $ver_info;
my $outfile;
my $filter_genes;
my $gene_file;
my $mois_only;
my $vartype;

GetOptions(
    'g|genes=s'      => \$filter_genes,
    'l|list=s'       => \$gene_file,
    'm|mois_only'    => \$mois_only,
    't|type=s'       => \$vartype,
    'o|outfile=s'    => \$outfile,
    'v|version'      => \$ver_info,
    'h|help'         => \$help
) or die "$usage";

sub help {
    print "$scriptname - $version\n$description\n$usage\n\n";
    exit;
}

sub version {
    print "$scriptname - $version\n\n";
    exit;
}

help() if $help;
version() if $ver_info;

my @moma_reports = @ARGV;
die "ERROR: You must input at least one MOMA report!\n" unless @moma_reports;

my @wanted_genes;
if ($filter_genes) {
    @wanted_genes = split(/,/, $filter_genes);
}
elsif ($gene_file) {
    open(my $fh, "<", $gene_file);
    @wanted_genes = map{ chomp; $_ } <$fh>;
    close $fh;
}

# Choose which variant types we want to report on, assuming we want to filter
# any specific types out.
my @vartypes;
if ($vartype and $vartype ne 'all') {
    my @temp = split(/,/, $vartype);
    for my $t (@temp) {
        unless (grep { $_ =~ /$t/i } qw(snv cnv fusion all)) {
            print "ERROR: '$t' is not valid. Choose from 'snv', 'cnv', " .
                 "'fusion', or 'all' only\n";
            exit 1;
        }
        push(@vartypes, ($t =~ /fusion/i) ? 'Fusion' : uc($t));
    }
} else {
    push(@vartypes, qw(SNV CNV Fusion));
}

##########----------  End Arg Parsing and Prog Setup  ----------##########

my $wanted = {
    'snv' => [qw(Ref Alt Gene Transcript Chr Pos CDS AA VAF Ref_Reads
        Alt_Reads Function Oncogenicity Effect)],
    'cnv' => [qw(Gene Chr CN Oncogenicity Effect)],
    'fusion' => [qw(Driver_Gene Fusion Read_Count Oncogenicity Effect)]
};

my %data;
for my $report (@moma_reports) {
    my ($sample_name) = $report =~ /^(.*?)_moma.*/;
    $data{$sample_name} = read_report($report);
}

print_results(\%data, $outfile);

##########----------  Methods and Subroutines  ----------##########
 
sub print_results {
    my ($data, $outfile) = @_;
    
    my $outfh;
    if ($outfile) {
        print "Writing results to $outfile.\n";
        open($outfh, ">", $outfile);
    } else {
        $outfh = \*STDOUT;
    }

    for my $sample (sort {versioncmp($a, $b)} keys %$data) {
        for my $type (@vartypes) {
            if (@{$data{$sample}->{$type}} < 1) {
                print {$outfh} join(',', $sample, $type, "None\n");
            }
            for my $var (@{$data{$sample}->{$type}}) {
                if ($mois_only) {
                    next unless $var->{'Oncogenicity'} ne '.';
                }
                if (@wanted_genes) {
                    next unless grep { $var->{'Gene'} eq $_ } @wanted_genes;
                }

                print {$outfh} join(',', $sample, $type, @$var{@{$wanted->{lc $type}}}), "\n";
            }
        }
        print {$outfh} "\n";
    }
}

sub read_report {
    my ($report, $results) = @_;

    my %parsed_data;

    local $/ = '';
    open(my $fh, "<", $report);

    my $snv_data = readline($fh);
    $parsed_data{'SNV'} = parse_block(\$snv_data, 'snv');

    my $cnv_data = readline($fh);
    $parsed_data{'CNV'} = parse_block(\$cnv_data, 'cnv');

    my $fus_data = readline($fh);
    $parsed_data{'Fusion'} = parse_block(\$fus_data, 'fusion');
    
    # dd \%parsed_data;
    # __exit__(__LINE__, "");
    return \%parsed_data;
}


sub parse_block {
    my ($block, $type) = @_;

    # If we didn't generate any data for a variant type, return an empty array.
    return [] unless $$block;

    my @ext_data;

    my $csv_reader = Text::CSV->new({ sep_char => ',' });
    my @var_data = split(/\n/, $$block);
    my @headers = splice(@var_data, 0, 2);
    my @header_elems = split(/,/, $headers[1]);

    for my $line (@var_data) {
        return [] if $line =~ /^No/;
        my %parsed;
        $csv_reader->parse($line);

        @parsed{@header_elems} = $csv_reader->fields();
        
        my %tmp;
        @tmp{@{$wanted->{$type}}} = @parsed{@{$wanted->{$type}}};

        # Add HGVS data for reporting if we have an SNV variant type.
        if ($type eq 'snv') {
            my ($hgvs_g, $hgvs_c, $hgvs_p) = __gen_hgvs(\%tmp);
            @tmp{qw(HGVS_g HGVS_c HGVS_p)} = ($hgvs_g, $hgvs_c, $hgvs_p);
        }

        if ($type eq 'cnv') {
            ($parsed{'CN'} >= 4)
                ? ($tmp{'function'} = 'Amplification')
                : ($tmp{'function'} = 'Copy Loss');
        }

        push(@ext_data, \%tmp);
    }
    # dd \@ext_data;
    # __exit__(__LINE__, '');
    return \@ext_data;
}

sub __gen_hgvs {
    my ($var_elems, $ret_type) = @_;

    local $SIG{__WARN__} = sub {
        my $msg = shift;
        print "ERROR: Can not generate a HGVS annotatino for the following:\n";
        dd $var_elems;
        print $msg;
        exit;
    };

    my %g_refseq = (
        chr1 => 'NC_000001.10', chr2 => 'NC_000002.11', chr3 => 'NC_000003.11',
        chr4 => 'NC_000004.11', chr5 => 'NC_000005.9', chr6 => 'NC_000006.11',
        chr7 => 'NC_000007.13', chr8 => 'NC_000008.10', chr9  => 'NC_000009.11',
        chr10 => 'NC_000010.10', chr11 => 'NC_000011.9', chr12 => 'NC_000012.11',
        chr13 => 'NC_000013.10', chr14 => 'NC_000014.8', chr15 => 'NC_000015.9',
        chr16 => 'NC_000016.9', chr17 => 'NC_000017.10', chr18 => 'NC_000018.9',
        chr19 => 'NC_000019.9', chr20 => 'NC_000020.10', chr21 => 'NC_000021.8',
        chr22 => 'NC_000022.10', chrX => 'NC_000023.10', chrY => 'NC_000024.9',
    );

    my @wanted_keys = qw(Gene Chr Pos Ref Alt CDS AA Transcript);

    my %data;
    @data{qw(gene chr start ref alt cds aa refseq)} = @$var_elems{@wanted_keys};

    $data{'aa'} = 'p.?' unless $data{'aa'};
    my %hgvs_annots = (
        'hgvs_g' => "$g_refseq{$data{'chr'}}:g.$data{'start'}$data{'ref'}>" .
            "$data{'alt'}",
        'hgvs_c' => "$data{'refseq'}($data{'gene'}):$data{'cds'}",
        'hgvs_p' => "$data{'refseq'}($data{'gene'}):$data{'aa'}",
    );

    ($ret_type) 
        ? return $hgvs_annots{$ret_type}
        : return @hgvs_annots{qw(hgvs_g hgvs_c hgvs_p)};
}

sub __exit__ {
    my ($lineno, $msg) = @_;
    print colored("Program exited at line $lineno with message: $msg\n",
        'bold green on_black');
    exit;
}

#!/usr/bin/env perl
# -*- coding: utf-8 -*-
# Read a hotspots BED file (must be Ion formatted) and a MAF file (from the 
# TSO500 ctDNA pipeline) and annotate variants as being MOIs or not based on the
# Hotspots BED file and the NCI-MATCH MOI Rules.
#
# 2018/07/03 - D Sims
################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;
use List::Util qw(min max);
use Sort::Versions;
use Term::ANSIColor; # Optional colorized output to terminal.
use Log::Log4perl qw(get_logger :levels);
use DateTime;
use Text::CSV;

# Import NonHotspotRules.pm for new rules engine.
use FindBin;
use lib "$FindBin::Bin/../lib/";
use NonHotspotRules;

use constant TRUE => 1;
use constant FALSE => 0;

my $DEBUG = 1;

my $version = "v0.36.092419";
my $scriptdir = dirname($0);

# Default lookup files.
my $tsg_file = "$scriptdir/../resource/gene_reference.csv";
my $hs_bed = "$scriptdir/../resource/mocha_tso500_ctdna_hotspots_v1.072018.bed";
my $oncokb_file = "$scriptdir/../resource/moma_hotspot_lookup.txt";
my $nhs_rules_json = "$scriptdir/../resource/non-hotspot_rules.json";

my $nhs_rules = NonHotspotRules->new($nhs_rules_json, $tsg_file);

for my $resource_file ($tsg_file, $hs_bed, $oncokb_file, $nhs_rules_json) {
    die "ERROR: Can not locate necessary resource file '$resource_file'! ",
        "Check that this file is in your package.\n" unless -e $resource_file;
}

my $scriptname = basename($0);
my $description = <<"EOT";
Read in a hotspots BED file (Ion Torrent formatted) or an OncoKB variants file 
(preferred!), and annotate a MAF file with the matching variant ID.  Also,
determine which variants are MOIs based on NCI-MATCH MOI rules.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] -a <annotation_method> <maf_file(s)>
    Filtering Options:
    -a, --annot        Method to use for annotation. Select from "hs_bed" or 
                       "oncokb".
    -b, --bed          Hotspots BED file to use for annotation. DEFAULT: $hs_bed
    -O, --Oncokb_file  OncoKB file to use for annotation. DEFAULT: $oncokb_file
    -p, --popfreq      Which population frequency data to use for filtering,
                       along with optional frequency threshold. Valid categories 
                       are "exac, gnomad, or 1000g.  If you want to include the 
                       threshold above which variants will be filtered, you can
                       add a colon and a float value for frequency.  For example, 
                       if you wanted to filter out any variants with a GnomAD 
                       score above 5%, you can use the option '-p gnomad:0.05' 
                       (this filter is the default).
    -n, --nhs_rules    Custom non-hotspot rules JSON file if not using the
                       default. Consult the documentation for the proper 
                       formatting of this file. DEFAULT: $nhs_rules_json

    Output Options
    -m, --mois_only    Only output variants that have passed our MOI rules.
    -t, --trim_file    Output a trimmed data file in addition to the MAF file.
                       Will only contain some of the most basic output data.
    -o, --outfile      Write output to custom filename rather than default 
                       "annotated.maf" file.
    -l, --logfile      Write log output to a custom filename rather than a file
                       called "moma_annotator_<today>.log". If you'd rather
                       not write a physical logfile, then pass `/dev/null` to 
                       this option.
    -V, --Verbose      Output messages to STDOUT as well as to a log file, along
                       with some additional loggin messages and info.

    Other Options
    -v, --version      Version information
    -h, --help         Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $annot_method;
my $mois_only = 0;
my $trim_file = 1;
my $verbose;
my $custom_logfile;
my $wanted_popfreq;
my $debugging = 0;
my $custom_nhs_json;

GetOptions( 
    "annot|a=s"     => \$annot_method,
    "mois_only|m"   => \$mois_only,
    "bed|b=s"       => \$hs_bed,
    "Oncokb|O=s"    => \$oncokb_file,
    "version|v"     => \$ver_info,
    "trim_file|t!"  => \$trim_file,
    "outfile|o=s"   => \$outfile,
    "Verbose|V"     => \$verbose,
    "help|h"        => \$help,
    "popfreq|p=s"   => \$wanted_popfreq,
    "logfile|l=s"   => \$custom_logfile,
    "debug|d"       => \$debugging,  # Undocumented.
    "nsh_rules|n=s" => \$custom_nhs_json,
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
 
# Add some color output info
my $err   = colored("ERROR:", 'bold red on_black');
my $warn  = colored("WARN:", 'bold yellow on_black');
my $info  = colored("INFO:", 'bold cyan on_black');
my $debug = colored("DEBUG:", 'bold magenta on_black');

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "$err No MAF files passed to script!\n";
    print "$usage\n";
    exit 1;
}

unless ($annot_method) {
    print "$err You must choose a method to use for annotation of these ",
        "data!\n\n";
    print "$usage\n";
    exit 1;
}

# Always going to want verbose output if we're doing debugging.
$DEBUG = TRUE if $debugging;
$verbose = TRUE if $DEBUG;

# Configure a logger.
my $logfile;
($custom_logfile)
    ? ($logfile = $custom_logfile)
    : ($logfile = 'mocha_oncokb_moi_annotator_' . now('short') . '.log');

my $loggers = 'DEBUG, Logfile';
$loggers .= ', Screen' if $verbose;

my $logger_conf = qq(
    log4perl.logger.main = $loggers

    log4perl.appender.Logfile = Log::Log4perl::Appender::File
    log4perl.appender.Logfile.filename = $logfile
    log4perl.appender.Logfile.layout = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Logfile.layout.message_chomp_before_newlist = 0
    log4perl.appender.Logfile.layout.ConversionPattern = %d [ %p ]: %m{indent=4}%n

);

my $extra_conf = qq(
    log4perl.appender.Screen = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr = 1
    log4perl.appender.Screen.layout = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Screen.layout.message_chomp_before_newlist = 0
    log4perl.appender.Screen.layout.ConversionPattern = %d [ %p ]: %m{indent=4}%n
);

$logger_conf .= $extra_conf if $verbose;

# my $logger_conf = "$scriptdir/../lib/log4perl.conf";
Log::Log4perl->init(\$logger_conf);
my $logger = get_logger();

# Levels: DEBUG, INFO, WARN, ERROR, FATAL
my $loglevel;
($DEBUG) ? ($loglevel = 'DEBUG') : ($loglevel = 'INFO');
$logger->level($loglevel);

my $intro_str = sprintf("\n%s\n\t\t    Starting MoCha OncoKB MOI Annotation " . 
    "Script\n%s\n", "="x80, "="x80);
$logger->info($intro_str);

# Load up annotation data from either hotspots BED or oncokb file.
my $annotation_data;
if ($annot_method eq 'hs_bed') {
    $logger->info("Using Hotspot BED file " . basename($hs_bed));
    $annotation_data = read_hs_bed($hs_bed);
} else {
    $logger->info("Using OncoKB file " . basename($oncokb_file));
    $annotation_data = read_oncokb_file($oncokb_file);
}

# Load up TSGs for the non-hs rules
# TODO: remove TSG reference here?
$logger->info("Using TSG file " . $nhs_rules->{'tsg_version'} . "\n");
# $logger->info("Using TSG file " . basename($tsg_file));
# my $tsgs = read_tsgs($tsg_file);

# Set up population frequency filter.
my $popfreq_category = 'gnomad';
my $popfreq_threshold = 0.05;
if ($wanted_popfreq) {
    ($popfreq_category, $popfreq_threshold) = split(/:/, $wanted_popfreq);
    $popfreq_threshold //= 0.05;
    if ( ! grep { $_ eq $popfreq_category } qw(gnomad exac 1000g) ) {
        print "$err Population frequency category not valid.\n";
        print "$usage\n";
        exit 1;
    }
}

$logger->info("Using a population frequency filter $popfreq_category -> " . 
    # "$popfreq_threshold");
    sprintf("%0.2f", $popfreq_threshold));

$nhs_rules_json = $custom_nhs_json if ($custom_nhs_json && -e $custom_nhs_json);

################------ END Arg Parsing and Script Setup ------#################

# Process each MAF file
for my $maf_file (@ARGV) {
    $logger->info( "Annotating '" . basename($maf_file) . "'..." );
    my $results;
    # TODO
    # $results = annotate_maf($maf_file, $annotation_data, $tsgs);
    $results = annotate_maf($maf_file, $annotation_data);

    # Print results.
    my $new_file;
    ($outfile)
        ? ($new_file = $outfile)
        : (($new_file = $maf_file) =~ s/\.truncmaf/.annotated.maf/);
    $logger->info( "Finished annotating. Printing results..." );
    print_results($results, $new_file, $mois_only, $trim_file);
    $logger->info("Done with $maf_file!\n\n");
}

sub print_results {
    # TODO: Do we want to use Text::CSV here to ensure good formatting, or would
    # it be better not to worry about it?  
    my ($data, $filename, $filter, $trim_file) = @_;
    
    my $filter_status;
    ($filter) ? ($filter_status = 'on') : ($filter_status = 'off');
    $logger->info( "Writing results to $filename (filter is $filter_status)");
    #print "Writing results to $filename (filter is $filter_status)\n";

    open(my $outfh, ">", $filename);

    my @trim_fields = qw(Hugo_Symbol Entrez_Gene_Id Chromosome Start_Position
         End_Position Strand Variant_Classification Variant_Type Reference_Allele
         Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode 
         Matched_Norm_Sample_Barcode HGVSc HGVSp HGVSp_Short Transcript_ID 
         t_depth t_ref_count t_alt_count n_depth n_ref_count n_alt_count 
         Existing_variation RefSeq SIFT PolyPhen ExAC_AF gnomAD_AF i_TumorVAF
         MOI_Type Oncogenicity Effect
    );

    my $trimfh;
    if ($trim_file) {
        $logger->info( "Also generating a trimmed output file.");
        (my $trim_outfile = $filename) =~ s/\.(.*?)$/_trimmed.tsv/;
        open($trimfh, ">", $trim_outfile);
        print {$trimfh} join("\t", @trim_fields), "\n";
    };

    # Print headers based on original order.
    my @header = sort{ 
        versioncmp($data->[0]{$a}, $data->[0]{$b}) }
    keys %{$data->[0]};
    shift @$data;
    print {$outfh} join("\t", @header), "\n";

    for my $var (@$data) {
        my @variant_data = @$var{@header};
        next if $variant_data[-1] eq '.' and $filter; 
        print {$outfh}  join("\t", @variant_data), "\n";
        print {$trimfh} join("\t", @$var{@trim_fields}), "\n" if ($trim_file);
    }
}

sub annotate_maf {
    my ($maf, $hotspots) = @_;
    my ($var_count, $filter_count) = (0, 0);
    my ($hsid, $category, $oncogenicity, $effect);
    my @results;

    open(my $fh, "<", $maf);
    readline($fh);  # Ditch the MAF version info
    my $csv = Text::CSV->new({ sep_char => "\t" });
    my $header = $csv->getline($fh);

    # Make the first and second elements of @results the headers (MAF version
    # and MAF elems) so that we can keep the same output order as the input.
    my $tot_header_elems = scalar(@$header);
    push(@results, {map { $header->[$_] => $_ } 0..$#{$header}});

    # Get a listing of all MOI categories for counts later.
    my $moi_summary = $nhs_rules->get_category_list();

    while (my $elems = $csv->getline($fh)) {
        $var_count++;

        my %var_data;
        @var_data{@$header} = @$elems;

        # DEBUG
        # next unless $var_data{'Hugo_Symbol'} eq 'TP53';

        # Filter out SNPs, Intronic Variants, etc.
        unless (filter_var(\%var_data, $popfreq_category, $popfreq_threshold)) {
            ($filter_count++);
            next;
        }

        my @wanted_terms = qw(Hugo_Symbol Chromosome Start_Position End_Position 
            Reference_Allele Tumor_Seq_Allele2 HGVSc HGVSp_Short Transcript_ID 
            Exon_Number Variant_Classification);
        my ($gene, $chr, $start, $end, $ref, $alt, $hgvs_c, $hgvs_p, $tscript_id, 
            $exon, $function) = @var_data{@wanted_terms};

        # If splice mutation, location may be intronic, and we need to have a
        # val for 'exon' anyway.
        # TODO
        # FIXME: Not sure if i really want it this way. Maybe shoudl keep it as
        # --- or use a '0' as an exon number. 
        $exon = 'splicesite' if ($exon eq '' and $function =~ /splice/);

        do { 
            print "$warn: Still no exon info for var: \n";
            dd @var_data{@wanted_terms} 
        } and exit if $exon eq '';

        # Annotate with a Hotspots BED file
        if ($annot_method eq 'hs_bed') {
            $logger->debug("testing $chr:$start-$end:$ref>$alt");
            $hsid = map_variant_hs($chr, $start, $end, $ref, $alt, $hotspots);
        } 

        # Annotate with an OncoKB Lookup file. 
        else {
            my $debug_header = sprintf("%s  testing: %s:%s  %s", '-'x10,
                $gene, $hgvs_p, '-'x10);
            $logger->debug($debug_header);
            ($hsid, $oncogenicity, $effect) = map_variant_oncokb($gene, $hgvs_p,
                $hotspots);
        }

        # DEBUG
        $logger->debug("> $tscript_id($gene):$hgvs_c:$hgvs_p maps to ==> $hsid");

        if ($hsid ne '.') {
            # We have a hotspot.
            $category = 'Hotspot';
            $moi_summary->{'Hotspots'}++;
        } else {
            # Check to see if captured by non-hotspot rule.
            ($category, $oncogenicity, $effect) = run_nonhs_rules($gene, $exon,
                $hgvs_p, $function);
        }

        @var_data{qw(variant_id MOI_Type Oncogenicity Effect)} = ($hsid, 
            $category, $oncogenicity, $effect);

        if ($DEBUG) {
            print "MOI category: $category\n";
            print "\t$_: $var_data{$_}\n" for qw(Hugo_Symbol MOI_Type Oncogenicity Effect);
            print "-"x75;
            print "\n\n";
        }
        push(@results, {%var_data});
    }

    # XXX
    dd \@results;
    __exit__(__LINE__, 'Stopping after running NHS rules. Just need to finish' .
        'the tally for the final output and log file.');

    # Add these new elements to our saved header in @results.
    # TODO: What am I doing here?
    for my $term (qw(MOI_Type Oncogenicity Effect)) {
        $results[0]->{$term} = $tot_header_elems++;
    }
        
    $logger->info(sprintf("\n>>>  Variant Summary  <<<\nTotal variants:\t%3s\nFiltered " .
        "out:\t%3s\nRetained:\t\t%3s\n", $var_count, $filter_count, 
        $var_count-$filter_count)
    );

    # TODO: FIXME: Need to pick this back up once we have a better count.
    # my $moi_count_string;
    # for (sort keys %moi_count) {
        # $moi_count_string .= sprintf("%-42s: %s\n", $_, $moi_count{$_});
    # }
    # $logger->info("\n>>>  MOI Counts  <<<\n$moi_count_string" );

    # dd \@results;
    exit;
    return \@results;
}

sub filter_var {
    # Remove any variants that are SNPs, Intronic, etc.  May have already been
    # filtered out prior to getting here, but nevertheless the buck stops here!
    my ($variant, $pop_filter_category, $pop_filter_threshold) = @_;

    # Functional annotation filter. 
    my $var_class = $variant->{'Variant_Classification'};
    if (grep { $var_class =~ /$_/} qw(Intron UTR Silent Flank IGR)) {
        $logger->debug(sprintf("Filtering out variant %s because it's '%s'", 
            __gen_hgvs($variant, 'hgvs_c')->[0], $var_class));
        return FALSE;
    }

    # Give some options to the population data you can use to filter. 
    my %pop_fields = (
        'gnomad'  => [qw(gnomAD_AF gnomAD_AFR_AF gnomAD_AMR_AF gnomAD_ASJ_AF
            gnomAD_EAS_AF gnomAD_FIN_AF gnomAD_NFE_AF gnomAD_OTH_AF 
            gnomAD_SAS_AF)],
        'exac'    => [qw(ExAC_AF ExAC_AF_Adj ExAC_AF_AFR ExAC_AF_AMR ExAC_AF_EAS
            ExAC_AF_FIN ExAC_AF_NFE ExAC_AF_OTH ExAC_AF_SAS)],
        '1000G'   => [qw(1000G_ALL 1000G_AFR 1000G_AMR 1000G_EAS 1000G_EUR
            1000G_SAS)],
    );

    # Get the population frequency data for the set requested, and sort them
    # largest to smallest for comparison.
    my @pop_freqs = sort {$b <=> $a } 
        map{ ($variant->{$_}) ? $variant->{$_} : 0 } @{$pop_fields{$pop_filter_category}};
    if ($pop_freqs[0] > $pop_filter_threshold) {
        $logger->debug(sprintf("Filtering out variant %s because %s frequency" .
            " is too high.", __gen_hgvs($variant, 'hgvs_c')->[0], 
            $pop_filter_category));
        return FALSE;
    }
    return TRUE;
}

sub map_variant_oncokb {
    my ($gene, $hgvsp, $lookup_data) = @_;
    my $oncogenicity = '.';
    my $effect = '.';
    my $varid = '.';
    if (exists $lookup_data->{$gene}) {
        if (exists $lookup_data->{$gene}{$hgvsp}) {
            ($oncogenicity, $effect, $varid) = @{$lookup_data->{$gene}{$hgvsp}};
        }
    }
    return ($varid, $oncogenicity, $effect);
}

sub run_nonhs_rules {
    # Look for non-hotspot MOIs in the data. Return after the first hit, even
    # though some variants may fit into more than one category.
    my ($gene, $location, $hgvs_p, $function) = @_;

    my $exon = (split(/\//, $location))[0];
    # if data coming from annovar, will not have the "exon#/total_exons" string.
    # instead will be more explicit (i.e. exon11). Need to strip the "exon"
    # string off to make it compatible.
    $exon =~ s/exon//; 

    my ($aa_start, $aa_end);
    if ($function =~ /splice/i) {
        $aa_start = $aa_end = -1 
    } else {
        ($aa_start, $aa_end) = $hgvs_p =~ /^p\.[\*A-Z]+(\d+)(?:_[A-Z]+(\d+))?.*/;
    }
    $aa_end //= $aa_start; # only get end if there is a range from indel.

    # DEBUG
    $logger->debug("incoming => gene: $gene, exon: $exon: " . 
        "aa_start: $aa_start, aa_end: $aa_end, function: $function");

    local $SIG{__WARN__} = sub {
        my $msg = shift;
        print "There is an issue with the following entry:\n";
        print "incoming =>\n\tgene: $gene\n\texon: $exon\n\tHGVS_p: $hgvs_p",
             "\n\tfunction: $function\n";
        print "Warning message: $msg\n";
    };

    # NHS Rules.
    # my ($category, $oncogenicity, $effect) = $nhs_rules->check_variant($gene, 
        # $aa_start, $aa_end, $function, $exon);
    return $nhs_rules->check_variant($gene, $aa_start, $aa_end, $function, $exon);
}

sub map_variant_hs {
    my ($chr, $maf_start, $maf_end, $ref, $alt, $hotspots) = @_;
    my $varid = '.';

    if ($hotspots->{$chr}) {
        for my $range (keys %{$hotspots->{$chr}}) {
            my ($r1, $r2) = split(/-/, $range);
            if ($maf_start >= $r1 and $maf_end <= $r2) {
                $varid = match_variant($hotspots->{$chr}{$range}, 
                    $maf_start, $maf_end, $ref, $alt);
                last;
            }
        }
    }
    return $varid;
}

sub match_variant {
    my ($vars, $maf_start, $maf_end, $ref, $alt) = @_;
    for my $var (@$vars) {
        my ($chr, $hs_start, $hs_end, $hs_ref, $hs_alt, $hsid, 
            $count) = split(/;/, $var);
        # MAF file is 1-based, while the HS BED file is 0-based. Also will have 
        # to use different position for mapping if indel compared to snv.
        my $pos;
        ($ref eq '-') ? ($pos = $maf_start) : ($pos = $maf_start-1);
        # TODO: May need to add an AA mapping algorithm here to match hotspots 
        # with different base changes.
        if ($pos == $hs_start and $ref eq $hs_ref and $alt eq $hs_alt) {
            return $hsid;
        }
    }
    return '.'
}

sub read_hs_bed {
    my $bedfile = shift;
    my @hotspots;
    my %positions;
    open(my $fh, '<', $bedfile);
    my $header = <$fh>;
    while (<$fh>) {
        my ($chr, $start, $end, $id, $allele, $amp) = split(/\t/);

        # If we have a count metric from Rajesh's PublicData file, then include it.
        my $count = '.';
        if ( $amp =~ /COUNT=(\d+)$/ ) {
            $count = $1;
        }

        my ($ref, $alt) = $allele =~ /REF=([ACTG]+)?;OBS=([ACTG]+)?(?:;.*)?/;
        $ref //= '-';
        $alt //= '-';
        
        push(@{$positions{$chr}}, $start, $end);
        push(@hotspots, [$chr, $start, $end, $ref, $alt, $id, $count]);
    }
    close $fh;
    return build_hs_table(\%positions, \@hotspots);
}

sub build_hs_table {
    my ($positions, $hotspots) = @_;
    my %ranged_hs_table;
    my $bin_width = 10000000;

    for my $chr (sort { versioncmp($a, $b) } keys %$positions) {
        my $min = min(@{$$positions{$chr}});
        my $max = max(@{$$positions{$chr}});
        my $floor = round($min, $bin_width, 'floor');

        while ($floor < $max) {
            my $range = sprintf( "%s-%s", $floor+1, $floor += $bin_width);
            $ranged_hs_table{$chr}->{$range} = [];
        }
    }

    # Insert hotspots where they belong in our hash.
    for my $var (@$hotspots) {
        if ($ranged_hs_table{$var->[0]}) {
            insert_var($ranged_hs_table{$var->[0]}, $var);
        }
    }
    return \%ranged_hs_table;
}

sub read_oncokb_file {
    my $oncokb_file = shift;
    my %data;
    open(my $fh, "<", $oncokb_file);
    my $okb_version = (split(' ', readline($fh)))[1];
    $logger->info("OncoKB lookup file version: v$okb_version.");
    my $header = <$fh>;
    while (<$fh>) {
        chomp(my @fields = split(/\t/));
        $data{$fields[0]}->{$fields[2]} = [@fields[3..5]];
    }
    return \%data;
}
    
sub insert_var {
    my ($regions, $var) = @_;
    for my $range (keys %$regions) {
        my ($start, $end) = split(/-/, $range);
        if ($var->[1] > $start and $var->[2] < $end) {
            push(@{$$regions{$range}}, join(';', @$var)) and return;
        }
    }
}

sub round {
    # Currently just a floor routine, but can add a ceil too if needed.
    my ($num, $bin, $op) = @_;
    if ($op eq 'floor') {
        return int($num/$bin) * $bin;
    }
}
    
sub now {
    my $format = shift;
    my $now = DateTime->now;
    ($format eq 'short') 
        ? return $now->strftime('%Y%m%d')
        : return $now->strftime('%Y%m%d %H:%m:%s');
}

sub __commify {
    # From Perl Cookbook.  Just for dev and testing to make visualization easier.
    my $number = reverse shift;
    $number =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $number;
}

sub __gen_hgvs {
    # Input the MAF variant string, and an optional hgvs output type, and 
    # retrieve a fully formatted HGVS annotation string. Not including the
    # output type will result in a list containg gHGVS, cHGVS, and pHGVS vals.
    # Valid return types are 'hgvs_g', 'hgvs_c', 'hgvs_p'.
    my ($var_elems, $ret_type) = @_;

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
    my @wanted_keys = qw(Hugo_Symbol Chromosome Start_Position Reference_Allele
        Tumor_Seq_Allele2 HGVSc HGVSp_Short Transcript_ID);
    my %data;

    @data{qw(gene chr start ref alt cds aa refseq)} = @$var_elems{@wanted_keys};

    # We have a splice variant or something that doesn't have a change.
    $data{'aa'} = 'p.?' unless ($data{'aa'});

    my %hgvs_annots = (
        'hgvs_g' => "$g_refseq{$data{'chr'}}:g.$data{'start'}$data{'ref'}>" .
           "$data{'alt'}",
        'hgvs_c' => "$data{'refseq'}($data{'gene'}):$data{'cds'}",
        'hgvs_p' => "$data{'refseq'}($data{'gene'}):$data{'aa'}",
    );

    ($ret_type) 
        ? return [$hgvs_annots{$ret_type}] 
        : return [@hgvs_annots{qw( hgvs_g hgvs_c hgvs_p )}];
}

sub __exit__ {
    # Better exit routine for dev and debugging purposes.
    my ($line, $msg) = @_;
    $msg //= '';
    print "\n\n";
    print colored("Got exit message at line $line with message: $msg", 
        'bold white on_green');
    print "\n";
    exit;
}

#!/usr/bin/env perl
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Sort::Versions;
use JSON -support_by_pp;
use Data::Dump;
use Term::ANSIColor;

use constant DEBUG => 0;

my $scriptname = basename($0);
my $version = "v1.0.070919-dev1";

# Remove when in prod.
#print "\n";
#print colored("*" x 75, 'bold yellow on_black'), "\n";
#print colored("      DEVELOPMENT VERSION OF $scriptname\n", 
#  'bold yellow on_black');
#print colored("*" x 75, 'bold yellow on_black');
#print "\n\n";

my $description = <<"EOT";
<<  Description  >>
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF_file(s)>
    Filter Options
    --cn  INT        Only report amplifications above this threshold (DEFAULT: 
                     CN >= 0)
    --cu  INT        Upper bound for amplifications based on 5% CI (DEFAULT: 
                     5% CI >= 4)
    --cl  INT        Lower bound for deletions based on 95% CI (DEFAULT: 
                     95% CI <= 1)

    Output Options
    -o, --outfile     Send output to custom file.  Default is STDOUT.
    -v, --version     Version information
    -h, --help        Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $copy_number = 0;
my $cu = 4;
my $cl = 1;

GetOptions( 
    "copy-number|cn=f"    => \$copy_number,
    "cu=i"                => \$cu,
    "cl=i"                => \$cl,
    "outfile|o=s"         => \$outfile,
    "version|v"           => \$ver_info,
    "help|h"              => \$help )
or die $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help() if $help;
version() if $ver_info;

# We don't need both cn and cu or cl
if ($cl or $cu) {
    undef $copy_number;
}

my %filters = (
    'cn'    => $copy_number,
    'cu'    => $cu,
    'cl'    => $cl,
);

if (DEBUG) {
    print '='x35, '  DEBUG  ', '='x35, "\n";
    print "Filters being employed\n";
    while (my ($keys, $values) = each %filters) {
        $values //= 'undef';
        if ($keys eq 'gene') {
            printf "\t%-7s => %s\n",$keys,join(',',@$values);
        } else {
            printf "\t%-7s => %s\n",$keys,$values;
        }
    }
    print '='x79, "\n";
}

# Make sure enough args passed to script
if ( @ARGV < 1 ) {
    print "ERROR: No VCF files passed to script!\n\n"; 
    print "$usage\n";
    exit 1;
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile );
} else {
	$out_fh = \*STDOUT;
}

#########--------------------- END ARG Parsing -----------------------#########
my $vcf = shift;
my %results;

my ($cnv_data, $sample_id)  = proc_vcf(\$vcf);
# dd $cnv_data;
# __exit__(__LINE__, "");

my %mapped_cnv_data;
for my $cnv (sort { versioncmp($a, $b) } keys %$cnv_data) {
    my @outfields = qw(END LEN NUMTILES RAW_CN REF_CN CN HS FUNC);
    my ($ci_5, $ci_95) = $cnv_data->{$cnv}->{'CI'} =~ /0\.05:(.*?),0\.95:(.*)$/;
    %mapped_cnv_data = map{ $_ => $cnv_data->{$cnv}->{$_} } @outfields;
    @mapped_cnv_data{qw(ci_05 ci_95)} = (sprintf("%.2f", $ci_5), 
        sprintf("%.2f", $ci_95));
    @mapped_cnv_data{qw(chr start gene)} = split(/:/, $cnv);

    $mapped_cnv_data{HS}    //= 'No'; 
    $mapped_cnv_data{ci_05} //= 0;
    $mapped_cnv_data{ci_95} //= 0;

    # Get OVAT Annot Data and add it.
    my ($gene_class, $variant_class);
    my $func = $mapped_cnv_data{FUNC};
    
    if ( $func && $func =~ /oncomine/ ) {
        my $json_annot = JSON->new->allow_singlequote->decode($func);
        my $parsed_annot = $$json_annot[0];
        $gene_class = $$parsed_annot{'oncomineGeneClass'};
        $variant_class = $$parsed_annot{'oncomineVariantClass'};
    } else {
        $gene_class = $variant_class = '---';
    }
    $mapped_cnv_data{GC} = $gene_class;
    $mapped_cnv_data{VC} = $variant_class;

    my @filtered_data = filter_results(\%mapped_cnv_data, \%filters);
    push(@{$results{$sample_id}}, \@filtered_data) if @filtered_data;
}

# dd \%results;;
# exit;
print_results(\%results);

sub print_results {
    my $data = shift;
    my $delimiter;
    my @header = qw( Sample Gender Cellularity MAPD Chr Gene Start End Length 
        Tiles Raw_CN Ref_CN CI_05 CI_95 CN Annot );

    select $out_fh;

    print join(',', @header), "\n";
    for my $sample (keys %$data) {
        my @elems = split(/:/, $sample);
        for my $cnv (@{$$data{$sample}}) {
            print join(',', @elems, @$cnv), "\n";
        }
    }
}

sub filter_results {
    # Filter out CNV data prior to printing it all out.
    my ($data, $filters) = @_;
    my @cn_thresholds = @$filters{qw(cn cu cl)};

    # OVAT Filter
    return if ($$filters{annot} and $$data{GC} eq '---');
    
    # Copy number filter
    (copy_number_filter($data, \@cn_thresholds)) 
        ? return return_data($data) 
        : return;
}

sub return_data {
    my $data = shift;
    my @fields = qw(chr gene start END LEN NUMTILES RAW_CN REF_CN ci_05 ci_95 
        CN GC);
    return @$data{@fields};
}

sub copy_number_filter {
    # if there is a 5% / 95% CI filter, use that, otherwise if there's a cn 
    # filter use that, and if there's nothing, just return it all.
    my ($data, $threshold) = @_;
    my ($cn, $cu, $cl) = @$threshold;

    if ($cn) {
        return 1 if ($$data{CN} >= $cn);
    }
    elsif ($cu) {
        return 1 if $cl and ($$data{ci_05} >= $cu || $$data{ci_95} <= $cl);
        return 1 if ($$data{ci_05} >= $cu);
    }
    elsif ($cl) {
        return 1 if $cu and ($$data{ci_05} >= $cu || $$data{ci_95} <= $cl);
        return 1 if ($$data{ci_95} <= $cl);
    } else {
        # Return everything if there are no filters.
        return 1;
    }
    return 0;
}

sub proc_vcf {
    my $vcf = shift;
    my ($sample_id, $gender, $mapd, $cellularity, $sample_name);
    my %results;

    open( my $vcf_fh, "<", $$vcf);
    while (<$vcf_fh>) {
        if ( /^##/ ) {
            if ( $_ =~ /sampleGender=(\w+)/ ) {
                $gender = $1 and next;
            }
            # Need to add to accomodate the new CNV plugin; may not have the 
            # same field as the normal IR data.
            if ($_ =~ /AssumedGender=([mf])/) {
                ($1 eq 'm') ? ($gender='Male') : ($gender='Female');
                next;
            }
            elsif ( $_ =~ /mapd=(\d\.\d+)/ ) {
                $mapd = $1 and next;
            }
            elsif ( $_ =~ /CellularityAsAFractionBetween0-1=(.*)$/ ) {
                $cellularity = $1 and next;
            }
        } 

        my @data = split;
        if ( $data[0] =~ /^#/ ) {
            $sample_name = $data[-1] and next;
        }
        $sample_id = join( ':', $sample_name, $gender, $cellularity, $mapd );
        next unless $data[4] eq '<CNV>';
        my $varid = join( ':', @data[0..2] );
        
        # Only get the "hotspot" entries since they are guaranteed to have
        # enough amplicons to get an accurate call.  Also, fix the field to
        # match the rest with the format of Key=Value.
        next unless $data[7] =~ /HS/;
        $data[7] =~ s/HS/HS=Yes/;

        # New v5.2 standard deviation value; sometimes data in this field and 
        # sometimes not. 
        $data[7] =~ s/SD;/SD=NA;/; 

        my @format = split( /;/, $data[7] );
        my ($cn) = $data[9] =~ /:([^:]+)$/;
        push( @format, "CN=$cn" );

        %{$results{$varid}} = map { split /=/ } @format;
    }
    if (DEBUG) {
        print "="x35, "  DEBUG  ", "="x35, "\n";
        print "\tSample Name:  $sample_name\n";
        print "\tCellularity:  $cellularity\n";
        print "\tGender:       $gender\n";
        print "\tMAPD:         $mapd\n";
        print "="x79, "\n";
    }
    return \%results, $sample_id;
}

sub __exit__ {
    my ($line, $msg) = @_;
    print "\n\n";
    print colored("Got exit message at line: $line with message:\n$msg", 
        'bold white on_green');
    print "\n";
    exit;
}

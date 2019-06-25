#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Data::Dump;
use Getopt::Long qw(:config bundling auto_abbrev no_ignore_case);
use File::Basename;

my $scriptname = basename($0);
my $version = 'v0.4.072718';
my $description = <<"EOT";
Read in an OncoKB allAnnotations file and output a lookup file that can be used
with the TSO500 annotator to annotate variants that are derived from the 
pipeline.  As a part of the development process, we can either output those that
would be included in the lookup file, or those that would be excluded, so that
we can see what should be considered for non-hotspot rules.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <OncoKB File>
    -e, --excluded   Print the excluded variants instead of the Hotspots. 
    -o, --outfile    Write output to file instead of STDOUT.
    -v, --version    Print the version information and exit.
    -h, --help       Print this help message and exit.
EOT

my ($help, $ver_info, $outfile);
my $excluded;

GetOptions(
    "excluded|e"     => \$excluded,
    "outfile|o=s"    => \$outfile,
    "version|v"      => \$version,
    "help|h"         => \$help,
) or die $usage;

sub help {
    print "$scriptname - $version\n\n$description\n\n$usage\n";
    exit;
}

sub version {
    print "$scriptname - $version\n";
    exit;
}

help if $help;
version if $ver_info;

if (@ARGV < 1) {
    print "ERROR: Not enough arguments passed to script!\n\n";
    print "$usage\n";
    exit 1;
}
my $oncokb_file = shift;

my $outfh;
if ($outfile) {
    print "Writing output to $outfile.\n";
    open($outfh, ">", $outfile);
} else {
    $outfh = \*STDOUT;
}

################------ END Arg Parsing and Script Setup ------#################

open(my $fh, "<", $oncokb_file);
chomp(my @header = split(/\t/, readline($fh)));

my @data;
my @excluded;
my @wanted_fields = qw(3 1 5 6 7);

while (<$fh>) {
    my @fields = split(/\t/);

    # Get rid of variants that aren't strong oncogenic (i.e. Neutral or
    # Inconclusive). There are still 3 "oncogenic" variants left with odd
    # mutations effects, but leave them in there for now.
    next if $fields[6] =~ /(neutral|inconclusive)/i;

    if ($fields[5] !~ /^[A-Z]\d+(:?[_\*A-Z]*|del|ins)?/) {
        push(@excluded, [@fields[@wanted_fields]]) and next;
    } else {
        # Still a few other things left to get rid of
        my @unwanted = qw(fusion trunc mut);
        
        # NOTE: MAF file has "X<AA>_splice" as notation for splice variants
        # in HGVSp_short column.
        next if grep { $fields[5] =~ /$_/i } @unwanted;

        $fields[5] = "p.$fields[5]";
        push(@data, [@fields[@wanted_fields]]);
    }
}
close $fh;

# Manually curated variants
# There are some that we want to add to the hotspots list that are not included
# in the OncoKB file by default.  Add them in.
my $manual_vars = [
    ['MET', 'NM_000245.2', 'p.X1010_splice', 'Likely Oncogenic', 
        'Likely Gain-of-Function'],
    ['RIT1', 'NM_006912.5', 'p.TA83del', 'Oncogenic', 'Gain-of-function']
];
push(@data, $_) for @$manual_vars;

# Print the header
print {$outfh} join("\t", @header[@wanted_fields], 'varid'), "\n";

# Sort it so that we can read it easier when we're looking for stuff.
my @sorted_data = sort { $a->[0] cmp $b->[0] || $a->[-1] cmp $b->[-1] } @data;

# Print the data.
if ($excluded) {
    print {$outfh} join("\t", @$_), "\n" for @excluded;
} else {
    print {$outfh} join("\t", @$_, "$_->[1]($_->[0]):$_->[2]"), "\n" for @sorted_data;
}

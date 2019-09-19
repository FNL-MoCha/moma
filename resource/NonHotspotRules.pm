# Non-hotspot rules object for MOMA. Will require a JSON blob of rules that can
# then be traversed to determine if an incoming variant is a non-hotspot rule
# oncogenic variant or not.
#
# 9/12/2019 - D Sims
################################################################################
package NonHotspotRules;

use strict;
use warnings;
use Carp;

use JSON;
use Data::Dump;  # TODO: Go away!

our $VERSION='0.2.091219';

sub new {
    my $class = shift;
    my $rule_set = shift;

    my $self = {};
    bless $self, $class;

    $self->__load_rules($rule_set);

    return $self;
}

sub __load_rules {
    my ($self, $rules_json) = @_;
    my $json;
    { 
        local $/;
        open(my $fh, "<", $rules_json);
        $json = <$fh>;
        close $fh;
    }
    my $data = JSON->new->utf8->decode($json);
    $self->{_version}        = $data->{"metadata"}{"version"};
    $self->{_oncokb_version} = $data->{"metadata"}{"oncokb_version"};
    $self->{rules}           = $data->{"rules"};
}

sub check_variant {
    # TODO: Do we want to use this as the main method entry point, which will
    # then branch out to check by aa, check by exon, and check by raw function?
    my ($self, $gene, $aa_start, $aa_end, $function, $exon) = @_;


    return $self->codon_in_range($gene, $aa_start, $aa_end, $function, $exon);
}

sub codon_in_range {
    ############################################################################
    # XXX TODO:
    #     - This method seems to work well.  See if we can simplify it further,
    #       along with cleaning up some redundancy (especially in the terms to be
    #       returned.
    #
    #     - Figure out where to consider the functional effect. We're already
    #       parsing rules here, so maybe I ougth to handle it here rather than
    #       waiting for later and having to re-run.
    #
    #     - What do we do with a null result? Need to figure out how we want to
    #       handle those cases.
    ############################################################################
    
    # Run the rules, returning either the category, oncogenicity, and effect, or
    # a list of empty strings if no match.
    # TODO: Figure out exactly what we want in this args list.
    my ($self, $gene, $aa_start, $aa_end, $function, $exon) = @_;
    my ($category, $oncogenicity, $effect) = ('', '', '');

    if (exists $self->{'rules'}{$gene}) {
        while (my ($index, $rule) = each @{$self->{'rules'}{$gene}}) {
            if ($rule->{'aa_start'} or $rule->{'aa_end'}) {
                # If only given a start bound or end bound, need to set the end
                # bound to a very high number or start bound to 0 (respectively)
                # just to have values for comparison.
                $rule->{'aa_start'} = '>0' unless $rule->{'aa_start'};
                $rule->{'aa_end'}   = '<100000' unless $rule->{'aa_end'};

                if ($aa_start == $aa_end) {
                    # We have a substitution call; only need to check once.
                    if ($self->__comp($aa_start, $rule->{'aa_start'}) 
                        && $self->__comp($aa_start, $rule->{'aa_end'})) {
                        # TODO: Check func before deciding what to report.
                        
                        if ($self->check_function($function, $rule->{'function'})) {
                            ($category, $oncogenicity, $effect) = 
                                @{$rule}{qw(category oncogenicity effect)};
                        }
                    }
                } else {
                    # We have an indel, and need to check the range.
                    my ($start_in_range, $end_in_range);
                    $start_in_range = $self->__comp($aa_start, $rule->{'aa_start'}) 
                        && $self->__comp($aa_start, $rule->{'aa_end'});

                    $end_in_range = $self->__comp($aa_end, $rule->{'aa_start'})
                        && $self->__comp($aa_end, $rule->{'aa_end'});

                    if ($start_in_range && $end_in_range) {
                        # The whole indel within one rule. Return match data.
                        # TODO: Check func before deciding what to report.
                        return @{$rule}{qw(category oncogenicity effect)};
                    }
                    elsif ($start_in_range || $end_in_range) {
                        # Only partial match to one rule. Check other rule and
                        # report either partial to initial rule or partial to
                        # both rules (assuming we only have 2 rules per gene).
                        $oncogenicity = "Possibly Oncogenic";

                        if (($index+1) < @{$self->{'rules'}{$gene}}) {
                            my $next_rule = $self->{'rules'}{$gene}[$index+1];
                            if ($start_in_range) {
                                if ($self->__comp($aa_end, $next_rule->{'aa_end'})) {
                                    $category = sprintf("(Partial) %s/%s",
                                        $rule->{'category'}, 
                                        $next_rule->{'category'});
                                    $effect = "Unknown";
                                    # TODO: check func before deciding what to
                                    # report.
                                    return $category, $oncogenicity, $effect;
                                }
                            }
                        }
                        # We are only partially within one of the two rules.
                        $category = "(Partial) $rule->{'category'}";
                        $effect = "(Possible) $rule->{'effect'}";
                        # TODO: Check func before deciding what to report.
                        return $category, $oncogenicity, $effect;
                    }
                }
            }
        }
    }
    return $category, $oncogenicity, $effect;
}

sub check_function {
    # Check to see if the functional annotation matches rule.
    my ($self, $function, $rule_func) = @_;
    (grep $function =~ /$_/, @$rule_func) ? return 1 : return 0;
    # print "function: $function\n";
    ## Check to see if function in the list as the last check.
    # if (grep $function =~ /$_/, @{$rule->{'function'}}) {
        # return ($rule->{'category'}, $rule->{'oncogenicity'},
            # $rule->{'effect'});
}

sub __comp {
    my ($self, $test_val, $rule) = @_;
    my ($sigil, $v2) = $rule =~ /^([<>=]+)(\d+)/;

    my %mapping = (
        "<"  => sub { $_[0] < $_[1] },
        ">"  => sub { $_[0] > $_[1] },
        "<=" => sub { $_[0] <= $_[1] },
        ">=" => sub { $_[0] >= $_[1] },
    );

    ($mapping{$sigil}->($test_val, $v2)) 
        ? return 1
        : return 0;
}

1;

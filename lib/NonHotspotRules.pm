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
use Data::Dump;
use File::Basename;

our $VERSION='0.7.092419';

sub new {
    my ($class, $rule_set, $tsg_file) = @_;

    # Set defaults if one is not passed.
    my $__default_rules = dirname(__FILE__) . "/../resource/non-hotspot_rules.json";
    my $__default_tsgs  = dirname(__FILE__) . "/../resource/gene_reference.csv";

    if (! defined $rule_set) {
        warn "WARN: Using default Non-Hotspot Rules JSON: " . 
            basename($__default_rules) . "\n";
        $rule_set = $__default_rules;
    }
    if (! defined $tsg_file) {
        warn "WARN: Using default TSGs File: " . basename($__default_tsgs) . "\n";
        $tsg_file = $__default_tsgs;
    }

    my $self = {};
    bless $self, $class;

    $self->__load_rules($rule_set);
    $self->__load_tsgs($tsg_file);

    return $self;
}

sub get_category_list {
    # Get a list of all non-hs rule categories for tallying later.
    my $self = shift;
    my %categories;

    for my $gene (keys %{$self->{'rules'}}) {
        for my $rule (@{$self->{'rules'}{$gene}}) {
            $categories{$rule->{'category'}} = 0;
        }
    }
    @categories{('Hotspots', 'Truncating in TSG')} = (0,0);
    return \%categories;
}

sub check_variant {
    # Main entry point for checking non-hotspot rules. Determine if there is a
    # rule for the incoming gene, and if so, determine whether it's an exon
    # location rule or an amino acid rule; it can only either be an exon or
    # amino acid rule, but not none or both.  If either one of those is a match,
    # return the matching rule, and then make sure the function matches.  If so,
    # report the category, oncogenicity, and effect of the variant.
    my ($self, $gene, $aa_start, $aa_end, $function, $exon) = @_;

    my @tsg_funcs = qw(nonsense stop frameshift splice);
    my ($category, $oncogenicity, $effect) = ('.', '.', '.');

    my $tsg_flag = 1 if ($self->__is_tsg($gene));

    if (exists $self->{'rules'}{$gene}) {
        my $selected_rule;
        my $gene_rules = $self->{'rules'}{$gene};
        # Special TSG or Oncogene
        # Check to see if we need to evaluate AA or Exon position.
        for (@$gene_rules) {
            if ($_->{'aa_start'} or $_->{'aa_end'}) {
                $selected_rule = $self->codon_in_range($aa_start, $aa_end, 
                    $gene_rules);
            }
            elsif ($_->{'exon'}) {
                $selected_rule = $_ if $self->check_exon($exon, $_);
            }
            last if $selected_rule;
        }

        # We have a TSG with dual gene specific and generic TSG rules, but
        # variant didn't fit gene specific rule; check functional rules for
        # generic TSG match.
        if (! $selected_rule && $tsg_flag) {
            if ($self->check_function($function, \@tsg_funcs)) {
                return 'Truncating in TSG', 'Oncogenic', 'Loss-of-Function';
            }
        }
        # If no rule, return the null strings, otherwise continue checking.
        $selected_rule or return $category, $oncogenicity, $effect;


        # Function Check
        if (@{$selected_rule->{'function'}}) {
            if (! $self->check_function($function, $selected_rule->{'function'})) {
                return $category, $oncogenicity, $effect;
            }
        } 
        # If we made it here, we selected a rule, but there are no
        # functional annotations to match. Return as is.
        $category = $selected_rule->{'category'};
        $oncogenicity = $selected_rule->{'oncogenicity'};
        $effect = $selected_rule->{'effect'};

        return $category, $oncogenicity, $effect;
    } 
    elsif ($tsg_flag) {
        # Just a plain TSG
        if ($self->check_function($function, \@tsg_funcs)) {
            return ('Truncating in TSG', 'Oncogenic', 'Loss-of-Function');
        }
    }
    return $category, $oncogenicity, $effect;
}

sub codon_in_range {
    # Check to see if the variant is within an amino acid range, and return the
    # rule that matches (if one does!), which we will then check to see if the
    # rest of the conditions match.
    my ($self, $aa_start, $aa_end, $gene_rules) = @_;

    for my $index (0..$#{$gene_rules}) {
        my $rule = $gene_rules->[$index];

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
                    return $rule;
                }
            } else {
                # We have an indel, and need to check the range.
                my ($start_in_range, $end_in_range);
                $start_in_range = $self->__comp($aa_start, $rule->{'aa_start'}) 
                    && $self->__comp($aa_start, $rule->{'aa_end'});

                $end_in_range = $self->__comp($aa_end, $rule->{'aa_start'})
                    && $self->__comp($aa_end, $rule->{'aa_end'});

                if ($start_in_range && $end_in_range) {
                    return $rule;
                }
                elsif ($start_in_range || $end_in_range) {
                    # Only partial match to one rule. Check other rule and
                    # report either partial to initial rule or partial to
                    # both rules depending on what matches.
                    if (($index+1) < @$gene_rules) {
                        my $next_rule = $gene_rules->[$index+1];
                        if ($start_in_range) {
                            if ($self->__comp($aa_end, $next_rule->{'aa_end'})) {
                                my $new_category = sprintf("(Partial) %s/%s",
                                    $rule->{'category'}, 
                                    $next_rule->{'category'});
                                # Functional rules should be the same, but just
                                # in case, get a union of the two and use that.
                                my %union = ();
                                for (
                                    @{$rule->{'function'}}, 
                                    @{$next_rule->{'function'}}
                                ) { $union{$_}++ };

                                # Create a new rule to parse downstream.  
                                my $new_rule = {
                                    "exon" => [],
                                    "aa_start" => "$rule->{'aa_start'}&$next_rule->{'aa_start'}",
                                    "aa_end" => "$rule->{'aa_end'}&$next_rule->{'aa_end'}",
                                    "function" => [keys %union],
                                    "category" => $new_category,
                                    "oncogenicity" => "Possibly Oncogenic",
                                    "effect" => "Unknown"
                                };
                                return $new_rule;
                            }
                        }
                    }

                    # We are only partially within one of the two rules.
                    $rule->{'category'} = "(Partial) $rule->{'category'}";
                    $rule->{'effect'} = "(Possible) $rule->{'effect'}";
                    return $rule;
                }
            }
        }
    }
    # If we got here, then no rules match, and we'll return False.
    return 0;
}

sub check_function {
    # Check to see if the functional annotation matches rule.
    my ($self, $function, $rule_func) = @_;
    (grep $function =~ /$_/i, @$rule_func) ? return 1 : return 0;
}

sub check_exon {
    # Check to see if the exon matches the rule.
    my ($self, $exon, $rule_exons) = @_;
    # Technically should be using numeric comparison op, but since we have TERT
    # promoter mutations and will get the '---' string in place of exon
    # location, use the string compare op instead. Should still work fine for
    # numbers since 1:1 comprison.
    (grep $exon eq $_, @{$rule_exons->{'exon'}}) ? return 1 : return 0;
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
    $self->{version}        = $data->{"metadata"}{"version"};
    $self->{oncokb_version} = $data->{"metadata"}{"oncokb_version"};
    $self->{rules}           = $data->{"rules"};
}

sub __load_tsgs {
    my ($self, $tsg_file) = @_;

    open(my $fh, "<", $tsg_file);
    $self->{'tsg_version'} = (split(' ', readline($fh)))[1];
    readline($fh); # throw out header...don't need it.
    $self->{'tsgs'}{'genes'} = [map { $_ =~ /Yes$/ ? ((split(/,/, $_))[3]) : () } <$fh>];
}

sub __is_tsg {
    my ($self, $gene) = @_;
    (grep { $gene eq $_ } @{$self->{'tsgs'}{'genes'}}) ? return 1 : return 0;
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

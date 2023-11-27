#!/usr/bin/perl 
use strict;

my $TE = shift @ARGV or die$!;

unless (-s "${TE}_20bp_CasOFFinder.output.txt"){    &Cas_OFFinder ($TE, 20); }

sub Cas_OFFinder {
    my ($title, $seed_length) = @_;
    system ("cas-offinder ${title}_${seed_length}bp_CRISPOR.txt G0 ${title}_${seed_length}bp_CasOFFinder.output.txt");
}


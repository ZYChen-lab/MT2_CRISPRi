#!/usr/bin/perl -w 

use strict;
use Data::Dumper;

my $CRISPOR_output = shift @ARGV or die$!;
my $TE = shift @ARGV or die$!;

#process CRISPOR output (20bp seed sequences) -> generate input for cas-offinder
#also return hash_ref (key = sgRNA sequence; value = sg_id)
my $sg20 = &process_CRISPOR_output4Cas_OFFinder ($TE, $CRISPOR_output, 20, 5); 

sub process_CRISPOR_output4Cas_OFFinder{

	#input as CRISPOR_output (csv)
	#output sgRNA fasta ${title}_${seedlength}_CRISPOR.txt (which serve as input for cas-offinder
    	#Also return hash_reference (key = sequence_without_PAM; value = sgRNA_id)
    	#ID (guideID,Doench_score,Moreno-Mateos_score,Graf_motif)
    	#options: 20bp guide

	my ($title, $file, $seed_length, $mismatch) = @_;
    	my %sgRNA;
	
	open(IN, $file) or die$!;
    	open(OUT, ">${title}_${seed_length}bp_CRISPOR.txt") or die$!;	

	#prepare header
	print OUT "/data/ZYChenlab/Zhiyuan/genomes_annotations/mouse/mm10/genomes\n";
	print OUT "N" x $seed_length . "NGG\t0\t0\n";

	while (my $line = <IN>){
        chomp $line;
        next if $line =~ /(^\#)|(^,)/;
	next if $line =~ /^</;

        my @split = split(/,/, $line);
        my $id = $split[0];

	#retrieve seed sequence & add "NGG" PAM
	my $seq = $split[1];
	my $gRNA = substr $seq, 20-$seed_length, $seed_length;

	my $Doench = $split[6];
        $Doench =~ s/NotEnoughFlankSeq/NA/g;

        my $Moreno = $split[13];
        $Moreno =~ s/NotEnoughFlankSeq/NA/g;

        my $GrafEtAlStatus = $split[18];
        next if $GrafEtAlStatus =~ "tt";
        $GrafEtAlStatus =~ s/GrafOK/1/g;

	$id = $id . "," . $Doench . "," . $Moreno . "," . $GrafEtAlStatus;
        $gRNA = $gRNA . "NNN";

        print OUT "$gRNA\t$mismatch\t$id\n";
        $sgRNA{$gRNA} = $id;
	}
	close IN;
    	close OUT;
    	return (\%sgRNA);
}
		

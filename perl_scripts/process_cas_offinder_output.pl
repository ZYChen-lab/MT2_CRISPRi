#!/usr/bin/perl 

use strict;

my $CRISPOR_output = shift @ARGV or die$!;
my $TE = shift @ARGV or die$!;

#simply to get sgRNA_id:
my $sg20 = &get_id  ($TE, $CRISPOR_output, 20, 3);;

#correct cas-offinder mismatch:
&recalculate_mismatch("${TE}_20bp_CasOFFinder.output.txt", "${TE}_20bp_CasOFFinder.output.corrected.txt");

#Convert cas-offinder output to BED format
&Cas_OFFinder_to_BED ("${TE}_20bp_CasOFFinder.output.corrected.txt", $TE, 20);

#IntersectBed to TE.GTF
unless (-s "${TE}_20bp_intersect_TE.txt" ) { &IntersectBed ( "${TE}_20bp_CasOFFinder.bed", $TE, 20); }

#Output sgRNA on-target&off-target information
&sgRNA_on_off_Target($TE, 20);

sub recalculate_mismatch{
	
	my $file = shift @_;
	my $out = shift @_;

	open (FILE, $file) or die$!;
	open (OUT, ">$out") or die$!;
 
	while (my $line = <FILE>){
		chomp $line;
		#TTCAAACCAAGATCCAGATANNN chr1 1  10982418        TTCAAACCAAGATCCAGATAAGGATCT     -       0
		my @split = split (/\s+/, $line);
		my $seq1 = substr $split[0], 0, 20;
		my $seq2 = substr $split[4], 0, 20;		
	
		$seq1  = uc $seq1;
		$seq2  = uc $seq2;
		my @seq1 = split(//, $seq1);
		my @seq2 = split(//, $seq2);
		my $mismatch = 0;
		foreach(my $i = 0; $i<=$#seq1; $i++){
			if($seq1[$i] ne $seq2[$i]){
				$mismatch++;
			}
		}
		my $tmp = pop @split;
		my $print = join("\t", @split);
		print OUT $print  . "\t" . $mismatch . "\n";
	}
	close FILE;
	close OUT;
}

sub sgRNA_on_off_Target{
    my ($title, $seed_length) = @_;
    
    my $all_TE = &TE_id_length($title);
    my ($onTarget_mis0, $onTarget_mis1, $onTarget_mis2, $onTarget_mis3, $offTarget_mis0, $offTarget_mis1, $offTarget_mis2, $offTarget_mis3) = &ProcessIntersectBed ( "${title}_${seed_length}bp_intersect_TE.txt", $TE, $seed_length);

    #ouput on target table
    open (OUT0, ">${TE}_${seed_length}bp_sgRNA_onTarget_mis0_count.txt") or die$!;
    open (OUT1, ">${TE}_${seed_length}bp_sgRNA_onTarget_mis1_count.txt") or die$!;
    open (OUT2, ">${TE}_${seed_length}bp_sgRNA_onTarget_mis2_count.txt") or die$!;
    open (OUT3, ">${TE}_${seed_length}bp_sgRNA_onTarget_mis3_count.txt") or die$!;
   
    #col = sgRNAs; row = all MT2_Mm copies 
    #print header;
    print OUT0 "LTR_id\tLength\t";    print OUT1 "LTR_id\tLength\t";    print OUT2 "LTR_id\tLength\t";    print OUT3 "LTR_id\tLength\t";
    foreach my $sg_id (sort keys %{$onTarget_mis0} ){
	print OUT0 $sg_id . "\t"; print OUT1 $sg_id . "\t"; print OUT2 $sg_id . "\t"; print OUT3 $sg_id . "\t";
    }
    print OUT0 "\n"; print OUT1 "\n"; print OUT2 "\n"; print OUT3 "\n";

    foreach my $te_id (sort keys %{$all_TE}){
	print OUT0 $te_id . "\t" . $all_TE->{$te_id} . "\t";
	print OUT1 $te_id . "\t" . $all_TE->{$te_id} . "\t";
	print OUT2 $te_id . "\t" . $all_TE->{$te_id} . "\t";
	print OUT3 $te_id . "\t" . $all_TE->{$te_id} . "\t";

	foreach my $sg_id (sort keys %{$onTarget_mis0}){
	    unless ($onTarget_mis0->{$sg_id}->{$te_id}) { $onTarget_mis0->{$sg_id}->{$te_id}=0; }
	    unless ($onTarget_mis1->{$sg_id}->{$te_id}) { $onTarget_mis1->{$sg_id}->{$te_id}=0; }
	    unless ($onTarget_mis2->{$sg_id}->{$te_id}) { $onTarget_mis2->{$sg_id}->{$te_id}=0; }
	    unless ($onTarget_mis3->{$sg_id}->{$te_id}) { $onTarget_mis3->{$sg_id}->{$te_id}=0; }

     	    #if sgRNA targets with 0 mismatch, it should also be counted into 1, 2, 3 mismatch
     	    if ($onTarget_mis0->{$sg_id}->{$te_id} == 1) { $onTarget_mis1->{$sg_id}->{$te_id}=1; $onTarget_mis2->{$sg_id}->{$te_id} = 1; $onTarget_mis3->{$sg_id}->{$te_id} = 1; }
     	    if ($onTarget_mis1->{$sg_id}->{$te_id} == 1) { $onTarget_mis2->{$sg_id}->{$te_id}=1; $onTarget_mis3->{$sg_id}->{$te_id} = 1; }
     	    if ($onTarget_mis2->{$sg_id}->{$te_id} == 1) { $onTarget_mis3->{$sg_id}->{$te_id}=1; }	
	    
	    print OUT0 $onTarget_mis0->{$sg_id}->{$te_id} . "\t";
	    print OUT1 $onTarget_mis1->{$sg_id}->{$te_id} . "\t";
	    print OUT2 $onTarget_mis2->{$sg_id}->{$te_id} . "\t";
	    print OUT3 $onTarget_mis3->{$sg_id}->{$te_id} . "\t";
	}
	print OUT0 "\n"; print OUT1 "\n"; print OUT2 "\n"; print OUT3 "\n";
    }  
    close OUT0; close OUT1; close OUT2; close OUT3;

    #output off target table -> summarize all TE types of off-target for each level of mismatches
    my (%offtarget_TE_mis0, %offtarget_TE_mis1, %offtarget_TE_mis2, %offtarget_TE_mis3);
    foreach my $sg_id (sort keys %{$offTarget_mis0} ){
	foreach my $te_name (sort keys %{ $offTarget_mis0->{$sg_id} } ){
	    $offtarget_TE_mis0{$te_name}->{$sg_id} = $offTarget_mis0->{$sg_id}->{$te_name};
	}
    }
    foreach my $sg_id (sort keys %{$offTarget_mis1} ){
        foreach my $te_name (sort keys %{ $offTarget_mis1->{$sg_id} } ){
            $offtarget_TE_mis1{$te_name}->{$sg_id} = $offTarget_mis1->{$sg_id}->{$te_name};
        }
    }
    foreach my $sg_id (sort keys %{$offTarget_mis2} ){
        foreach my $te_name (sort keys %{ $offTarget_mis2->{$sg_id} } ){
            $offtarget_TE_mis2{$te_name}->{$sg_id} = $offTarget_mis2->{$sg_id}->{$te_name};
        }
    }
    foreach my $sg_id (sort keys %{$offTarget_mis3} ){
        foreach my $te_name (sort keys %{ $offTarget_mis3->{$sg_id} } ){
            $offtarget_TE_mis3{$te_name}->{$sg_id} = $offTarget_mis3->{$sg_id}->{$te_name};
        }
    }

    open (OUF0, ">${TE}_${seed_length}bp_sgRNA_offTarget_mis0_count.txt") or die$!;
    open (OUF1, ">${TE}_${seed_length}bp_sgRNA_offTarget_mis1_count.txt") or die$!;
    open (OUF2, ">${TE}_${seed_length}bp_sgRNA_offTarget_mis2_count.txt") or die$!;
    open (OUF3, ">${TE}_${seed_length}bp_sgRNA_offTarget_mis3_count.txt") or die$!;

    #print header;
    print OUF0 "sg_id\t"; print OUF1 "sg_id\t"; print OUF2 "sg_id\t"; print OUF3 "sg_id\t";
    foreach my $te_name (sort keys %offtarget_TE_mis0 ){ print OUF0 "$te_name\t"; }
    foreach my $te_name (sort keys %offtarget_TE_mis1 ){ print OUF1 "$te_name\t"; }
    foreach my $te_name (sort keys %offtarget_TE_mis2 ){ print OUF2 "$te_name\t"; }
    foreach my $te_name (sort keys %offtarget_TE_mis3 ){ print OUF3 "$te_name\t"; }
    print OUF0 "\n"; print OUF1 "\n"; print OUF2 "\n"; print OUF3 "\n";

    foreach my $sg_id (sort keys %{$onTarget_mis0} ){
	print OUF0 $sg_id . "\t"; print OUF1 $sg_id . "\t"; print OUF2 $sg_id . "\t"; print OUF3 $sg_id . "\t";
	foreach my $te_name (sort keys %offtarget_TE_mis0){
	    unless ( $offtarget_TE_mis0{$te_name}->{$sg_id} ) { $offtarget_TE_mis0{$te_name}->{$sg_id} = 0; } 
	    print OUF0 $offtarget_TE_mis0{$te_name}->{$sg_id} . "\t"; }
	
	foreach my $te_name (sort keys %offtarget_TE_mis1){
            unless ( $offtarget_TE_mis1{$te_name}->{$sg_id} ) { $offtarget_TE_mis1{$te_name}->{$sg_id} = 0; } 
	    print OUF1 $offtarget_TE_mis1{$te_name}->{$sg_id} . "\t"; }

	foreach my $te_name (sort keys %offtarget_TE_mis2){
            unless ( $offtarget_TE_mis2{$te_name}->{$sg_id} ) { $offtarget_TE_mis2{$te_name}->{$sg_id} = 0; } 
	    print OUF2 $offtarget_TE_mis2{$te_name}->{$sg_id} . "\t"; }

	foreach my $te_name (sort keys %offtarget_TE_mis3){
            unless ( $offtarget_TE_mis3{$te_name}->{$sg_id} ) { $offtarget_TE_mis3{$te_name}->{$sg_id} = 0; } 
	    print OUF3 $offtarget_TE_mis3{$te_name}->{$sg_id} . "\t";
	}
	print OUF0 "\n"; print OUF1 "\n"; print OUF2 "\n"; print OUF3 "\n";
    }
    close OUF0; close OUF1; close OUF2; close OUF3;
}

sub ProcessIntersectBed {
    my ($file, $title, $seed_length) = @_;
    
    #make two output: 1) On-target of sgRNA on each $TE; 2) Off-target counts for non-$TE
    open (FILE, "$file") or die$!;

    my (%onTarget_mis0, %onTarget_mis1, %onTarget_mis2, %onTarget_mis3);
    my (%offTarget_mis0, %offTarget_mis1, %offTarget_mis2, %offTarget_mis3);

    while (my $line = <FILE>){
	chomp $line;
	my @split = split ("\t", $line);	
	
	#sgRNA id
	my $sg_id = $split[3];
	my @tmp = split (/,/, $sg_id); pop(@tmp); $tmp[-1] =~ s/N//g;
	$sg_id = join (",", @tmp);

	#TE id #chr:pos1-pos2
	my $te_id = $split[6] . ":" . $split[9] . "-" . $split[10];

	#TE name
	my @split2 = split (/\s+/, $split[-1]);
	my $te_name  = $split2[1];
	$te_name =~ s/[";]//g;

	#on target record
	if ( $te_name eq $title ){ 
	    if ( $split[4] == 0 ){ $onTarget_mis0{$sg_id}->{$te_id} = 1;  }
	    if ( $split[4] == 1 ){ $onTarget_mis1{$sg_id}->{$te_id} = 1;  }
	    if ( $split[4] == 2 ){ $onTarget_mis2{$sg_id}->{$te_id} = 1;  }
	    if ( $split[4] == 3 ){ $onTarget_mis3{$sg_id}->{$te_id} = 1;  }
	}
	
	#off target record
	if ( $te_name ne $title ){ 
	    if ( $split[4] == 0 ){ $offTarget_mis0{$sg_id}->{$te_name}++; }
	    if ( $split[4] == 1 ){ $offTarget_mis1{$sg_id}->{$te_name}++; }
	    if ( $split[4] == 2 ){ $offTarget_mis2{$sg_id}->{$te_name}++; }
	    if ( $split[4] == 3 ){ $offTarget_mis3{$sg_id}->{$te_name}++; }
	}	
    }   
    return (\%onTarget_mis0, \%onTarget_mis1, \%onTarget_mis2, \%onTarget_mis3, \%offTarget_mis0, \%offTarget_mis1, \%offTarget_mis2, \%offTarget_mis3);
}
	
sub TE_id_length {
    #return a hash reference: key = TE_id (chr_pos); value = length (bp)
    my $title = shift @_;
    my $GTF = "/data/ZYChenlab/Zhiyuan/genomes_annotations/mouse/mm10/annotations/mm10_rmsk_TE.gtf";
    my $grep = "\'gene_id " . '"' . $title . '"\'';

    system("grep $grep $GTF > mm10_${title}.gtf");    
    open (GTF, "mm10_${title}.gtf") or die$!;
    my %hash;
    while (my $line = <GTF>){
	chomp $line;
	my @split = split (/\s+/, $line);
	my $id = $split[0] . ":" . $split[3] . "-" . $split[4];
	my $length = $split[4] - $split[3] + 1;
	$hash{$id} = $length;
    }
    return (\%hash);
}

sub IntersectBed {
    my ($bed, $title, $seed_length) =  @_;
    my $TE_GTF = "/data/ZYChenlab/Zhiyuan/genomes_annotations/mouse/mm10/annotations/mm10_rmsk_TE.gtf";
    system("bedtools intersect -a $bed -b $TE_GTF -wa -wb > ${title}_${seed_length}bp_intersect_TE.txt");
}


sub Cas_OFFinder_to_BED {
    #convert cas-offinder output to bed
    #chr  start  end  Name(sgRNAid) Score(#_of_mismatches)  Strand
          
    my ($inFile, $title, $seed_length) = @_;

    open (BED, ">${title}_${seed_length}bp_CasOFFinder.bed") or die$!;
    open (INFILE, "$inFile") or die$!;

    while (my $line = <INFILE>){
    	chomp $line;
    	my ($sgRNA, $chr, $chr2, $start, $ref, $strand, $mismatch) = split (/\s+/, $line);
    	my $end;
    	
	if ($strand eq '+'){ $end = $start + $seed_length - 1; }
    	if ($strand eq '-'){ $start = $start + 7; $end = $start + $seed_length - 1; }
    	if ($seed_length == 20){    print BED "$chr\t$start\t$end\t$sg20->{$sgRNA},$sgRNA,$ref\t$mismatch\t$strand\n"; }        
    }
    close INFILE;
    close BED;
}

sub get_id{

	#input as CRISPOR_output (csv)
	#Also return hash_reference (key = sequence_without_PAM; value = sgRNA_id)
	#ID (guideID,Doench_score,Moreno-Mateos_score,Graf_motif)
	#options: 20bp guide

	my ($title, $file, $seed_length, $mismatch) = @_;
    	my %sgRNA;
	
	open(IN, $file) or die$!;
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

        	$sgRNA{$gRNA} = $id;
	}
	close IN;
    	return (\%sgRNA);
}
	

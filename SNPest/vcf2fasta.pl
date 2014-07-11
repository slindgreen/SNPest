#!/usr/bin/perl -w use strict;
use warnings;

my $reffile=$ARGV[0];
my $insertions=0;
my $deletions=0;
my $substitutions=0;
my $seqpositions=0;
my $refgenomelength=-1;

#print STDERR $reffile."\n";

my $INFILE;
my $temp;
open(INFILE,$reffile) or die "Cannot read file $reffile\n";
# First line is FASTA header
$temp=<INFILE>;
#    print $temp;
$temp="";
$refgenomelength=0;
while(<INFILE>){
    $temp=$temp.$_;
    chomp($temp);
    $refgenomelength=$refgenomelength + length($temp);
}
my @refgenome=split(//,$temp);


#print "Min depth=".$mindepth."\nMin qual=".$minqual."\nFrac deletion=".$minfracdel."\nRef file=".$reffile."\n";

my $refpos=0;
my $id;
my $pos;
my $dummy;
my $ref;
my $geno;
my $qual;
my $info;
my $depth;
my $avmapq;
my $del;
my $genomeseq="";

while(<STDIN>){
    $temp=$_;
    chomp($temp);
    if( $temp !~ m/^\#/ ){
#	print $temp."\n";
	($id,$pos,$dummy,$ref,$geno,$qual,$dummy,$info)=split(/\t/,$temp);
	$refpos++;

	# First, fill in with reference if needed
	while($refpos < $pos){
	    $genomeseq=$genomeseq.$refgenome[$refpos-1];
	    $refpos++;
	}

	# Then print the new genotype (either a substitution or an INDEL containing the ref nucleotide where needed
	$genomeseq=$genomeseq.$geno;

	if( $info =~ m/FRACINS=/ ){
	    $refpos=$pos;
	    $insertions=$insertions+length($geno)-1;
	}
	elsif( $info =~ m/DEL/ ){
	    $refpos=$pos+length($ref)-1;
	    $deletions=$deletions+length($ref)-1;
	}
	else{
	    $refpos=$pos;
	    $substitutions++;
	}
    }
}
while($refpos < $refgenomelength){
    $genomeseq=$genomeseq.$refgenome[$refpos-1];
    $refpos++;
}


print ">My sequenced genome: ".$reffile." Number of insertions: ".$insertions." Number of deletions: ".$deletions." Number of SNPs: ".$substitutions."\n".$genomeseq."\n"

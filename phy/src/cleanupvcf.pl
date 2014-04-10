#!/usr/bin/perl -w use strict;
use Getopt::Long;
use warnings;

my $mindepth=10;
my $minqual=30;
my $minfracdel=0.9;
my $minfracins=0.9;

# Read VCF on STDIN. 
GetOptions( "mindepth:i" => \$mindepth,
            "minqual:i" => \$minqual,
            "minfracdel:f" => \$minfracdel,
            "minfracins:f" => \$minfracins)
    or die("Unrecognized arguments.\n");

#print STDERR "Min depth=".$mindepth."\nMin qual=".$minqual."\nMin frac deletion=".$minfracdel."\nMin frac insertion=".$minfracins."\n";

my $id;
my $pos;
my $dummy;
my $ref;
my $geno;
my $qual;
my $info;
my $depth;
my $avmapq;

my $numins=0;
my $insdepthline="";
my $insfracline="";
my $insgenoline="";
my $fracins;

my $numdel=0;
my $deldepthline="";
my $delfracline="";
my $delgenoline="";
my $fracdel;
my $deletedreads;

my $lastpos;
my $lastref;
my $lastgeno;
my $lastqual;
my $lastinfo;
my $lastdepth=-1;

my $substitution=0;

while(<STDIN>){
    $temp=$_;
    chomp($temp);
    if( $temp =~ m/^\#/ ){
	print $temp."\n";
    }
    else{	
	($id,$pos,$dummy,$ref,$geno,$qual,$dummy,$info)=split(/\t/,$temp);
	($depth = $info )=~ s/DP=([0-9]+);.*/$1/;
	($avmapq = $info )=~ s/.*AVMQ=([0-9]+)/$1/;
	
	if( $info =~ m/DEL/ ){
	    ($fracdel = $info ) =~ s/.*FRACDEL=([0-9]\.[0-9]+);.*/$1/;
	    ($deletedreads = $info ) =~ s/.*\DEL=([0-9]+);.*/$1/;
	    $depth = $depth + $deletedreads;
	    if($fracdel >= $minfracdel && $depth >= $mindepth ){
		if($numdel==0){
		    $deldepthline="DEL=".$depth;
		    $delfracline="FRACDEL=".$fracdel;
		}
		else{
		    $deldepthline=$deldepthline.",".$depth;
		    $delfracline=$delfracline.",".$fracdel;
		}
		$delgenoline=$delgenoline.$ref;
		$numdel++;
	    }
	    else{
		$info =~ s/DEL.*PP/PP/;
	    }
	}

	if( $info =~ m/INS=/){
	    $fracins=$depth/$lastdepth;
	    if( $fracins >= $minfracins && $depth >= $mindepth ){
		if($numins==0){
		    $insdepthline="INS=".$depth;
		    $insfracline="FRACINS=".$fracins; 
		}
		else{
		    $insdepthline=$insdepthline.",".$depth;
		    $insfracline=$insfracline.",".$fracins; 
		}
		$insgenoline=$insgenoline.$geno;
		$numins++;
	    }
	}

	if( $info !~ m/INS=/ && $info !~ m/DEL/ ){
	    # Output a line if:
	    # It is a high quality substitution
	    # Or we have a high quality insertion
	    # Or we have a high quality deletion

	    # Check if we have an insertion, deletion or substitution that needs to be written
	    if($numins > 0){
		print "$id\t$lastpos\t.\t$lastref\t$lastgeno$insgenoline\t$lastqual\t.\t$lastinfo;$insdepthline;$insfracline\n";
		$numins=0;
		$substitution=0;
		$insdepthline="";
		$insfracline="";
		$insgenoline="";
	    }
	    
	    if($numdel > 0){
		print "$id\t$lastpos\t.\t$lastref$delgenoline\t$lastgeno\t$lastqual\t.\t$lastinfo;$deldepthline;$delfracline\n";
		$numdel=0;
		$substitution=0;
		$deldepthline="";
		$delfracline="";
		$delgenoline="";
	    }

	    if($substitution == 1){
		print "$id\t$lastpos\t.\t$lastref\t$lastgeno\t$lastqual\t.\t$lastinfo\n";
		$substitution=0;
	    }
	    
	    if( $depth>=$mindepth && $qual >= $minqual && $avmapq >= $minqual && $geno ne "."){
		$substitution=1;
	    }
	    else{
		$geno=$ref;
	    }
	    $lastpos=$pos;
	    $lastref=$ref;
	    $lastdepth=$depth;
	    $lastgeno=$geno;
	    $lastqual=$qual;
	    $lastinfo=$info;
	}
    }
}
    
if($numins > 0){
    print "$id\t$lastpos\t.\t$lastref\t$lastgeno$insgenoline\t$lastqual\t.\t$lastinfo;$insdepthline;$insfracline\n";
}
if($numdel > 0){
    print "$id\t$lastpos\t.\t$lastref$delgenoline\t$lastgeno\t$lastqual\t.\t$lastinfo;$deldepthline;$delfracline\n";
}
if($substitution == 1){
    print "$id\t$lastpos\t.\t$lastref\t$lastgeno\t$lastqual\t.\t$lastinfo\n";
}

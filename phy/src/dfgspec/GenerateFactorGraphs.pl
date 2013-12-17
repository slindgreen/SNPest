#!/usr/bin/perl -w
use strict;

# Script to autogenerate factor graphs.
# The files are called depthN_variables.txt and depthN_factorGraph.txt


my $filename;

my $min=$ARGV[0];
my $max=$ARGV[1];
print $min."\n".$max."\n";

for(my $depth=$min;$depth<=$max;$depth++){
    # First, generate variables.txt
    $filename="depth".$depth."_variables.txt";
    open FILE, ">", $filename or die $!;
    print FILE "STATE_MAP_NAME:	nucleotideMap\nVAR_NAMES:\tC ";
    # Loop to generate state A for each read
    for(my $i=1;$i<=$depth;$i++){
	print FILE "A$i ";
    }
    print FILE "\nSTATE_MAP_NAME:\tgenotypeMap\nVAR_NAMES:\tG\n\nSTATE_MAP_NAME:\tinputMap\nVAR_NAMES:\t";
    # Loop to generate state O for each read
    for(my $i=1;$i<=$depth;$i++){
	print FILE "O$i ";
    }
    print FILE "\n\n";
    close(FILE);

    # Now, generate the file factorGraph.txt
    $filename="depth".$depth."_factorGraph.txt";
    open FILE, ">", $filename or die $!;
    print FILE "NAME:\t\tC.prior\nNB1:\t\tC\nPOT:\t\tprior\n\nNAME:\t\tC.G\nNB1:\t\tC\nNB2:\t\tG\nPOT:\t\tgenotype\n\n";
    # Loop to generate state factors for each read
    for(my $i=1;$i<=$depth;$i++){
	print FILE "# Read depth $i\nNAME:\t\tG.A$i\nNB1:\t\tG\nNB2:\t\tA$i\nPOT:\t\toriginal\n\nNAME:\t\tA$i.O$i\nNB1:\t\tA$i\nNB2:\t\tO$i\nPOT:\t\tobservation\n\n";
    }
    close(FILE);
}

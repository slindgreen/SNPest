To install SNPest, you first need to install the phy library as described here:
http://github.com/jakob-skou-pedersen/phy/

For SNPest to work, you have to place the 'dfgEval_SNPest.cpp' file in the '/phy/src' folder along with the modified version of the 'Makefile.am' provided with SNPest.

The folder 'dfgspec' contains all the model specifications and should be placed in the '/phy/src' folder as well.

SNPest.pl reads input on STDIN and outputs the genotype data on STDOUT.

You should use the provided script 'cleanupvcf.pl' to generate a high quality set of SNPs and indels from the output. The default is to use a minimum read depth of 10X, a minimum phred scaled quality of 30, and - for insertions and deletions - a minimum of 90% of reads agreeing with the indel.

Run SNPest.pl -h to see the possible parameters.

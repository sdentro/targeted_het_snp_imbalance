# Mouse heterozygous SNP imbalance pipeline

For targeted sequencing data

## Running

Run in this order:
 * allelecount
 * segmentation
 * combine
 * plotting

## Dependencies
 * alleleCounter
 * R
 * optparse
 * Battenberg
 * ggplot2
 * gridExtra
 * gtools
 * GenomicRanges 

## Generating reference data

First download the SNP reference data from the mouse genomes project here: (ftp://ftp-mouse.sanger.ac.uk/current_snps)[ftp://ftp-mouse.sanger.ac.uk/current_snps]

Then `cd reference` and check the correct filename is referenced in split.sh.

Finally run `./split.sh`, which takes a while to generate the alleles

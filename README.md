# Mouse heterozygous SNP imbalance pipeline

For targeted sequencing data

## Running

First make sure to place all BAM files in the `raw_data` directory

Then add samples to the `sample_mapping.txt` file. This file contains two columns, separated by a tab: samplename and controlsamplename

Please note that samplenames must be unique. The name provided in `sample_mapping.txt` is used to pattern match BAM files and output of the pipeline

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

First download the SNP reference data from the mouse genomes project here: [ftp://ftp-mouse.sanger.ac.uk/current_snps](ftp://ftp-mouse.sanger.ac.uk/current\_snps)

Then `cd reference` and check the correct filename is referenced in split.sh.

Finally run `./split.sh`, which takes a while to generate the alleles

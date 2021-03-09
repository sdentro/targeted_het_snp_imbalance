# Mouse heterozygous SNP imbalance pipeline

For targeted sequencing data

## Running

First make sure to place all BAM files in the `raw_data` directory

Then add samples to the `sample_mapping.txt` file. This file contains two columns, separated by a tab: samplename and controlsamplename

Please note that samplenames must be unique. The name provided in `sample_mapping.txt` is used to pattern match BAM files and output of the pipeline

### hg19
Run in this order:
 * allelecount
 * phasing
 * combine
 * plotting

### mm10
Run in this order:
 * allelecount
 * segmentation
 * combine
 * plotting

## Dependencies - overall
 * alleleCounter
 * R
 * optparse
 * Battenberg
 * ggplot2
 * gridExtra
 * gtools
 * GenomicRanges 

### hg19 - additional dependencies
 * Impute2

## Reference data

### hg19
Reference data to count alleles and perform haplotype reconstruction is available on the Battenberg github page
 * Battenberg GRCh37 reference files from [here](https://github.com/Wedge-lab/battenberg)

### mm10
First download the SNP reference data from the mouse genomes project here: [ftp://ftp-mouse.sanger.ac.uk/current_snps](ftp://ftp-mouse.sanger.ac.uk/current_snps)

Then `cd reference` and check the correct filename is referenced in split.sh.

Finally run `./split.sh`, which takes a while to generate the alleles

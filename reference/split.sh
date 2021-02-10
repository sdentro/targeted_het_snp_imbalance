#!/bin/bash
zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz | grep -v ^\# | cut -f 1,2,4,5 | grep -v "," | awk '{print > ("mm10_SNP_alleles_" $1 ".txt")}'
mv mm10_SNP_alleles_X.txt mm10_SNP_alleles_20.txt
mv mm10_SNP_alleles_Y.txt mm10_SNP_alleles_21.txt
rm mm10_SNP_alleles_MT.txt

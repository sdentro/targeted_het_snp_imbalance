#!/bin/bash
alleleprefix="../../reference/mm10_SNP_alleles_"
mem=15000
for item in `find ../../raw_data | grep "\.bam" | grep -v bai`; do
	sam=`basename ${item} | cut -f 1 -d .`
	for i in $(seq 1 21); do
		bsub -R"select[mem>${mem}] rusage[mem=${mem}]" -M${mem} -o logs/${sam}_${i}_%J.out "alleleCounter -d -l ${alleleprefix}${i}.txt -b ${item} -o output/${sam}_alleleFrequencies_chr${i}.txt"
	done
done

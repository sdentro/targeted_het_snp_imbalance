#!/bin/bash
alleleprefix="../../reference/mm10_SNP_alleles_"
mem=15000
find output/ | grep txt > all_done.lst
for item in `find ../../raw_data | grep "\.bam" | grep -v bai`; do
	sam=`basename ${item} | cut -f 1 -d .`
	for i in $(seq 1 21); do
	if grep ${sam}_alleleFrequencies_chr${i}.txt all_done.lst > /dev/null; then 
		continue; 
	else
		echo -e ${sam}"\t"${i}
	fi
	done
done > todo.txt

#!/bin/bash
alleleprefix=`cat ../config.txt | grep human_g1000lociprefix | cut -f 2`
mem=15000
nlines=`cat todo.txt | wc -l`
#for item in `find ../../raw_data | grep "\.bam" | grep -v bai | tail -n 684`; do
for j in $(seq 1 ${nlines}); do
	sam=`cut -f 1 todo.txt | head -n ${j} | tail -n 1`
	i=`cut -f 2 todo.txt | head -n ${j} | tail -n 1`
	if [[ "${i}" == 23 ]]; then
		chrom="X"
	else
		chrom=${i}
	fi
	item=`find ../../raw_data | grep "\.bam" | grep -v bai | grep "${sam}\."`
	bsub -R"select[mem>${mem}] rusage[mem=${mem}]" -M${mem} -o logs/${sam}_${i}_%J.out "alleleCounter -l ${alleleprefix}${i}.txt -b ${item} -o output/${sam}_alleleFrequencies_chr${chrom}.txt"
done

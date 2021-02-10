#!/bin/bash
mem=10000
nthreads=1
tumbam="NA"
normbam="NA"
nsamples=`cat ../../sample_mapping.txt | wc -l`
#nsamples=1
for i in $(seq 1 ${nsamples}); do
	tumour=`head -n ${i} ../../sample_mapping.txt | cut -f 1 | tail -n 1`;
	normal=`head -n ${i} ../../sample_mapping.txt | cut -f 2 | tail -n 1`;
	if ls output/ | grep "${tumour}.BAF" > /dev/null; then 
		# Skip when already done
		continue
	else
		bsub -n ${nthreads} -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" -M${mem} -o logs/${tumour}_%J.out Rscript pipeline.R ${tumour} ${normal} "mouse"
	fi
done

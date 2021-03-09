#!/bin/bash
mem=70000
nthreads=3
tumbam="NA"
normbam="NA"
nsamples=`cat ../../sample_mapping.lst | wc -l`
for i in $(seq 1 ${nsamples}); do
	tumour=`head -n ${i} ../../sample_mapping.lst | cut -f 1 | tail -n 1`;
	normal=`head -n ${i} ../../sample_mapping.lst | cut -f 2 | tail -n 1`;
	bsub -n ${nthreads} -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" -M${mem} -o logs/${item}_%J.out Rscript pipeline.R -t ${tumour} -n ${normal} -o output --tb ${tumbam} --nb ${normbam} --sex male --cpu ${nthreads}
done

#!/bin/bash
mem=5000
mapping="../../sample_mapping.txt"
num_samples=`cat ${mapping} | wc -l`
for i in $(seq 1 ${num_samples}); do
  item=`head -n ${i} ${mapping} | tail -n 1 | cut -f 1`
  controlsample=`head -n ${i} ${mapping} | tail -n 1 | cut -f 2`
  bsub -R"select[mem>${mem}] rusage[mem=${mem}]" -M${mem} -o logs/${item}_%J.out Rscript pipeline.R -s ${item} -c ${controlsample} -o output/ --species mouse --min_baf_divergence 0.1 --reads_het_snp_allele 10 --min_cov 25 --seg ../segmentation/output/${item}.BAFsegmented.txt
done

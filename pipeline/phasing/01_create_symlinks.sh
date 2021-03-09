#!/bin/bash
cd output/
for samplename in `cat ../../sample_mapping.lst | cut -f 1`; do
	mkdir -p ${samplename};
	cd ${samplename};
	normalname=`grep ${samplename} ../../../sample_mapping.lst | cut -f 2`;
	find ../../allelecount/output | grep ${samplename} | grep -v logs | xargs -i ln -s {}
	find ../../allelecount/output | grep ${normalname} | grep -v logs | xargs -i ln -s {}
	cd ../
done

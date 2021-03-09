## Phasing

Use this step to reconstruct haplotypes for human samples. If you're analysing mouse data, please run the segmentation step instead, as there are no haplotype blocks for mm10. This step requires Impute2 to be available in $PATH and that the `config.txt` file points to Battenberg reference files for hg19/GRCh37. Note: The impute_info.txt file also contains file paths that need to be correct!

### Running
 * First run `01_create_symlinks.sh` to create symlinks for the allelecounts files generated in the previous step in the output directory
 * Then run `02_submit.sh` to submit all jobs into the queue

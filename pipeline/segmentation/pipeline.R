#' pipeline for when phasing is not available

args = commandArgs(T)
samplename = args[1]
controlname = args[2]
species = args[3]

# samplename = "MD6356d"
# controlname = "MD6356b"
# species = "mouse"

if (species=="human") {
  # 1000G reference files
  loci_dir_path = "../../reference/snptest_1000genomesloci_phase3/"
  loci_file_prefix = "1000genomesloci_phase3_chr"
  chroms = 1:23
} else {
  loci_dir_path = "../../reference/"
  loci_file_prefix = "mm10_SNP_alleles_"
  chroms = 1:21
}

segmentation_gamma = 6
segmentation_kmin = 3
min_het_vaf_control = 0.1 #0.4
max_het_vaf_control = 0.9 #0.6

counts_dir = "../allelecount/output"
control_counts_file_prefix = paste0(controlname, "_alleleFrequencies_chr")
sample_counts_file_prefix = paste0(samplename, "_alleleFrequencies_chr")

min_reads_allele_het_snp = 30
min_coverage = 75
nucleotides = 3:6
names(nucleotides) = c("A", "C", "G", "T")

selections = list()
allele_data = list()
control_data = list()
for (chrom in as.character(chroms)) {
  dat = readr::read_tsv(file.path(counts_dir, paste0(control_counts_file_prefix, chrom, ".txt")))
  sel = dat$Good_depth >= min_coverage
  selections[[chrom]] = sel
  control_data[[chrom]] = dat[sel,]
  rm(dat)
  dat = readr::read_tsv(file.path(loci_dir_path, paste0(loci_file_prefix, chrom, ".txt")), col_names=F)
  allele_data[[chrom]] = dat[sel,]
}

output = list()
for (chrom in as.character(names(allele_data))) {
  if (sum(selections[[chrom]])==0) {
    output[[chrom]] = NULL
  } else {
    dat = readr::read_tsv(file.path(counts_dir, paste0(sample_counts_file_prefix, chrom, ".txt")))
    dat = dat[selections[[chrom]],]
    
    n_s_a0 = as.numeric(as.data.frame(dat)[cbind(1:nrow(dat), nucleotides[allele_data[[chrom]]$X3])])
    n_s_a1 = as.numeric(as.data.frame(dat)[cbind(1:nrow(dat), nucleotides[allele_data[[chrom]]$X4])])
    n_c_a0 = as.numeric(as.data.frame(control_data[[chrom]])[cbind(1:nrow(dat), nucleotides[allele_data[[chrom]]$X3])])
    n_c_a1 = as.numeric(as.data.frame(control_data[[chrom]])[cbind(1:nrow(dat), nucleotides[allele_data[[chrom]]$X4])])
    output[[chrom]] = data.frame(chrom=dat[[1]], pos=dat[[2]], a0=allele_data[[chrom]][[3]], a1=allele_data[[chrom]][[4]], n_a0_sample=n_s_a0, n_a1_sample=n_s_a1, n_a0_control=n_c_a0, n_a1_control=n_c_a1, stringsAsFactors=F)
  }
}
rm(selections, allele_data, control_data, dat)
output = do.call(rbind, output)
print(head(output))
print(dim(output))
# remove any double alleles
output = output[output$a1 %in% c("A","C","G","T"),]
output$is_het = F
output$vaf_control = output$n_a0_control / (output$n_a0_control+output$n_a1_control)
output$vaf_sample = output$n_a0_sample / (output$n_a0_sample+output$n_a1_sample)
print(sum(output$is_het))
if (!is.na(min_het_vaf_control) & !is.na(max_het_vaf_control)) {
	# new approach that only takes bonafide hets
	output$is_het = output$n_a0_control > min_reads_allele_het_snp & output$n_a1_control > min_reads_allele_het_snp & (output$vaf_control < 0.6 & output$vaf_control > 0.4)
} else {
	# old approach that selects any het SNPs
	output$is_het = output$n_a0_control > min_reads_allele_het_snp & output$n_a1_control > min_reads_allele_het_snp
}


if (sum(output$is_het) == 0) {
	output$baf_segment = NA
} else {

	# dat = read.table("output/MD6356d_snp_inventory.txt.gz", header=T, stringsAsFactors=F)
	het_snps = output$is_het
	alt_vaf = ifelse(output$vaf_sample[het_snps] < 0.5, 1-output$vaf_sample[het_snps], output$vaf_sample[het_snps])
	
	sdev <- Battenberg:::getMad(alt_vaf, k=25)
	if (sdev<0.09 | is.na(sdev)) { sdev = 0.09 }
	res = Battenberg:::selectFastPcf(alt_vaf, segmentation_kmin, segmentation_gamma * sdev, T)
	
	output$baf_segment = NA
	output$baf_segment[het_snps] = res$yhat
}

output = output[het_snps,]
output = output[, c("chrom", "pos", "vaf_sample", "baf_segment", "baf_segment")]
colnames(output) = c("Chromosome", "Position", "BAF", "BAFphased", "BAFseg")

write.table(output, file=file.path("output", paste0(samplename, ".BAFsegmented.txt")), row.names=F, sep="\t", quote=F)




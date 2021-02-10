# pipeline to test for significant deviation from 0.5 of heterozygous SNPs - human

library(optparse)
library(ggplot2)

option_list = list(
  make_option(c("-s", "--samplename"), type="character", default=NULL, help="Samplename of the tumour", metavar="character"),
  make_option(c("-c", "--controlname"), type="character", default=NULL, help="Samplename of the control", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, help="Output directory", metavar="character"),
  make_option(c("--seg"), type="character", default=NULL, help="Segmented file as output from phasing", metavar="character"),
  make_option(c("--species"), type="character", default="human", help="Species", metavar="character"),
  make_option(c("--min_baf_divergence"), type="numeric", default=0.01, help="Minimum BAF divergence required before testing", metavar="character"),
  make_option(c("--reads_het_snp_allele"), type="numeric", default=30, help="Minimum number of reads per allele required to call a SNP heterozygous", metavar="character"),
  make_option(c("--min_cov"), type="numeric", default=75, help="Minimum coverage required to keep a SNP", metavar="character"),
  make_option(c("--snps_are_phased"), type="logical", action="store_true", default=FALSE, help="Whether SNPs have been phased", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

samplename = opt$samplename
controlname = opt$controlname
outdir = opt$outdir
if ("seg" %in% names(opt)) {
  segmented_file = opt$seg
} else {
  segmented_file = NULL
}
min_baf_divergence = opt$min_baf_divergence
min_reads_allele_het_snp = opt$reads_het_snp_allele
min_coverage = opt$min_cov
species = opt$species
snps_are_phased = opt$snps_are_phased

#args = commandArgs(T)
#samplename = args[1]
#controlname = args[2]
#outdir = args[3]
#if (length(args) > 3) {
#  segmented_file = args[4]
#} else {
#  segmented_file = NULL
#}

# parameter of how far the BAF has to change before considering it for testing
# So far used: 0.01 for exome/genome and 0.10 for targeted
# (for segmentation only so far)
#if (length(args) > 4) {
#  min_baf_divergence = as.numeric(args[5])
#} else {
#  min_baf_divergence = 0.01
#}

# outdir = "output_exome"
# controlname = "PD31008b"
# samplename = "PD31008g"
# segmented_file = "phasing/exome/output/PD31008g/PD31008g.BAFsegmented.txt"

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

min_het_vaf_control = 0.4
max_het_vaf_control = 0.6

# path to where allelecounts have been stored
counts_dir = "../allelecount/output"
control_counts_file_prefix = paste0(controlname, "_alleleFrequencies_chr")
sample_counts_file_prefix = paste0(samplename, "_alleleFrequencies_chr")

#min_reads_allele_het_snp = 30
#min_coverage = 75
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
# remove any double alleles
output = output[output$a1 %in% c("A","C","G","T"),]
output$is_het = F
output$vaf_control = output$n_a0_control / (output$n_a0_control+output$n_a1_control)
output$vaf_sample = output$n_a0_sample / (output$n_a0_sample+output$n_a1_sample)

if (!is.na(min_het_vaf_control) & !is.na(max_het_vaf_control)) {
        # new approach that only takes bonafide hets
        output$is_het = output$n_a0_control > min_reads_allele_het_snp & output$n_a1_control > min_reads_allele_het_snp & (output$vaf_control < max_het_vaf_control & output$vaf_control > min_het_vaf_control)
} else {
        # old approach that selects any het SNPs
        output$is_het = output$n_a0_control > min_reads_allele_het_snp & output$n_a1_control > min_reads_allele_het_snp
}

#save.image("checkpoint0.RData")
print("Checking segmentation")
# consider switching alleles when there is phasing information
if (!is.null(segmented_file)) {
  print("Getting to switch alleles")
  # to_switch = c()
  to_switch_chrom = list()
  output$vaf_segment_sample = NA
output$segmented_p_value_balance = NA
output$segmented_p_value_balance_adj = NA

  skipped = c() # keep track of which are skipped, to check
  # add in phasing info, if available
  segmented = read.table(segmented_file, header=T, stringsAsFactors=F)
  segmented_split = split(segmented, segmented$Chromosome)
  #segmented_chrpos = paste0(segmented$Chromosome, "_", segmented$Position)
  #output_chrpos = paste0(output$chrom, "_", output$pos)
  output_split = split(output, output$chrom)
  # n_iter = nrow(segmented)

  print("Starting main for loop")
  for (chrom in names(segmented_split)) {
    print(chrom)
    segmented_chrom = segmented_split[[chrom]]
    output_chrom = output_split[[chrom]]
    n_iter = nrow(segmented_chrom)
    for (i in 1:nrow(segmented_chrom)) {
      if (i %% 100 == 0) { print(paste("Iter", i, "out of", n_iter)) }
      seg_entry = segmented_chrom[i,]
      #output_index = which(output$chrom==seg_entry$Chromosome & output$pos==seg_entry$Position)
      #output_index = which(output_chrpos == paste0(seg_entry$Chromosome, "_", seg_entry$Position) & output$is_het)
      output_index = which(output_chrom$pos==seg_entry$Position)
      
    
      if (length(output_index) == 0) { skipped=c(skipped, i); next() }
      else if (length(output_index) > 1) { 
        
        if (any(output_chrom$is_het[output_index])) {
          #output_index = which(output$chrom==seg_entry$Chromosome & output$pos==seg_entry$Position & output$is_het)
      	#output_index = which(output_chrpos == paste0(seg_entry$Chromosome, "_", seg_entry$Position) & output$is_het)
          output_index = which(output_chrom$pos==seg_entry$Position & output_chrom$is_het)
        } else {
          # none are hets anyway, so skip - this should not occur
          next()
        }
      }
      # entry = output[output$chrom=="1" & output$pos==1001177,]
      # seg_entry = segmented[segmented$Chromosome=="1" & segmented$Position==1001177,]
      entry = output_chrom[output_index,]
      
      # check if af different and more than machine precision
      if (seg_entry$BAFphased != entry$vaf_sample & abs(seg_entry$BAFphased - entry$vaf_sample) > (3*.Machine$double.eps)) {
        #to_switch = c(to_switch, output_index)
        to_switch_chrom[[chrom]] = c(to_switch_chrom[[chrom]], output_index)
      } 
      # annotate segment BAF
      output_split[[chrom]]$vaf_segment_sample[output_index] = seg_entry$BAFseg
    }
  }
  
  print("Switching alleles")
  #if (length(to_switch) > 0) {
    for (chrom in names(output_split)) {
    print(chrom)
    if (length(to_switch_chrom[[chrom]])==0) { next(); }
      entry_copy = output_split[[chrom]][to_switch_chrom[[chrom]],,drop=F]
      output_split[[chrom]]$a0[to_switch_chrom[[chrom]]] = entry_copy$a1
      output_split[[chrom]]$a1[to_switch_chrom[[chrom]]] = entry_copy$a0
      output_split[[chrom]]$n_a0_sample[to_switch_chrom[[chrom]]] = entry_copy$n_a1_sample
      output_split[[chrom]]$n_a1_sample[to_switch_chrom[[chrom]]] = entry_copy$n_a0_sample
      output_split[[chrom]]$n_a0_control[to_switch_chrom[[chrom]]] = entry_copy$n_a1_control
      output_split[[chrom]]$n_a1_control[to_switch_chrom[[chrom]]] = entry_copy$n_a0_control
      output_split[[chrom]]$vaf_control[to_switch_chrom[[chrom]]] = output_split[[chrom]]$n_a0_control[to_switch_chrom[[chrom]]] / (output_split[[chrom]]$n_a0_control[to_switch_chrom[[chrom]]]+output_split[[chrom]]$n_a1_control[to_switch_chrom[[chrom]]])
      output_split[[chrom]]$vaf_sample[to_switch_chrom[[chrom]]] = output_split[[chrom]]$n_a0_sample[to_switch_chrom[[chrom]]] / (output_split[[chrom]]$n_a0_sample[to_switch_chrom[[chrom]]]+output_split[[chrom]]$n_a1_sample[to_switch_chrom[[chrom]]])
    }
    output = do.call(rbind, output_split)
  #}
}

if (!snps_are_phased) {
  output$vaf_control = ifelse(output$vaf_control < 0.5, 1-output$vaf_control, output$vaf_control)
  output$vaf_sample = ifelse(output$vaf_sample < 0.5, 1-output$vaf_sample, output$vaf_sample)
}


print("Testing SNPs")
# test SNPs individually
output$vaf_p_value = NA
output$vaf_p_value[output$is_het] = unlist(lapply(which(output$is_het), function(i) {
  x = output[i,]
  if ((x["n_a0_sample"]+x["n_a1_sample"]) == 0) {
	NA
  } else {
  prop.test(x=as.numeric(x["n_a0_sample"]), 
            n=as.numeric(x["n_a0_sample"]+x["n_a1_sample"]),
            p=as.numeric(x["vaf_control"]),
            alternative="two.sided")$p.value 
  }}))

output$vaf_p_value_adj = p.adjust(output$vaf_p_value, method="bonferroni")

# test VAF is different from 0.5
output$exphetvaf_p_value = NA
output$exphetvaf_p_value[output$is_het] = unlist(lapply(which(output$is_het), function(i) {
  x = output[i,]
  if ((x["n_a0_sample"]+x["n_a1_sample"]) == 0) {
        NA
  } else {
  prop.test(x=as.numeric(x["n_a0_sample"]),
            n=as.numeric(x["n_a0_sample"]+x["n_a1_sample"]),
            p=0.5,
            alternative="two.sided")$p.value
  }}))
output$exphetvaf_p_value_adj = p.adjust(output$exphetvaf_p_value, method="bonferroni")

# if both tests are significant then mark as signif
output$is_signif = output$vaf_p_value_adj < 0.05 & output$exphetvaf_p_value_adj < 0.05
output$is_signif[is.na(output$is_signif)] = F

print("Adding control segmentation values and performing segment test")
output$vaf_segment_control = NA
output$segment_p_value = NA
output$segment_p_value_adj = NA
output$segment_p_value_balanced = NA
output$segment_p_value_balanced_adj = NA
output$segment_is_signif = NA
output$segment_is_signif_adj = NA

#save.image("checkpoint1.RData")

if (!is.null(segmented_file)) {
  output$segment_is_signif = FALSE
  # Add af of segment for control and test per segment
  # is_phased_index = which(!is.na(output$vaf_segment_sample))
  # seg_breaks = cumsum(rle(output$vaf_segment_sample[is_phased_index])$lengths)
  
  num_segs = 0
  chrom_order = unique(output$chrom)
  for (j in 1:length(chrom_order)) {
    # for (i in 1:2) {
    output_chrom = output[output$chrom==chrom_order[j],]
    is_phased_index = which(!is.na(output_chrom$vaf_segment_sample) & output_chrom$is_het)
    if (length(is_phased_index) > 1) {
      seg_breaks = cumsum(rle(output_chrom$vaf_segment_sample[is_phased_index])$lengths)
      # keep track of number of segments (number of tests done) for multiple testing correction
      num_segs = num_segs + length(seg_breaks)
      
      for (i in 1:length(seg_breaks)) {
        if (i==1) {
          startpoint = 1
        } else {
          startpoint = seg_breaks[i-1]+1
        }
        endpoint = seg_breaks[i]
        if (snps_are_phased) {
          output_chrom$vaf_segment_control[is_phased_index[startpoint:endpoint]] = median(output_chrom$vaf_control[is_phased_index[startpoint:endpoint]])
        } else {
          selected_snps = output_chrom$vaf_control[is_phased_index[startpoint:endpoint]]
          selected_snps = ifelse(selected_snps < 0.5, 1-selected_snps, selected_snps)
          output_chrom$vaf_segment_control[is_phased_index[startpoint:endpoint]] = median(selected_snps)
        }
        
        # do not select too many SNPs, otherwise the test will always be significant
        if ((endpoint-startpoint) > 500) {
          selected_snps = sample(startpoint:endpoint, 500)
        } else {
          selected_snps = startpoint:endpoint
        }
        
        # set a minimim deviation required for calling a divergence
        do_test_control = abs(output_chrom$vaf_segment_sample[is_phased_index[selected_snps]][1] - output_chrom$vaf_segment_control[is_phased_index[selected_snps]][1]) > min_baf_divergence
        do_test_balanced = abs(output_chrom$vaf_segment_sample[is_phased_index[selected_snps]][1] - 0.5) > min_baf_divergence
        #if (abs(output_chrom$vaf_segment_sample[is_phased_index[selected_snps]][1] - output_chrom$vaf_segment_control[is_phased_index[selected_snps]][1]) > min_baf_divergence & abs(output_chrom$vaf_segment_sample[is_phased_index[selected_snps]][1] - 0.5) > min_baf_divergence) {
        if (do_test_control & do_test_balanced & length(selected_snps) > 1) {
          # do a simple two sided t-test    
          output_chrom$segment_p_value[is_phased_index[startpoint:endpoint]] = t.test(x=output_chrom$vaf_sample[is_phased_index[selected_snps]],
                                                                                      y=output_chrom$vaf_control[is_phased_index[selected_snps]],
                                                                                      alternative="two.sided")$p.value
        	# adjust control to exact 0.5 mean to test against diff from 0.5
        	mean_control = mean(output_chrom$vaf_control[is_phased_index[selected_snps]])
        	if (mean_control < 0.5) {
        		balanced_control = output_chrom$vaf_control[is_phased_index[selected_snps]] + abs(0.5-mean_control)
        	} else {
        		balanced_control = output_chrom$vaf_control[is_phased_index[selected_snps]] - abs(0.5-mean_control)
        	}
		if (length(selected_snps) > 1) {
        	output_chrom$segment_p_value_balanced[is_phased_index[startpoint:endpoint]] = t.test(x=output_chrom$vaf_sample[is_phased_index[selected_snps]],
        											    y=balanced_control,
        											    alternative="two.sided")$p.value				
		}
        } else {
          output_chrom$segment_p_value[is_phased_index[startpoint:endpoint]] = 1
  	      output_chrom$segment_p_value_balanced[is_phased_index[startpoint:endpoint]] = 1
        }
      }
      output[output$chrom==chrom_order[j],] = output_chrom
    }
  }
  print("Adjusting for multiple testing") 
  # adjust - bonferroni
  chrom_order = unique(output$chrom)
  for (j in 1:length(chrom_order)) {
    # for (i in 1:2) {
    output_chrom = output[output$chrom==chrom_order[j],]
    is_phased_index = which(!is.na(output_chrom$vaf_segment_sample) & output_chrom$is_het)
    if (length(is_phased_index) > 0) {
      seg_breaks = cumsum(rle(output_chrom$vaf_segment_sample[is_phased_index])$lengths)
    
      for (i in 1:length(seg_breaks)) {
        if (i==1) {
          startpoint = 1
        } else {
          startpoint = seg_breaks[i-1]+1
        }
        endpoint = seg_breaks[i]
        
        #output_chrom$segment_p_value_adj[is_phased_index[startpoint:endpoint]] = pmin(1, num_segs * output_chrom$segment_p_value[is_phased_index[startpoint]])
        #output_chrom$segment_p_value_balanced_adj[is_phased_index[startpoint:endpoint]] = pmin(1, num_segs * output_chrom$segment_p_value_balanced[is_phased_index[startpoint]])
  
        output_chrom$segment_p_value_adj[is_phased_index[startpoint:endpoint]] = p.adjust(output_chrom$segment_p_value[is_phased_index[startpoint]], n=num_segs)
        output_chrom$segment_p_value_balanced_adj[is_phased_index[startpoint:endpoint]] = p.adjust(output_chrom$segment_p_value_balanced[is_phased_index[startpoint]], n=num_segs)
        output_chrom$segment_is_signif_adj[is_phased_index[startpoint:endpoint]] = (output_chrom$segment_p_value_adj[is_phased_index[startpoint]] < 0.05) & (output_chrom$segment_p_value_balanced_adj[is_phased_index[startpoint]] < 0.05)
        output_chrom$segment_is_signif[is_phased_index[startpoint:endpoint]] = (output_chrom$segment_p_value[is_phased_index[startpoint]] < 0.05) & (output_chrom$segment_p_value_balanced[is_phased_index[startpoint]] < 0.05)
        #output_chrom$segment_is_signif[is_phased_index[startpoint:endpoint]] = (output_chrom$segment_p_value[is_phased_index[startpoint]] < 0.05) & (output_chrom$segment_p_value_balanced[is_phased_index[startpoint]] < 0.05)
      }
      output[output$chrom==chrom_order[j],] = output_chrom
    }
  }
  
  
  # is_phased_index = which(!is.na(output$vaf_segment_sample))
  
  # for (i in 1:length(seg_breaks)) {
  #   if (i==1) {
  #     startpoint = 1
  #   } else {
  #     startpoint = seg_breaks[i-1]+1
  #   }
  #   endpoint = seg_breaks[i]
  #   output$vaf_segment_control[is_phased_index[startpoint:endpoint]] = median(output$vaf_control[is_phased_index[startpoint:endpoint]])
  #   
  #   # do a simple two sided t-test    
  #   output$segment_p_value[is_phased_index[startpoint:endpoint]] = t.test(x=output$vaf_sample[is_phased_index[startpoint:endpoint]],
  #                                                                         y=output$vaf_control[is_phased_index[startpoint:endpoint]],
  #                                                                         alternative="two.sided")$p.value
  #   # adjust - bonferroni
  #   output$segment_p_value_adj[is_phased_index[startpoint:endpoint]] = pmin(1, length(seg_breaks) * output$segment_p_value[is_phased_index[startpoint]])
  #   output$segment_is_signif[is_phased_index[startpoint:endpoint]] = output$segment_p_value_adj[is_phased_index[startpoint]] < 0.05
  # }
}

outfile = file.path(outdir, paste0(samplename, "_snp_inventory.txt"))
write.table(output, file=file.path(outdir, paste0(samplename, "_snp_inventory.txt")), row.names=F, sep="\t", quote=F)
system(paste("gzip", file.path(outdir, paste0(samplename, "_snp_inventory.txt"))))

output$depth_sample = output$n_a0_sample+output$n_a1_sample
output$depth_control = output$n_a0_control+output$n_a1_control
p = ggplot(output) + geom_density(mapping=aes(x=depth_control, y=..density..), colour="blue") + geom_density(mapping=aes(x=depth_sample, y=..density..), colour="red")
png(file.path(outdir, paste0(samplename, "_depth.png")), height=500, width=800)
print(p)
dev.off()



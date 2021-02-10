#
# How to run:
# for item in `find ../../human/snp_inventory/ | grep "inventory.txt.gz"`; do Rscript ~/Documents/repo/cnloh_analysis/human_het_snp_analysis.R ${item}; done
#
#
args = commandArgs(T)
# samplename = args[1]
infile = args[1]
outdir = args[2]
species = args[3]

dir.create(file.path(outdir, "plot_all_data", "1"), recursive=T, showWarnings=F)
dir.create(file.path(outdir, "plot_all_data", "2"), recursive=T, showWarnings=F)
dir.create(file.path(outdir, "plot_all_data", "3"), recursive=T, showWarnings=F)
dir.create(file.path(outdir, "plot_detail_notch1"), recursive=T, showWarnings=F)
dir.create(file.path(outdir, "summary_data"), recursive=T, showWarnings=F)
if (species=="human") {
  dir.create(file.path(outdir, "plot_chr9"), recursive=T, showWarnings=F)
} else {
  dir.create(file.path(outdir, "plot_chr2"), recursive=T, showWarnings=F)
}

# infile = "/Users/sdentro/Documents/projects/notch1_cnloh/human/analysis_with_phasing/input_data/testing/PD31008g_snp_inventory.txt.gz"
#outdir = "/Users/sdentro/Documents/projects/notch1_cnloh/human/analysis_with_phasing/testing/"

# dat = read.table("../combine/output/MD6356d_snp_inventory.txt.gz", header=T, stringsAsFactors=F)
# samplename = "MD6356d"

library(ggplot2)
library(gridExtra)
library(gtools)

print("Read in data...")
# dat = read.table(file.path(outdir, "snp_inventory", paste0(samplename, "_snp_inventory.txt.gz")), header=T, stringsAsFactors=F)
dat = read.table(infile, header=T, stringsAsFactors=F)
samplename = unlist(strsplit(basename(infile), "_"))[1]

print(samplename)

print("Checking whether all columns are there...")
# Deal with older run that did not yet have these data added
if (!"is_signif" %in% colnames(dat)) {
  min_reads_allele_het_snp = 20
  dat$is_het = dat$n_a0_control > min_reads_allele_het_snp & dat$n_a1_control > min_reads_allele_het_snp
  dat$vaf_control = dat$n_a0_control / (dat$n_a0_control+dat$n_a1_control)
  dat$vaf_sample = dat$n_a0_sample / (dat$n_a0_sample+dat$n_a1_sample)
  dat$vaf_p_value = NA
  
  dat$vaf_p_value[dat$is_het] = unlist(lapply(which(dat$is_het), function(i) {
              x = dat[i,]
              prop.test(x=as.numeric(x["n_a0_sample"]),
                        n=as.numeric(x["n_a0_sample"]+x["n_a1_sample"]),
                        p=as.numeric(x["vaf_control"]),
                        alternative="two.sided")$p.value }))
  dat$vaf_p_value_adj = p.adjust(dat$vaf_p_value, method="bonferroni")
  dat$is_signif = dat$vaf_p_value_adj < 0.05
  dat$is_signif[is.na(dat$is_signif)] = F
  # write.table(dat, file=file.path(outdir, paste0(samplename, "_snp_inventory_anno.txt")), quote=F, sep="\t", row.names=F)
}

print("Make genome wide plot...")
dat$is_signif = factor(dat$is_signif, levels=c(TRUE, FALSE))
dat$chrom = factor(dat$chrom, levels=gtools::mixedsort(unique(dat$chrom)))
selected = dat$vaf_control < 0.999 & dat$vaf_control > 0.001 & !is.na(dat$vaf_control)
dat_hom = dat[selected & !dat$is_het,]
dat_het = dat[selected & dat$is_het,]
mycolours = c("red", "black")
p = ggplot() +
  scale_colour_manual(values=mycolours, drop=F) + 
  # geom_segment(data=anno, mapping=aes(x=start, xend=end, y=0, yend=1), colour="darkgrey") + # this was used to draw targeted regions
  geom_hline(yintercept=0.5, linetype=2) + 
  geom_point(data=dat_hom, mapping=aes(x=pos, y=vaf_sample), colour="grey", size=0.5) + 
  geom_point(data=dat_het, mapping=aes(x=pos, y=vaf_sample, colour=is_signif), size=0.5) + 
  # scale_x_continuous(breaks=seq(min(anno$start), max(anno$end), 10000000), labels=paste0(round(seq(min(anno$start), max(anno$end), 10000000)/1000000), "Mb")) +
  theme_bw() + xlab("Position") + ylab("BAF") + ylim(0,1) +
  scale_x_continuous(expand = c(0, 0)) +
  #annotate(geom="label", x=90000000, y=1, label=samplename) + 
  # facet_wrap(~samplename, ncol=2, strip.position="right") +
  theme(axis.text.x = element_text(colour="black",size=14,face="plain"),
        axis.title.x = element_text(colour="black",size=16,face="plain"),
        axis.text.y = element_text(colour="black",size=14,face="plain"),
        axis.title.y = element_text(colour="black",size=16,face="plain"),
        strip.text.x = element_text(colour="black",size=16,face="plain"),
        strip.text.y = element_text(colour="black",size=15,face="plain"),
        legend.position = "none") + 
  facet_wrap(~chrom, ncol=2, scales="free_x", strip.position="right")

# collect and summarise segment data
if (!all(is.na(dat$segment_is_signif))) {
  segment_data = data.frame()
  # is_phased_index = which(!is.na(dat$vaf_segment_sample))
  # seg_breaks = cumsum(rle(dat$vaf_segment_sample[is_phased_index])$lengths)
  
  print("Fetching segmentation data...")
  chrom_order = unique(dat$chrom)
  chrom_char = as.character(dat$chrom)
  for (j in 1:length(chrom_order)) {
    dat_chrom = dat[chrom_char==chrom_order[j],]
    is_phased_index = which(!is.na(dat_chrom$vaf_segment_sample) & dat_chrom$is_het)
    if (length(is_phased_index) > 1) {
    seg_breaks = cumsum(rle(dat_chrom$vaf_segment_sample[is_phased_index])$lengths)
    for (i in 1:length(seg_breaks)) {
      if (i==1) {
        startpoint = 1
      } else {
        startpoint = seg_breaks[i-1]+1
      }
      endpoint = seg_breaks[i]
      segment_data = rbind(segment_data, data.frame(chrom=dat_chrom$chrom[is_phased_index[startpoint]], 
                                                    start=dat_chrom$pos[is_phased_index[startpoint]],
                                                    end=dat_chrom$pos[is_phased_index[endpoint]],
                                                    value=dat_chrom$vaf_segment_sample[is_phased_index[startpoint]],
                                                    value_control=dat_chrom$vaf_segment_control[is_phased_index[startpoint]], 
                                                    is_signif=dat_chrom$segment_is_signif[is_phased_index[startpoint]],
                                                    is_signif_adj=dat_chrom$segment_is_signif_adj[is_phased_index[startpoint]], 
						    num_het_snps=sum(dat_chrom$pos >= dat_chrom$pos[is_phased_index[startpoint]] & dat_chrom$pos <= dat_chrom$pos[is_phased_index[endpoint]] & dat_chrom$is_het),
						    stringsAsFactors=F))
      # add a second segment to make data look symmetrical
      # it overplots when there is no significance
      segment_data = rbind(segment_data, data.frame(chrom=dat_chrom$chrom[is_phased_index[startpoint]], 
                                                    start=dat_chrom$pos[is_phased_index[startpoint]],
                                                    end=dat_chrom$pos[is_phased_index[endpoint]],
                                                    value=1-dat_chrom$vaf_segment_sample[is_phased_index[startpoint]], 
                                                    value_control=1-dat_chrom$vaf_segment_control[is_phased_index[startpoint]], 
                                                    is_signif=dat_chrom$segment_is_signif[is_phased_index[startpoint]],
                                                    is_signif_adj=dat_chrom$segment_is_signif_adj[is_phased_index[startpoint]], 
						    num_het_snps=sum(dat_chrom$pos >= dat_chrom$pos[is_phased_index[startpoint]] & dat_chrom$pos <= dat_chrom$pos[is_phased_index[endpoint]] & dat_chrom$is_het),
						    stringsAsFactors=F))
    }
    }
  }
  print(head(segment_data))
  segment_data = segment_data[!is.na(segment_data$is_signif), ]
  print(head(segment_data))
  p = p + geom_rect(data=segment_data[!segment_data$is_signif,], mapping=aes(xmin=start, xmax=end, ymin=value-0.01, ymax=value+0.01), fill="black")
  
  if (any(segment_data$is_signif)) {
    p = p + geom_rect(data=segment_data[segment_data$is_signif,], mapping=aes(xmin=start, xmax=end, ymin=ifelse((value-0.01) > 0, value-0.01, 0), ymax=ifelse((value+0.01) < 1, value+0.01, 1)), fill="orange")
  }
  if (any(segment_data$is_signif_adj)) {
    p = p + geom_rect(data=segment_data[segment_data$is_signif_adj,], mapping=aes(xmin=start, xmax=end, ymin=ifelse((value-0.01) > 0, value-0.01, 0), ymax=ifelse((value+0.01) < 1, value+0.01, 1)), fill="red")
  }
  
  write.table(segment_data, file=file.path(outdir, "summary_data", paste0(samplename, "_segmentation.txt")), quote=F, sep="\t", row.names=F)
  
} else {
  segment_data = NULL
}
png(file.path(outdir, "plot_all_data", "2", paste0(samplename, "_all_data_2.png")), height=1500, width=2000)
print(p)
dev.off()

print("Make condensed genome wide plot...")
p = ggplot() +
  #scale_colour_manual(values=mycolours, drop=F) + 
  # geom_segment(data=anno, mapping=aes(x=start, xend=end, y=0, yend=1), colour="darkgrey") + # this was used to draw targeted regions
  geom_hline(yintercept=0.5, linetype=2) + 
  # geom_point(data=dat_hom, mapping=aes(x=pos, y=vaf_sample), colour="grey", size=0.5) + 
  geom_point(data=dat_het, mapping=aes(x=pos, y=vaf_sample), colour="darkgrey", size=0.5) + 
  scale_x_continuous(expand = c(0, 0)) +
  # scale_x_continuous(breaks=seq(min(anno$start), max(anno$end), 10000000), labels=paste0(round(seq(min(anno$start), max(anno$end), 10000000)/1000000), "Mb")) +
  theme_bw() + xlab("Position") + ylab("BAF") + ylim(0,1) +
  #annotate(geom="label", x=90000000, y=1, label=samplename) + 
  # facet_wrap(~samplename, ncol=2, strip.position="right") +
  theme(axis.text.x = element_blank(), #element_text(colour="black",size=14,face="plain"),
        axis.title.x = element_text(colour="black",size=16,face="plain"),
        axis.text.y = element_text(colour="black",size=14,face="plain"),
        axis.title.y = element_text(colour="black",size=16,face="plain"),
        strip.text.x = element_text(colour="black",size=16,face="plain"),
        strip.text.y = element_text(colour="black",size=15,face="plain"),
        legend.position = "none") + 
  facet_grid(~chrom, scales="free_x", space = "free")

if (!is.null(segment_data)) {
  p = p + geom_rect(data=segment_data[!segment_data$is_signif,], mapping=aes(xmin=start, xmax=end, ymin=value-0.01, ymax=value+0.01), fill="black")
  
  if (any(segment_data$is_signif)) {
    p = p + geom_rect(data=segment_data[segment_data$is_signif,], mapping=aes(xmin=start, xmax=end, ymin=value-0.01, ymax=value+0.01), fill="orange")
  }
  if (any(segment_data$is_signif_adj)) {
    p = p + geom_rect(data=segment_data[segment_data$is_signif_adj,], mapping=aes(xmin=start, xmax=end, ymin=value-0.01, ymax=value+0.01), fill="red")
  }
}
png(file.path(outdir, "plot_all_data", "1", paste0(samplename, "_all_data_1.png")), height=300, width=1500)
# png(file.path("~/Desktop/", paste0(samplename, "_all_data.png")), height=500, width=2000)
print(p)
dev.off()

p2 = ggplot() +
  #scale_colour_manual(values=mycolours, drop=F) + 
  # geom_segment(data=anno, mapping=aes(x=start, xend=end, y=0, yend=1), colour="darkgrey") + # this was used to draw targeted regions
  geom_hline(yintercept=0.5, linetype=2) + 
  # geom_point(data=dat_hom, mapping=aes(x=pos, y=vaf_sample), colour="grey", size=0.5) + 
  geom_point(data=dat_het, mapping=aes(x=pos, y=vaf_control), colour="darkgrey", size=0.5) + 
  scale_x_continuous(expand = c(0, 0)) +
  # scale_x_continuous(breaks=seq(min(anno$start), max(anno$end), 10000000), labels=paste0(round(seq(min(anno$start), max(anno$end), 10000000)/1000000), "Mb")) +
  theme_bw() + xlab("Position") + ylab("BAF") + ylim(0,1) +
  #annotate(geom="label", x=90000000, y=1, label=samplename) + 
  # facet_wrap(~samplename, ncol=2, strip.position="right") +
  theme(axis.text.x = element_blank(), #element_text(colour="black",size=14,face="plain"),
        axis.title.x = element_text(colour="black",size=16,face="plain"),
        axis.text.y = element_text(colour="black",size=14,face="plain"),
        axis.title.y = element_text(colour="black",size=16,face="plain"),
        strip.text.x = element_text(colour="black",size=16,face="plain"),
        strip.text.y = element_text(colour="black",size=15,face="plain"),
        legend.position = "none") + 
  facet_grid(~chrom, scales="free_x", space = "free")

if (!is.null(segment_data)) {
  p2 = p2 + geom_rect(data=segment_data, mapping=aes(xmin=start, xmax=end, ymin=value_control-0.01, ymax=value_control+0.01), fill="black")
}
p = p + ggtitle(samplename)
p2 = p2 + ggtitle("Control")
png(file.path(outdir, "plot_all_data", "3", paste0(samplename, "_all_data_3.png")), height=650, width=2000)
grid.arrange(p, p2, ncol=1)
dev.off()

# panel.border = element_blank(),
# axis.line = element_line(colour = "black"))

print("Make NOTCH1 plot...")
if (species=="human") {
  # detailed plot of NOTCH1 - 9	139388896	139440314
  dat_s = dat[as.character(dat$chrom)=="9" & dat$pos > 139388896 & dat$pos < 139440314,]
} else {
  # detailed plot of NOTCH1 - 2	26457903 26516663 
  dat_s = dat[as.character(dat$chrom)=="2" & dat$pos > 26457903 & dat$pos < 26516663 ,]
}
write.table(dat_s, file=file.path(outdir, "summary_data", paste0(samplename, "_data_notch1.txt")), quote=F, sep="\t", row.names=F)

p1 = ggplot() +
  scale_colour_manual(values=mycolours, drop=F) + 
  # geom_segment(data=anno, mapping=aes(x=start, xend=end, y=0, yend=1), colour="darkgrey") + # this can be used to force drawing of whole chromosome
  geom_hline(yintercept=0.5, linetype=2) + 
  geom_point(data=dat_s[!dat_s$is_het,], mapping=aes(x=pos, y=vaf_sample), colour="grey", size=2) + 
  geom_point(data=dat_s[dat_s$is_het,], mapping=aes(x=pos, y=vaf_sample, colour=is_signif), size=2) + 
  scale_x_continuous(expand = c(0, 0)) +
  # scale_x_continuous(breaks=seq(min(anno$start), max(anno$end), 10000000), labels=paste0(round(seq(min(anno$start), max(anno$end), 10000000)/1000000), "Mb")) +
  theme_bw() + xlab("Position") + ylab("BAF") + ylim(-0.01,1.01) +
  #annotate(geom="label", x=90000000, y=1, label=samplename) + 
  # facet_wrap(~samplename, ncol=2, strip.position="right") +
  ggtitle("BAF position") +
  theme(axis.text.x = element_text(colour="black",size=14,face="plain"),
        axis.title.x = element_text(colour="black",size=16,face="plain"),
        axis.text.y = element_text(colour="black",size=14,face="plain"),
        axis.title.y = element_text(colour="black",size=16,face="plain"),
        strip.text.x = element_text(colour="black",size=16,face="plain"),
        strip.text.y = element_text(colour="black",size=15,face="plain"),
        plot.title = element_text(colour="black",size=16,face="plain", hjust=0.5),
        legend.position = "none")

if (!is.null(segment_data)) {
  if (species=="human") {
    segment_data_s = segment_data[segment_data$chrom=="9" & (segment_data$start > 139388896 & segment_data$end < 139440314) | 
                                    segment_data$chrom=="9" & (segment_data$start <= 139388896 & segment_data$end > 139388896) |
                                    segment_data$chrom=="9" & (segment_data$start <= 139440314 & segment_data$end > 139440314),,drop=F]
    # remap start/end of most first and last segment, as we only want to see NOTCH1
    segment_data_s$start = ifelse(segment_data_s$start < 139388896, 139388896, segment_data_s$start)
    segment_data_s$end = ifelse(segment_data_s$end > 139440314, 139440314, segment_data_s$end)
    segment_data_s$num_het_snps = sum(dat$chrom=="9" & dat$pos >= 139388896 & dat$pos <= 139440314 & dat$is_het)
  } else {
      # detailed plot of NOTCH1 - 2	26457903 26516663 
      segment_data_s = segment_data[segment_data$chrom=="2" & (segment_data$start > 26457903 & segment_data$end < 26516663) | 
                                      segment_data$chrom=="2" & (segment_data$start <= 26457903 & segment_data$end > 26457903) |
                                      segment_data$chrom=="2" & (segment_data$start <= 26516663 & segment_data$end > 26516663),,drop=F]
      # remap start/end of most first and last segment, as we only want to see NOTCH1
      segment_data_s$start = ifelse(segment_data_s$start < 26457903, 26457903, segment_data_s$start)
      segment_data_s$end = ifelse(segment_data_s$end > 26516663, 26516663, segment_data_s$end)

      if (nrow(segment_data_s) > 0) {
          segment_data_s$num_het_snps = sum(dat$chrom=="2" & dat$pos >= 26457903 & dat$pos <= 26516663 & dat$is_het, na.rm=T)
      }
  }
  write.table(segment_data_s, file=file.path(outdir, "summary_data", paste0(samplename, "_segmentation_notch1.txt")), quote=F, sep="\t", row.names=F)
  
  
  print(segment_data_s)
  if (any(segment_data_s$is_signif)) {
    p1 = p1 + geom_rect(data=segment_data_s[segment_data_s$is_signif,], mapping=aes(xmin=start, xmax=end, ymin=value-0.005, ymax=value+0.005), fill="orange")
  }
  if (any(segment_data_s$is_signif_adj)) {
    p1 = p1 + geom_rect(data=segment_data_s[segment_data_s$is_signif_adj,], mapping=aes(xmin=start, xmax=end, ymin=value-0.005, ymax=value+0.005), fill="red")
  }
  if (any(!segment_data_s$is_signif)) {
    p1 = p1 + geom_rect(data=segment_data_s[!segment_data_s$is_signif,], mapping=aes(xmin=start, xmax=end, ymin=value-0.005, ymax=value+0.005), fill="black")
  }
}

p2 = ggplot() +
  scale_colour_manual(values=mycolours, drop=F) + 
  geom_abline(intercept=0, slope=1, linetype=2) +
  geom_point(data=dat_s[!dat_s$is_het,], mapping=aes(x=vaf_control, y=vaf_sample), colour="grey", size=2) + 
  geom_point(data=dat_s[dat_s$is_het,], mapping=aes(x=vaf_control, y=vaf_sample, colour=is_signif), size=2) + 
  ggtitle("BAF control") +
  theme_bw() + xlab("BAF - control") + ylab("BAF - biopsy") + ylim(0,1) + xlim(0,1) +
  theme(axis.text.x = element_text(colour="black",size=14,face="plain"),
        axis.title.x = element_text(colour="black",size=16,face="plain"),
        axis.text.y = element_text(colour="black",size=14,face="plain"),
        axis.title.y = element_text(colour="black",size=16,face="plain"),
        strip.text.x = element_text(colour="black",size=16,face="plain"),
        strip.text.y = element_text(colour="black",size=15,face="plain"),
        plot.title = element_text(colour="black",size=16,face="plain", hjust=0.5),
        legend.position = "none")
# png(file.path(outdir, "plot_controlVSsample_notch1", paste0(samplename, "_sampleVScontrol_notch1.png")), width=557, height=557)
# print(p)
# dev.off()

dat$depth_sample = dat$n_a0_sample+dat$n_a1_sample
dat_s$depth_sample = dat_s$n_a0_sample+dat_s$n_a1_sample
p3 = ggplot() + 
  geom_density(data=dat, mapping=aes(x=depth_sample, y=..count..), colour="blue") + 
  geom_histogram(data=dat_s, mapping=aes(x=depth_sample, y=..count..), binwidth=2, fill="grey", colour="red") + 
  geom_vline(xintercept=mean(dat$depth_sample), colour="blue", linetype=2, size=1.2) +
  geom_vline(xintercept=mean(dat_s$depth_sample), colour="red", linetype=2, size=1.2) +
  xlab("Depth biopsy") + ylab("Count") + xlim(0, max(quantile(dat$depth_sample, 0.9), quantile(dat_s$depth_sample, 0.9))) +
  ggtitle(paste0("Depth - Blue: all SNPs, Red: NOTCH1 - Means as line")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black",size=12,face="plain"),
        axis.title.x = element_text(colour="black",size=14,face="plain"),
        axis.text.y = element_text(colour="black",size=12,face="plain"),
        axis.title.y = element_text(colour="black",size=14,face="plain"),
        strip.text.x = element_text(colour="black",size=16,face="plain"),
        strip.text.y = element_text(colour="black",size=15,face="plain"),
        plot.title = element_text(colour="black",size=16,face="plain", hjust=0.5),
        legend.position = "none")
# png(file.path(outdir, "plot_depth", paste0(samplename, "_depth.png")), width=557, height=557)
# print(p)
# dev.off()

png(file.path(outdir, "plot_detail_notch1", paste0(samplename, "_detail_notch1.png")), width=1671, height=557)
grid.arrange(p1, p2, p3, ncol=3, top=samplename)
dev.off()

print("Make chrom plot...")
if (species=="human") {
  dat_s = dat[as.character(dat$chrom)=="9",]
} else {
  dat_s = dat[as.character(dat$chrom)=="2",]
}
p1 = ggplot() +
  scale_colour_manual(values=mycolours, drop=F) + 
  # geom_segment(data=anno, mapping=aes(x=start, xend=end, y=0, yend=1), colour="darkgrey") + # this was used to draw targeted regions
  geom_hline(yintercept=0.5, linetype=2) + 
  geom_point(data=dat_s[!dat_s$is_het,], mapping=aes(x=pos, y=vaf_sample), colour="grey", size=0.7) + 
  geom_point(data=dat_s[dat_s$is_het,], mapping=aes(x=pos, y=vaf_sample, colour=is_signif), size=0.7) + 
  scale_x_continuous(expand = c(0, 0)) +
  # scale_x_continuous(breaks=seq(min(anno$start), max(anno$end), 10000000), labels=paste0(round(seq(min(anno$start), max(anno$end), 10000000)/1000000), "Mb")) +
  theme_bw() + xlab("Position") + ylab("BAF") + ylim(0,1) +
  #annotate(geom="label", x=90000000, y=1, label=samplename) + 
  # facet_wrap(~samplename, ncol=2, strip.position="right") +
  theme(axis.text.x = element_text(colour="black",size=14,face="plain"),
        axis.title.x = element_text(colour="black",size=16,face="plain"),
        axis.text.y = element_text(colour="black",size=14,face="plain"),
        axis.title.y = element_text(colour="black",size=16,face="plain"),
        strip.text.x = element_text(colour="black",size=16,face="plain"),
        strip.text.y = element_text(colour="black",size=15,face="plain"),
        legend.position = "none",
        plot.title = element_text(colour="black",size=16,face="plain", hjust=0.5))

if (!is.null(segment_data)) {
  segment_data_s = segment_data[segment_data$chrom=="9",]
  p1 = p1 + geom_rect(data=segment_data_s[!segment_data_s$is_signif,], mapping=aes(xmin=start, xmax=end, ymin=value-0.005, ymax=value+0.02), fill="black") + 
    geom_rect(data=segment_data_s[segment_data_s$is_signif,], mapping=aes(xmin=start, xmax=end, ymin=value-0.005, ymax=value+0.02), fill="red")
}

dat_s$depth_sample = dat_s$n_a0_sample+dat_s$n_a1_sample
dat_s$depth_control = dat_s$n_a0_control+dat_s$n_a1_control

p2 = ggplot() +
  scale_colour_manual(values=mycolours, drop=F) + 
  # geom_segment(data=anno, mapping=aes(x=start, xend=end, y=0, yend=1), colour="darkgrey") + # this was used to draw targeted regions
  # geom_hline(yintercept=0.5, linetype=2) + 
  geom_point(data=dat_s, mapping=aes(x=pos, y=depth_sample), colour="black", size=0.7) + 
  # geom_point(data=dat_s[dat_s$is_het,], mapping=aes(x=pos, y=vaf_sample, colour=is_signif), size=0.7) + 
  # scale_x_continuous(breaks=seq(min(anno$start), max(anno$end), 10000000), labels=paste0(round(seq(min(anno$start), max(anno$end), 10000000)/1000000), "Mb")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw() + xlab("Position") + ylab("Depth") + #ylim(0,1) +
  #annotate(geom="label", x=90000000, y=1, label=samplename) + 
  # facet_wrap(~samplename, ncol=2, strip.position="right") +
  theme(axis.text.x = element_text(colour="black",size=14,face="plain"),
        axis.title.x = element_text(colour="black",size=16,face="plain"),
        axis.text.y = element_text(colour="black",size=14,face="plain"),
        axis.title.y = element_text(colour="black",size=16,face="plain"),
        strip.text.x = element_text(colour="black",size=16,face="plain"),
        strip.text.y = element_text(colour="black",size=15,face="plain"),
        legend.position = "none",
        plot.title = element_text(colour="black",size=16,face="plain", hjust=0.5))

if (species=="human") {
  png(file=file.path(outdir, "plot_chr9", paste0(samplename, "_chr9.png")), height=250, width=1000)
  print(p1 + ggtitle(paste0(samplename, " - chromosome 9")))
  dev.off()
} else {
  png(file=file.path(outdir, "plot_chr2", paste0(samplename, "_chr2.png")), height=250, width=1000)
  print(p1 + ggtitle(paste0(samplename, " - chromosome 2")))
  dev.off()
}



library(GenomicRanges)
infiles = list.files("output/", recursive=T, full.names=T, pattern="segmentation_notch1")

dat = lapply(infiles, function(x) { 
temp=read.table(x, header=T, stringsAsFactors=F); 
if (nrow(temp)==0) { return(temp) };
temp=temp[seq(1, nrow(temp), 2),,drop=F]; 

fulldat = unique(makeGRangesFromDataFrame(read.table(gsub("_notch1", "", x), header=T, stringsAsFactors=F)[, c("chrom", "start", "end", "num_het_snps")], keep.extra.columns=T))
temp_gr = makeGRangesFromDataFrame(temp)
overlap = findOverlaps(temp_gr, fulldat)
#if (length(overlap)==1) {
temp$start = start(fulldat)[subjectHits(overlap)]
temp$end = end(fulldat)[subjectHits(overlap)]
temp$num_het_snps_notch1 = temp$num_het_snps
temp$num_het_snps = fulldat$num_het_snps[subjectHits(overlap)]
#} else  {
#print(x)
#}


sam=unlist(strsplit(basename(x), "_"))[1];
temp$samplename=sam;

if (length(overlap)>1) {
print(overlap)
print(temp)
}


temp
})


dat = do.call(rbind, dat)
write.table(dat, file="oll_utput_combined.txt", quote=F, sep="\t", row.names=F)

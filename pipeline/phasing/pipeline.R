library(Battenberg)
library(optparse)
library(parallel)
library(foreach)
library(doParallel)

option_list = list(
  make_option(c("-t", "--tumourname"), type="character", default=NULL, help="Samplename of the tumour", metavar="character"),
  make_option(c("-n", "--normalname"), type="character", default=NULL, help="Samplename of the normal", metavar="character"),
  make_option(c("--tb"), type="character", default=NULL, help="Tumour BAM file", metavar="character"),
  make_option(c("--nb"), type="character", default=NULL, help="Normal BAM file", metavar="character"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Directory where output will be written", metavar="character"),
  make_option(c("--skip_allelecount"), type="logical", default=FALSE, action="store_true", help="Provide when alleles don't have to be counted. This expects allelecount files on disk", metavar="character"),
  make_option(c("--cpu"), type="numeric", default=8, help="The number of CPU cores to be used by the pipeline (Default: 8)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

tumourbam.file = opt$tb
normalbam.file = opt$nb
tumourname = opt$tumourname
normalname = opt$normalname
is.male = opt$sex=="male" | opt$sex=="Male"
skip_allelecount = opt$skip_allelecount
run_dir = opt$output
nthreads = opt$cpu

config = read.table("../config.txt", header=T, stringsAsFactors=F)
g1000allelesprefix = config$value[config$item=="human_g1000allelesprefix"]
imputeinfofile = config$value[config$item=="human_imputeinfofile"]
problemloci = config$value[config$item=="human_problemloci"]
#g1000allelesprefix = "/gpfs/nobackup/gerstung/sdentro/reference/human/battenberg/battenberg_1000genomesloci_phase3/1000genomesloci2012_chr"
#imputeinfofile = "/gpfs/nobackup/gerstung/sdentro/reference/human/battenberg/battenberg_impute_phase3/impute_info.txt"
#problemloci = "/gpfs/nobackup/gerstung/sdentro/reference/human/battenberg/battenberg_probloci/probloci_270415.txt.gz"

usebeagle = TRUE
           beaglemaxmem=10
           beaglenthreads=1
           beaglewindow=40
           beagleoverlap=4
        GENOME_VERSION = "b37"
        GENOMEBUILD = "hg19"
        BEAGLE_BASEDIR = "/hps/research/gerstung/sdentro/reference/human/battenberg/battenberg_beagle"
        beaglejar = file.path(BEAGLE_BASEDIR, "beagle.24Aug19.3e8.jar")
        beagleref.template = file.path(BEAGLE_BASEDIR, GENOME_VERSION, "chrCHROMNAME.1kg.phase3.v5a.b37.bref3")
        beagleplink.template = file.path(BEAGLE_BASEDIR, GENOME_VERSION, "plink.chrCHROMNAME.GRCh37.map")


min_normal_depth = 15
segmentation_gamma = 5
phasing_gamma = 1
segmentation_kmin = 3
phasing_kmin = 1
calc_seg_baf_option = 1

if (!dir.exists(file.path(run_dir, tumourname))) {
	dir.create(file.path(run_dir, tumourname), recursive=T)
}
setwd(file.path(run_dir, tumourname))
print(getwd())
chrom_names = Battenberg:::get.chrom.names(imputeinfofile=imputeinfofile, is.male=is.male)
print(chrom_names)
#if (F) {
clp = parallel::makeCluster(nthreads)
doParallel::registerDoParallel(clp)

foreach::foreach(i=1:length(chrom_names)) %dopar% {
#for (i in 1:length(chrom_names)) {
  chrom = chrom_names[i]
  print(chrom)
  print(file.exists(paste(tumourname, "_alleleFrequencies_chr", i, ".txt", sep="")))
  if (!file.exists(paste(tumourname, "_alleleFrequencies_chr", i, ".txt", sep=""))) {
    Battenberg:::getAlleleCounts(tumourbam.file, 
                    output.file=paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep=""), 
                    g1000.loci=paste(g1000allelesprefix, i, ".txt", sep=""), 
                    min.base.qual=20, 
                    min.map.qual=35, 
                    allelecounter.exe="alleleCounter")
  } else { print("Skipping counting tumour reads") }
  if (!file.exists(paste(normalname, "_alleleFrequencies_chr", i, ".txt", sep=""))) {
    Battenberg:::getAlleleCounts(normalbam.file, 
                    output.file=paste(normalname, "_alleleFrequencies_chr", i, ".txt", sep=""), 
                    g1000.loci=paste(g1000allelesprefix, i, ".txt", sep=""), 
                    min.base.qual=20, 
                    min.map.qual=35, 
                    allelecounter.exe="alleleCounter")
  } else { print("Skipping counting normal reads") }
 print("starting haplotyping") 
  Battenberg:::run_haplotyping(chrom=chrom, 
                   tumourname=tumourname, 
                   normalname=normalname, 
                   ismale=is.male, 
                   imputeinfofile=imputeinfofile, 
                   problemloci=problemloci, 
                   impute_exe="impute2", 
                   min_normal_depth=min_normal_depth, 
                   chrom_names=chrom_names,
                   snp6_reference_info_file=NA, heterozygousFilter=NA,
		   usebeagle=usebeagle,
                        beaglejar=beaglejar,
                        beagleref=gsub("CHROMNAME", chrom, beagleref.template),
                        beagleplink=gsub("CHROMNAME", chrom, beagleplink.template),
                        beaglemaxmem=beaglemaxmem,
                        beaglenthreads=beaglenthreads,
                        beaglewindow=beaglewindow,
                        beagleoverlap=beagleoverlap,
                        externalhaplotypeprefix=NA,
                        use_previous_imputation=F)

}

parallel::stopCluster(clp)
#}
Battenberg:::combine.baf.files(inputfile.prefix=paste(tumourname, "_chr", sep=""), 
                  inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt", 
                  outputfile=paste(tumourname, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
		  chr_names=chrom_names)
#chr_names=chrom_names)

#Battenberg:::segment.baf.phased(samplename=tumourname,
#                   inputfile=paste(tumourname, "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
#                   outputfile=paste(tumourname, ".BAFsegmented.txt", sep=""),
#                   gamma=segmentation_gamma,
#                   phasegamma=phasing_gamma,
#                   kmin=segmentation_kmin,
#                   phasekmin=phasing_kmin,
#                   calc_seg_baf_option=calc_seg_baf_option)

Battenberg:::segment.baf.phased.legacy(samplename=tumourname, 
	inputfile=paste(tumourname, "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
	outputfile=paste(tumourname, ".BAFsegmented.txt", sep=""), 
	gamma=segmentation_gamma, 
	phasegamma=phasing_gamma,
	kmin=segmentation_kmin, 
	phasekmin=phasing_kmin)

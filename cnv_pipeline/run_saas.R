options <- commandArgs(trailingOnly = T)
sample_id <- options[1]
sample_dir <- options[2]
seq_path <- options[3]
min.snps <- as.numeric(options[4])
max.chpts <- as.numeric(options[5])
use.null.data <- as.logical(options[6])
merge.pvalue.cutoff <- as.numeric(options[7])
cnvcall.pvalue.cutoff <- as.numeric(options[8])

source('NGS.CNV_mod.R')

# dirname = paste0('saas_minsnps',min.snps, '_chpts', max.chpts, '_null', use.null.data, '_mergep', merge.pvalue.cutoff, '_callp', cnvcall.pvalue.cutoff)
dirname = 'saasCNV_results'
# outputLoc = options[2]
# thresh = as.numeric(options[3])
# chrom = unlist(strsplit(options[5],","))

library(feather)
setwd(sample_dir)
# vcf_feather <- paste0(sample_id, '_bac.feather')  # BAF feather path
# seq_path <- 'seq.feather'  # BAF feather path
seq_data <- read_feather(seq_path)
# vcf_table$CHROM <- paste0('chr', vcf_table$CHROM)  # Add chr to chrom column

seq_data$pos <- seq_data$position

## NGS pipeline analysis
output.dir <- file.path(sample_dir, dirname)
library(saasCNV)
NGS.CNV_sgg(seq.data=seq_data, output.dir=output.dir, sample.id=sample_id,
        min.chr.probe=100,
        min.snps=min.snps,
        joint.segmentation.pvalue.cutoff=1e-4,
        max.chpts=max.chpts,
        do.merge=TRUE, use.null.data=use.null.data, num.perm=1000, #maxL=2000,
        merge.pvalue.cutoff=merge.pvalue.cutoff,
        do.cnvcall.on.merge=TRUE,
        cnvcall.pvalue.cutoff=cnvcall.pvalue.cutoff,
        do.plot=TRUE, cex=0.3, ref.num.probe=1000,
	do.gene.anno=FALSE, seed=123456789,
	verbose=TRUE)

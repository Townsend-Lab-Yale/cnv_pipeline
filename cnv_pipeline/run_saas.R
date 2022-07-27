options <- commandArgs(trailingOnly = T)
sample_id <- options[1]
sample_dir <- options[2]
seq_path <- options[3]
min.snps <- as.numeric(options[4])
max.chpts <- as.numeric(options[5])
use.null.data <- as.logical(options[6])
merge.pvalue.cutoff <- as.numeric(options[7])
cnvcall.pvalue.cutoff <- as.numeric(options[8])

dirname <- 'saasCNV_results'

library(arrow)
setwd(sample_dir)
vcf.data <- read_parquet(seq_path)
vcf.data$CHROM <- paste0('chr', vcf.data$CHROM)  # Add chr to chrom column

## NGS pipeline analysis
output.dir <- file.path(sample_dir, dirname)
library(saasCNV)
NGS.CNV(vcf=vcf.data, output.dir=output.dir, sample.id=sample_id,
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

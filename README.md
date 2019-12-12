# cnv_pipeline

Easily run ADTEx and saasCNV analyses on sequence data. Works with multi-sample GATK vcf files.

## Installation

Install into your python distribution (tested on python3).

```bash
pip install git+https://github.com/sggaffney/cnv_pipeline.git
```

Install required R packages, including saasCNV.

```R
source("https://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
install.packages("RANN")
install.packages("saasCNV")
install.packages("feather")
```

Install [ADTEx](https://sourceforge.net/projects/adtex/).

## Example command line usage

The only required arguments are tumor and normal bam files, corresponding sample names, a vcf file that includes those sample names, and an output directory.

```bash
run_cnv -t bams/23986.bam -n bams/23986N.bam \
  -v exome_calls.vcf \
  -s cnv_results -tid 23986 -nid 23986N
```

Creates the following output:
```
cnv_results
├── adtex_output/
├── baf.txt
├── genome.txt
├── normal_cov.bed
├── saasCNV_results/
├── saas.feather
├── snps_trimmed.vcf
└── tumor_cov.bed
```

Command line help
```
$ run_cnv -h
usage: CNV PIPELINE [-h] [-v VCF] -s SAMPLE_DIR -t TUMOR -n NORMAL
                    [-a ADTEX_DIR] -tid TUMOR_ID [-nid NORMAL_ID]
                    [-rmin RATIO_MIN] [-rmax RATIO_MAX] [-tmin MIN_TUMOR]
                    [-nmin MIN_NORMAL] [-gq MIN_GQ] [--ploidy PLOIDY]
                    [--minReadDepth MINREADDEPTH] [-ao ADTEX_STDOUT]

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     VCF file for sample pair [REQUIRED]
  -s SAMPLE_DIR, --sample_dir SAMPLE_DIR
                        Sample-specific intermediate output dir [OPTIONAL]
  -t TUMOR, --tumor TUMOR
                        Tumor BAM [REQUIRED]
  -n NORMAL, --normal NORMAL
                        Normal BAM [REQUIRED]
  -a ADTEX_DIR, --adtex_dir ADTEX_DIR
                        ADTEx output dir
  -tid TUMOR_ID, --tumor_id TUMOR_ID
                        Tumor name, for vcf extraction
  -nid NORMAL_ID, --normal_id NORMAL_ID
                        Normal name, for vcf extraction
  -rmin RATIO_MIN, --ratio_min RATIO_MIN
                        Min alt ratio in normal sample [0.4]
  -rmax RATIO_MAX, --ratio_max RATIO_MAX
                        Max alt ratio in normal sample [0.6]
  -tmin MIN_TUMOR, --min_tumor MIN_TUMOR
                        Min reads from tumor sample [20]
  -nmin MIN_NORMAL, --min_normal MIN_NORMAL
                        Min reads from normal sample [10]
  -gq MIN_GQ, --min_gq MIN_GQ
                        Genotype quality cutoff for normal sample [90]
  --ploidy PLOIDY       Most common ploidy in the tumour sample
  --minReadDepth MINREADDEPTH
                        The threshold for minimum read depth for each exon
                        [10]
  -ao ADTEX_STDOUT, --adtex_stdout ADTEX_STDOUT
                        ADTEx stdout path

```

## Configuration

You can modify the config.ini file in the python site-packages directory, which is set up for Yale's Ruddle cluster
(October 2016):
e.g. path_to_python_distribution/lib/python3.5/site-packages/cnv_pipeline/config.ini

```
[paths]
CODING_REGIONS = /ycga-ba/home/bioinfo/software/knightlab/genomes/hs37d5/bed_files/hs37d5_refgene_coding_Nov2015.bed
FASTA = /ycga-ba/home/bioinfo/software/knightlab/genomes/hs37d5/human_g1k_v37_decoy.fasta
PYTHON2 = /ycga-gpfs/apps/hpc/Langs/Python/2.7.11/bin/python
SAMTOOLS = /ycga-gpfs/apps/hpc/Apps/SAMtools/1.3/bin/samtools
ADTEX = /ycga-gpfs/project/fas/townsend/software/ADTEx.v.2.0/ADTEx_sgg.py
BEDTOOLS = /ycga-gpfs/project/fas/townsend/software/bedtools-2.26.0/bedtools

[gatk]
cmd = java -Xmx12g -XX:ParallelGCThreads=8 -jar /home/ky89/ngstools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar
```

# cnv_pipeline

Easily run ADTEx and saasCNV analyses on sequence data. Works with multi-sample GATK vcf files.

## Installation

Install into your python distribution (tested on python3).

```bash
pip install git+https://github.com/sggaffney/cnv_pipeline.git
```

### R setup

- ADTEx requires `wmtsa`, which was removed from CRAN in 2020. The following 
    instructions assume use of R 4.0.3, and specify the latest MRAN snapshot that
    includes `wmtsa`. If you don't want to run ADTEx, later versions of R and 
    CRAN could be used.
- `saasCNV` is a CRAN package, still in CRAN as of 2022-07-28, and required 
    dependencies are `DNAcopy` (from Bioconductor), and `RANN`.


```R
install.packages('BiocManager', repos='https://mran.microsoft.com/snapshot/2020-06-08')
r <- BiocManager::repositories()
r['CRAN'] = 'https://mran.microsoft.com/snapshot/2020-06-08'
options(repos=r)
install.packages(c('wmtsa', 'DNAcopy', 'RANN', 'saasCNV', 'arrow'), dependencies=TRUE)
```

- Manually install saasCNV using Rscript.
```shell
 Rscript -e 'install.packages("saasCNV", repos = "http://cran.us.r-project.org")'
```

### ADTEx installation

The 2016 v2.0.0 version of [ADTEx](https://sourceforge.net/projects/adtex/) requires 
python2. A python3-compatible version is available at https://github.com/Townsend-Lab-Yale/ADTEx.  


## Example command line usage

The only required arguments are tumor and normal bam files, corresponding sample 
names, a vcf file that includes those sample names, and an output directory.

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

You can modify the config.ini file in the python site-packages directory, which 
is set up for Yale's Ruddle / Farnam / McCleary clusters (October 2022):
e.g. path_to_python_distribution/lib/python3.5/site-packages/cnv_pipeline/config.ini

```
[paths]
CODING_REGIONS = /gpfs/gibbs/pi/ycga/mane/jk2269/knightlab/genomes/hs37d5/bed_files/hs37d5_refgene_coding_Nov2015.bed
FASTA = /gpfs/gibbs/pi/ycga/mane/jk2269/knightlab/genomes/hs37d5/human_g1k_v37_decoy.fasta
ADTEX = /gpfs/gibbs/pi/townsend/software/cnv_pipeline/ADTEx/ADTEx.py
SAMTOOLS = samtools
BEDTOOLS = bedtools
RSCRIPT = Rscript

[gatk]
cmd = java -Xmx12g -XX:ParallelGCThreads=8 -jar /home/ky89/ngstools/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar
```

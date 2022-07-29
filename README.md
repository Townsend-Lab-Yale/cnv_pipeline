# cnv_pipeline

Easily run ADTEx and saasCNV analyses on sequence data. Works with multi-sample 
VCF files.

## Installation

### Installing all dependencies other than GATK

#### The easy way: `conda`/`mamba` 

The `cnv.environment.yaml` file specifies all dependencies (other than GATK) 
needed to run this software. SAMtools, BEDTools, Python, R, and the necessary 
packages for Python and R will all be installed in an isolated environment that's
easy to activate and deactivate.

1. If you don't have `conda` on your system, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).  
2. Clone or [download](https://github.com/Townsend-Lab-Yale/cnv_pipeline/archive/main.zip) 
    the cnv_pipeline repository, and unzip it somewhere convenient for long-term storage
    —your Python environment will store a link to the folder as part of installation.
3. In your terminal, navigate to the folder and run the following to create a 
    conda virtual environment with all dependencies:
    ```bash
    conda env create -n cnv -f cnv.environment.yaml
    ```
    - or to create an environment precisely mimicking a tested version, specify 
    `cnv.environment.freeze.yaml` instead of `cnv.environment.yaml` above. 
4. Activate the new environment with `conda activate cnv`. (This environment
    will need to be activated each time you want to run the pipeline.)
5. Install the CRAN package `saasCNV` into the environment's R library: 
    ```bash
    Rscript -e 'install.packages("saasCNV", repos = "http://cran.us.r-project.org")'
    ```
6. Install the cnv_pipeline package by running `python setup.py develop`. This 
    will add a `run_cnv` CLI executable to your PATH.
7. Clone or [download and unzip](https://github.com/Townsend-Lab-Yale/ADTEx/archive/master.zip) 
    the [ADTEx repository](https://github.com/Townsend-Lab-Yale/ADTEx) into the top
    level of your cnv_pipeline folder. (Or place it elsewhere and add a symlink 
    called `ADTEx` pointing to the ADTEx directory.)


#### The harder way: manual setup

1. Make sure `samtools`, `bedtools`, `python` (≥3.6), and `Rscript` are on your path.
2. Install `cnv_pipeline` into your python environment:
```bash
pip install -e git+https://github.com/sggaffney/cnv_pipeline.git
```
3. Clone or [download and unzip](https://github.com/Townsend-Lab-Yale/ADTEx/archive/master.zip) 
    the [ADTEx repository](https://github.com/Townsend-Lab-Yale/ADTEx) into the top
    level of your cnv_pipeline folder. (Or place it elsewhere and add a symlink 
    called `ADTEx` pointing to the ADTEx directory.)
4. Now set up the R dependencies as outlined below. 

##### R dependencies

- ADTEx requires `wmtsa`, which was removed from CRAN in 2020. The following 
    instructions assume use of R 4.0, and specify the latest MRAN snapshot that
    includes `wmtsa`. If you don't want to run ADTEx, later versions of R and 
    CRAN could be used.
- `saasCNV` is a CRAN package, still in CRAN as of 2022-07-28, and other required 
    saasCNV dependencies are `DNAcopy` (from Bioconductor), and `RANN`.

```R
install.packages('BiocManager', repos='https://mran.microsoft.com/snapshot/2020-06-08')
r <- BiocManager::repositories()
r['CRAN'] = 'https://mran.microsoft.com/snapshot/2020-06-08'
options(repos=r)
install.packages(c('wmtsa', 'DNAcopy', 'RANN', 'saasCNV', 'arrow'), dependencies=TRUE)
```

- Manually install `saasCNV` using Rscript.
```shell
 Rscript -e 'install.packages("saasCNV", repos = "http://cran.us.r-project.org")'
```

### GATK

- The app uses the GATK tool `SelectVariants` to subset the input VCF file, extracting 
sites that are high quality, well-covered, and heterozygous in the normal sample.
A future version may switch to BCFtools to make life easier. In the meantime, install 
a copy of [GATK v4](https://github.com/broadinstitute/gatk) (tested on 4.2.6.1) 
or obtain a container image from 
[DockerHub](https://hub.docker.com/r/broadinstitute/gatk/).
- The app looks for an environment variable called `GATK_ALIAS` which contains a 
  command for running gatk, for use in `<GATK_ALIAS> SelectVariants [arguments]`.
  If the variable is not present, the app will assume your path contains an 
  executable/alias called `gatk`.


## Example command line usage

The required arguments are:
- tumor and normal bam files
- tumor and normal sample names
- a vcf file that includes the above sample names
    - NOTE: the too-simplistic VCF parser assumes FORMAT data starts with GT then AD
- a genome reference FASTA file
- a BED file with capture targets
- an output directory (which will be automatically created if necessary)

```bash
# EXAMPLE GATK ALIAS, using gatk within a singularity container
export GATK_ALIAS='singularity exec /path/to/gatk.sif gatk'

run_cnv -v variants.vcf  \
-t tumor1.bam -tid tumor1 \
-n normal1.bam -nid normal1 \
--sample_dir cnv_tumor1 \
--ref_fasta /path/to/genome.fasta \
--bed /path/to/targets.bed
```

Creates the following output:
```
cnv_tumor1
├── adtex_output/
├── baf.txt
├── genome.txt
├── normal_cov.bed
├── saasCNV_results/
├── saas.parquet
├── snps_trimmed.vcf
├── snps_trimmed.vcf.idx
└── tumor_cov.bed
```

Further arguments are listed in the help documentation of the run_cnv CLI:
```
$ run_cnv -h
usage: CNV PIPELINE [-h] -v VCF -s SAMPLE_DIR -t TUMOR_BAM -n NORMAL_BAM -tid TUMOR_ID
                    -nid NORMAL_ID -R REF_FASTA [--saas_only | --adtex_only]
                    [-rmin RATIO_MIN] [-rmax RATIO_MAX] [-tmin MIN_TUMOR]
                    [-nmin MIN_NORMAL] [-gq MIN_GQ] [-a ADTEX_DIR] [-b BED]
                    [--ploidy PLOIDY] [--minReadDepth MINREADDEPTH] [-ao ADTEX_STDOUT]

options:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     VCF file for sample pair
  -s SAMPLE_DIR, --sample_dir SAMPLE_DIR
                        Sample-specific output dir
  -t TUMOR_BAM, --tumor_bam TUMOR_BAM
                        Tumor BAM
  -n NORMAL_BAM, --normal_bam NORMAL_BAM
                        Normal BAM
  -tid TUMOR_ID, --tumor_id TUMOR_ID
                        Tumor name, for vcf extraction
  -nid NORMAL_ID, --normal_id NORMAL_ID
                        Normal name, for vcf extraction
  -R REF_FASTA, --ref_fasta REF_FASTA
                        Reference genome fasta path
  --saas_only           Only run saasCNV, not ADTEx.
  --adtex_only          Only run ADTEx, not saasCNV.
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
  -a ADTEX_DIR, --adtex_dir ADTEX_DIR
                        ADTEx: output dir
  -b BED, --bed BED     ADTEx: BED file for targeted regions [REQUIRED FOR ADTEx]
  --ploidy PLOIDY       ADTEx: most common ploidy in the tumour sample
  --minReadDepth MINREADDEPTH
                        The ADTEx threshold for minimum read depth for each exon [10]
  -ao ADTEX_STDOUT, --adtex_stdout ADTEX_STDOUT
                        ADTEx stdout path if overriding STDOUT
```

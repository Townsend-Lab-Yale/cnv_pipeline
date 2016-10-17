# WRITES saasCNV seq_data table as feather
library(feather)


options <- commandArgs(trailingOnly = T)
vcf_path <- options[1]  # input vcf path
out_feather <- options[2]  # output filename
tumor.col <- as.numeric(options[3])  # e.g. 10
normal.col <- as.numeric(options[4])  # e.g. 11
format.ad.indx <- as.numeric(options[5])  # e.g. 2
if (length(args) > 5) {
    chroms <- unlist(strsplit(options[6],","))
    } else {
    chroms <- c(1:22, "X", "Y")
}
if (length(args) > 6) {
    MQ.cutoff <- as.numeric(options[7])
    } else {
    MQ.cutoff <- 30
}


vcf2txt_sgg <- function (vcf.file, normal.col=11, tumor.col=10, MQ.cutoff=30,
    format.ad.indx=2, format.gt.indx=1, chrs=c(1:22, "X", "Y"))
{
  
  inputFile <- file(vcf.file, "r")
  while (length(lines <- readLines(inputFile, n = 1, warn = FALSE)) > 0) {
    if (length(grep("\\#CHROM", lines)) == 1) 
      break
  }
  header <- sub("\\#", "", lines)
  header <- strsplit(header, "\t")[[1]]
  header[c(normal.col, tumor.col)] <- c("Normal", "Tumor")
  vcf <- read.table(inputFile, header = FALSE, sep = "\t", 
                    as.is = TRUE, comment.char = "")
  names(vcf) <- header
  vcf <- vcf[, c(1:9, normal.col, tumor.col)]
  close(inputFile)
  # chrs <- paste0("", c(1:22, "X", "Y"))  # SGG: prev used chr prefix
  idx <- with(vcf, CHROM %in% chrs & FILTER == "PASS")
  vcf <- vcf[idx, ]
  idx <- grep(",", vcf$ALT)
  if (length(idx) >= 1) 
    vcf <- vcf[-idx, ]
  normal.RD <- strsplit(vcf$Normal, ":")
  normal.RD <- do.call(rbind, normal.RD)
  vcf$Normal.GT <- normal.RD[, format.gt.indx]
  normal.AD <- do.call(rbind, strsplit(normal.RD[, format.ad.indx], ","))  # SGG: orig indx 2
  vcf$Normal.REF.DP <- as.integer(normal.AD[, 1])
  vcf$Normal.ALT.DP <- as.integer(normal.AD[, 2])
  tumor.RD <- strsplit(vcf$Tumor, ":")
  tumor.RD <- do.call(rbind, tumor.RD)
  vcf$Tumor.GT <- tumor.RD[, format.gt.indx]
  tumor.AD <- do.call(rbind, strsplit(tumor.RD[, format.ad.indx], ","))
  vcf$Tumor.REF.DP <- as.integer(tumor.AD[, 1])
  vcf$Tumor.ALT.DP <- as.integer(tumor.AD[, 2])
  vcf <- vcf[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                 "INFO", "Normal.GT", "Normal.REF.DP", "Normal.ALT.DP", 
                 "Tumor.GT", "Tumor.REF.DP", "Tumor.ALT.DP")]
  idx <- which(vcf$Normal.GT == "./." | vcf$Tumor.GT == "./.")
  if (length(idx) >= 1) 
    vcf <- vcf[-idx, ]
  info <- strsplit(vcf$INFO, ";")
  extract.MQ <- function(x) {
    idx.MQ <- grep("^MQ\\=", x)
    if (length(idx.MQ) == 1) 
      return(x[idx.MQ])
    else return(NA)
  }
  info1 <- sapply(info, FUN = extract.MQ)
  vcf$MQ <- as.numeric(sub("^MQ\\=", "", info1))
  vcf <- vcf[vcf$MQ >= MQ.cutoff, ]
  vcf <- vcf[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                 "MQ", "Normal.GT", "Normal.REF.DP", "Normal.ALT.DP", 
                 "Tumor.GT", "Tumor.REF.DP", "Tumor.ALT.DP")]
  a <- vcf[, c("CHROM", "POS")]
  idx <- which(duplicated(a))
  if (length(idx) >= 1) 
    vcf <- vcf[-idx, ]
  return(vcf)
}

#chrs=c(1:22, "X", "Y")

# BUILD DATA TABLE
vcf_table <- vcf2txt_sgg(vcf_path, normal.col=normal.col, tumor.col=tumor.col,
    format.ad.indx=format.ad.indx, format.gt.indx=1, MQ.cutoff=MQ.cutoff,
    chrs=chroms)

# STRIP OUT ROWS WITH '.' GENOTYPE
vcf_table <- vcf_table[(vcf_table$Normal.GT!='.') & (vcf_table$Tumor.GT!='.'),]

# print(head(vcf_table))
write_feather(vcf_table, out_feather)

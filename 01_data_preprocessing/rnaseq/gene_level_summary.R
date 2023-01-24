library(tximport)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly=TRUE)
target_path <- args[1]

# List of quant.sf files for each sample
fpaths <- file.path(list.dirs(target_path, full.names = TRUE, recursive = FALSE),
                    "quant.sf")
names(fpaths) <- list.dirs(target_path, full.names = FALSE, recursive = FALSE)

# List the transcript to gene map
t2g <- read_csv(paste0(target_path, "/tx2gene.csv"), 
                col_names = c("txid", "geneid", "symbol")) %>% 
  select(txid, symbol)

# Perform the gene-level summarisation
txi <- tximport(fpaths, type = "salmon", tx2gene = t2g)

# Create a SummarizedExperiment object
library(SummarizedExperiment)
libSize <- colSums(txi$counts)
fpkm <- (txi$counts / (txi$length * matrix(libSize, 
                                          byrow = TRUE, 
                                          nrow(txi$counts), 
                                          ncol(txi$counts)))) * 1e9

se <- SummarizedExperiment(assays = list(counts = txi$counts,
                                         tpm = txi$abundance,
                                         fpkm = fpkm,
                                         length = txi$length))

saveRDS(se, file = paste0(target_path, "/se_fpkm_tpm.RDS"))


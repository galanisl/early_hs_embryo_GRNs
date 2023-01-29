
library(SingleCellExperiment)

load("final_data/sce_integrated_lateBlast_clean.RData")

epi <- sce[, sce$cell_type == "Epiblast"]
te <- sce[, sce$cell_type == "Trophectoderm"]
pe <- sce[, sce$cell_type == "Primitive endoderm"]

# TPM
epi_tpm <- assay(epi, "tpm")
te_tpm <- assay(te, "tpm")
pe_tpm <- assay(pe, "tpm")

saveRDS(epi_tpm, file = "final_data/epi_tpm.rds")
saveRDS(te_tpm, file = "final_data/te_tpm.rds")
saveRDS(pe_tpm, file = "final_data/pe_tpm.rds")

# logcounts
epi_lgc <- assay(epi, "logcounts")
te_lgc <- assay(te, "logcounts")
pe_lgc <- assay(pe, "logcounts")

saveRDS(epi_lgc, file = "final_data/epi_logcounts.rds")
saveRDS(te_lgc, file = "final_data/te_logcounts.rds")
saveRDS(pe_lgc, file = "final_data/pe_logcounts.rds")

# batch-corrected counts
epi_btc <- assay(epi, "batch_corrected")
te_btc <- assay(te, "batch_corrected")
pe_btc <- assay(pe, "batch_corrected")

saveRDS(epi_btc, file = "final_data/epi_batchcorr.rds")
saveRDS(te_btc, file = "final_data/te_batchcorr.rds")
saveRDS(pe_btc, file = "final_data/pe_batchcorr.rds")


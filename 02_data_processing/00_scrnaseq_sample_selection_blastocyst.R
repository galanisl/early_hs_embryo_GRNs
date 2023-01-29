library(dplyr)
library(ggplot2)
library(SummarizedExperiment)

# Filtering based on samples used for Wamaitha *et al.*, Nature Communications 2020
sissy <- readr::read_csv("gene_expr/samples_Wamaitha_NatComm2020.csv") %>% 
  filter(exp_group %in% c("Epi", "PE", "TE"))

meta_blk <- readr::read_csv("gene_expr/blakeley_sample_details.csv") %>%
  select(sample_accession, sample_title) %>%
  filter(!duplicated(sample_accession)) %>%
  mutate(cell_type = case_when(
    grepl("TE", sample_title) ~ "Trophectoderm",
    grepl("PE", sample_title) ~ "Primitive endoderm",
    grepl("EPI", sample_title) ~ "Epiblast"
  ), batch = "Blk")

yan_alias <- readr::read_tsv("gene_expr/yan_sample_alias.tsv") %>% 
  select(sample_accession, sample_alias)

yan_stirp <- readr::read_csv("gene_expr/stirparo_revlineages_lateBlast.csv") %>% 
  filter(grepl("Yan", Study)) %>% 
  filter(RevisedLineage %in% c("primitive_endoderm", "trophectoderm")) %>% 
  select(Cell, RevisedLineage)

meta_yan <- readr::read_csv("gene_expr/yan_sample_details.csv") %>%
  select(sample_accession, run_accession, simple_title) %>%
  filter(!duplicated(run_accession)) %>%
  mutate(original_ctype = case_when(
    grepl("LateBlast", simple_title) ~ "Blastocyst",
    TRUE ~ "Other"
  ), batch = "Yan") %>%
  filter(original_ctype != "Other") %>% 
  left_join(yan_alias, by = "sample_accession") %>%
  left_join(select(sissy, sample_name, exp_group), 
            by = c("sample_alias" = "sample_name")) %>% 
  left_join(yan_stirp, by = c("simple_title" = "Cell")) %>% 
  filter(!(is.na(exp_group) & is.na(RevisedLineage))) %>% 
  mutate(exp_group = paste(exp_group, RevisedLineage)) %>% 
  mutate(cell_type = case_when(
    grepl("trophectoderm", exp_group) ~ "Trophectoderm",
    grepl("primitive", exp_group) ~ "Primitive endoderm",
    grepl("Epi", exp_group) ~ "Epiblast"
  )) %>% 
  select(run_accession, simple_title, cell_type, batch) %>% 
  dplyr::rename(sample_accession = run_accession, 
                sample_title = simple_title)

pet_stirp <- readr::read_csv("gene_expr/stirparo_revlineages_lateBlast.csv") %>% 
  filter(grepl("Petro", Study)) %>% 
  filter(RevisedLineage == "trophectoderm") %>% 
  select(Cell, RevisedLineage)
pet_stirp <- pet_stirp[sample(1:nrow(pet_stirp), 20),]

pet_bad_epi <- readr::read_csv("pet_epi_bad.csv")

meta_pet <- readr::read_csv("gene_expr/petro_sample_details.csv") %>%
  select(run_accession, sample_title, fastq_ftp) %>%
  mutate(sample_name = sapply(strsplit(fastq_ftp, "/"),
                              function(x) x[length(x)])) %>%
  mutate(sample_name = sapply(strsplit(sample_name, "[.]"), `[`, 1)) %>%
  mutate(batch = "Pet") %>%
  select(sample_name, run_accession, batch) %>%
  dplyr::rename(sample_accession = sample_name, sample_title = run_accession) %>%
  left_join(select(sissy, sample_name, exp_group), 
            by = c("sample_accession" = "sample_name")) %>% 
  left_join(pet_stirp, by = c("sample_title" = "Cell")) %>% 
  filter(!(is.na(exp_group) & is.na(RevisedLineage))) %>% 
  mutate(exp_group = paste(exp_group, RevisedLineage)) %>% 
  mutate(cell_type = case_when(
    grepl("trophectoderm", exp_group) ~ "Trophectoderm",
    grepl("PE", exp_group) ~ "Primitive endoderm",
    grepl("Epi", exp_group) ~ "Epiblast"
  )) %>%
  select(sample_accession, sample_title, cell_type, batch) %>% 
  filter(!sample_accession %in% pet_bad_epi$sample_accession)


# Load SummarizedExperiment objects and keep samples from the meta tables only
blk <- readRDS("gene_expr/blakeley_salmon_fpkm_tpm.RDS")
blk <- blk[, meta_blk$sample_accession]

yan <- readRDS("gene_expr/yan_salmon_fpkm_tpm.RDS")
yan <- yan[, meta_yan$sample_accession]

pet <- readRDS("gene_expr/petro_salmon_fpkm_tpm.RDS")
pet <- pet[, meta_pet$sample_accession]

# Form a single metadata table and a single expression matrix
meta <- bind_rows(meta_blk, meta_yan, meta_pet)

cts <- as.matrix(cbind(assays(blk)$counts,
                       assays(yan)$counts, 
                       assays(pet)$counts))
tpm <- as.matrix(cbind(assays(blk)$tpm,
                       assays(yan)$tpm, 
                       assays(pet)$tpm))
fpkm <- as.matrix(cbind(assays(blk)$fpkm,
                       assays(yan)$fpkm, 
                       assays(pet)$fpkm))
lg <- as.matrix(cbind(assays(blk)$length,
                      assays(yan)$length, 
                      assays(pet)$length))

# Prepare a tibble of gene details
feats <- tibble(symbol = rownames(cts))
feats <- feats %>% 
  left_join(select(annotables::grch38, symbol, ensgene, entrez, chr, biotype)) %>% 
  filter(!duplicated(symbol))

# If samples in the metadata and the matrices are in the same order, create
# a SingleCellExperiment object
if(all(meta$sample_accession == colnames(cts)) & 
   all(feats$symbol == rownames(cts))){
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cts, 
                                                                  tpm = tpm, 
                                                                  fpkm = fpkm, 
                                                                  length = lg),
                                                    colData = meta, 
                                                    rowData = feats)
}else{
  stop("Fix the sample or gene order!!!") 
}

save(sce, file = "final_data/sce_integrated_lateBlast.RData")

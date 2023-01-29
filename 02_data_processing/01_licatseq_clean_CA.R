
# Clean TF-gene associations from ATAC-seq

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# Change to ICM or TE accordingly
ctype <- "licat_peaks/ann_ICM_binding_sites.csv"

tb <- read_csv(ctype) %>% 
  select(motif, nearest_prom, dist2tss)

tb <- tb %>% 
  mutate(tf = map_chr(str_split(motif, "_HUMAN"), `[`, 1)) %>% 
  mutate(tf = map_chr(str_split(tf, "\\.\\d\\."), 
                      function(x) x[length(x)])) %>% 
  mutate(tf = map_chr(str_split(tf, "\\(var"), `[`, 1))

# Duplicate rows in which more than one TF is part of a motif
tb_expanded <- tb %>% 
  separate_rows(tf, sep = "::")

dup_idx <- tb_expanded %>% 
  unite(col = "regInt", tf, nearest_prom, sep = "-")

tb_expanded <- tb_expanded[!duplicated(dup_idx$regInt),]

# Save the processed table
saveRDS(tb_expanded, file = "final_data/icm_tf_binding.rds")



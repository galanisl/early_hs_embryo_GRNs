
library(dplyr)
library(readxl)
library(tidyr)
library(purrr)
library(ggplot2)
library(viridis)
library(patchwork)

# Function definitions ----------------------------------------------------

read_organise <- function(fpath){
  normal <- sapply(strsplit(fpath, "V"), `[`, 2)
  normal <- sapply(strsplit(normal, "[.]"), `[`, 1)
  
  xl_sheets <- excel_sheets(fpath)
  ctypes <- sapply(strsplit(xl_sheets, "_TF"), `[`, 1)
  pred <- sapply(strsplit(xl_sheets, "_TF"), `[`, 2)
  xl_tb <- vector("list", length(pred))
  names(xl_tb) <- pred
  for(i in seq_along(xl_sheets)){
    xl_tb[[i]] <- read_excel(fpath, sheet = xl_sheets[i]) %>% 
      select(ends_with("p-val")) %>% 
      pivot_longer(cols = everything(), 
                   names_to = "gene", values_to = "pval") %>% 
      mutate(gene = sapply(strsplit(gene, " "), `[`, 1), ctype = ctypes[i], 
             method = pred[i], norm_method = normal)
  }
  xl_tb <- bind_rows(xl_tb)
  return(xl_tb)
}

# Load data into memory ---------------------------------------------------

# Path to the CVS files with the results
dpath <- "path/to/GRN/inferences"
fnames <- list.files(dpath, pattern = "^TFnetV")

Vstat <- map_dfr(paste0(dpath, fnames), read_organise)
Vstat <- Vstat %>% 
  mutate(method = case_when(
    method == "netL0L2" ~ "L0L2+CA",
    method == "net0L0L2" ~ "L0L2",
    method == "netGENIE3" ~ "GENIE3+CA",
    method == "net0GENIE3" ~ "GENIE3",
    method == "netMI" ~ "MI+CA",
    method == "net0MI" ~ "MI",
    method == "netCor" ~ "Spearman+CA",
    method == "net0Cor" ~ "Spearman",
    is.na(method) ~ "Rnd."
  ),
  method = factor(method, levels = c("Rnd.", "L0L2", "GENIE3", "Spearman", "MI", 
                                     "L0L2+CA", "GENIE3+CA", "Spearman+CA", "MI+CA"), 
                  ordered = TRUE),
  ctype = ifelse(grepl("random", ctype), 
                 sapply(strsplit(ctype, "_"), `[`, 1), 
                 ctype),
  ctype = ifelse(ctype == "Epiblast", "EPI", "TE"))

# Plot V-stat distributions -----------------------------------------------

# This is how the V-stat is calculated for one method
Vstat_example <- Vstat %>% 
  filter(ctype == "EPI" & norm_method == "bc" & method == "L0L2") %>% 
  group_by(gene) %>% 
  summarise(V = sum(pval < 0.1)) %>% 
  pull(V) %>% 
  mean()

# Now we compute V-stats for all methods
Vstat_full <- Vstat %>% 
  group_by(ctype, norm_method, method, gene) %>% 
  summarise(V0 = sum(pval < 0.1)) %>% 
  group_by(ctype, norm_method, method) %>% 
  summarise(V = mean(V0))

# Finally, we generate bar plots for each norm-ctype combo
col_cell <- c("EPI" = "#009c00", 
              "TE" = "#287dff")

norm_methods <- c("tpm", "fpkm", "lc", "bc")
ctypes <- c("EPI", "TE")
lims <- c(max(Vstat_full$V[Vstat_full$ctype == "EPI"]),
          max(Vstat_full$V[Vstat_full$ctype == "TE"]))
vstat_plots <- list()
k <- 1
for(i in seq_along(norm_methods)){
  for(j in seq_along(ctypes)){
    vstat_plots[[k]] <- Vstat_full %>% 
      filter(ctype == ctypes[j] & norm_method == norm_methods[i]) %>% 
      ggplot(aes(x = method, y = V)) + 
      geom_col(fill = col_cell[ctypes[j]]) +
      coord_cartesian(ylim = c(0, lims[j])) +
      labs(x = "", y = "V", title = paste0(ctypes[j], "-", norm_methods[i])) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
    k <- k + 1
  }
}
vstat_grid <- wrap_plots(vstat_plots, nrow = 4, ncol = 2)

# Output based on Nature's standard figure sizes
# 89 mm (single column) and 183 mm (double column),full depth of page is 247mm
ggsave(filename = paste0("figs/vstat_barplot.svg"), 
       plot = vstat_grid, width = 150, height = 300, units = "mm")

# Heatmap for main figure
vstat_hm <- list()
for(j in seq_along(ctypes)){
  vstat_hm[[j]] <- Vstat_full %>% 
    filter(ctype == ctypes[j]) %>% 
    mutate(V_norm = V/max(V)) %>% 
    ggplot(aes(method, norm_method, fill = V_norm)) +
    geom_tile() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_fill_viridis(limits = c(0, 1))+
    labs(x = "Inference method", y = "Normalisation", fill = "norm. V",
         title = ctypes[j])
}
vstat_hm <- wrap_plots(vstat_hm, nrow = 2, ncol = 1)

ggsave(filename = paste0("figs/Vstat_heatmap.svg"), 
       plot = vstat_hm, width = 100, height = 130, units = "mm")

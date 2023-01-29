library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(viridis)
library(readxl)
library(patchwork)

dpath <- "~/Dropbox (The Francis Crick)/net_inference/Tom_outputs/"

normali <- c("bc", "lc", "fpkm", "tpm")
ctypes <- c("Epiblast", "PE", "TE")
pred_methods <- c("L0L2", "GENIE3", "Cor", "MI")
ctype_method <- expand.grid(ctypes, pred_methods, stringsAsFactors = FALSE)

sheet_list <- vector("list", length(ctypes)*length(pred_methods))
names(sheet_list) <- paste0(ctype_method$Var1, "_", ctype_method$Var2)

norm_list <- vector("list", length(normali))

for(i in seq_along(normali)){
  for(j in seq_along(sheet_list)){
    sheet_list[[j]] <- read_excel(paste0(dpath, "TFnetCV", normali[i], ".xlsx"),
                                  sheet = names(sheet_list)[j]) %>% 
      mutate(ctype = ctype_method$Var1[j], method = ctype_method$Var2[j]) %>% 
      select(-`Rnd w/ chrom`) %>% 
      pivot_longer(cols = c("Actual", "no chrom"), 
                   names_to = "method_detail", 
                   values_to = "L2loss") %>% 
      mutate(method_detail = case_when(
        method_detail == "Actual" ~ paste0(method, "+CA"),
        method_detail == "no chrom" ~ method
      )) %>% 
      mutate(method_detail = factor(method_detail, 
                                    levels = c("L0L2", "GENIE3", "Cor", "MI",
                                               "L0L2+CA","GENIE3+CA", "Cor+CA", "MI+CA"),
                                    ordered = TRUE)) %>% 
      mutate(normali = normali[i])
  }
  norm_list[[i]] <- bind_rows(sheet_list)
}

all_l2_data <- bind_rows(norm_list)
max_perf <- 1/min(all_l2_data$L2loss)
min_perf <- 1/max(all_l2_data$L2loss)

p_l2 <- vector("list", length(ctypes))
names(p_l2) <- ctypes
for(i in seq_along(ctypes)){
  p_l2[[i]] <- all_l2_data %>% 
    filter(ctype == ctypes[i]) %>% 
    group_by(method_detail, normali) %>% 
    summarise(mu = 1/median(L2loss)) %>% 
    ggplot(aes(method_detail, normali, fill = mu)) +
    geom_tile() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_viridis(limits = c(min_perf, max_perf))+
    labs(x = "Inference method", y = "Normalisation", fill = "Inv. Norm. L2 loss",
         title = ctypes[i])
}
p_l2_grid <- wrap_plots(p_l2, nrow = 3, ncol = 1) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "L2")

ggsave(filename = paste0("../figs/l2_heatmap.svg"), 
       plot = p_l2_grid, width = 100, height = 160, units = "mm")

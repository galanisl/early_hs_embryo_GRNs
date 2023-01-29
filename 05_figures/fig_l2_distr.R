library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(viridis)
library(readxl)

col_pal <- c("Epiblast" = "#009c00", "TE" = "#287dff", "PE" = "#ff284b")

dpath <- "~/Dropbox (The Francis Crick)/net_inference/Tom_outputs/"

normali <- c("bc", "lc", "fpkm", "tpm")
ctypes <- c("Epiblast", "PE", "TE")
pred_methods <- c("L0L2", "GENIE3", "Cor", "MI")
ctype_method <- expand.grid(ctypes, pred_methods, stringsAsFactors = FALSE)

sheet_list <- vector("list", length(ctypes)*length(pred_methods))
names(sheet_list) <- paste0(ctype_method$Var1, "_", ctype_method$Var2)

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
                                    levels = c("L0L2", "L0L2+CA", "GENIE3", "GENIE3+CA", 
                                               "Cor", "Cor+CA","MI", "MI+CA"),
                                    ordered = TRUE))
  }
  tmp <- bind_rows(sheet_list)
  
  p_L2 <- ggplot(tmp, aes(method_detail, L2loss, colour = ctype)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = "width") +
    scale_colour_manual(values = col_pal) +
    coord_cartesian(ylim=c(1, 1.45))+
    labs(x = "", y = "Norm. L2 loss") +
    theme_bw(base_size = 13) +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~ctype)
  
  ggsave(filename = paste0("../figs/l2_violin_plots/", normali[i], "_L2.svg"), 
         plot = p_L2, width = 180, height = 90, units = "mm")
  
}

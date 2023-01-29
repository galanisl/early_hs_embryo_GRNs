
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

path <- "../../result_analysis/Zimeng_results/"

mtd_names <- c("L0L2", "GENIE3", "Spearman", "MI")

plot_tf <- vector("list", length(mtd_names))
for(i in seq_along(mtd_names)){
  load(paste0(path, "TFnetSim", mtd_names[i], "res.rd"))
  
  res_list <- vector("list", length = length(aucValidSymAll))
  for(j in seq_along(aucValidSymAll)){
    res_list[[j]] <- as_tibble(aucValidSymAll[[j]], rownames = "n") %>% 
      pivot_longer(cols = `10`:`100`, names_to = "p", values_to = "aucValid") %>% 
      mutate(n = as.integer(n), p = as.integer(p))
  }
  
  res_tb <- bind_rows(res_list) %>% 
    group_by(n, p) %>% 
    summarise(avgAUC = mean(aucValid), std = sd(aucValid), num = n(), se = std/sqrt(num))
  
  plot_tf[[i]] <- ggplot(res_tb, aes(n, avgAUC, colour = factor(p))) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = avgAUC - se, ymax = avgAUC + se)) +
    scale_color_brewer(palette = "Paired") +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Single cells", y = "AUROC", colour = "TFs", title = mtd_names[i]) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  if(i > 1){
    plot_tf[[i]] <- plot_tf[[i]] +
      theme(axis.title.y = element_blank())
  }
}

pA <- wrap_plots(plot_tf, nrow = 1, ncol = 4) + 
  plot_layout(guides = 'collect')

ggsave(filename = paste0("../figs/artif_grns_updated_cellxTF_100epochs.svg"), 
       plot = pA, width = 280, height = 70, units = "mm")

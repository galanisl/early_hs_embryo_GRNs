library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(viridis)
library(patchwork)

get_perf_table <- function(norm_method){
  # Path to the RDS files with the results
  path <- paste0("../results/results_cv_merged/", norm_method, "/")
  
  # Get files in the path
  fnames <- list.files(path)
  
  # Read the files and form a tibble with all the results
  res <- lapply(paste0(path, fnames), readRDS)
  ctypes <- toupper(sapply(strsplit(fnames, "_"), `[`, 1))
  res <- lapply(seq_along(ctypes), 
                function(i) mutate(res[[i]], ctype = ctypes[i], 
                                   normalisation = norm_method))
  
  tb_res <- bind_rows(res) %>% 
    mutate(method = factor(method, 
                           levels = c("l0l2", "genie3", "spearman", "mi", "rnd",
                                      "l0l2CA", "genie3CA", "spearmanCA", "miCA", "rndCA"), 
                           ordered = TRUE),
           normalisation = factor(normalisation, levels = c("bcr", "lgc", "fpkm", "tpm"), 
                                  ordered = TRUE))
  return(tb_res)
}

get_diff_perf <- function(tb){
  tb <- tb %>% 
    group_by(method, normalisation) %>% 
    summarise(mu = mean(prc))
  # Extract the performances associated with the rnd and rndCA predictors
  rnd <- tb %>% 
    filter(method == "rnd")
  rndCA <- tb %>% 
    filter(method == "rndCA")
  tb <- tb %>% 
    filter(!(method %in% c("rnd", "rndCA")))
  for(i in seq_along(rnd$normalisation)){
    idx_rndCA <- grepl("CA", tb$method) & tb$normalisation == rndCA$normalisation[i]
    idx_rnd <- !grepl("CA", tb$method) & (tb$normalisation == rndCA$normalisation[i])
    tb$mu[idx_rndCA] <- tb$mu[idx_rndCA] - rndCA$mu[i]
    tb$mu[idx_rnd] <- tb$mu[idx_rnd] - rnd$mu[i]
  }
  return(tb)
}

norm_methods <- c("bcr", "lgc", "tpm", "fpkm")
perf_tb <- map_df(norm_methods, get_perf_table) %>% 
  nest(perf = c(roc, prc, method, normalisation))

# Difference heatmaps ----------------------------------------------------

perf_diff <- perf_tb %>% 
  mutate(diffs = purrr::map(perf, get_diff_perf ))

max_perf <- perf_diff$diffs %>% 
  map_dbl(function(x) max(x$mu)) %>% 
  max()

p_diff <- vector("list", nrow(perf_diff))
names(p_diff) <- perf_diff$ctype
for(i in seq_along(perf_diff$ctype)){
  p_diff[[i]] <- ggplot(perf_diff$diffs[[i]], aes(method, normalisation, fill = mu)) +
    geom_tile() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_viridis(limits = c(0, max_perf))+
    labs(x = "Inference method", y = "Normalisation", fill = "Difference wrt rnd.",
         title = perf_diff$ctype[i])
}
p_diff_grid <- wrap_plots(p_diff, nrow = 3, ncol = 1) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Cross-validation")

ggsave(filename = paste0("../figs/cv_heatmap.svg"), 
       plot = p_diff_grid, width = 100, height = 160, units = "mm")






# get_avg_perf <- function(tb){
#   tb <- tb %>% 
#     group_by(method, normalisation) %>% 
#     summarise(mu = mean(prc))
#   return(tb)
# }
# 
# get_fold_change_perf <- function(tb){
#   tb <- tb %>% 
#     group_by(method, normalisation) %>% 
#     summarise(mu = mean(prc))
#   # Extract the performances associated with the rnd and rndCA predictors
#   rnd <- tb %>% 
#     filter(method == "rnd")
#   rndCA <- tb %>% 
#     filter(method == "rndCA")
#   tb <- tb %>% 
#     filter(!(method %in% c("rnd", "rndCA")))
#   for(i in seq_along(rnd$normalisation)){
#     idx_rndCA <- grepl("CA", tb$method) & tb$normalisation == rndCA$normalisation[i]
#     idx_rnd <- !grepl("CA", tb$method) & (tb$normalisation == rndCA$normalisation[i])
#     tb$mu[idx_rndCA] <- tb$mu[idx_rndCA]/rndCA$mu[i]
#     tb$mu[idx_rnd] <- tb$mu[idx_rnd]/rnd$mu[i]
#   }
#   return(tb)
# }
# 
# # Avg. AUPR heatmaps ------------------------------------------------------
# 
# perf_aupr <- perf_tb %>% 
#   mutate(avg_perf = purrr::map(perf, get_avg_perf ))
# 
# p_aupr <- vector("list", nrow(perf_aupr))
# names(p_aupr) <- perf_aupr$ctype
# for(i in seq_along(perf_aupr$ctype)){
#   p_aupr[[i]] <- ggplot(perf_aupr$perf[[i]], aes(method, normalisation, fill = prc)) +
#     geom_tile() +
#     scale_fill_viridis()+
#     labs(x = "Inference method", y = "Normalisation", fill = "Avg. AUPR",
#          title = perf_aupr$ctype[i])
# }
# 
# 
# 
# # Fold change heatmaps ----------------------------------------------------
# 
# perf_fc <- perf_tb %>% 
#   mutate(fold_changes = purrr::map(perf, get_fold_change_perf ))
# 
# p_fc <- vector("list", nrow(perf_fc))
# names(p_fc) <- perf_fc$ctype
# for(i in seq_along(perf_fc$ctype)){
#   p_fc[[i]] <- ggplot(perf_fc$fold_changes[[i]], aes(method, normalisation, fill = mu)) +
#     geom_tile() +
#     scale_fill_viridis()+
#     labs(x = "Inference method", y = "Normalisation", fill = "Fold change wrt rnd.",
#          title = perf_fc$ctype[i])
# }

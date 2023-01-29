library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(viridis)
library(patchwork)


# Function definitions ----------------------------------------------------

get_tibble_for_R_violin <- function(tbs_list, norm_method, rnd_list){
  tb4viol <- bind_rows(tbs) %>% 
    filter(normali == norm_method) %>% 
    bind_rows(tibble(SCORE = 0, 
                     R = rnd_list$Epiblast$samp0R, 
                     method = "Rnd.", 
                     normali = norm_method, 
                     ctype = "EPI"),
              tibble(SCORE = 0, 
                     R = rnd_list$Epiblast$sampR, 
                     method = "Rnd.+CA", 
                     normali = norm_method, 
                     ctype = "EPI"),
              tibble(SCORE = 0, 
                     R = rnd_list$PE$samp0R, 
                     method = "Rnd.", 
                     normali = norm_method, 
                     ctype = "PE"),
              tibble(SCORE = 0, 
                     R = rnd_list$PE$sampR, 
                     method = "Rnd.+CA", 
                     normali = norm_method, 
                     ctype = "PE"),
              tibble(SCORE = 0, 
                     R = rnd_list$TE$samp0R, 
                     method = "Rnd.", 
                     normali = norm_method, 
                     ctype = "TE"),
              tibble(SCORE = 0, 
                     R = rnd_list$TE$sampR, 
                     method = "Rnd.+CA", 
                     normali = norm_method, 
                     ctype = "TE")) %>% 
    mutate(method = factor(method, levels = c("Rnd.", "L0L2", "GENIE3", "Spearman", 
                                              "MI", "Rnd.+CA", "L0L2+CA", "GENIE3+CA", 
                                              "Spearman+CA", "MI+CA"), ordered = TRUE),
           ctype = factor(ctype, levels = c("EPI", "PE", "TE"), ordered = TRUE))
  return(tb4viol)
}

get_diff_perf <- function(tb){
  tb <- tb %>% 
    group_by(method, normali) %>% 
    summarise(mu = median(R))
  # Extract the performances associated with the rnd and rndCA predictors
  rnd <- tb %>% 
    filter(method == "Rnd.")
  rndCA <- tb %>% 
    filter(method == "Rnd.+CA")
  tb <- tb %>% 
    filter(!(method %in% c("Rnd.", "Rnd.+CA")))
  for(i in seq_along(rnd$normali)){
    idx_rndCA <- grepl("CA", tb$method) & tb$normali == rndCA$normali[i]
    idx_rnd <- !grepl("CA", tb$method) & (tb$normali == rndCA$normali[i])
    tb$mu[idx_rndCA] <- tb$mu[idx_rndCA] - rndCA$mu[i]
    tb$mu[idx_rnd] <- tb$mu[idx_rnd] - rnd$mu[i]
  }
  return(tb)
}



# Load data into memory ---------------------------------------------------

# Path to the CVS files with the results
dpath <- "~/Dropbox (The Francis Crick)/net_inference/Tom_outputs/"
fnames <- list.files(dpath, pattern = ".csv")
fnames <- fnames[!grepl("Imp", fnames)]

fname_simple <- str_replace(fnames, "TFnet", "")
fname_simple <- str_replace(fname_simple, "[R|r]es", "")
fname_simple <- str_replace(fname_simple, "[.]csv", "")

# Reading in the tables
tbs <- list()
for(i in seq_along(fnames)){
  tbs[[i]] <- readr::read_csv(paste0(dpath, fnames[i])) %>% 
    select(-FROM, -TO) %>% 
    mutate(method = toupper(sapply(strsplit(fname_simple[i], "_"), `[`, 1)),
           normali = sapply(strsplit(fname_simple[i], "_"), `[`, 2),
           ctype = sapply(strsplit(fname_simple[i], "_"), `[`, 3)) %>% 
    mutate(
      method = case_when(
        method == "0COR" ~ "Spearman",
        method == "0L0L2" ~ "L0L2",
        method == "0GENIE3" ~ "GENIE3",
        method == "0MI" ~ "MI",
        method == "COR" ~ "Spearman+CA",
        method == "L0L2" ~ "L0L2+CA",
        method == "GENIE3" ~ "GENIE3+CA",
        method == "MI" ~ "MI+CA"),
      ctype = case_when(
        ctype == "Epiblast" ~ "EPI",
        TRUE ~ ctype
      ))
}

# Random R values
rnd_list <- vector("list", 4)
load(paste0(dpath, "TFnetRrnd_bc.rd"))
rnd_list[[1]] <- allSampR
load(paste0(dpath, "TFnetRrnd_lc.rd"))
rnd_list[[2]]  <- allSampR
load(paste0(dpath, "TFnetRrnd_fpkm.rd"))
rnd_list[[3]]  <- allSampR
load(paste0(dpath, "TFnetRrnd_tpm.rd"))
rnd_list[[4]]  <- allSampR


norm_method <- c("bc", "lc", "fpkm", "tpm")
tmp_list <- list()

for(i in seq_along(norm_method)){
  tmp_list[[i]] <- get_tibble_for_R_violin(tbs, norm_method[i], rnd_list[[i]])
}

all_R_data <- bind_rows(tmp_list) %>% 
  mutate(normali = factor(normali, levels = c("bc", "lc", "fpkm", "tpm"), 
                          ordered = TRUE))

all_R_data <- all_R_data %>% 
  nest(perf = c(SCORE, R, method, normali))

# Difference heatmaps ----------------------------------------------------

r_diff <- all_R_data %>% 
  mutate(diffs = purrr::map(perf, get_diff_perf))

max_perf <- r_diff$diffs %>% 
  map_dbl(function(x) max(x$mu)) %>% 
  max()

p_diff <- vector("list", nrow(r_diff))
names(p_diff) <- r_diff$ctype
for(i in seq_along(r_diff$ctype)){
  p_diff[[i]] <- ggplot(r_diff$diffs[[i]], aes(method, normali, fill = mu)) +
    geom_tile() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_viridis(limits = c(0, max_perf))+
    labs(x = "Inference method", y = "Normalisation", fill = "Difference wrt rnd.",
         title = r_diff$ctype[i])
}
p_diff_grid <- wrap_plots(p_diff, nrow = 3, ncol = 1) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "R")

ggsave(filename = paste0("../figs/r_heatmap.svg"), 
       plot = p_diff_grid, width = 100, height = 160, units = "mm")

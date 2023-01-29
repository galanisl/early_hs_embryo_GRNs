
library(readr)
library(ggplot2)
library(patchwork)
library(purrr)
library(dplyr)
library(tidyr)
library(igraph)
library(stringr)
library(viridis)

get_common_spec_nets <- function(tb){
  ctype_edges <- map(tb$fname, read_csv)
  ctype_nets <- map(ctype_edges, graph_from_data_frame, directed = TRUE)
  names(ctype_nets) <- tb$ctype
  
  net_inter <- igraph::intersection(ctype_nets$Epiblast, 
                                    ctype_nets$PE, 
                                    ctype_nets$TE, 
                                    keep.all.vertices = FALSE)
  
  epi_specific <- difference(ctype_nets$Epiblast, ctype_nets$PE) %>% 
    difference(ctype_nets$TE) %>% 
    delete_vertices(v = degree(.) == 0)
  
  pe_specific <- difference(ctype_nets$PE, ctype_nets$Epiblast) %>% 
    difference(ctype_nets$TE) %>% 
    delete_vertices(v = degree(.) == 0)
  
  te_specific <- difference(ctype_nets$TE, ctype_nets$Epiblast) %>% 
    difference(ctype_nets$PE) %>% 
    delete_vertices(v = degree(.) == 0)
  return(list(common = net_inter, 
              epi = epi_specific, pe = pe_specific, te = te_specific))
}

get_mrk_activity <- function(net, cell_type, mrks, method, normalisation){
  net_edg <- igraph::as_data_frame(net) %>% 
    group_by(from) %>% 
    summarise(out_deg = n()) %>%
    # mutate(norm_outdeg = ifelse(max(out_deg) == 0, 0, out_deg/max(out_deg))) %>% 
    mutate(norm_outdeg = out_deg/max(degree(net, mode = "out"))) %>% 
    filter(from %in% mrks) %>% 
    select(from, out_deg, norm_outdeg) %>% 
    mutate(ctype = cell_type, method = method, normalisation = normalisation)
  missing_mrks <- setdiff(mrks, net_edg$from)
  if(length(missing_mrks) > 0){
    net_edg <- bind_rows(net_edg, 
                         tibble(from = missing_mrks, out_deg = 0, norm_outdeg = 0,
                                ctype = cell_type, method = method, 
                                normalisation = normalisation))
  }
  
  return(net_edg)
}

dpath <- "~/Dropbox (The Francis Crick)/net_inference/Tom_outputs/"

tb_files <- list.files(dpath, pattern = "*.csv")

tb_details <- sapply(strsplit(tb_files, split = "[.]"), `[`, 1)
tb_details <- strsplit(tb_details, split = "_")

tb_details <- map_dfr(tb_details, 
                      function(x) tibble(method = x[1],
                                         normal = x[2],
                                         ctype = x[3])) %>% 
  mutate(fname = paste0(dpath, tb_files),
         method = str_replace(method, "TFnet", ""),
         method = str_replace(method, "res", ""),
         method = str_replace(method, "Res", "")) %>% 
  nest(data = c(ctype, fname)) %>% 
  filter(!grepl("Imp", method))

mrk_of_interest <- c("NANOG", "GATA3", "SOX17")

activity_list <- vector("list", nrow(tb_details))

for(i in 1:nrow(tb_details)){
  tmp <- get_common_spec_nets(tb_details$data[[i]])
  tmp_activities <- map2(tmp, names(tmp), get_mrk_activity, 
                         mrks = mrk_of_interest, method = tb_details$method[i],
                         normalisation = tb_details$normal[i])
  activity_list[[i]] <- bind_rows(tmp_activities) 
    # mutate(norm_odeg = norm_outdeg/max(norm_outdeg))
    # mutate(sqrt_odeg = sqrt(out_deg)) %>% 
    # mutate(norm_odeg = sqrt_odeg/max(sqrt_odeg))
}

activity_results <- bind_rows(activity_list) %>% 
  mutate(method = case_when(
    method == "0cor" ~ "Spearman",
    method == "0GENIE3" ~ "GENIE3",
    method == "0L0L2" ~ "L0L2",
    method == "0MI" ~ "MI",
    method == "Cor" ~ "Spearman+CA",
    method == "GENIE3" ~ "GENIE3+CA",
    method == "L0L2" ~ "L0L2+CA",
    method == "MI" ~ "MI+CA"
  ),
  ctype = str_to_upper(ctype),
  normalisation = case_when(
    normalisation == "bc" ~ "Batch-corrected",
    normalisation == "fpkm" ~ "FPKM",
    normalisation == "lc" ~ "log-counts",
    normalisation == "tpm" ~ "TPM")
  ) %>% 
  mutate(method = factor(method, 
                         levels = c("L0L2", "GENIE3", "Spearman", "MI",
                                    "L0L2+CA", "GENIE3+CA", "Spearman+CA", "MI+CA"),
                         ordered = TRUE),
         normalisation = factor(normalisation, 
                                levels = c("log-counts", "Batch-corrected", 
                                           "FPKM", "TPM"),
                                ordered = TRUE),
         actv = norm_outdeg >= (quantile(norm_outdeg)["50%"]-0.005))


ctypes <- unique(activity_results$ctype)
norms <- levels(activity_results$normalisation)

# Build heatmap for each Cell type-Normalisation group

p_activ <- vector("list", length(ctypes)*length(norms))
k <- 1
for(i in seq_along(ctypes)){
  for(j in seq_along(norms)){
    p_activ[[k]] <- activity_results %>% 
      filter(ctype == ctypes[i] & normalisation == norms[j]) %>% 
      ggplot(aes(method, from, fill = as.numeric(actv))) +
      geom_tile() +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      scale_fill_viridis(limits = c(0, 1)) +
      labs(x = "Inference method", y = "", fill = "Active or Inactive",
           title = paste0(ctypes[i], " ", norms[j]))
    k <- k + 1
  }
}

p_activ_grid <- wrap_plots(p_activ, nrow = 4, ncol = 4) + 
  plot_layout(guides = "collect") 

ggsave(filename = paste0("figs/mrk_heatmap_YESNO_v2.svg"), 
       plot = p_activ_grid, width = 220, height = 190, units = "mm")

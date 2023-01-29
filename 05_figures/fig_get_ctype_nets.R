library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(igraph)
library(stringr)

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
  filter(!grepl("Imp", method)) %>% 
  mutate(method = case_when(
    method == "0cor" ~ "Spearman",
    method == "0GENIE3" ~ "GENIE3",
    method == "0L0L2" ~ "L0L2",
    method == "0MI" ~ "MI",
    method == "Cor" ~ "Spearman+CA",
    method == "GENIE3" ~ "GENIE3+CA",
    method == "L0L2" ~ "L0L2+CA",
    method == "MI" ~ "MI+CA"
  ))

mrk_of_interest <- c("NANOG", "SOX17", "GATA3")

# Focus on a Method+Normalisation combination
mtd <- "MI+CA"
nrm <- "tpm"
tb_flt <- tb_details %>% 
  filter(method == mtd & normal == nrm)

# Get the regulatory networks
nets <- get_common_spec_nets(tb_flt$data[[1]])
nets <- nets[-1]

# Write the networks in table format
for(i in seq_along(nets)){
  mrk_neighs <- ego(nets[[i]], order = 1, nodes = mrk_of_interest[i], mode = "all")[[1]]
  mrk_net <- induced_subgraph(nets[[i]], vids = names(mrk_neighs))
  
  net_tb <- igraph::as_data_frame(mrk_net, what = "edges") 
  write_csv(net_tb, paste0("../figs/ctype_nets/", names(nets)[i], "_", 
                           mrk_of_interest[i], "_",
                           mtd, "_", nrm, ".csv"))
}

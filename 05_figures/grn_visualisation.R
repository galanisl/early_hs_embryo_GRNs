library(readr)
library(dplyr)
library(igraph)
library(visNetwork)

ctype <- "Epiblast"
dpath <- paste0("path_to_GRN_inferences/TFnetMIres_tpm_", ctype, ".csv")
tf_colour <- c("Epiblast" = "#009c00", "TE" = "#287dff", "PE" = "#ff284b")
gn_colour <- c("Epiblast" = "#a1d99b", "TE" = "#9ecae1", "PE" = "#fc9272")

mrk <- c("Epiblast" = "NANOG", "TE" = "GATA3", "PE" = "SOX17")

tfs <- read_tsv("../02_data_processing/final_data/Homo_sapiens_TF.txt")

# GRN around gene ---------------------------------------------------------

grn_tb <- read_csv(dpath) %>% 
  filter(FROM %in% mrk[ctype] | TO %in% mrk[ctype]) %>% 
  arrange(desc(SCORE))

# Take top 25 interactions
grn_tb <- grn_tb[1:25,]

gene_grn <- graph_from_data_frame(grn_tb, directed = TRUE)

V(gene_grn)$isTF <- V(gene_grn)$name %in% tfs$Symbol
V(gene_grn)$shape <- ifelse(V(gene_grn)$isTF, "diamond", "circle")
V(gene_grn)$color <- ifelse(V(gene_grn)$isTF, tf_colour[ctype], gn_colour[ctype])

E(gene_grn)$value <- E(gene_grn)$SCORE

visIgraph(gene_grn)

# Top GRN -----------------------------------------------------------------

grn_tb <- read_csv(dpath) %>% 
  arrange(desc(SCORE))

# Take top 25 interactions
grn_tb <- grn_tb[1:25,]

top_grn <- graph_from_data_frame(grn_tb, directed = TRUE)

V(top_grn)$isTF <- V(top_grn)$name %in% tfs$Symbol
V(top_grn)$shape <- ifelse(V(top_grn)$isTF, "diamond", "circle")
V(top_grn)$color <- ifelse(V(top_grn)$isTF, tf_colour[ctype], gn_colour[ctype])

E(top_grn)$value <- E(top_grn)$SCORE

visIgraph(top_grn)


# TF-TF -------------------------------------------------------------------

grn_tb <- read_csv(dpath) %>% 
  filter(FROM %in% mrk[ctype] | TO %in% mrk[ctype]) %>% 
  arrange(desc(SCORE))

g <- graph_from_data_frame(grn_tb, directed = TRUE)

tf_tf <- delete_vertices(g, v = !(V(g)$name %in% tfs$Symbol))
# Uncomment to focus on top TF-TF interactions
# tf_tf <- delete_edges(tf_tf, edges = which(E(tf_tf)$SCORE < sort(E(tf_tf)$SCORE, decreasing = TRUE)[25]))
# tf_tf <- delete_vertices(tf_tf, v = degree(tf_tf) == 0)

V(tf_tf)$isTF <- V(tf_tf)$name %in% tfs$Symbol
V(tf_tf)$shape <- ifelse(V(tf_tf)$isTF, "diamond", "circle")
V(tf_tf)$color <- ifelse(V(tf_tf)$isTF, tf_colour[ctype], gn_colour[ctype])

E(tf_tf)$value <- E(tf_tf)$SCORE

visIgraph(tf_tf)


# Save GRNs as tab-separated files ----------------------------------------

# GRN around genes
tb <- igraph::as_data_frame(gene_grn)
write_tsv(tb, paste0("grn_vis/", mrk[ctype], "_", ctype, "_TPM_MI.tsv"))

# Top interactions
tb <- igraph::as_data_frame(top_grn)
write_tsv(tb, paste0("grn_vis/TOP_", ctype, "_TPM_MI.tsv"))

# TF-TF
tb <- igraph::as_data_frame(tf_tf)
write_tsv(tb, paste0("grn_vis/TF-TF_", mrk[ctype], "_", ctype, "_TPM_MI.tsv"))

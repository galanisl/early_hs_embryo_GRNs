library(scater)
library(patchwork)

load("data/sce_integrated_lateBlast_clean.RData")

col_cell <- c("Epiblast" = "#009c00", 
              "Trophectoderm" = "#287dff", 
              "Primitive endoderm" = "#ff284b")

mrk_genes <- c("Ref", "NANOG", "GATA3", "SOX17")

p <- vector("list", length(mrk_genes))
names(p) <- mrk_genes

for(i in seq_along(mrk_genes)){
  if(mrk_genes[i] == "Ref"){
    p[[i]] <- plotReducedDim(sce, dimred = "UMAP", 
                             colour_by = "cell_type", add_legend = FALSE) +
      scale_fill_manual(values = col_cell)
  }else{
    p[[i]] <- plotReducedDim(sce, dimred = "UMAP", 
                             colour_by = mrk_genes[i]) 
  }
}

p_grid <- wrap_plots(p, ncol = 2, nrow = 2)

# Output based on Nature's standard figure sizes
# 89 mm (single column) and 183 mm (double column),full depth of page is 247mm
ggsave(filename = paste0("figs/marker_umap.pdf"), 
       plot = p_grid, width = 183, height = 160, units = "mm")

save(p, file = "analysis_outputs/mrk_umap_plots.rds")

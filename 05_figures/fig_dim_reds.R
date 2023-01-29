
library(scater)
library(patchwork)

load("../data/sce_integrated_lateBlast_clean.RData")

col_cell <- c("Epiblast" = "#009c00", 
             "Trophectoderm" = "#287dff", 
             "Primitive endoderm" = "#ff284b")

col_batch <- c("Blk" = "#1b9e77", 
             "Pet" = "#d95f02", 
             "Yan" = "#7570b3")

norm_methods <- c("batch_corrected", "logcounts", "tpm", "fpkm")
p <- vector("list", length(norm_methods))
names(p) <- norm_methods

for(i in seq_along(norm_methods)){
  if(norm_methods[i] == "batch_corrected"){
    panelA <- plotReducedDim(sce, dimred = "UMAP", 
                             colour_by = "cell_type", add_legend = FALSE) +
      scale_fill_manual(values = col_cell)
    panelB <- plotReducedDim(sce, dimred = "UMAP", 
                             colour_by = "batch", add_legend = FALSE) +
      scale_fill_manual(values = col_batch)
    p[[i]] <- panelA | panelB
  }else{
    tmp <- runUMAP(sce, exprs_values = norm_methods[i], name = "UMAP")
    panelA <- plotReducedDim(tmp, dimred = "UMAP", 
                             colour_by = "cell_type", add_legend = FALSE) +
      scale_fill_manual(values = col_cell)
    panelB <- plotReducedDim(tmp, dimred = "UMAP", 
                             colour_by = "batch", add_legend = FALSE) +
      scale_fill_manual(values = col_batch)
    p[[i]] <- panelA | panelB
  }
  # Output based on Nature's standard figure sizes
  # 89 mm (single column) and 183 mm (double column),full depth of page is 247mm
  ggsave(filename = paste0("figs/", norm_methods[i], "_umap.pdf"), 
         plot = p[[i]], width = 183, height = 80, units = "mm")
}

save(p, file = "analysis_outputs/umap_plots.rds")

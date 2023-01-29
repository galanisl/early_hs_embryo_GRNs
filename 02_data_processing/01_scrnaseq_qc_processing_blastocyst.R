library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scran)
library(batchelor)

load("final_data/sce_integrated_lateBlast.RData")

# Bioconductor scRNA-seq processing pipeline ------------------------------
sce <- sce[!is.na(rowData(sce)$chr), ]

# Cell QC metrics
is_mito <- rowData(sce)$chr == "MT"
sce <- addPerCellQC(sce, subsets = list(Mito = is_mito))

qc <- tibble(tot_counts = sce$total / 1e6,
             expr_genes = sce$detected,
             pct_mt = sce$subsets_Mito_percent)

lbl <- c("Library size (millions)", "Number of expressed genes", 
         "Counts mapped to mito-genes (%)")
p <- list()
for(i in 1:ncol(qc)){
  p[[i]] <- ggplot(qc, aes_string(x = colnames(qc)[i])) + 
    geom_histogram(bins = 15) +
    geom_vline(xintercept = median(qc[[i]]), linetype = 2, colour = "blue") +
    labs(x = lbl[i], y = "Number of cells") + theme_bw()
}
wrap_plots(p)

# Adaptive thresholds on cells
libsize_drop <- isOutlier(sce$total, nmads = 3, type = "lower", log = TRUE)
feature_drop <- isOutlier(sce$detected, nmads = 3, type = "lower", log = TRUE)
mito_drop <- isOutlier(sce$subsets_Mito_percent , nmads = 3, type = "higher")

sce <- sce[, !(libsize_drop | feature_drop | mito_drop)]

# Remove PE sample that clusters with TE
sce <- sce[, colnames(sce) != "ERR1042037"]

tibble(Stat = c("Removed due to Lib. Size", "Removed due to Features", 
                "Removed due to Mito.", "Remaining"), 
       `Number of cells` = c(sum(libsize_drop), sum(feature_drop), 
                             sum(mito_drop), ncol(sce)))

# Get rid of no-show genes
sce <- sce[rowSums(counts(sce)) != 0, ]

# Only keep genes with avg. expression across cells >= 1
ave.counts <- rowMeans(counts(sce))
keep <- ave.counts >= 1
sum(keep)

sce <- sce[keep, ]

# Size-factor normalisation
clust_sce <- quickCluster(sce) 
sce <- computeSumFactors(sce, cluster = clust_sce)
sce <- logNormCounts(sce)

# HVG selection.
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop = 0.1)

rowData(sce)$is_hvg <- rownames(sce) %in% hvg

# Dimensionality reduction.
set.seed(1234)
sce <- runPCA(sce, ncomponents = 25, subset_row = hvg)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors = TRUE)

plotPCA(sce, colour_by = "batch") + plotPCA(sce, colour_by = "cell_type")
plotUMAP(sce, colour_by = "batch") + plotUMAP(sce, colour_by = "cell_type")

# Remove batch effect for better visualisation and clustering
sce_resc <- rescaleBatches(sce, batch = factor(sce$batch))

assay(sce, "batch_corrected") <- assays(sce_resc)$corrected

# Redo dimensionality reduction using batch-corrected gene expression
set.seed(1234)
sce <- runPCA(sce, ncomponents = 25, 
              subset_row = hvg, 
              exprs_values="batch_corrected")
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors = TRUE)

plotPCA(sce, colour_by = "batch") + plotPCA(sce, colour_by = "cell_type")
plotUMAP(sce, colour_by = "batch") + plotUMAP(sce, colour_by = "cell_type")

# Clustering
g <- buildSNNGraph(sce, use.dimred = 'PCA', k = 10)
sce$louvain_clusters <- factor(igraph::cluster_louvain(g)$membership)

plotUMAP(sce, colour_by = "louvain_clusters")

# Look at some markers
plotReducedDim(sce, "UMAP", colour_by = "GATA3", by_exprs_values = "batch_corrected")
plotReducedDim(sce, "UMAP", colour_by = "NANOG", by_exprs_values = "batch_corrected")
plotReducedDim(sce, "UMAP", colour_by = "SOX17", by_exprs_values = "batch_corrected")

# griph exploration -------------------------------------------------------

library(griph)

res_griph <- griph_cluster(counts(sce), 
                           BatchAssignment = factor(sce$batch),
                           ClassAssignment = factor(sce$cell_type),
                           use.par = FALSE, plot = FALSE)

g.true <- plotGraph(res_griph, fill.type = "true", mark.type = "predicted")

emb <- tibble(x = g.true$y[,1], y = g.true$y[,2], 
              cluster = factor(res_griph$MEMB), 
              btch = factor(sce$batch),
              cell_type = factor(sce$cell_type))

# Add griph results to the SCE object
reducedDim(sce, type = "griph") <- as.data.frame(select(emb, x, y))
sce$griph_clusters <- emb$cluster

# Make a summary plot
umap_batch <- plotUMAP(sce, colour_by = "batch", point_size = 2) + 
  theme_bw(base_size = 12)
umap_cell <- plotUMAP(sce, colour_by = "cell_type", point_size = 2) + 
  theme_bw(base_size = 12)
griph_cell <- plotReducedDim(sce, "griph", colour_by = "cell_type", 
                             point_size = 2) + theme_bw(base_size = 12)
umap_epi <- plotReducedDim(sce, "UMAP", colour_by = "NANOG", point_size = 2, 
                     by_exprs_values = "batch_corrected") + 
  theme_bw(base_size = 12)
umap_pe <- plotReducedDim(sce, "UMAP", colour_by = "SOX17", point_size = 2, 
                          by_exprs_values = "batch_corrected") + 
  theme_bw(base_size = 12)
umap_te <- plotReducedDim(sce, "UMAP", colour_by = "GATA3", point_size = 2, 
                          by_exprs_values = "batch_corrected") + 
  theme_bw(base_size = 12)

(umap_batch | umap_cell | griph_cell)/(umap_epi | umap_pe | umap_te) +
  plot_layout(guides = 'collect')


save(res_griph, emb, file = "final_data/griph_lateBlast_clean.RData")
save(sce, file = "final_data/sce_integrated_lateBlast_clean.RData")


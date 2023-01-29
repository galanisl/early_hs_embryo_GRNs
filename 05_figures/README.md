# Scripts to construct the figures in the manuscript

The scripts provided in this folder require large data outputs to generate the figure in the manuscript. Some of these files are provided in the `../02_data_processing/final_data` folder and the rest can be downloaded from `URL`.

- `grn_visualisation.R`: Script for interactive visualisation of the predicted GRNs.
- `fig_simulations.R`: GRN inference on simulated scRNA-seq data (Fig. 2).
- `fig_R_heatmap.R`: Summary of reproducibility score (R) GRN evaluation as a heatmap (Fig. 3A).
- `fig_cv_heatmap.R`: Summary of cross-validation of inferred GRNs as a heatmap (Fig. 3B).
- `fig_l2_heatmap.R`: Summary of L2 loss evaluation of inferred GRNs as a heatmap (Fig. 3C).
- `fig_l2_distr.R`: Exploration of the L2 loss distributions per method and inferred GRN (linked to Fig. 3C).
- `fig_dim_reds.R`: Representation of the scRNA-seq data in UMAP space (Fig. 4A).
- `fig_markers.R`: Representation of gene expression levels for selected markers in UMAP space (Fig. 4A).
- `fig_mrk_heatmap.R`: Summary of marker activity analysis as a heatmap (Fig. 4B).
- `fig_Vstat.R`: Summary of V-statistic analysis as a barplot and heatmap (Fig. 4C).
- `grn_visualisation_v4.0_nb_commented.Rmd`: R Markdown notebook for GRN visualisation with directionality (Fig. 5 and 6A)






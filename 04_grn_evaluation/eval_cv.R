
# Evaluate all methods for cell type X

source("l0l2.R")
source("genie3.R")
source("mutual_info.R")
source("spearman.R")
source("rnd_predictor.R")

source("kfold_cv_evaluation.R")

library(rsample)
library(purrr)

# Reading necessary data --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

ctype <- c("epi", "te", "pe")
norm_method <- c("tpm", "bcr", "lgc", "fpkm")
inf_method <- c("genie3", "genie3CA", "l0l2", "l0l2CA", 
                "mi", "miCA", "spearman", "spearmanCA", "rnd", "rndCA")

analysis_opts <- expand.grid(ctype, norm_method, inf_method, 
                             stringsAsFactors = FALSE) %>% 
  as_tibble() %>% 
  rename(ctype = Var1, norm_method = Var2, inf_method = Var3) %>% 
  arrange(ctype, norm_method, inf_method)

sel_ctype <- analysis_opts$ctype[as.integer(args[1])]
sel_norm <- analysis_opts$norm_method[as.integer(args[1])]
sel_inf <- analysis_opts$inf_method[as.integer(args[1])]

exprMatr <- readRDS(paste0("data/", sel_ctype, "_", sel_norm, ".rds"))
exprMatr <- exprMatr[rowSums(exprMatr) > 0, ]

if(sel_ctype == "te"){
  atac <- readRDS("data/te_tf_binding.rds")
}else{
  atac <- readRDS("data/icm_tf_binding.rds")  
}
atac <- atac %>% 
  select(tf, nearest_prom) %>% 
  rename(regulatoryGene = tf, targetGene = nearest_prom)

# Only keep genes from the expression matrix for which there is chromatin 
# accessibility information
exprMatr <- exprMatr[intersect(rownames(exprMatr), 
                               union(atac$regulatoryGene, atac$targetGene)), ]

regulators <- readRDS("data/regulators.rds")
regulators <- regulators[regulators %in% rownames(exprMatr)]

if(grepl("CA", sel_inf)){
  fn_name <- paste0("refined_", stringr::str_sub(sel_inf, end = -3))
}else{
  fn_name <- paste0("refined_", sel_inf)
  atac <- NULL
}

# Setting up folds and gold standard --------------------------------------

# Interactions to be considered gold standard in test set
top_interactions <- round((length(regulators) * nrow(exprMatr) -
                             length(regulators)) * 0.01)

folds <- vfold_cv(tibble(sam = colnames(exprMatr)), v = 2, repeats = 10)

# Evaluation of inference method ------------------------------------------

print(paste0("Message from job ", args[1], ": "))
print(paste0("   Working on ", sel_ctype, "_", sel_norm, "_", sel_inf))

train_set <- map(folds$splits, apply_method_to_fold, fld_type = "analysis",
                 netinf_method = fn_name,
                 exprMtx = exprMatr, tfs = regulators, tf_binding = atac)

test_set <- map(folds$splits, apply_method_to_fold, fld_type = "assessment",
                netinf_method = fn_name,
                exprMtx = exprMatr, tfs = regulators,
                topInts = top_interactions, tf_binding = atac)

ev_result <- map2(train_set, test_set, eval_method) %>%
  bind_rows() %>%
  mutate(method = sel_inf)

print(paste0("   Done with ", sel_ctype, "_", sel_norm, "_", sel_inf))

saveRDS(ev_result, file = paste0("results_cv/", sel_ctype, "_", 
                                  sel_norm, "_", sel_inf, ".rds"))

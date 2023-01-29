
library(dplyr)
library(precrec)

# Function definitions ----------------------------------------------------

apply_method_to_fold <- function(fld, fld_type = "analysis", netinf_method,
                                 exprMtx, tfs, tf_binding = NULL, topInts = NA){
  sample_names <- as.data.frame(fld, data = fld_type)$sam
  exprMtx <- exprMtx[, sample_names]
  
  links <- do.call(what = netinf_method, 
                   args = list(exprMatr = exprMtx, regulators = tfs, 
                               tf_binding = tf_binding, 
                               nCores = detectCores()-1))
  
  if(!is.na(topInts) && fld_type == "assessment"){
    links <- links[1:topInts, ]
  }
  
  return(links)
}

eval_method <- function(train_int, test_int){
  lbl_training <- left_join(train_int, test_int, 
                            by = c("regulatoryGene", "targetGene")) %>% 
    rename(lbl = weight.y) %>% 
    mutate(lbl = ifelse(is.na(lbl), FALSE, TRUE))
  
  ev <- evalmod(scores = lbl_training$weight.x, labels = lbl_training$lbl)
  ev <- auc(ev)
  ev_tidy <- tibble(roc = ev$aucs[ev$curvetypes == "ROC"],
                    prc = ev$aucs[ev$curvetypes == "PRC"])
  
  return(ev_tidy)
}

library(GENIE3)
library(dplyr)

#' GENIE3 regulatory network inference
#' 
#' Given an n x m gene expression (n genes and m samples) and a list of 
#' regulators assign likelihood scores to all possible regulator x gene 
#' interactions using GENIE3.
#' 
#' @param exprMatr matrix; A n x m expression matrix (n genes and m samples).
#' Rows must be named with gene symbols or IDs and columns with sample IDs.
#' @param regulators character; A vector of regulator genes (e.g. transcription 
#' factors). If NULL, all genes in the expression matrix are used as regulators.
#' Regulators that are not listed the expression matrix are discarded.
#' @param nCores integer; Number of cores to run the network inference in 
#' parallel.
#' 
#' @return A tibble with regulatory links regulatoryGene-->targetGene sorted
#' by the likelihood score weight.
#' 
#' @author Gregorio Alanis-Lobato \email{gregorio.alanis@crick.ac.uk}
#' 
#' @examples
#' # Generate a synthetic dataset
#' exprMatr <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
#' rownames(exprMatr) <- paste("Gene", 1:20, sep="")
#' colnames(exprMatr) <- paste("Sample", 1:5, sep="")
#' regulators <- c("Gene2", "Gene4", "Gene7")
#' 
#' # Infer a regulatory network with L0L2 regression
#' net <- GENIE3(exprMatr, regulators, 1)
#' 
#' @export
#' @import GENIE3
#' @importFrom dplyr tibble as_tibble %>% arrange desc
#'
genie3_net <- function(exprMatr, regulators = NULL, nCores = 1){
  if(is.null(regulators)){
    regulators <- rownames(exprMatr)
  }else{
    regulators <- regulators[regulators %in% rownames(exprMatr)]
  }
  weightMat <- GENIE3(exprMatr, regulators=regulators, nCores = nCores)
  linkList <- getLinkList(weightMat) %>% 
    as_tibble() %>% 
    arrange(desc(weight))
  return(linkList)
}

# Refine GENIE3 predictions with TF-gene binding data from ChIP-seq or ATAC-seq.
refined_genie3 <- function(exprMatr, regulators = NULL, tf_binding = NULL, 
                         nCores = 1){
  linkList <- genie3_net(exprMatr, regulators, nCores)
  if(!is.null(tf_binding)){
    linkList <- left_join(tf_binding, linkList,
                          by = c("regulatoryGene", "targetGene")) %>% 
      filter(!is.na(weight)) %>% 
      arrange(desc(weight))
  }
  return(linkList)
}
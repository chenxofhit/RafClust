#' perform the RafClust dissimilarity learning and clustering algorithm
#'
#' @title RafClust
#
#' @param data an (n x m) data matrix of gene expression measurements of single cells
#' @param NumC number of clusters to be extracted over data
#' @param gene_filter logical should gene filtering be performed?
#' @param minClusterSize minimize cluster size 
#' @param deepSplit deep split 
#' @param verbose should the help information be output?
#' @return list of 3 elements describing the results:
#'		SF:  Spearman Feature Space,
#'  	D:   RafClust dissimilarity matrix,
#'    lab: inferred cluster labels
#'  	hc:  hierarchical clustering solution
#  
#' @importFrom e1071 kurtosis
#' @importFrom dclone make.symmetric
#' @import stats
#' @import fastcluster
#' @import cutreeDynamic
#' @import dendextend
#' @export 
#'
RafClust<-function(data,NumC=NULL,gene_filter=TRUE,minClusterSize=30,deepSplit=1,verbose=FALSE)
{
  if(verbose){
    message("Calculating feature space of the cells...")
  }
  
  #- feature construction
  FE <- RafClustFE(data, gene_filter=gene_filter, frq=0.06, verbose=verbose)
  
  if(verbose){
    message("Learning similarities  of the cells...")
  }

  #- random forest learning
  rfdis <- RafClustRF(FE$F)

  dis <- rfdis

  #- "normalized" distance
  X            <- rfdis
  X            <- t(t(X) - colMeans(X))
  X            <- X / max(abs(X))
  rfdis2       <- make.symmetric(X)
  rfdis2       <- rfdis2 + abs(min(rfdis2))
    
  #- normalize?
  t1 <- abs(kurtosis(as.vector(rfdis)))
  t2 <- abs(kurtosis(as.vector(rfdis2)))
  if(t1 > 5*t2){
    diag(rfdis2) <- 0
    dis <- rfdis2
  }else{
    dis <- rfdis
  }

  hc = NULL
  lab = NULL
  if(!is.null(NumC)){
    if(verbose){
      message("hierarchical clustering the cells...")
    }
	  hc  <- hclust(as.dist(dis), method="average")
  	lab <- cutree(hc, NumC)

    if(verbose){
      library(dendextend)

      nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
      dend <- as.dendrogram(hc)
      plot(dend, xlab = "", sub="", ylab = "Distance",
         main = "Dendrogram", nodePar = nodePar)
    
      rect.dendrogram(dend , k = length(unique(lab)), border = "red")
    }
  }else{
    library(dynamicTreeCut)
    library(fastcluster)

    if(verbose){
      message("hierarchical clustering the cells with k automatic selection")
    }
    d <-  as.dist(dis)
    hc <- fastcluster::hclust(d,method = "average")
    lab <- cutreeDynamic(dendro = hc, cutHeight = NULL,
                           minClusterSize = minClusterSize,
                           method = "hybrid", deepSplit = deepSplit,
                           pamStage = TRUE,  distM = as.matrix(d), maxPamDist = 0,
                           verbose = 0)
    
    outliers_ids <- which(lab==0)
    hc_labs_clean  = lab[-outliers_ids]
    
    if(verbose){
      library(dendextend)

      nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
      dend <- as.dendrogram(hc)
      plot(dend, xlab = "", sub="", ylab = "Distance",
         main = "Dendrogram", nodePar = nodePar)
      
      rect.dendrogram(dend , k = length(unique(lab)), border = 8, lty = 5, lwd = 2)
    }
  }

  res = list(SF = FE,
	             D = rfdis, #- never normalize for now
	             lab = lab,
	             hc  = hc)

   return(res)
}

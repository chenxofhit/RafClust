#' SIMILARITY LEARNING FOR RafClust
#'
#' @param F matrix Feature matrix cells x genes
#' @param oob logical Should out of bag proximity be used?
#' @return matrix of dissimilarities between cells
#' @importFrom randomForest randomForest
#' @importFrom fpc pamk
RafClustRF <- function(F, oob=TRUE){
  featFu <- function(i){
    if(dim(F)[1] >= 2000){
      cl   = pamk(F[,i],usepam = FALSE)
    }
    else{
      cl   = pamk(F[,i],usepam = TRUE)
    }
    labs = as.factor(cl$pamobject$clustering)
    return(randomForest(x=F[,-i], y=labs, proximity=TRUE, oob.prox=oob)$proximity)
  }

  sims = vapply(1:ncol(F), featFu, array(0,dim=c(nrow(F),nrow(F))))
  sim  = apply(sims, c(1,2), mean)
  dis  = (1 - sim)

  return(dis)

}

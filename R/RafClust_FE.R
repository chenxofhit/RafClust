#' FEATURE ENGINEERING FOR RafClust
#'
#' @param texpr matrix Expression values, cells x genes
#' @param frq numeric Frequency filtering fraction
#' @param gene_filter logical shold gene filtering be performed?
#' @return list with two components: feature matrix F and the result of gene clustering cl.gene
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom fpc pamk
RafClustFE <- function(texpr, gene_filter=TRUE, frq = 0.06, verbose = FALSE){
  
  #- GENE FILTERING
  if (gene_filter == TRUE){
    texpr  = texpr[, colSums(texpr>0) > floor(frq*nrow(texpr)) ]
  }
  
  #- PCA EMBEDDING OF GENES
  # pc.gen  = rpca(t(texpr), retx=TRUE, center=TRUE, scale=TRUE)
  # nc      = max(which(cumsum(pc.gen$sdev)/sum(pc.gen$sdev)< .9))
  # nc      = length(elbow_detection((pc.gen$sdev[1:nc])))+1 #- could use log
  # gen.dat = scale(pc.gen$x[,1:nc])
  if(verbose){
    message("  Calculating cells  embedding  expression...")
  }
  gen.dat <- fast.pca(scale(texpr, scale=TRUE, center=TRUE))

  # #- K-MEANS CLUSTERING
  # clres   = sapply(2:10,  function(nc)
  #                             sum(KMeans_rcpp(gen.dat, max_iters = 50,
  #                             clusters = nc,initializer="kmeans++")$WCSS_per_cluster))
  # nc      = max(elbow_detection(log(clres))) + 1
  # cl.gene = KMeans_rcpp(gen.dat, clusters = nc, num_init = 5, max_iters = 100, initializer = 'kmeans++')


  # #- SPEARMAN FEATURE CONSTRUCTION
  # sfc <- function(mat){

  # }

  # NULLing the variables to avoid notes in R CMD CHECK
  i <- NULL
    
  distances <- c("euclidean", "pearson", "spearman")
  if(verbose){
    message("  Calculating distances between the cells...")
  }
  
  # Calculate the number of cores
  library(parallel)
  n_cores <- min(detectCores() - 1, 3)
  
  cl <- makeCluster(n_cores, type="FORK")
  registerDoParallel(cl)
    
  # calculate distances in parallel
  dists <- foreach(i = distances, .export=c("calculate_distance")) %dorng% {
    try({
        fast.pca(calculate_distance(texpr, i))
      })
  }
    
  # stop local cluster
  stopCluster(cl)

  names(dists) <- distances
  
  #dists.reduce <- lapply(1:length(distances),function(i) (dists[i]))
  #F    = sapply(1:nc,function(cind) sfc(texpr[, cl.gene$clusters == cind]))
  dists.mat = matrix(unlist(dists), nrow=nrow(texpr))
  
  features <- cbind(dists.mat,gen.dat)
  
  Fmat = scale(features, scale=TRUE, center=TRUE)

  # #- components per cluster in names
  # ncomps         = unlist(lapply(F,ncol))
  # names          = sapply(1:nc, function(ci) paste("C",rep(ci, ncomps[ci]), "-F", 1:ncomps[ci], sep=""))
  # names          = unlist(names)
  # colnames(Fmat) = names
  rownames(Fmat) = rownames(texpr)
  # #- return FEATURE MATRIX
  return(list(
    F = Fmat
    #,
    #geneClusters = cl.gene
  ))
}


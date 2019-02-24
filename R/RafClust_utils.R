normal <- function(x) {
 	x = x - min(x, na.rm=TRUE)
 	x = x/max(x,na.rm=TRUE)
 	return(x)
}

#' Calculate a distance matrix
#'
#' Distance between the cells, i.e. columns, in the input expression matrix are
#' calculated using the Euclidean, Pearson and Spearman metrics to construct
#' distance matrices.
#'
#' @param data expression matrix
#' @param method one of the distance metrics: 'spearman', 'pearson', 'euclidean'
#' @return distance matrix
#'
#' @importFrom stats cor dist
#' 
#' @useDynLib SC3
#' @importFrom Rcpp sourceCpp
#'
calculate_distance <- function(data, method) {
  dist <- NULL
  if (method == "spearman") {
    dist  <- as.matrix(1 - cor(t(data), method = "spearman"))
  } else if (method == "pearson") {
    dist <- as.matrix(1 - cor(t(data), method = "pearson"))
  } else {
    dist <- as.matrix(dist(data))
  }
  
  dist.norm <- scale(dist,center = TRUE, scale = TRUE)
  return(dist.norm)
}

#' ELBOW DETECTION
#'
#' @param scores vector
elbow_detection <-function(scores, if_plot = FALSE) {
 # We included this function from uSORT package, version 1.6.0 .
 #=====================
  num_scores <- length(scores)
  if (num_scores < 2) {
    stop("Input scores must be a vector with length more than 1!")
  }
  scores <- data.frame(id = seq_len(num_scores), value = scores)
  sorted_scores <- scores[order(scores$value, decreasing = TRUE), 
                          ]
  xy_coordinates <- cbind(x = seq_len(num_scores), y = sorted_scores$value)
  start_point <- xy_coordinates[1, ]
  end_point <- xy_coordinates[num_scores, ]
  x1 <- start_point[1]
  x2 <- end_point[1]
  y1 <- start_point[2]
  y2 <- end_point[2]
  a <- y1 - y2
  b <- x2 - x1
  c <- x1 * y2 - x2 * y1
  dist_to_line <- abs(a * xy_coordinates[, "x"] + b * xy_coordinates[, 
                                                                     "y"] + c)/sqrt(a^2 + b^2)
  best_point_id <- which.max(dist_to_line)
  score_cutoff <- xy_coordinates[best_point_id, "y"]
  select_ID <- scores$id[which(scores$value >= score_cutoff)]
  if (if_plot) {
    plot(seq(nrow(scores)), sorted_scores$value, col = ifelse(sorted_scores$value >= 
                                                                score_cutoff, "red", "black"), xlab = "ID", ylab = "Score", 
         main = paste0("Optimal number = ", length(select_ID), 
                       " with cutoff value = ", round(score_cutoff, 
                                                      digits = 4)))
  }
  return(select_ID)
}


##RSV

#' RANDOM PROJECTION SVD
#'
#' @param A matrix
#' @param K rank
#' @return list with components U S and V
#' @importFrom pracma orth
fast.rsvd<-function( A, K ) {
#============================

  M = dim(A)[1]
  N = dim(A)[2]
  P = min(2*K,N)
  X = matrix(rnorm(N*P),nrow=N,ncol=P)
  Y = A%*%X
  W1 = orth(Y)
  B = t(W1)%*%A
  res = svd(B,nu=min(dim(B)),nv=min(dim(B)))
  W2 = res$u
  tmp_S = res$d
  S = array(0,c(length(tmp_S),length(tmp_S)))
  diag(S) = tmp_S
  V = res$v
  U = W1%*%W2
  K = min(K,dim(U)[2])
  U = U[,1:K]
  S = S[1:K,1:K]
  V = V[,1:K]

  return(list(U=U,S=S,V=V))
}



#' PCA built on fast.rsvd
#'
#' @param X matrix
#' @return projection of X on principal top principal components
#' @importFrom matrixStats colVars
fast.pca<-function(X){
#=====================

  K<-400
  #X = t(X)
  tmp_val = as.vector(colSums(X)/nrow(X))
  X = X - t(apply(array(0,dim(X)),MARGIN=1,FUN=function(x) {x=tmp_val}))
  res = fast.rsvd(X,K)
  U = res$U
  S = res$S
  K = min(dim(S)[2],K)
  diag_val = sqrt(diag(S[1:K,1:K]))
  diag_mat = array(0,c(length(diag_val),length(diag_val)))
  diag(diag_mat) = diag_val
  X = U[,1:K]%*%diag_mat
  normalization_val = sqrt(rowSums(X*X))
  X = X / apply(array(0,c(length(normalization_val),K)),MARGIN=2,FUN=function(x) {x=normalization_val})
  pcgeneres<-X
  varlist<-colVars(pcgeneres)
  ordered_varlist<-order(varlist,decreasing = TRUE)
  LM1<-pcgeneres[,ordered_varlist]
  varf<-colVars(LM1)

  num<-length(elbow_detection(varf))
  pcgene<-LM1[,c(1:num)]

  return(pcgene)
}

#' Plot a (dis)similarity matrix using the heatmap function
#'
#' @param x matrix (Dis)Similarities
#' @param labels integer Labels for classes/categories
#' @param col character "qnt" for quantile-based color breaks, "lin" for linear color breaks
#' @return list The list of components heatmap returns.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @export
plotSIM = function(x,labels,col="qnt") {
 #--------------------------------

   diag(x) = NA
   if(col=="lin"){
     x = x-min(x,na.rm=TRUE)
     x = x/max(x,na.rm=TRUE)
     brks  = seq(0,1,len=65)
     col  = rev(colorRampPalette(brewer.pal(9,"PuBuGn"))(64))
   } else {
     brks = quantile(as.vector(x),seq(0,1,len=65),na.rm=TRUE)
     col  = rev(colorRampPalette(brewer.pal(9,"PuBuGn"))(64))
   }
   nc  = length(unique(labels))
   cl = brewer.pal(nc,"Set1")

   #col = paste(col,"77",sep="")
   hm = heatmap(x,Rowv=NA,Colv=NA,scale="none",labRow="",labCol="",
                margins=c(.7,.7),col=col, breaks=brks,ColSideColors=cl[labels],
                RowSideColors=cl[labels],revC=TRUE)
   return(invisible(hm))
}

#' Plot a tSNE plot of a (dis)similarity matrix 
#'
#' @param x matrix (Dis)Similarities
#' @param labels integer Labels for classes/categories
#' @param seed integer seed for RNG (->tSNE)
#' @param ... further arguments to Rtsne
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics plot
#' @importFrom Rtsne Rtsne
#' @export
plotTSNE <- function(x,labels,seed=3749829,...) {
  #-------------------------------
  set.seed(seed)
  a   = Rtsne(x,...)
  nc  = length(unique(labels))
  col = brewer.pal(nc,"Set1")
  col = paste(col,"77",sep="")
  plot(a$Y,col=col[labels],pch=19,cex=1,xaxt="n",yaxt="n",xlab="",ylab="")
}


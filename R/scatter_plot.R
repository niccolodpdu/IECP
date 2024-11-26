#' @name scatter_plot
#' @title A function for illustrating data points in scatter plot with optional marker color-coding
#' @param data With samples on the columns and features on the rows.
#' @param mg.info Do you have the info for marker genes?
#' @param mg A list of marker genes
#' @param mg.col A vector of colors
#' @export

scatter_plot<-function(data, mg.info = FALSE, mg = NULL, 
                       mg.col = c('red','orange2','dodgerblue','green4','purple','pink')){
  
  library(corpcor)
  
  X <- data
  Xproj <- X
  
  Xproj <- apply(Xproj, 1, function(x) x / sum(x))
  
  A <- diag(dim(Xproj)[1])
  
  K = dim(A)[2]
  
  PS <- t(matrix(c(cos((seq_len(K) - 1) * 2 * pi / K),
                   sin((seq_len(K) - 1) * 2 * pi / K)), K))
  
  if (is.null(colnames(data))){
    name<- paste("Column", 1:K, sep="")
  } else {
    name<- colnames(data)
  }
  
  colnames(PS)<-name
  row.names(PS)<-c('X','Y')
  message("Coordinates of the vertices and their meaning:")
  print(PS)
  
  tmp <- PS %*% pseudoinverse(A)
  Xp <- tmp %*% Xproj
  
  plot(t(PS), xlab = NA, ylab = NA, asp = 1, pch = '.')
  points(Xp[1,], Xp[2,], col = rgb(128,128,128, max = 255, alpha = 125),
         xlab = NA, ylab = NA, asp = 1)
  
  #### If you want to add some marker genes, here is your place ####
  
  if (mg.info == TRUE) {
    for(i in seq_along(mg)){
      points(Xp[1,mg[[i]]], Xp[2,mg[[i]]], col=mg.col[i],pch=16)
    }
  }
  
  
  points(t(PS), pch = 17,cex=0.5)
}
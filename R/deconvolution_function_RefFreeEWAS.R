#' RefFreeEwas package is no longer available online. Here we provide the R scripts for users' reference.
#################################################################################
# RefFreeEwasModel: Reference-free cell-mixture-adjusted EWAS
# Author & paper: Houseman, E.A., Kile, M.L., Christiani, D.C. et al. Reference-free deconvolution of DNA methylation data and mediation by cell composition effects. BMC Bioinformatics 17, 259 (2016). https://doi.org/10.1186/s12859-016-1140-4
#################################################################################
RefFreeEwasModel <- function(
    Y,   
    X,   
    K,   
    smallOutput=FALSE  
){
  n1 <- dim(Y)[1]
  n2 <- dim(X)[1]
  pdim <- dim(X)[2]
  
  noMiss <- !apply(is.na(Y),1,any)
  allObs <- all(noMiss)
  
  HX <- solve(t(X)%*%X)
  PX <- X %*% HX
  
  if(allObs){ # If no missing value
    Bstar <- Y %*% PX
    
    muStar <- Bstar %*% t(X)
    Estar <- Y - muStar
    
    sstar <- sqrt(apply(Estar*Estar,1,sum)/(n2-pdim))
    
    BetaE <- cbind(Bstar, Estar)
    svdStar <- svdSafe(BetaE)
    
    Lambda <- t(svdStar$d[1:K]*t(svdStar$u[,1:K]))
    U <- svdStar$v[,1:K]
  }
  else{ #If missing values, do as much as possible on one fell swoop,
    # and for the rest, do on a CpG-by-CpG basis
    Bstar <- matrix(NA, n1, pdim)
    degfree <- rep(n2-pdim, n1)
    Bstar[noMiss,] <- Y[noMiss,] %*% PX
    whichMiss <- which(!noMiss)
    nMiss <- length(whichMiss)
    
    for(j in 1:nMiss){
      jj <- whichMiss[j]
      mflag <- !is.na(Y[jj,])
      Xj <- X[mflag,,drop=FALSE]
      HXj <- solve(t(Xj)%*%Xj)
      PXj <- Xj %*% HXj
      Bstar[jj,] <- Y[jj,mflag,drop=FALSE]%*%PXj
      degfree[jj] <- sum(mflag)-pdim
    }
    
    muStar <- Bstar %*% t(X)
    Estar <- Y - muStar
    
    sstar <- sqrt(apply(Estar*Estar,1,sum,na.rm=TRUE)/degfree)
    BetaE <- cbind(Bstar, Estar)
    svdStar <- svdSafe(BetaE[noMiss,])
    
    Lambda <- matrix(NA, n1, K)
    Lambda[noMiss,] <- t(svdStar$d[1:K]*t(svdStar$u[,1:K]))
    U <- svdStar$v[,1:K]
    for(j in 1:nMiss){
      jj <- whichMiss[j]
      mflag <- c(rep(TRUE,pdim), !is.na(Y[jj,]))
      Uj <- U[mflag,,drop=FALSE]
      Lambda[jj,] <- solve(t(Uj)%*%Uj,t(Uj)%*%BetaE[jj,mflag])
    }
  }
  
  LambdaProjBstar <- solve(t(Lambda)%*%Lambda, t(Lambda)%*%Bstar)
  Beta <- Bstar - Lambda %*% LambdaProjBstar
  
  out <- list(Bstar=Bstar, Beta=Beta, sstar=sstar, Lambda=Lambda, U=U, d=svdStar$d)
  
  if(!smallOutput) {
    muStar <- ifelse(muStar<0.00001,0.00001,muStar)
    muStar <- ifelse(muStar>0.99999,0.99999,muStar)
    out$dispersion <- sqrt(muStar*(1-muStar))
    out$E <- Estar/out$dispersion
    out$X <- X
  }
  out
  class(out) <- "RefFreeEwasModel"
  out
}

#' @export
print.RefFreeEwasModel <- function(x,...){
  cat("Reference Free EWAS Model\n\n")
  cat("Assay matrix: ", dim(x$Beta)[1], " features\n")
  if(!is.null(x$X)) cat("Design matrix: ", dim(x$X)[1], 
                        " subjects x ", dim(x$X)[2]," covariates\n\n", sep="")
  else cat("(small output version)\n\n")
  
  cat(dim(x$Lambda)[2]," latent variables\n\n")
}


#################################################################################
# BootOneRefFreeEwasModel: Create *one( bootstrap sample from
#     reference-free cell-mixture-adjusted EWAS
#################################################################################
BootOneRefFreeEwasModel <- function(mod){
  n2 <- dim(mod$X)[1]
  iboot <- sample(1:n2, n2, replace=TRUE)
  
  mu <- mod$Bstar %*% t(mod$X)
  
  return(mu + mod$dispersion*mod$E[,iboot])
}



#################################################################################
# BootRefFreeEwasModel: Create several bootstrap samples from
#     reference-free cell-mixture-adjusted EWAS
#     and return Beta estimates
#################################################################################
BootRefFreeEwasModel <- function(
    mod,   
    nboot  
){
  BetaBoot <- array(NA, dim=c(dim(mod$Beta),2,nboot))
  dimnames(BetaBoot)[1:2] <- dimnames(mod$Beta) 
  dimnames(BetaBoot)[[3]] <- c("B","B*")
  dimnames(BetaBoot)[[4]] <- 1:nboot
  attr(BetaBoot,"nSample") <- dim(mod$X)[1]
  
  for(r in 1:nboot){
    isError <- TRUE
    while(isError){
      catchError <- try({
        Yboot <- BootOneRefFreeEwasModel(mod)
        bootFit <- RefFreeEwasModel(Yboot, mod$X, 
                                    dim(mod$Lambda)[2], smallOutput=TRUE)
        BetaBoot[,,1,r] <- bootFit$Beta
        BetaBoot[,,2,r] <- bootFit$Bstar
      })
      isError <- inherits(catchError,"try-error")
    }
    if(r%%10==0) cat(r,"\n")
  }
  class(BetaBoot) <- "BootRefFreeEwasModel"
  BetaBoot
}

#' @export
summary.BootRefFreeEwasModel <- function(object,...){
  
  x <- object
  
  out <- array(NA, dim=c(dim(x)[1:3],2))
  dimnames(out)[1:3] <- dimnames(x)[1:3]
  dimnames(out)[[4]] <- c("mean","sd")
  
  out[,,,1] <- apply(x,c(1:3),mean)
  out[,,,2] <- apply(x,c(1:3),sd)
  
  class(out) <- "summaryBootRefFreeEwasModel"
  attr(out,"nBoot") <- dim(x)[4]
  attr(out,"nSample") <- attr(x,"nSample")
  
  out
}

#' @export
print.summaryBootRefFreeEwasModel <- function(x,...){
  cat(attr(x,"nBoot"),"bootstrap samples, n =", attr(x,"nSample"),"subjects\n\nBeta Mean\n")
  print(x[1:6,,1,1],...)
  cat("\nBeta Standard Deviation\n")
  print(x[1:6,,1,2],...)
}

#' @export
print.BootRefFreeEwasModel <- function(x,...){
  print(summary(x),...)
}

#################################################################################
# svdSafe: svd that traps errors and switches to QR when necessary
#################################################################################
svdSafe <- function(X){
  sv <- try(svd(X), silent=TRUE)
  if(inherits(sv,"try-error")){
    warning("SVD algorithm failed, using QR-decomposition instead")
    QR <- qr(X)
    sv <- list(d=rep(1,dim(X)[2]))
    sv$u <- qr.Q(QR)
    sv$v <- t(qr.R(QR))
  }
  sv
}

#################################################################################
# EstDimIC: Estimate dimension using AIC and BIC
#################################################################################
EstDimIC <- function (
    Rmat,             
    Krange=0:25
){
  N1 <- dim(Rmat)[1]
  N2 <- dim(Rmat)[2]
  svdRmat <- svdSafe(Rmat)
  nK = length(Krange)
  tmpAIC <- tmpBIC <- rep(NA,nK)
  for(Ktest in Krange){
    if(Ktest==0) tmpRminLU <- Rmat
    else tmpRminLU <- (Rmat -
                         svdRmat$u[,1:Ktest] %*% (svdRmat$d[1:Ktest] * t(svdRmat$v[,1:Ktest])))
    tmpSigSq <- apply(tmpRminLU*tmpRminLU,1,sum)/N2
    tmpAIC[Ktest+1] <- 2*(N1+Ktest*(N1+N2)) + N1*N2 + N2*sum(log(tmpSigSq))
    tmpBIC[Ktest+1] <- log(N2)*(N1+Ktest*(N1+N2)) + N1*N2 + N2*sum(log(tmpSigSq))
  }
  list(icTable=cbind(K=Krange,AIC=tmpAIC,BIC=tmpBIC), 
       best=Krange[c(AIC=which.min(tmpAIC),BIC=which.min(tmpBIC))])
}

#################################################################################
# omnibusBoot: Omnibus test of association using bootstraps
#################################################################################
omnibusBoot <- function(est, boots, denDegFree){
  nFeature <- length(est)
  se <- apply(boots,1,sd)
  pv <- 2*pt(-abs(est)/se,denDegFree)
  pvNull <- 2*pt(-(1/se)*abs(sweep(boots,1,apply(boots,1,mean),"-")),denDegFree+1)
  
  ks <- max(abs(sort(pv)-(1:nFeature-0.5)/nFeature))
  ksNull <-apply(pvNull,2,function(u)max(abs(sort(u)-(1:nFeature-0.5)/nFeature)))
  mean(ks<=ksNull)
}

#################################################################################
# EstDimRMT: Estimate dimension using Random Matrix Theory
#  Note:  this function was originally authored by A. Teschendorff in the
#         package isva.
#         Previous versions of RefFreeEWAS used the isva version of the function.
#         However, because of dependency issues in that package, this version
#         simply reproduces the function found in version 1.9 of isva and
#         removes the dependency on the isva package
#         Plotting functionality also removed from original source.
#################################################################################
EstDimRMT <- function (Rmat) 
{
  data.m <- Rmat
  M <- apply(data.m, 2, function(X) {
    (X - mean(X))/sqrt(var(X))
  })
  sigma2 <- var(as.vector(M))
  Q <- nrow(data.m)/ncol(data.m)
  ns <- ncol(data.m)
  lambdaMAX <- sigma2 * (1 + 1/Q + 2 * sqrt(1/Q))
  lambdaMIN <- sigma2 * (1 + 1/Q - 2 * sqrt(1/Q))
  delta <- lambdaMAX - lambdaMIN
  roundN <- 3
  step <- round(delta/ns, roundN)
  while (step == 0) {
    roundN <- roundN + 1
    step <- round(delta/ns, roundN)
  }
  lambda.v <- seq(lambdaMIN, lambdaMAX, by = step)
  dens.v <- vector()
  ii <- 1
  for (i in lambda.v) {
    dens.v[ii] <- (Q/(2 * pi * sigma2)) * sqrt((lambdaMAX - 
                                                  i) * (i - lambdaMIN))/i
    ii <- ii + 1
  }
  thdens.o <- list(min = lambdaMIN, max = lambdaMAX, step = step, 
                   lambda = lambda.v, dens = dens.v)
  C <- 1/nrow(M) * t(M) %*% M
  eigen.o <- eigen(C, symmetric = TRUE)
  estdens.o <- density(eigen.o$values, from = min(eigen.o$values), 
                       to = max(eigen.o$values), cut = 0)
  intdim <- length(which(eigen.o$values > thdens.o$max))
  evalues.v <- eigen.o$values
  return(list(cor = C, dim = intdim, estdens = estdens.o, thdens = thdens.o, 
              evals = eigen.o$values))
}


projectMix <- function(Y, Xmat, nonnegative=TRUE, sumLessThanOne=TRUE, lessThanOne=!sumLessThanOne){
  
  nCol = dim(Xmat)[2]
  nSubj = dim(Y)[2]
  
  mixCoef = matrix(0, nSubj, nCol)
  rownames(mixCoef) = colnames(Y)
  colnames(mixCoef) = colnames(Xmat)
  
  if(nonnegative){
    
    if(sumLessThanOne){
      Amat = cbind(rep(-1,nCol), diag(nCol))
      b0vec = c(-1,rep(0,nCol))
    }
    else if(lessThanOne){
      Amat = cbind(-diag(nCol), diag(nCol))
      b0vec = c(rep(-1,nCol),rep(0,nCol))
    }
    else{
      Amat = diag(nCol)
      b0vec = rep(0,nCol)
    }
    
    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec)$sol
    }
  }
  else{
    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
    }
  }
  
  return(mixCoef)
}


RefFreeCellMix <- function(Y,mu0=NULL,K=NULL,iters=10,Yfinal=NULL,verbose=TRUE){
  if(is.null(mu0)){
    if(K==1) {
      if(!is.null(Yfinal)) Y <- Yfinal
      n <- dim(Y)[2]
      
      mu <- matrix(apply(Y,1,mean,na.rm=TRUE),ncol=1)
      omega <- matrix(1, n, 1)
      o <- list(Mu=mu, Omega=omega)
      class(o) <- "RefFreeCellMix"
      return(o)
    }
    else mu0 <- RefFreeCellMixInitialize(Y,K=K,method="ward")
  }
  
  incrementalChangeSummary <- list()
  for(i in 1:iters){
    flag <- !apply(is.na(mu0),1,any)
    omega <- projectMix(Y[flag,],mu0[flag,])
    mu <- projectMix(t(Y),omega,sumLessThanOne=FALSE)
    incrementalChangeSummary[[i]] <- summary(abs(as.vector(mu-mu0)))  
    if(verbose) print(incrementalChangeSummary[[i]])
    mu0 <- mu
  }
  if(!is.null(Yfinal)){
    mu <- projectMix(t(Yfinal),omega,sumLessThanOne=FALSE)
  }
  
  o <- list(Mu=mu, Omega=omega, incrementalChangeSummary=incrementalChangeSummary)
  class(o) <- "RefFreeCellMix"
  o
}

#' @export
print.RefFreeCellMix <- function(x,...){
  cat("Reference Free Deconvolution\n\n")
  cat("Mu: ",dim(x$Mu)[1]," cpgs x ",dim(x$Mu)[2],"cell types\n")
  cat("Omega :",dim(x$Omega)[1]," subjects x ",dim(x$Omega)[2],"cell types\n")
}

#' @export
summary.RefFreeCellMix <- function(object,...){
  list(Mu=apply(object$Mu,2,summary), Omega=apply(object$Omega,2,summary), MuCorr=cor(object$Mu))
}

RefFreeCellMixArray <- function(Y,Klist=1:5,iters=10,Yfinal=NULL,verbose=FALSE,dist.method = "euclidean",...){
  D <- dist(t(Y),method=dist.method)
  hc <- hclust(D,...)
  
  rfcmArray <- list()
  nK <- length(Klist)
  for(r in 1:nK){
    cat("Fitting K =",Klist[r],"\n")
    if(Klist[r]==1){
      rfcmArray[[r]] <- RefFreeCellMix(Y,K=1,Yfinal=Yfinal,iters=iters)
    }
    else{
      rfcmArray[[r]] <- RefFreeCellMix(Y,mu0=RefFreeCellMixInitialize(Y,K=Klist[r],Y.Cluster=hc),
                                       Yfinal=Yfinal, verbose=verbose,iters=iters)
    }
  }
  names(rfcmArray) <- Klist
  rfcmArray
}

#' @export
deviance.RefFreeCellMix <- function(object, Y, Y.oob=NULL, EPSILON=1E-9, 
                                    bootstrapIterations=0, bootstrapIndices=NULL, ...){
  
  N <- dim(Y)[2]
  if(bootstrapIterations>0){# Do the bootstrap and replace x (but initialize with x$Mu)
    if(is.null(bootstrapIndices)){
      boots <- sample(1:N, N, replace=TRUE)
    }  
    else {
      boots <- bootstrapIndices
    }
    Y.oob <- Y[,-unique(boots)]
    Y <- Y[,boots]
    if(dim(object$Mu)[2]==1) {
      object <- RefFreeCellMix(Y,K=1,iters=bootstrapIterations,verbose=FALSE)
    }
    else{
      object <- RefFreeCellMix(Y,mu0=object$Mu,iters=bootstrapIterations,verbose=FALSE)
    }
  }
  
  Y.mu <- object$Mu %*% t(object$Omega)
  R <- Y-Y.mu
  Y.n <- apply(!is.na(Y),1,sum)
  Y.SSQ <- apply(R*R,1,sum,na.rm=TRUE)
  logSigma2 <- log(pmax(EPSILON,Y.SSQ)) - log(Y.n)
  
  if(!is.null(Y.oob)){
    Omega.oob <- projectMix(Y.oob,object$Mu)
    Y.mu <- object$Mu %*% t(Omega.oob)
    R.oob <- Y.oob-Y.mu
    n.oob <- apply(!is.na(Y.oob),1,sum)
    SSQ.oob <- apply(R.oob*R.oob,1,sum,na.rm=TRUE)
    N <- dim(Y.oob)[2]
  }
  else{
    SSQ.oob <- Y.SSQ
    n.oob <- Y.n
  }
  
  sum( n.oob*log(2*pi)+n.oob*logSigma2+SSQ.oob/exp(logSigma2))/N
}

RefFreeCellMixInitialize <- function(Y,K=2,Y.Distance=NULL, Y.Cluster=NULL, 
                                     largeOK=FALSE, dist.method = "euclidean", ...){
  
  if(!is.matrix(Y) | !is.numeric(Y)){
    stop("Y is not a numeric matrix\n")
  }
  n <- dim(Y)[2]
  
  if(is.null(Y.Cluster)){
    if(is.null(Y.Distance)){
      if(n>2500 & !largeOK){
        stop("Y has a large number of subjects!  If this is what you really want, change 'largeOK' to TRUE\n")
      }
      Y.Distance <- dist(t(Y),method=dist.method)
    }
    Y.Cluster <- hclust(Y.Distance,...)
  } 
  
  classes <- cutree(Y.Cluster, K)
  s <- split(1:n,classes)
  
  sapply(s, function(u) apply(Y[,u,drop=FALSE],1,mean,na.rm=TRUE))
}

RefFreeCellMixArrayDevianceBoot <- function(rfArray, Y, EPSILON=1E-9, bootstrapIterations=5){
  N <- dim(Y)[2]
  boots <- sample(1:N, N, replace=TRUE)
  sapply(rfArray, deviance.RefFreeCellMix, Y=Y, EPSILON=EPSILON, 
         bootstrapIterations=bootstrapIterations,bootstrapIndices=boots)
} 

RefFreeCellMixArrayDevianceBoots <- function(rfArray, Y, R=5, EPSILON=1E-9, bootstrapIterations=5){
  dv <- sapply(rfArray, deviance, Y=Y)
  nK <- length(dv)
  devs <- matrix(NA,R,nK)
  for(r in 1:R){
    if(r %% 10==0) cat("Bootstrap",r,"\n")
    devs[r,] <- RefFreeCellMixArrayDevianceBoot(rfArray, Y, 
                                                EPSILON=EPSILON, bootstrapIterations=bootstrapIterations)
  }
  out <- rbind(dv,devs)
  rownames(out) <- 0:R
  out
} 

ImputeByMean = function(Y){
  Yimpute=Y
  whichMiss = apply(is.na(Y),1,which)
  nmiss = sapply(whichMiss,sum)
  if(any(nmiss>0)) {
    for(i in which(nmiss>0)){
      Yimpute[i,whichMiss[[i]]] = 
        mean(Y[i,-whichMiss[[i]]])
    }
  }
  Yimpute
}

SVDwithMissing = function(Y){
  Yimpute = ImputeByMean(Y)
  svd(Yimpute)
}

RefFreeCellMixArrayWithCustomStart = function(Y, mu.start, 
                                              Klist = 1:5, iters = 10, Yfinal = NULL,  verbose = FALSE){
  
  rfcmArray <- list()
  nK <- length(Klist)
  for (r in 1:nK) {
    cat("Fitting K =", Klist[r], "\n")
    if (Klist[r] == 1) {
      rfcmArray[[r]] <- RefFreeCellMix(Y, K = 1, Yfinal = Yfinal, iters=iters)
    }
    else {
      rfcmArray[[r]] <- RefFreeCellMix(Y, mu0 = mu.start[,1:Klist[r]], 
                                       Yfinal = Yfinal, verbose = verbose, iters=iters)
    }
  }
  names(rfcmArray) <- Klist
  rfcmArray
  
}

RefFreeCellMixInitializeBySVD = function(Y, type=1){
  
  Y.svd = SVDwithMissing(Y)
  
  nn = ncol(Y.svd$u)
  Y.svd.sign = sapply(1:nn, function(i)sign(mean(sign(Y.svd$v[,i]))))
  Y.svd.sign[Y.svd.sign==0] = 1
  Y.svd.u = Y.svd.sign*t(Y.svd$u)
  if(type==1){
    mu.start = apply(Y.svd.u,1,function(x)(sign(x-median(x))+1)/2)
  }
  else {
    mu.start = apply(Y.svd.u, 1, function(x) rank(x)-0.5)/ncol(Y.svd.u)
  }
  mu.start
}


bootstrapPairs <- function(obs, pairID){
  
  pairIx <- split(obs, pairID)
  nobsTab <- table(sapply(pairIx, length))
  if(length(nobsTab)>1) {
    stop("All clusters must have the same number of observations.\n")
  }
  
  n <- length(pairIx)
  pairBoot <- pairIx[sample(1:n, n, replace=TRUE)]
  
  obsBoot <- unlist(lapply(pairBoot, function(u) sample(u,replace=FALSE)))
  names(obsBoot) <- NULL
  obsBoot
}

PairsBootRefFreeEwasModel <- function (mod, nboot, pairID) 
{
  BetaBoot <- array(NA, dim = c(dim(mod$Beta), 2, nboot))
  dimnames(BetaBoot)[1:2] <- dimnames(mod$Beta)
  dimnames(BetaBoot)[[3]] <- c("B", "B*")
  dimnames(BetaBoot)[[4]] <- 1:nboot
  attr(BetaBoot, "nSample") <- dim(mod$X)[1]
  for (r in 1:nboot) {
    isError <- TRUE
    while (isError) {
      catchError <- try({
        Yboot <- PairsBootOneRefFreeEwasModel(mod, pairID)
        bootFit <- RefFreeEwasModel(Yboot, mod$X, dim(mod$Lambda)[2], 
                                    smallOutput = TRUE)
        BetaBoot[, , 1, r] <- bootFit$Beta
        BetaBoot[, , 2, r] <- bootFit$Bstar
      })
      isError <- inherits(catchError, "try-error")
    }
    if (r%%10 == 0) 
      cat(r, "\n")
  }
  class(BetaBoot) <- "BootRefFreeEwasModel"
  BetaBoot
}

PairsBootOneRefFreeEwasModel <- function (mod, pairID) 
{
  n2 <- dim(mod$X)[1]
  iboot <- bootstrapPairs(1:n2, pairID)
  mu <- mod$Bstar %*% t(mod$X)
  return(mu + mod$dispersion * mod$E[, iboot])
}

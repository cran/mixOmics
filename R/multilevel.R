# Copyright (C) 2012 
# Benoit Liquet, Université de Bordeaux, France
# Kim-Anh Lê Cao, Queensland Facility for Advanced Bioinformatics, University of Queensland, Brisbane, Australia
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# ----------------------------------------
# generic function for multilevel analysis 
# ----------------------------------------

multilevel <-  function (X, Y = NULL, 
                         cond = NULL,                  
                         sample = NULL,
                         ncomp = 1,
                         keepX = rep(ncol(X), ncomp),
                         keepY = NULL,  # for spls
                         method = NULL,
                         tab.prob.gene=NULL, # remove?
                         max.iter = 500, 
                         tol = 1e-06,...) {
  
  # check parameter input
  check.one.level(X, Y, 
                  cond ,                  
                  sample,
                  ncomp,
                  keepX,
                  keepY,  # for spls
                  method ,
                  tab.prob.gene, # remove?
                  max.iter, 
                  tol,...)
  
  #   # checking general input parameters
  #   X = as.matrix(X)
  #   # X input
  #   if (length(dim(X)) != 2 || !is.numeric(X)) 
  #     stop("'X' must be a numeric matrix.")
  #   
  #   # Testing the cond vector
  #   if(is.null(cond)) stop('Vector cond is missing', call. = FALSE)
  #     
  #   # ncomp
  #   if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0) stop("invalid number of components, 'ncomp'.")
  #   
  #   #checking sample
  #   if(is.null(sample)) stop('Vector sample is missing', call. = FALSE)
  #   if(!is.null(dim(sample))) stop('sample should be a vector indicating the repeated measurements')
  #   if(length(sample) != nrow(X)) stop('X and the vector sample should have the same number of subjects')
  #   # check that the sample numbers are repeated
  #   if(length(summary(as.factor(sample))) == nrow(X)) stop('Check that the vector sample reflects the repeated measurements')
  
  if(is.factor(sample)){
    sample = as.numeric(sample)
    warning('the vector sample was converted into a numeric vector', call. = FALSE)
  }
  
  #   #check that the sample are numbered from 1
  #   if(!any(names(summary(as.factor(sample))) == '1')) {
  #     cat('The vector sample includes the values: ', as.vector(names(summary(as.factor(sample)))), '\n')
  #     stop('sample vector', call. =FALSE)
  #   }
  #   
  #   # method
  #   if(is.null(method)) stop('Input method missing, should be set to splsda or spls', call. = FALSE)
  
  
  # call the multilevel approach
  if(method == 'splsda'){
    # apply one or two-factor analysis with splsda
    result = multilevel.splsda(X = X, cond = cond, sample = sample, ncomp=ncomp, keepX = keepX, tab.prob.gene = tab.prob.gene)
  }else{
    # spls
    result = multilevel.spls(X = X, Y=Y, cond = cond, sample = sample, ncomp = ncomp, keepX=keepX, keepY = keepY, tab.prob.gene = tab.prob.gene)
  }
  
  return(result)
}


# ----------------------------------------
# multilevel analysis with sPLS-DA (1 data set X)
# ----------------------------------------
multilevel.splsda <-  function (X, 
                                cond = NULL,
                                sample = NULL,
                                ncomp = 2,
                                keepX = rep(ncol(X), ncomp),
                                tab.prob.gene=NULL, # remove?
                                max.iter = 500, 
                                tol = 1e-06,...) {
  
  #internal booleans to choose between 1 or 2 factors analysis
  factor1 = FALSE
  factor2 = FALSE
  
  #-- check input parameters --#
  
  # Testing the input cond and setting the 1 or 2 factor analysis   
  if (is.null(dim(cond))) {  # if cond is a vector: 1 factor analysis
    factor1 = TRUE
    if(!is.factor(cond)){
      cond = as.factor(cond)
      warning('cond was set as a factor', call. = FALSE)
    }
  }else{              #else if cond is a matrix: 2 factor analysis
    if(ncol(cond) == 2){
      factor2 = TRUE
      if(!is.factor(cond[,1])){ 
        warning('First cond response was set as a factor', call. = FALSE)
      }
      if(!is.factor(cond[,2])){ 
        warning('Second cond response was set as a factor', call. = FALSE)
      }
      cond1 = as.factor(cond[,1])
      cond2 = as.factor(cond[,2])
    }else {
      stop("'cond' must be a matrix with max. 2 columns for the 2 factor analysis.")
    } # end if ncol(cond)
  }
  
  
  if(factor1 == TRUE){
    if(nrow(X) != length(cond)) stop('X and cond should have the same number of subjects')
  }
  if(factor2 == TRUE){
    if(nrow(X) != nrow(cond)) stop('X and cond should have the same number of subjects')
  }
  
  
  
  # --multilevel analysis  
  # Variance decomposition for 1 or 2 factors
  if(factor1 == TRUE){
    Xw <- Split.variation.one.level(X, cond, sample)$Xw
  }else{ #if factor2 = TRUE 
    Xw <- Split.variation.two.level(X, cond1, cond2, sample)$Xw
    cond <- as.factor(paste(cond1, cond2, sep="."))
  }
  
  #Apply sPLS-DA on the within matrix
  res <- splsda(Xw, cond, ncomp = ncomp, keepX, max.iter, tol, ...)
  class(res) <- "list"
  if(factor1){
    result <- c(res,list(Xw=Xw,sample=sample,name.condition=factor(cond),tab.prob.gene=tab.prob.gene))
    class(result) <- c("splsda1fact","splsda")
  }else{ #if factor2 = TRUE
    result <- c(res,list(Xw=Xw,sample=sample, name.condition=cond1,tab.prob.gene=tab.prob.gene, name.time=cond2))
    class(result) <- c("splsda2fact","splsda")}
  return(invisible(result))
}


# ----------------------------------------
# multilevel analysis with sPLS (2 data sets X and Y)
# ----------------------------------------
multilevel.spls <-  function (X, Y, 
                              cond, 
                              sample,
                              ncomp = 2,
                              keepX = rep(ncol(X), ncomp),
                              keepY = rep(ncol(Y), ncomp),
                              tab.prob.gene=NULL,
                              max.iter = 500, tol = 1e-06,...) {
  
  #-- check input parameters --#
  
  # Y input
  Y = as.matrix(Y)
  if (length(dim(Y)) != 2 || !is.numeric(Y)) 
    stop("'Y' must be a numeric matrix.")
  
  
  # set the mode 
  mode="canonical"
  
  # Decomposition of the variance
  Xw <- Split.variation.one.level(X,Y=cond,sample)$Xw
  Yw <- Split.variation.one.level(X=Y,Y=cond,sample)$Xw
  
  # call sPLS
  res <- spls(Xw,Yw, mode=mode, ncomp = ncomp, keepY = keepY, keepX= keepX) 
  result <- c(res,list(Xw=Xw,Yw=Yw,sample=sample,name.condition=factor(cond),tab.prob.gene=tab.prob.gene))
  
  if(!is.null(tab.prob.gene)){ 
    #############match probes and gene for the graph #######
    probeX <- result$names$X
    geneX <- result$tab.prob.gene[match(probeX,result$tab.prob.gene[,1]),2]
    result$names$X <- as.character(geneX)
    ###for the network and CIM representation
    colnames(result$X) <- as.character(geneX)}
  class(result) <- c("splslevel","spls") #,"pls")
  return(invisible(result))
}



# -----------------------------------------------------------------------------------------
#                               internal functions
# ----------------------------------------------------------------------------------------

# ----------------
# check.one.level
# -----------------
check.one.level = function(X, Y, 
                           cond ,                  
                           sample,
                           ncomp,
                           keepX,
                           keepY,  # for spls
                           method ,
                           tab.prob.gene, # remove?
                           max.iter, 
                           tol,...){
  
  # checking general input parameters
  X = as.matrix(X)
  # X input
  if (length(dim(X)) != 2 || !is.numeric(X)) 
    stop("'X' must be a numeric matrix.")
  
  # Testing the cond vector
  if(is.null(cond)) stop('Vector cond is missing', call. = FALSE)
  
  # ncomp
  if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0) stop("invalid number of components, 'ncomp'.")
  
  #checking sample
  if(is.null(sample)) stop('Vector sample is missing', call. = FALSE)
  if(!is.null(dim(sample))) stop('sample should be a vector indicating the repeated measurements')
  if(length(sample) != nrow(X)) stop('X and the vector sample should have the same number of subjects')
  # check that the sample numbers are repeated
  if(length(summary(as.factor(sample))) == nrow(X)) stop('Check that the vector sample reflects the repeated measurements')
  
  
  #check that the sample are numbered from 1
  if(!any(names(summary(as.factor(sample))) == '1')) {
    cat('The vector sample includes the values: ', as.vector(names(summary(as.factor(sample)))), '\n')
    stop('sample vector', call. =FALSE)
  }
  
  # method
  if(is.null(method)) stop('Input method missing, should be set to splsda or spls', call. = FALSE)
  
}

# --------------------------------------
### Split.variation.one.level: decomposition of the variance (within matrix and between matrix) for one level
# --------------------------------------

Split.variation.one.level <- function(X,Y,sample){
  
  X = as.matrix(X)
  
  # check sample is numeric
  if(is.factor(sample)){
    sample = as.numeric(sample)
    warning('the vector sample was converted into a numeric vector', call. = FALSE)
  }
  
  Xmi <- colMeans(X)
  Xm <- matrix(Xmi,nrow=nrow(X),ncol=ncol(X),byrow=T)
  indX <- cbind(sample,X)
  
  indsample <- unique(sample)
  n.sample <- length(indsample)
  Xbi <- t(apply(matrix(indsample,ncol=1,nrow=n.sample),MARGIN=1,FUN=function(x,indX){indice<-which(indX[,1]==x[1]);res <- colMeans(indX[indice,])[-1];return(c(x,res))},indX=indX))
  
  Xb <- apply(matrix(sample,ncol=1,nrow=length(sample)),MARGIN=1,FUN=function(x,Xbi){Xbi[which(Xbi[,1]==x),-1]},Xbi=Xbi)
  Xb <- t(Xb)-Xm
  
  Xw <- X-Xm-Xb
  res <- list(Xw=Xw,Xb=Xb,Xm=Xm)
}


###--------------------------------------
### Split.variation.two.level: decomposition of the variance (within matrix and between matrix) for two level
###--------------------------------------


Split.variation.two.level <- function(X, factor1, factor2, sample){
  
  if(is.factor(sample)){
    sample = as.numeric(sample)
    warning('the vector sample was converted into a numeric vector', call. = FALSE)
  }
  
  ######## off set term
  Xmi <- colMeans(X)
  Xm <- matrix(Xmi,nrow=nrow(X),ncol=ncol(X),byrow=T)
  #######
  
  ###### compute Xb and Xs ######## Compute the between subject variation 
  indX <- cbind(sample,X)
  Xb <- apply(indX,MARGIN=1,FUN=function(x,indX){indice<-which(indX[,1]==x[1]);res <- colMeans(indX[indice,]);return(res[-1])},indX=indX)
  Xs <- t(Xb)
  Xb <- t(Xb)-Xm
  #################################
  
  
  xbfactor1 <- X
  for (i in levels(factor(factor1))){
    indice <- which(factor1==i)
    indXX <- indX[indice,] 
    res1 <- apply(indXX,MARGIN=1,FUN=function(x,indXX){indice<-which(indXX[,1]==x[1]);if(length(indice)==1){res <- colMeans(matrix(indXX[indice,],nrow=1,ncol=dim(indXX)[2]))}else{res <- colMeans(indXX[indice,])};return(res[-1])},indXX=indXX)
    xbfactor1[indice,] <- t(res1)
  }
  
  xbfactor2 <- X
  for (i in levels(factor(factor2))){
    indice <- which(factor2==i)
    indXX <- indX[indice,] 
    res1 <- apply(indXX,MARGIN=1,FUN=function(x,indXX){indice<-which(indXX[,1]==x[1]);if(length(indice)==1){res <- colMeans(matrix(indXX[indice,],nrow=1,ncol=dim(indXX)[2]))}else{res <- colMeans(indXX[indice,])};return(res[-1])},indXX=indXX)
    xbfactor2[indice,] <- t(res1)
  }
  
  ###fixed effect###
  matfactor1 <- matrix(factor1,nrow=1,ncol=length(factor1))
  XFACTOR1 <- apply(matfactor1,MARGIN=2,FUN=function(x,matfactor1){indice<-which(matfactor1==x[1]);res <- colMeans(X[indice,]);return(res)},matfactor1=matfactor1)
  
  matfactor2 <- matrix(factor2,nrow=1,ncol=length(factor2))
  XFACTOR2 <- apply(matfactor2,MARGIN=2,FUN=function(x,matfactor2){indice<-which(matfactor2==x[1]);res <- colMeans(X[indice,]);return(res)},matfactor2=matfactor2)
  #########
  
  XCS <- xbfactor1-Xs+Xm-t(XFACTOR1)
  XTS <- xbfactor2-Xs+Xm-t(XFACTOR2)
  Xw <- X-Xb-Xm-XCS-XTS
  
  res <- list(Xw=Xw,Xb=Xb,Xm=Xm,XCS=XCS,XTS=XTS)
}



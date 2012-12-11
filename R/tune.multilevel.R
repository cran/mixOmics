# Copyright (C) 2012 
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


# ------------------------------------
# tune  function
# ------------------------------------

tune.multilevel <- function(X, 
                            Y = NULL, 
                            cond = NULL,
                            sample, 
                            ncomp=1,
                            test.keepX=c(5, 10, 15), 
                            test.keepY=NULL, 
                            already.tested.X = NULL, 
                            already.tested.Y = NULL, 
                            method = NULL, 
                            dist,
                            validation = c("Mfold", "loo"),
                            folds = 10){  
  
  # check the input parameters that have not already been checked in valid.splsdalevel1 and valid.splsdalevel2
  
  # X input
  X = as.matrix(X)
  if (length(dim(X)) != 2 || !is.numeric(X)) stop("'X' must be a numeric matrix.")
  
  # Testing the cond vector
  if(is.null(cond)) stop('Vector cond is missing', call. = FALSE)
  
  # Checking sample
  if(!is.null(dim(sample))) stop('The sample  vector should indicate the repeated measurements')
  if(length(sample) != nrow(X)) stop('X and the vector sample should have the same number of subjects')
  # check that the sample numbers are repeated
  if(length(summary(as.factor(sample))) == nrow(X)) stop('Check that the vector sample reflects the repeated measurements')
  if(is.factor(sample)){
    sample = as.numeric(sample)
    warning('the vector sample was converted into a numeric vector', call. = FALSE)
  }
  #check that the sample are numbered from 1
  if(!any(names(summary(as.factor(sample))) == '1')) {
    cat('The vector sample includes the values: ', as.vector(names(summary(as.factor(sample)))), '\n')
    stop('sample vector', call. =FALSE)
  }
  
  # ncomp
  if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0) stop("invalid number of variates, 'ncomp'.")
  
  # method
  if(is.null(method)) stop('Input method missing, should be set to splsda or spls', call. = FALSE)
  
  
  #method and already.tested.Y
  if((method == 'splsda')&& (!is.null(already.tested.Y))) warning('Unecessary already.tested.Y parameter')
  
  if(method == 'splsda'){
    # apply one or two-factor analysis
    if (is.null(dim(cond))) {
      # one factor analysis
      result = tune.splsdalevel1(X = X,cond = cond, sample, ncomp=ncomp,test.keepX=test.keepX, dist=dist, already.tested.X = already.tested.X,
                                   validation, folds)
    }else{
      # two-factor analysis
      result = tune.splsdalevel2(X = X,cond = cond,sample, ncomp=ncomp,test.keepX=test.keepX, already.tested.X = already.tested.X)
    }
  }else{
    # apply spls multlevel
    result = tune.splslevel(X = X,Y = Y,cond = cond, sample = sample, ncomp=ncomp, 
                            test.keepX=test.keepX, test.keepY=test.keepY, 
                            already.tested.X = already.tested.X, already.tested.Y = already.tested.Y)    
  } # end if method
  
  
  return(result)
}




# ----------------------------------------
# tune for splsda with one factor: use cross validation
# ----------------------------------------
tune.splsdalevel1 <- function(X, cond, ncomp=1,test.keepX,sample,dist= NULL, already.tested.X = NULL,
                                validation,folds){  
  
  # put a warning message explain the tuning parameter
  cat(paste('For a one-factor analysis, the tuning criterion is based on', folds, 'cross validation'), '\n')
  
  # Check input cond
  if(!is.factor(cond)){
    cond = as.factor(cond)
    warning('cond was set as a factor', call. = FALSE)
  }
  
  
  # check dist
  if (is.null(dist)) stop('Input dist is missing')
  
  # check the already.tested input
  if(length(already.tested.X) != (ncomp-1)) stop('The number of already tested parameters should be ', ncomp-1, ' since you set ncomp = ', ncomp)
  if((!is.null(already.tested.X)) && (!is.numeric(already.tested.X))) stop('Expecting a numerical value in already.tested.X', call. = FALSE)
  
  # put a warning message to make things clear
  if(!is.null(already.tested.X)) cat('Number of variables selected on the first ', ncomp -1, 'component(s) was ', already.tested.X)
  
  # compute error rate
  vect.error = vector(length  = length(test.keepX))
  names(vect.error) = paste('var', test.keepX, sep = '')
  error.sw = matrix(nrow = length(unique(sample)), ncol = length(test.keepX), data = 0) #collects the errors
  
  
  #-- M fold or loo cross validation --#
  # we need to sample the SAME inviduals across the folds
  k=0
  n = length(sample)
  
  #- define the folds  # !!! to redefine
  if(validation == "Mfold"){
    if (is.list(folds)){
      if (length(folds) < 2 | length(folds) > n)
        stop("Invalid number of folds.")
      if (length(unique(unlist(folds))) != n)
        stop("Invalid folds.")
      
      M = length(folds)
    }
    else{
      if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
        stop("Invalid number of folds.")
      else{
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n)) 
      }
    }
  } 
  else{ 
    folds = split(1:n, rep(1:n, length = n)) 
    M = n
  }
  
  # start loop on the fold cross validation
  for(j in 1:M){
    k = k+1     # to fill the error matrix
    
    #define training and testing sets
    # !!! change here
    omit = folds[[j]]
    X.train = X[-c(omit), ]
    cond.train = cond[-c(omit)]
    X.test = matrix(X[omit, ], nrow = length(omit))
    cond.test =   cond[omit]  #matrix(cond[omit], nrow = length(omit))
    
    #renormalise the test set
    xwtest <- X.test - matrix(colMeans(X.test),ncol=ncol(X.test),nrow=length(omit),byrow=TRUE)
    
    #decompose matrix
    samplew <- sample[-omit]
    X.train.decomp <- Split.variation.one.level(X.train,cond.train,sample=samplew)
    trainxw <- X.train.decomp$Xw 
    
    # remove variables which have zero values accross most samples using the nearZeroVar function
    remove.zero = nearZeroVar(trainxw)$Position  # might be long to compute
    if(length(remove.zero) !=0){
      trainxw = trainxw[, -c(remove.zero)]
      # also remove from the test set
      xwtest = xwtest[, -c(remove.zero)]
    }
    
    #splsda prediction on the within training matrix
    # set new argument near.zero.var to FALSE since the test has been performed before (more computationally efficient)
    for(i in 1: length(test.keepX)){
      if(ncomp ==1){
        result.sw <- splsda(trainxw, cond.train ,keepX = test.keepX[i], ncomp = ncomp) #,near.zero.var = FALSE) 
      }
      else{
        result.sw <- splsda(trainxw, cond.train ,keepX = c(already.tested.X, test.keepX[i]), 
                            ncomp = ncomp) #, near.zero.var = FALSE)
      }
      test.predict.sw <- predict(result.sw, xwtest, dist = dist)
      Prediction.sw <- levels(cond)[test.predict.sw$class[[dist]][, ncomp]]
      error.sw[k, i] <- sum(as.character(cond.test)!=Prediction.sw)
    }    
  }	
  result <- apply(error.sw, 2, sum)/length(cond)
  names(result) = paste('var', test.keepX, sep = '')  
  return(list(error = result))
}


# ----------------------------------------
# tune for splsda with 2 factors: maximise the correlation between latent variables
# ----------------------------------------

tune.splsdalevel2 <- function(X,
                              cond, 
                              sample, 
                              ncomp=1,
                              test.keepX=c(5, 10, 15), 
                              already.tested.X = NULL){  
  
  # put a warning message explain the tuning parameter
  cat('For a two-factor analysis, the tuning criterion is based on the maximisation of the correlation between the components on the whole data set', '\n')
  
  # put a warning message to make things clear
  if(!is.null(already.tested.X)) cat('Number of variables selected on the first ', ncomp -1, 'component(s) was ', already.tested.X)
  
  # check the already.tested input
  if(length(already.tested.X) != (ncomp-1)) stop('The number of already tested parameters should be ', ncomp-1, ' since you set ncomp = ', ncomp)
  if((!is.null(already.tested.X)) && (!is.numeric(already.tested.X))) stop('Expecting a numerical value in already.tested.X', call. = FALSE)
  
  # cond input
  if((is.matrix(cond)) && ncol(cond) >2) stop("'cond' must be a matrix with 2 columns for the 2 factor analysis.")
  
  # check it is two factor
  if(ncol(cond) == 2){
    if(!is.factor(cond[,1])){ 
      warning('First cond response was set as a factor', call. = FALSE)
    }
    if(!is.factor(cond[,2])){ 
      warning('Second cond response was set as a factor', call. = FALSE)
    }
    cond1 = as.factor(cond[,1])
    cond2 = as.factor(cond[,2])
  }
  cond.fact = as.factor(paste(cond1, cond2, sep="."))
  
  # decompose the variance
  Xw <- Split.variation.two.level(X, cond1, cond2, sample)$Xw
  #cond <- as.factor(paste(cond1, cond2, sep="."))
  
  # compute the correlation
  cor.value = vector(length  = length(test.keepX))
  names(cor.value) = paste('var', test.keepX, sep = '')
  
  #for each number of variables selected:
  for(i in 1: length(test.keepX)){
    if(ncomp ==1){spls.train = splsda(Xw, cond.fact, ncomp = ncomp, keepX = test.keepX[i])}
    else{spls.train = splsda(Xw, cond.fact, ncomp = ncomp, keepX = c(already.tested.X, test.keepX[i]))}
    for(h in 1: ncomp){
      if(h==1) {X.deflated = Xw; cond.deflated = unmap(as.numeric(cond.fact))}
      #deflate for h >1 such that X.temp = X.temp - t %*% t(c) and Y.temp = Y.temp - u %*% t(d) for regression mode
      #here mat.t = spls.train$variates$X
      else{
        X.deflated =  X.deflated - spls.train$variates$X[,h-1] %*%  t(spls.train$mat.c[,h-1])
        cond.deflated = cond.deflated - spls.train$variates$Y[,h -1]%*% t(spls.train$mat.d[, h-1])
      }
      cor.value[i] = cor(as.matrix(X.deflated) %*%spls.train$loadings$X[,h], as.matrix(cond.deflated) %*% spls.train$loadings$Y[,h])
    } # end h
  } # end i    
  
  return(list(cor.value = cor.value)
  )
}

# ----------------------------------------
# tune for spls (unsupervised, with one factor): maximise the correlation between latent variables
# ----------------------------------------
tune.splslevel <- function(X, Y, cond  = NULL, sample = NULL, ncomp=NULL, test.keepX=rep(ncol(X), ncomp), test.keepY=rep(ncol(Y), ncomp), already.tested.X = NULL, already.tested.Y = NULL){  
  
  # general check on inputs
  
  # Y input
  Y = as.matrix(Y)
  if (length(dim(Y)) != 2 || !is.numeric(Y)) 
    stop("'Y' must be a numeric matrix.")
  
  
  
  # --- specific check for this function  
  # put a warning message explain the tuning parameter
  cat('For a multilevel spls analysis, the tuning criterion is based on the maximisation of the correlation between the components from both data sets', '\n')
  
  # put a warning message to make things clear
  if(!is.null(already.tested.X)) cat('Number of X variables selected on the first ', ncomp -1, 'component(s) was ', already.tested.X, '\n')
  if(!is.null(already.tested.Y)) cat('Number of Y variables selected on the first ', ncomp -1, 'component(s) was ', already.tested.Y, '\n')
  
  # chcekcing input parameter already.tested
  if((!is.null(already.tested.X)) && is.null(already.tested.Y)) stop('Input already.tested.Y is missing')
  if((!is.null(already.tested.Y)) && is.null(already.tested.X)) stop('Input already.tested.X is missing')  
  
  if(length(already.tested.X) != (ncomp-1)) stop('The number of already.tested.X parameters should be ', ncomp-1, ' since you set ncomp = ', ncomp)
  if(length(already.tested.Y) != (ncomp-1)) stop('The number of already.tested.Y parameters should be ', ncomp-1, ' since you set ncomp = ', ncomp)
  
  if((!is.null(already.tested.X)) && (!is.numeric(already.tested.X))) stop('Expecting a numerical value in already.tested.X', call. = FALSE)
  if((!is.null(already.tested.Y)) && (!is.numeric(already.tested.Y))) stop('Expecting a numerical value in already.tested.X', call. = FALSE)
  
  
  
  # decompose the variance in both data sets
  Xw <- Split.variation.one.level(X,Y=cond,sample)$Xw
  Yw <- Split.variation.one.level(X=Y,Y=cond,sample)$Xw
  
  # create the correlation matrix
  cor.value = matrix(nrow = length(test.keepX), ncol = length(test.keepY))
  rownames(cor.value) = paste('varX', test.keepX, sep = '')
  colnames(cor.value) = paste('varY', test.keepY, sep = '')
  
  # set the mode in spls
  mode = 'canonical'
  
  #for each number of variables selected:
  for(i in 1: length(test.keepX)){
    for(j in 1:length(test.keepY)){
      if(ncomp ==1){spls.train = spls(Xw, Yw, ncomp = ncomp, keepX = test.keepX[i], keepY = test.keepY[j], mode = mode)}
      else{spls.train = spls(Xw, Yw, ncomp = ncomp, keepX = c(already.tested.X, test.keepX[i]), keepY = c(already.tested.Y, test.keepY[j]), mode = mode)}
      for(h in 1: ncomp){
        if(h==1) {X.deflated = Xw; Y.deflated = Yw}
        #deflate for h >1 such that X.temp = X.temp - t %*% t(c) and Y.temp = Y.temp - u %*% t(e) for canonical mode
        #here mat.t = spls.train$variates$X
        else{
          X.deflated =  X.deflated - spls.train$variates$X[,h-1] %*%  t(spls.train$mat.c[,h-1])
          Y.deflated = Y.deflated - spls.train$variates$Y[,h -1]%*% t(spls.train$mat.e[, h-1])
        }
        cor.value[i,j] = cor(as.matrix(X.deflated) %*%spls.train$loadings$X[,h], as.matrix(Y.deflated) %*% spls.train$loadings$Y[,h])
      } # end h
    } # end j
  } # end i    
  
  return(list(cor.value = cor.value)
  )
}


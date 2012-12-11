# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
#  Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence in Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# Amrit Singh, University of British Columbia, Vancouver.
#
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


# --------------------------
# declare the S3 function:
# -------------------------
perf <- function(object, ...) UseMethod("perf")


#------------------------------------------------------#
#-- Includes perf for PLS, sPLS, PLS-DA and sPLS-DA --#
#------------------------------------------------------#

# ---------------------------------------------------
# perf for pls object
# ---------------------------------------------------
perf.pls <-
  function(object,
           criterion = c("all", "MSEP", "R2", "Q2"), 
           validation = c("Mfold", "loo"),
           folds = 10,
           #max.iter = 500, 
           #tol = 1e-06, 
           progressBar = TRUE,
           ...)
  {
    
    #-- validation des arguments --#
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    
    if (length(dim(X)) != 2) 
      stop("'X' must be a numeric matrix for validation.")
    
    if(object$mode == 'canonical') stop('PLS mode should be set to regression, invariant or classic')
    
    # this will be used only for loocv computation (see below)
    means.Y = attr(Y, "scaled:center")
    sigma.Y = attr(Y, "scaled:scale")
    
    validation = match.arg(validation)
    
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    if (any(criterion == "all" | criterion == "Q2") & ncomp == 1)
      stop("'ncomp' must be > 1 for Q2 criterion.")
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    
    #-- define  M fold or loo cross validation --------------------#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    #-- compute MSEP and/or R2 --#
    
    ## add Q2
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = array(0, c(n, q, ncomp+1))
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    #press = matrix(nrow = n, ncol = q)
    
    #KA: all criteria included in the computation
    if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))){
      
      press.mat = Ypred = array(0, c(n, q, ncomp))
      MSEP = R2 = matrix(0, nrow = q, ncol = ncomp)
      
      # set up dimnames
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
      dimnames(press.mat)[[2]] = colnames(Y)
      
      # in case the test set only includes one sample, it is better to advise the user to perform loocv
      stop.user = FALSE
      
      if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
      
      for (i in 1:M) {
        if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
        
        omit = folds[[i]]
        
        # see below, we stop the user if there is only one sample drawn on the test set using MFold
        if(length(omit) == 1) stop.user = TRUE
        
        # the training set is scaled
        X.train = scale(X[-omit, ], center = TRUE, scale = TRUE)
        Y.train = scale(Y[-omit, ], center = TRUE, scale = TRUE)
        
        # the test set is scaled either in the predict function directly (for X.test)
        # or below for Y.test
        X.test = matrix(X[omit, ], nrow = length(omit))
        Y.test = matrix(Y[omit, ], nrow = length(omit))
        
        
        #-- pls --#
        pls.res = pls(X = X.train, Y = Y.train, ncomp = ncomp, 
                      mode = mode, max.iter = max.iter, tol = tol)
        
        if (!is.null( pls.res$nzv$Position)) X.test = X.test[, - pls.res$nzv$Position]
        # in the predict function, X.test is already normalised w.r.t to training set X.train, so no need to do it here
        Y.hat = predict( pls.res, X.test)$predict
        
        for (h in 1:ncomp) {
          Ypred[omit, , h] = Y.hat[, , h]
          
          # compute the press and the RSS
          # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
          if(validation == 'Mfold'){
            # for Mfold, Y.test is simply scaled (seemed to be ok when there are enough samples per fold)
            press.mat[omit, , h] = (scale(Y.test) - Y.hat[, , h])^2
            RSS.indiv[omit, ,h+1] = (scale(Y.test) - Y.hat[, , h])^2
          } else{ 
            # in the case of loo we need to scale w.r.t the parameters in Y.train
            Y.test = sweep(Y.test, 2, means.Y, FUN = "+")
            Y.test = sweep(Y.test, 2, sigma.Y, FUN = "*")
            press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
            RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
          }
        } # end h
      } #end i (cross validation)
      
      # warn the user that at least test set had a length of 1
      if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
      
      # these criteria are computed across all folds.
      for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(scale(Y), Ypred[, , h], use = "pairwise")))^2
        
        # on en profite pour calculer le PRESS par composante et le Q2
        if(q>1){
          RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
          PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
          RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
          PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
      }
      
      colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
      
      if (q == 1){
        rownames(MSEP) = rownames(R2) = ""
      }
      
      if(ncomp>1){
        # compute Q2 total
        if(q>1){
          Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
        }else{ # for q == 1
          Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
          #           colnames(Q2.total) = paste('comp', 1:ncomp, sep = " ")
          #           rownames(Q2.total) = ""
        }
      }else{
        Q2.total = NA
      }
      
      names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
      
    } # end all, MSEP, R2, Q2
    
    
    if (progressBar == TRUE) cat('\n')
    
    # ----------  final outputs:
    if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
    if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    if (any(criterion %in% c("all", "Q2"))) res$Q2 = t(Q2.inside)
    if (any(criterion %in% c("all", "Q2"))) res$Q2.total = Q2.total
    
    
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
  }


# ===========================================================================================
# ---------------------------------------------------
# perf for spls object
# ---------------------------------------------------
perf.spls <-
  function(object,
           criterion = c("all","MSEP", "R2", "Q2"), 
           validation = c("Mfold", "loo"),
           folds = 10,
           #max.iter = 500, 
           #tol = 1e-06, 
           progressBar = TRUE,
           ...)
  {
    
    #-- validation des arguments --#
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    
    # tells which variables are selected in X and in Y:
    keepX = object$keepX   #####(object$loadings$X != 0) 
    keepY = object$keepY   ######(object$loadings$Y != 0)
    
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    # these attributes are saved for leave on out normalisation
    means.Y = attr(Y, "scaled:center")
    sigma.Y = attr(Y, "scaled:scale")
    
    validation = match.arg(validation)
    
    # initialize new objects:= to record feature stability
    featuresX  = featuresY =  list()
    for(k in 1:ncomp){
      featuresX[[k]] = featuresY[[k]] = NA
    }
    
    
    if (length(dim(X)) != 2) 
      stop("'X' must be a numeric matrix for validation.")
    
    if(object$mode == 'canonical') stop('sPLS mode should be set to regression, invariant or classic')
    
    if (any(criterion == "Q2") & ncomp == 1)
      stop("'ncomp' must be > 1 for Q2 criterion.")
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    
    #-- M fold or loo cross validation --#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } 
    else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    
    #-- compute criteria ------------ --#
    ## KA: add Q2
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = array(NA, c(n, q, ncomp+1))
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    
    #KA: all criteria were included.
    if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))) {
      press.mat = Ypred = array(NA, c(n, q, ncomp))
      MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
      
      # set up dimnames
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
      dimnames(press.mat)[[2]] = colnames(Y)
      
      # in case the test set only includes one sample, it is better to advise the user to
      # perform loocv
      stop.user = FALSE
      if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
      
      for (i in 1:M) {
        if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
        
        omit = folds[[i]]
        # see below, we stop the user if there is only one sample drawn on the test set using MFold
        if(length(omit) == 1) stop.user = TRUE
        
        # the training set is scaled
        X.train = scale(X[-omit, ], center = TRUE, scale = TRUE)
        Y.train = scale(Y[-omit, ], center = TRUE, scale = TRUE)
        # the test set is scaled either in the predict function directly (for X.test)
        # or below for Y.test
        X.test = matrix(X[omit, ], nrow = length(omit))
        Y.test = matrix(Y[omit, ], nrow = length(omit))
        
        
        #-- spls --#
        spls.res = spls(X.train, Y.train, ncomp, mode, max.iter, tol, keepX=keepX, keepY=keepY)     ## change
        
        # added: record selected features in each set
        for(k in 1:ncomp){
          featuresX[[k]] = c(unlist(featuresX[[k]]), select.var(spls.res, comp = k)$name.X)
          featuresY[[k]] = c(unlist(featuresY[[k]]), select.var(spls.res, comp = k)$name.Y)
        }
        
        
        if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
        Y.hat = predict(spls.res, X.test)$predict
        
        #compute prediction of Y
        for (h in 1:ncomp) {
          Ypred[omit, , h] = Y.hat[, , h]
          
          # KA: this bunch was added:
          # compute the press and the RSS
          # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
          if(validation == 'Mfold'){
            # for Mfold, Y.test is simply scaled (seemed to be ok when there are enough samples per fold)
            press.mat[omit, , h] = (scale(Y.test) - Y.hat[, , h])^2
            RSS.indiv[omit, ,h+1] = (scale(Y.test) - Y.hat[, , h])^2
          }else{ 
            # in the case of loo we need to scale w.r.t the parameters in Y.train
            Y.test = sweep(Y.test, 2, means.Y, FUN = "-")
            Y.test = sweep(Y.test, 2, sigma.Y, FUN = "/")
            press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
            RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
          }
        } # end h
      } #end i (cross validation)
      
      # KA added
      # warn the user that at least test set had a length of 1
      if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
      
      # these criteria are computed across all folds.
      for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[,,h], use = "pairwise")))^2
        #KA:  PRESS is also computed as well as Q2 inside this procedure 
        if(q>1){
          RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
          PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
          RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
          PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
      }
      
      colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside)  = colnames(Y)
      
      if (q == 1) rownames(MSEP) = rownames(R2) = ""
      
      if(ncomp>1){
        # compute Q2 total
        if(q>1){
          Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
        }else{ # for q == 1
          Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
          #           colnames(Q2.total) = paste('comp', 1:ncomp, sep = " ")
          #           rownames(Q2.total) = ""
        }
      }else{
        Q2.total = NA
      }
      
      names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
      
    } # end all, MSEP, R2
    
    if (progressBar == TRUE) cat('\n')
    
    
    # ---- extract stability of features ----- # NEW
    list.featuresX = list.featuresY =list()
    for(k in 1:ncomp){
      #remove the NA value that was added for initialisation
      remove.naX = which(is.na(featuresX[[k]]))
      remove.naY = which(is.na(featuresY[[k]]))
      # then summarise as a factor and output the percentage of appearance
      list.featuresX[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
      list.featuresY[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
    }
    
    # extract features selected from the full model ---------
    features.finalX = features.finalY =list()
    for(k in 1:ncomp){
      features.finalX[[k]] = select.var(object, comp = k)$value.X
      features.finalY[[k]] = select.var(object, comp = k)$value.Y
    }
    
    names(features.finalX)  = names(features.finalY) = names(list.featuresX) = names(list.featuresX) = paste('comp', 1:ncomp)
    
    
    # ----------  final outputs:
    if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
    if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    if (any(criterion %in% c("all", "Q2"))) res$Q2 = t(Q2.inside)
    if (any(criterion %in% c("all", "Q2"))) res$Q2.total = Q2.total
    
    
    #features
    res$features$stable.X = list.featuresX
    res$features$stable.Y = list.featuresY
    res$features$final.X = features.finalX
    res$features$final.Y = features.finalY
    
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
  }


# ---------------------------------------------------
# perf for plsda object
# ---------------------------------------------------
perf.plsda <-
  function(object,
           method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
           validation = c("Mfold", "loo"),
           folds = 10,
           #max.iter = 500, 
           #tol = 1e-06,
           progressBar = TRUE, ...) 
  {
    
    #-- validation des arguments --#
    X = object$X
    level.Y = object$names$Y  # to make sure the levels are ordered
    Y = object$ind.mat
    Y = map(Y)
    Y = as.factor(level.Y[Y])
    ncomp = object$ncomp
    n = nrow(X)
    
    tol = object$tol
    max.iter = object$max.iter
    
    method.predict = match.arg(method.predict, several.ok = TRUE)
    if (any(method.predict == "all")) nmthdd = 3 
    else nmthdd = length(method.predict)
    
    error.fun = function(x, y) {
      error.vec = sweep(x, 1, y, FUN = "-")
      error.vec = (error.vec != 0)
      error.vec = apply(error.vec, 2, sum) / length(y)
      return(error.vec)
    }
    
    #-- define the folds --#
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } 
    else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    error.mat = array(0, dim = c(ncomp, nmthdd, M))
    
    # in case the test set only includes one sample, it is better to advise the user to
    # perform loocv
    stop.user = FALSE
    # set up a progress bar
    if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
    
    for (i in 1:M) {
      if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
      
      #set up leave out samples.
      omit = folds[[i]]
      
      # see below, we stop the user if there is only one sample drawn on the test set using MFold
      if(length(omit) == 1) stop.user = TRUE
      
      # the training set is scaled
      X.train = scale(X[-omit, ], center = TRUE, scale = TRUE)
      ##Y.train = scale(Y[-omit], center = TRUE, scale = TRUE)
      
      Y.train = Y[-omit]      
      X.test = matrix(X[omit, ], nrow = length(omit))
      
      
      plsda.res = plsda(X = X.train, Y = Y.train, ncomp = ncomp, 
                        max.iter = max.iter, tol = tol)
      
      if (!is.null(plsda.res$nzv$Position)) X.test = X.test[, -plsda.res$nzv$Position]
      Y.predict = predict(plsda.res, X.test, method = method.predict)$class
      error.mat[, , i] = sapply(Y.predict, error.fun, y = as.numeric(Y[omit]))
    }
    
    # warn the user that at least test set had a length of 1
    if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
    
    if (progressBar == TRUE) cat('\n')
    
    #-- compute the error --#
    error.rate = apply(error.mat, 1:2, mean)
    rownames(error.rate) = paste('ncomp', 1:ncomp, sep = " ")
    colnames(error.rate) = names(Y.predict)
    
    result = list()
    result$error.rate = error.rate
    
    method = "plsda.mthd"
    class(result) = c("perf", method)
    #updated outputs
    return(invisible(result))
  }

# ---------------------------------------------------
# perf for splsda object
# ---------------------------------------------------
perf.splsda <- function(object,                                               
                        method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),     
                        validation = c("Mfold", "loo"),                                          
                        folds = 10,                                                              
                        #max.iter = 500,                                                          
                        #tol = 1e-06,
                        progressBar = TRUE, ...)                                                        
{
  
  #-- initialising arguments --#
  X = object$X
  level.Y = object$names$Y  #to make sure the levels are ordered
  Y = object$ind.mat
  Y = map(Y)
  Y = as.factor(level.Y[Y])
  ncomp = object$ncomp
  n = nrow(X)
  keepX = object$keepX  
  
  tol = object$tol
  max.iter = object$max.iter
  
  # initialize new objects:
  features <- list()
  for(k in 1:ncomp){
    features[[k]] = NA
  }
  
  method.predict = match.arg(method.predict, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
  if (any(method.predict == "all")) nmthdd = 3 
  else nmthdd = length(method.predict)  
  
  error.fun = function(x, y) {
    error.vec = sweep(x, 1, y, FUN = "-")
    error.vec = (error.vec != 0)
    error.vec = apply(error.vec, 2, sum) / length(y)
    return(error.vec)
  }
  
  #-- define the folds --#
  if (validation == "Mfold") {
    if (is.list(folds)) {
      if (length(folds) < 2 | length(folds) > n)
        stop("Invalid number of folds.")
      if (length(unique(unlist(folds))) != n)
        stop("Invalid folds.")
      
      M = length(folds)
    }
    else {
      if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
        stop("Invalid number of folds.")
      else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n)) 
      }
    }
  } 
  else { 
    folds = split(1:n, rep(1:n, length = n)) 
    M = n
  }
  
  
  error.mat = array(0, dim = c(ncomp, nmthdd, M))
  
  # in case the test set only includes one sample, it is better to advise the user to perform loocv
  stop.user = FALSE
  # set up a progress bar
  if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
  
  for (i in 1:M) {
    if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
    
    #set up leave out samples.
    omit = folds[[i]]
    
    # see below, we stop the user if there is only one sample drawn on the test set using MFold
    if(length(omit) == 1) stop.user = TRUE
    
    # the training set is scaled
    X.train = scale(X[-omit, ], center = TRUE, scale = TRUE)
    ##Y.train = scale(Y[-omit], center = TRUE, scale = TRUE)
    
    #X.train = X[-omit, ]
    Y.train = Y[-omit]
    X.test = matrix(X[omit, ], nrow = length(omit))
    
    spls.res = splsda(X.train, Y.train, ncomp, max.iter, tol, keepX=keepX)     ## change
    # added: record selected features
    for(k in 1:ncomp){
      features[[k]] = c(unlist(features[[k]]), select.var(spls.res, comp = k)$name)
    }
    
    if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
    Y.predict = predict(spls.res, X.test, method = method.predict)$class
    error.mat[, , i] = sapply(Y.predict, error.fun, y = as.numeric(Y[omit]))
    
    
  } # end loop on i
  
  # warn the user that at least test set had a length of 1
  if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
  
  if (progressBar == TRUE) cat('\n')
  
  #-- compute the error --#
  res = apply(error.mat, 1:2, mean)
  
  rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
  colnames(res) = names(Y.predict)
  
  # ---- extract stability of features ----- # NEW
  list.features = list()
  for(k in 1:ncomp){
    #remove the NA value that was added for initialisation
    remove.na = which(is.na(features[[k]]))
    # then summarise as a factor and output the percentage of appearance
    list.features[[k]] = sort(summary(as.factor(features[[k]][-remove.na]))/M, decreasing = TRUE)
  }
  
  # extract features selected from the full model ---------
  features.final = list()
  for(k in 1:ncomp){
    features.final[[k]] = select.var(object, comp = k)$value
  }
  
  names(features.final)  = names(list.features) = paste('comp', 1:ncomp)
  
  result = list()
  result$error.rate = res
  result$features$stable = list.features
  result$features$final = features.final
  
  
  method = "plsda.mthd"
  result$meth = "splsda.mthd"
  class(result) = c("perf", method)
  #updated outputs
  return(invisible(result))
}



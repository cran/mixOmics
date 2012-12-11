# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
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


valid <-
function(object, ...) UseMethod("valid")


#------------------------------------------------------#
#-- Includes valid for PLS, sPLS, PLS-DA and sPLS-DA --#
#------------------------------------------------------#

# ---------------------------------------------------
# valid for pls object
# ---------------------------------------------------
valid.pls <-
function(object,
         criterion = c("all", "MSEP", "R2", "Q2"), 
         validation = c("Mfold", "loo"),
         folds = 10,
         max.iter = 500, 
         tol = 1e-06, ...)
{
	
    #-- validation des arguments --#
    X = object$X
    Y = object$Y
	
    if (length(dim(X)) != 2) 
        stop("'X' must be a numeric matrix for validation.")
		
    means.X = attr(X, "scaled:center")
    means.Y = attr(Y, "scaled:center")
    sigma.X = attr(X, "scaled:scale")
    sigma.Y = attr(Y, "scaled:scale")
	
    X = sweep(X, 2, sigma.X, FUN = "*")
    X = sweep(X, 2, means.X, FUN = "+")
    Y = sweep(Y, 2, sigma.Y, FUN = "*")
    Y = sweep(Y, 2, means.Y, FUN = "+")
	
    validation = match.arg(validation)
	
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
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
	
    #-- compute MSEP and/or R2 --#
    if (any(criterion %in% c("all", "MSEP", "R2"))) {
        press.mat = Ypred = array(0, c(n, q, ncomp))
        MSEP = R2 = matrix(0, nrow = q, ncol = ncomp)
		
        for (i in 1:M) {
            omit = folds[[i]]
            X.train = X[-omit, ]
            Y.train = Y[-omit, ]
            X.test = matrix(X[omit, ], nrow = length(omit))
            Y.test = matrix(Y[omit, ], nrow = length(omit))
			
            X.train = scale(X.train, center = TRUE, scale = FALSE)
            xmns = attr(X.train, "scaled:center")
			
            Y.train = scale(Y.train, center = TRUE, scale = FALSE)
            ymns = attr(Y.train, "scaled:center")
			
            X.test = scale(X.test, center = xmns, scale = FALSE)
			
            #-- pls --#
            result = pls(X = X.train, Y = Y.train, ncomp = ncomp, 
                         mode = mode, max.iter = max.iter, tol = tol)
			
            if (!is.null(result$nzv$Position)) X.test = X.test[, -result$nzv$Position]
            Y.hat = predict(result, X.test)$predict
			
            for (h in 1:ncomp) {
                Y.mat = matrix(Y.hat[, , h], nrow = dim(Y.hat)[1], ncol= dim(Y.hat)[2])
                Y.hat[, , h] = sweep(Y.mat, 2, ymns, FUN = "+")
                press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
                Ypred[omit, , h] = Y.hat[, , h]
            }
        } #end i
		
        for (h in 1:ncomp) { 
            MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
            R2[, h] = (diag(cor(Y, Ypred[, , h], use = "pairwise")))^2
        }
		
        colnames(MSEP) = colnames(R2) = paste('ncomp', c(1:ncomp), sep = " ")
        rownames(MSEP) = rownames(R2) = colnames(Y)
		
        if (q == 1) rownames(MSEP) = rownames(R2) = ""
		
        #-- valeurs sortantes --#
        if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
        if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    }
	
    #-- compute Q2 --#
    if (any(criterion %in% c("all", "Q2"))) {
        Q2 = q2.pls(X, Y, ncomp, mode, M, folds, max.iter, tol)
        Y.names = dimnames(Y)[[2]]
		
        if (is.null(Y.names)) Y.names = paste("Y", 1:q, sep = "")
		
        if (q > 1) {
            res$Q2$variables = t(Q2[, 1:q])
            res$Q2$total = Q2[, q + 1]
            rownames(res$Q2$variables) = Y.names
            colnames(res$Q2$variables) = paste('comp', 1:ncomp, sep = " ")
            names(res$Q2$total) = paste('comp', 1:ncomp, sep = " ")    
        }
        else {
            colnames(Q2) = ""
            rownames(Q2) = paste('comp', 1:ncomp, sep = " ")
            res$Q2 = t(Q2)
        }
    }
	
    method = "pls.mthd"
    class(res) = c("valid", method)
    return(invisible(res))
}


# ---------------------------------------------------
# valid for spls object
# ---------------------------------------------------
valid.spls <-
function(object,
         criterion = c("all", "MSEP", "R2", "Q2"), 
         validation = c("Mfold", "loo"),
         folds = 10,
         max.iter = 500, 
         tol = 1e-06, ...)
{
	
    #-- validation des arguments --#
    X = object$X
    Y = object$Y
    keepX = (object$loadings$X != 0)
    keepY = (object$loadings$Y != 0)
	
    if (length(dim(X)) != 2) 
        stop("'X' must be a numeric matrix for validation.")
		
    means.X = attr(X, "scaled:center")
    means.Y = attr(Y, "scaled:center")
    sigma.X = attr(X, "scaled:scale")
    sigma.Y = attr(Y, "scaled:scale")
	
    X = sweep(X, 2, sigma.X, FUN = "*")
    X = sweep(X, 2, means.X, FUN = "+")
    Y = sweep(Y, 2, sigma.Y, FUN = "*")
    Y = sweep(Y, 2, means.Y, FUN = "+")
	
    validation = match.arg(validation)
	
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
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
	
    #-- compute MSEP and/or R2 --#
    if (any(criterion %in% c("all", "MSEP", "R2"))) {
        press.mat = Ypred = array(0, c(n, q, ncomp))
        MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
		
        for (i in 1:M) {
            omit = folds[[i]]
            X.train = X[-omit, ]
            Y.train = Y[-omit, ]
            X.test = matrix(X[omit, ], nrow = length(omit))
            Y.test = matrix(Y[omit, ], nrow = length(omit))
			
            X.train = scale(X.train, center = TRUE, scale = FALSE)
            xmns = attr(X.train, "scaled:center")
			
            Y.train = scale(Y.train, center = TRUE, scale = FALSE)
            ymns = attr(Y.train, "scaled:center")
			
            X.test = scale(X.test, center = xmns, scale = FALSE)
			
            #-- spls --#
            result = spls.model(X.train, Y.train, ncomp, mode, 
                                max.iter, tol, keepX, keepY)
			
            if (!is.null(result$nzv$Position)) X.test = X.test[, -result$nzv$Position]
            Y.hat = predict(result, X.test)$predict
			
            for (h in 1:ncomp) {
                Y.mat = matrix(Y.hat[, , h], nrow = dim(Y.hat)[1], ncol= dim(Y.hat)[2])
                Y.hat[, , h] = sweep(Y.mat, 2, ymns, FUN = "+")
                press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
                Ypred[omit, , h] = Y.hat[, , h]
            }
        } #end i
		
        for (h in 1:ncomp) { 
            MSEP[keepY[, h], h] = apply(as.matrix(press.mat[, keepY[, h], h]), 2, mean, na.rm = TRUE)
            R2[keepY[, h], h] = (diag(cor(Y[, keepY[, h]], Ypred[, keepY[, h], h], use = "pairwise")))^2
        }
		
        colnames(MSEP) = colnames(R2) = paste('ncomp', c(1:ncomp), sep = " ")
        rownames(MSEP) = rownames(R2) = colnames(Y)
		
        if (q == 1) rownames(MSEP) = rownames(R2) = ""
		
        #-- valeurs sortantes --#
        if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
        if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    }
	
    #-- compute Q2 --#
    if (any(criterion %in% c("all", "Q2"))) {
        Q2 = q2.spls(X, Y, ncomp, mode, M, folds, max.iter, tol, keepX, keepY)
        Y.names = dimnames(Y)[[2]]
		
        if (is.null(Y.names)) Y.names = paste("Y", 1:q, sep = "")
		
        if (q > 1) {
            res$Q2$variables = t(Q2[, 1:q])
            res$Q2$total = Q2[, q + 1]
            rownames(res$Q2$variables) = Y.names
            colnames(res$Q2$variables) = paste('comp', 1:ncomp, sep = " ")
            names(res$Q2$total) = paste('comp', 1:ncomp, sep = " ")    
        }
        else {
            colnames(Q2) = ""
            rownames(Q2) = paste('comp', 1:ncomp, sep = " ")
            res$Q2 = t(Q2)
        }
    }
	
    method = "pls.mthd"
    class(res) = c("valid", method)
    return(invisible(res))
}


# ---------------------------------------------------
# valid for plsda object
# ---------------------------------------------------
valid.plsda <-
function(object,
         method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
         validation = c("Mfold", "loo"),
         folds = 10,
         max.iter = 500, 
         tol = 1e-06, ...)
{
	
    #-- validation des arguments --#
    X = object$X
    lev = object$names$Y
    Y = object$ind.mat
    Y = map(Y)
    Y = as.factor(lev[Y])
    ncomp = object$ncomp
    n = nrow(X)

    means.X = attr(X, "scaled:center")
    sigma.X = attr(X, "scaled:scale")
    X = sweep(X, 2, sigma.X, FUN = "*")
    X = sweep(X, 2, means.X, FUN = "+")
	
    method = match.arg(method, several.ok = TRUE)
    if (any(method == "all")) nmthd = 3 
    else nmthd = length(method)
	
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

    error.mat = array(0, dim = c(ncomp, nmthd, M))
	
    for (i in 1:M) {
        omit = folds[[i]]
        X.train = X[-omit, ]
        Y.train = Y[-omit]
        X.test = matrix(X[omit, ], nrow = length(omit))
		
        X.train = scale(X.train, center = TRUE, scale = FALSE)
        xmns = attr(X.train, "scaled:center")
		
        X.test = scale(X.test, center = xmns, scale = FALSE)
		
        result = plsda(X = X.train, Y = Y.train, ncomp = ncomp, 
                       max.iter = max.iter, tol = tol)
		
        if (!is.null(result$nzv$Position)) X.test = X.test[, -result$nzv$Position]
        Y.hat = predict(result, X.test, method = method)$class
        error.mat[, , i] = sapply(Y.hat, error.fun, y = as.numeric(Y[omit]))
    }
	
    #-- compute the error --#
    res = apply(error.mat, 1:2, mean)
    rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
    colnames(res) = names(Y.hat)
	
    method = "plsda.mthd"
    class(res) = c("valid", method)
    return(invisible(res))
}


# ---------------------------------------------------
# valid for splsda object
# ---------------------------------------------------
valid.splsda <-
function(object,
         method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
         validation = c("Mfold", "loo"),
         folds = 10,
         max.iter = 500, 
         tol = 1e-06, ...)
{
	
    #-- validation des arguments --#
    X = object$X
    lev = object$names$Y
    Y = object$ind.mat
    Y = map(Y)
    Y = as.factor(lev[Y])
    ncomp = object$ncomp
    n = nrow(X)
    keepX = (object$loadings$X != 0)
    keepY = (object$loadings$Y != 0)

    means.X = attr(X, "scaled:center")
    sigma.X = attr(X, "scaled:scale")
    X = sweep(X, 2, sigma.X, FUN = "*")
    X = sweep(X, 2, means.X, FUN = "+")
	
    method = match.arg(method, several.ok = TRUE)
    if (any(method == "all")) nmthd = 3 
    else nmthd = length(method)	
    
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

    error.mat = array(0, dim = c(ncomp, nmthd, M))
	
    for (i in 1:M) {
        omit = folds[[i]]
        X.train = X[-omit, ]
        Y.train = Y[-omit]
        X.test = matrix(X[omit, ], nrow = length(omit))
		
        X.train = scale(X.train, center = TRUE, scale = FALSE)
        xmns = attr(X.train, "scaled:center")
		
        X.test = scale(X.test, center = xmns, scale = FALSE)
		
        result = splsda.model(X.train, Y.train, ncomp, 
                              max.iter, tol, keepX, keepY)
        
        if (!is.null(result$nzv$Position)) X.test = X.test[, -result$nzv$Position]
        Y.hat = predict(result, X.test, method = method)$class
        error.mat[, , i] = sapply(Y.hat, error.fun, y = as.numeric(Y[omit]))
    }

    #-- compute the error --#
    res = apply(error.mat, 1:2, mean)
    rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
    colnames(res) = names(Y.hat)
	
    method = "plsda.mthd"
    class(res) = c("valid", method)
    return(invisible(res))
}

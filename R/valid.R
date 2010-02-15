# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
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
function(X, 
         Y, 
         ncomp = min(6, ncol(X)), 
         mode = c("regression", "invariant", "classic"),
         method = c("pls", "spls"),
         keepX = if(method == "pls") NULL else c(rep(ncol(X), ncomp)),
         keepY = if(method == "pls") NULL else c(rep(ncol(Y), ncomp)),
         validation = c("loo", "Mfold"),
         M = if(validation == "Mfold") 10 else nrow(X),
         max.iter = 500, 
         tol = 1e-06,
         na.action = c("omit", "predict"),
         predict.par = NULL)
{

    #-- validation des arguments --#
    #-- do warning for mode + other warnings --#
    if (missing(validation)) 
        stop("Choose a cross-validation method: 'Mfold' or 'loo'.")
    if (missing(method)) 
        stop("Choose a method: 'pls' or 'spls'.")
    if (missing(mode)) 
        stop("Choose a mode: 'regression', 'invariant' or 'classic'.")
    if (mode == 'canonical') 
        stop("Only 'regression', 'classic' or 'invariant' mode are allowed !")
    if ((method == 'spls') & (mode == 'invariant')) 
        stop("No 'invariant' mode with sPLS.")
    if (length(dim(X)) != 2) 
        stop("'X' must be a numeric matrix.")
	
    X = as.matrix(X)
    Y = as.matrix(Y)
	
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)

    if (!is.numeric(X) || !is.numeric(Y)) 
    stop("'X' and/or 'Y' must be a numeric matrix.")
	
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || ncomp > p)
        stop("invalid number of components, 'ncomp'.")
		
    if (any(is.na(X))) {
        if (missing(na.action)) 
            stop("Missing data in 'X'. Choose an action for dealing with NAs: 'omit' or 'predict'.")
		
        na.action = match.arg(na.action)		
        if (na.action == "omit") {
            X = na.omit(X)
			idx.na = attr(X, "na.action")
			Y = as.matrix(Y[-idx.na, ])
            n = nrow(X)
        }
		else {
        idx.na = is.na(X)
        if (is.null(predict.par$ncomp)) predict.par$ncomp = ncol(X)
        if (is.null(predict.par$max.iter)) predict.par$max.iter = 500
        if (is.null(predict.par$tol)) predict.par$tol = 1e-09
        X[idx.na] = nipals(X, ncomp = predict.par$ncomp, reconst = TRUE,  
                    max.iter = predict.par$max.iter, tol = predict.par$tol)$rec[idx.na]
        }
    }   
     
    #-- criteria --#
    press.mat = Ypred = array(0, c(n, q, ncomp))
    MSEP = R2 = matrix(0, nrow = q, ncol = ncomp) 
     
    #-- M fold or loo cross validation --#
    ##- define the folds
    if (validation == "Mfold") { 
        fold = split(sample(1:n), rep(1:M, length = n)) 
    } 
    else { 
        fold = split(1:n, rep(1:n, length = n)) 
        M = n
    }
     
    for (i in 1:M) {
        omit = fold[[i]]
        X.train = X[-omit, ]
        Y.train = Y[-omit, ]
        X.test = matrix(X[omit, ], nrow = length(omit))
		Y.test = matrix(Y[omit, ], nrow = length(omit))
		
        X.train = scale(X.train, center = TRUE, scale = FALSE)
        xmns = attr(X.train, "scaled:center")

        Y.train = scale(Y.train, center = TRUE, scale = FALSE)
        ymns = attr(Y.train, "scaled:center")

        X.test = scale(X.test, center = xmns, scale = FALSE)
		
        #-- pls or spls --#
        if (method == "pls") {
            object = pls(X = X.train , Y = Y.train, ncomp = ncomp, 
                         mode = mode, max.iter = max.iter, tol = tol)
        } 
        else {
            object = spls(X = X.train , Y = Y.train, ncomp = ncomp, 
                          mode = mode, max.iter = max.iter, tol = tol, 
                          keepX = keepX, keepY = keepY)
        }
		
        Y.hat = predict(object, X.test)$predict

        for (h in 1:ncomp) {
			Y.mat = matrix(Y.hat[, , h], nrow = dim(Y.hat)[1], ncol= dim(Y.hat)[2])		
            Y.hat[, , h] = sweep(Y.mat, 2, ymns, FUN = "+")		
            press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
            Ypred[omit, , h] = Y.hat[, , h]
        }
        
    }  #end i	
	
    #-- compute MSEP and/or RMSEP --#
    for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = diag(cor(Y, Ypred[, , h], use = "pairwise"))		
    }
	
	RMSEP = sqrt(MSEP)
	
    colnames(MSEP) = colnames(RMSEP)= colnames(R2) = paste('ncomp', c(1:ncomp), sep = " ")
    rownames(MSEP) = rownames(RMSEP) = rownames(R2) = colnames(Y)        		

    #-- valeurs sortantes --#
    res = list(msep = MSEP, rmsep = RMSEP, r2 = R2)
	
    class(res) = "valid"
    return(invisible(res))
}

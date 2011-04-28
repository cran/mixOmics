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


#------------------------------------------------------#
#-- Includes valid for PLS, sPLS, PLS-DA and sPLS-DA --#
#------------------------------------------------------#

valid <-
function(X, 
         Y, 
         ncomp = min(6, ncol(X)), 
         method = c("pls", "spls", "plsda", "splsda"),
         mode = c("regression", "invariant", "classic"),
         pred.method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
         criterion = c("all", "MSEP", "R2", "Q2"),
         keepX = NULL, keepY = NULL, 
         validation = c("Mfold", "loo"),
         M = 10,
         max.iter = 500, 
         tol = 1e-06, ...)
{

    method = match.arg(method)
		
    #--------------------- PLS and sPLS ---------------------#
    if (any(c("pls", "spls") == method)) {
    	
        #-- validation des arguments --#
        #-- do warning for mode + other warnings --#
        if (length(dim(X)) != 2) 
            stop("'X' must be a numeric matrix.")
			
        mode = match.arg(mode)			
        if ((method == 'spls') & (mode == 'invariant')) 
            stop("No 'invariant' mode with sPLS.")
        if ((method == 'spls') & (mode == 'classic')) 
            stop("No 'classic' mode with sPLS.")			
			
        validation = match.arg(validation)
         
        X = as.matrix(X)
        Y = as.matrix(Y)
         	
        n = nrow(X)
        p = ncol(X)
        q = ncol(Y)
        res = list()
		 
        if (!is.numeric(X) || !is.numeric(Y)) 
            stop("'X' and/or 'Y' must be a numeric matrix.")
			
        if ((n != nrow(Y))) 
            stop("unequal number of rows in 'X' and 'Y'.")
			
        if (any(is.na(X)) || any(is.na(Y))) 
            stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
         
		nzv = nearZeroVar(X, ...)
        if (length(nzv$Position > 0)) {
            warning("Zero- or near-zero variance predictors. 
  Reset predictors matrix to not near-zero variance predictors.
  See $nzv for problematic predictors.")
            X = X[, -nzv$Position]
		    res$nzv = nzv
        }
	    p = ncol(X)     
        	
        if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || ncomp > p)
            stop("Invalid number of components, 'ncomp'.")
        ncomp = round(ncomp)
		  
        if (method == "spls") {		
            if(is.null(keepX)) keepX = c(rep(ncol(X), ncomp))
            if(is.null(keepY)) keepY = c(rep(ncol(Y), ncomp))
             
            if (length(keepX) != ncomp) 
                stop("length of 'keepX' must be equal to ", ncomp, ".")
                 
            if (length(keepY) != ncomp) 
                stop("length of 'keepY' must be equal to ", ncomp, ".")
                 
            if (any(keepX > p)) 
                stop("each component of 'keepX' must be lower or equal than ", p, ".")
                 
            if (any(keepY > q)) 
                stop("each component of 'keepY' must be lower or equal than ", q, ".")	
        }			
         
        #-- M fold or loo cross validation --#
        ##- define the folds
        if (validation == "Mfold") { 
            if (is.null(M) || !is.numeric(M) || M < 2 || M > n)
                stop("Invalid number of folds, 'M'.")
            M = round(M)
            fold = split(sample(1:n), rep(1:M, length = n)) 
        } 
        else { 
            fold = split(1:n, rep(1:n, length = n)) 
            M = n
        }
		 
        #-- compute MSEP and/or R2 --#
        if (any(criterion %in% c("all", "MSEP", "R2"))) {		
            press.mat = Ypred = array(0, c(n, q, ncomp))
            MSEP = R2 = matrix(0, nrow = q, ncol = ncomp)
         
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
            } #end i	
         	 
            for (h in 1:ncomp) { 
                MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
                R2[, h] = diag(cor(Y, Ypred[, , h], use = "pairwise"))		
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
            if (method == "pls") {
                Q2 = q2.pls(X, Y, ncomp, mode, M, fold, max.iter, tol)
            }
            else {
                Q2 = q2.spls(X, Y, ncomp, mode, keepX, keepY, M, fold, max.iter, tol)
            }
			    
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
    }

    #-------------------------- for plsda and splsda ------------------------#
    if (any(c("plsda", "splsda") == method)) {
    	
        #-- validation des arguments --#
        #-- do warning for mode + other warnings --#
        X = as.matrix(X)
		 
        if (length(dim(X)) != 2 || !is.numeric(X)) 
            stop("'X' must be a numeric matrix.")
	    	
        if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
            stop("invalid number of variates, 'ncomp'.")
			
        if(is.null(keepX) && method == "splsda") keepX = c(rep(ncol(X), ncomp))
	     	
        if (is.null(dim(Y))) {
            Y = as.factor(Y)					
        }
        else {
            stop("'Y' should be a factor or a class vector.")						
        }		
         
        pred.method = match.arg(pred.method, several.ok = TRUE)		
        if (any(pred.method == "all")) nmthd = 3 
        else nmthd = length(pred.method)
		
        cl = split(1:length(Y), Y)
        M = min(c(M, sapply(cl, length)))
        cl.fold = fold = omit = list()
     	
        error.mat = array(0, dim = c(ncomp, nmthd, M))
        error.fun = function(x, y) {
            error.vec = sweep(x, 1, y, FUN = "-")
            error.vec = (error.vec != 0)		
            error.vec = apply(error.vec, 2, sum)/length(y)
            return(error.vec)
        }
         
        for (k in 1:length(cl)) {
            fold[[k]] = split(sample(cl[[k]]), rep(1:M, length = length(cl[[k]]))) 
        }
         	
        for (k in 1:length(cl)) {
            id = sample(1:M)
            fold[[k]] = fold[[k]][id]  
        }
          
        for (k in 1:M) {
            for (i in 1:length(cl)) {
                cl.fold[[i]] = fold[[i]][[k]]
            }
            omit[[k]] = unlist(cl.fold)
        }
         
        for (i in 1:M) {	
            X.train = X[-omit[[i]], ]
            Y.train = Y[-omit[[i]]]
            X.test = matrix(X[omit[[i]], ], nrow = length(omit[[i]]))
     	    	
            X.train = scale(X.train, center = TRUE, scale = FALSE)
            xmns = attr(X.train, "scaled:center")
             
            X.test = scale(X.test, center = xmns, scale = FALSE)
	        	
            #-- plsda or splsda --#
            if (method == "plsda") {
                object = plsda(X = X.train, Y = Y.train, ncomp = ncomp, 
                               max.iter = max.iter, tol = tol)
            } 
            else {
                object = splsda(X = X.train, Y = Y.train, ncomp = ncomp, 
                                max.iter = max.iter, tol = tol, 
                                keepX = keepX)
            }
      	    	
            Y.hat = predict(object, X.test, method = pred.method)$class		
            error.mat[, , i] = sapply(Y.hat, error.fun, y = as.numeric(Y[omit[[i]]]))
        }  #end i	
    	 
        #-- compute the error --#
        res = apply(error.mat, 1:2, mean)
        rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
        colnames(res) = names(Y.hat)		
    } 
     
    method = paste(method, "mthd", sep = ".")	
    class(res) = c("valid", method)
    return(invisible(res))	
}

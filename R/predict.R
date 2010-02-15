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


predict.pls <- predict.spls <-
function(object, newdata, ...)
{

    #-- validation des arguments --#
    if (missing(newdata))
    stop("No new data available.")
     
    X = object$X
    Y = object$Y
    q = ncol(Y)
    p = ncol(X)
    mode = object$mode
     
    if (length(dim(newdata)) == 2) {
        if (ncol(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ", p, 
            " or a vector of length = ", p, ".")
    }
     
    if (length(dim(newdata)) == 0) {
        if (length(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ", p, 
            " or a vector of length = ", p, ".")
        dim(newdata) = c(1, p) 
    }
     
    if (mode == 'canonical') stop("Only 'regression', 'classic' or 'invariant' mode are allowed !")
     
    #-- initialisation des matrices --#	
    ncomp = object$ncomp
    a = object$loadings$X
    b = object$loadings$Y
    c = object$mat.c
     
    means.X = attr(X, "scaled:center")
    means.Y = attr(Y, "scaled:center")
    sigma.X = attr(X, "scaled:scale")
    sigma.Y = attr(Y, "scaled:scale")
     
    newdata = as.matrix(newdata)
    ones = matrix(rep(1, nrow(newdata)), ncol = 1)
    ##- coeff de regression 
    B.hat = array(0, dim = c(p, q, ncomp))
    ##- prediction
    Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
    ##- variates
    t.pred = array(0, dim = c(nrow(newdata), ncomp))
      
    #-- calcul de la prediction --# 
    for(h in 1:ncomp){
	    W = a[, 1:h] %*% solve(t(c[, 1:h]) %*% a[, 1:h]) 	
	    B = W %*% drop(t(b[, 1:h]))   
	    B = scale(B, center = FALSE, scale = 1 / sigma.Y)
	    B = as.matrix(scale(t(B), center = FALSE, scale = sigma.X))		
	    intercept = -scale(B, center = FALSE, scale = 1 / means.X)
	    intercept = matrix(apply(intercept, 1, sum) + means.Y, nrow = 1)
	    Y.hat[, , h] = newdata %*% t(B) + ones %*% intercept
	    t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[, h]
		B.hat[, , h] = B
    }  #end h
     
    #-- valeurs sortantes --#
    rownames(t.pred) = rownames(newdata)
    colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
    rownames(Y.hat) = rownames(newdata)
    colnames(Y.hat) = colnames(Y)
     
    return(invisible(list(predict = Y.hat, variates = t.pred, B.hat = B.hat)))
}


# -------------------------- for plsda and splsda ---------------------------------------
predict.plsda <- predict.splsda <-
function(object, newdata, method = c("class.dist", "centroids.dist", "Sr.dist", "max.dist"), ...)  
{
	#-- validation des arguments --#
    if (missing(newdata))
    stop("No new data available.")
     
    X = object$X
    Y = object$Y 
	Yprim = object$Yprim   
	q = ncol(Yprim)          
    p = ncol(X)
    mode = object$mode
	
    if (length(dim(newdata)) == 2) {
        if (ncol(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ", p, 
            " or a vector of length = ", p, ".")
    }
     
    if (length(dim(newdata)) == 0) {
        if (length(newdata) != p)
            stop("'newdata' must be a numeric matrix with ncol = ", p, 
            " or a vector of length = ", p, ".")
        dim(newdata) = c(1, p) 
    }
     
    if (mode == 'canonical') stop("Only 'regression', 'classic' or 'invariant' mode are allowed !")
     
    #-- initialisation des matrices --#	
    ncomp = object$ncomp
    a = object$loadings$X
    b = object$loadings$Y
    c = object$mat.c
     
    means.X = attr(X, "scaled:center")
    means.Y = attr(Y, "scaled:center")
    sigma.X = attr(X, "scaled:scale")
    sigma.Y = attr(Y, "scaled:scale")
     
    newdata = as.matrix(newdata)
    ones = matrix(rep(1, nrow(newdata)), ncol = 1)
    ##- coeff de regression 
    B.hat = array(0, dim = c(p, q, ncomp))
    ##- prediction
    Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
    ##- variates
    t.pred = array(0, dim = c(nrow(newdata), ncomp))
      
    #-- calcul de la prediction --# 
    for(h in 1:ncomp){
	    W = a[, 1:h] %*% solve(t(c[, 1:h]) %*% a[, 1:h]) 	
	    B = W %*% drop(t(b[, 1:h]))   
	    B = scale(B, center = FALSE, scale = 1 / sigma.Y)
	    B = as.matrix(scale(t(B), center = FALSE, scale = sigma.X))		
	    intercept = -scale(B, center = FALSE, scale = 1 / means.X)
	    intercept = matrix(apply(intercept, 1, sum) + means.Y, nrow = 1)
	    Y.hat[, , h] = newdata %*% t(B) + ones %*% intercept
	    t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[, h]
        B.hat[, , h] = B
    }  #end h
	
	G = matrix(0, nrow = q, ncol = ncomp)
	obsLevels = c(1:q)   
##	cl = list()
	cl = matrix(nrow = nrow(newdata), ncol = ncomp)
	
	for (i in 1:q) {
		if(ncomp > 1) {
			G[i, ] = apply(object$variates$X[Yprim[, i] == 1, ], 2, mean)
		}
		else {
			G[i, ] = mean(object$variates$X[Yprim[, i] == 1, ])
		}
	}	
			
	# ----    max distance -----------------
	if (method == "max.dist") {

	function.pred = function(x){
		tmp = numeric(nrow(x))
		for(j in 1:nrow(x)){
			tmp[j] = (which(x[j,] == max(x[j,]))[1])
		}
		return(tmp)
	}
	cl = apply(Y.hat, 3, function.pred)
	} #end method
	
	# ----    class distance -----------------

	if (method == "class.dist") {
	
		class.fun = function(x, q) {
			x = matrix(x, q, q, byrow = TRUE)
			diag(x) = diag(x) - 1
			d = apply(x^2, 1, sum)
			cl.id = which.min(d)
		}
	
		for (h in 1:ncomp) {
			cl.id = apply(Y.hat[, , h], 1, class.fun, q = q)
			cl[, h]  = as.factor(obsLevels[cl.id])
		}
	}	

	# ----    centroids distance -----------------

	if (method == "centroids.dist") {
	
		centroids.fun = function(x, G, h) {
			q = nrow(G)
			x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
			if (h > 1) {
				d = apply((x - G[, 1:h])^2, 1, sum)
			}
			else {
				d = (x - G[, 1])^2
			}
			cl.id = which.min(d)
		}
		
		for (h in 1:ncomp) {
			cl.id = apply(as.matrix(t.pred[, 1:h]), 1, centroids.fun, G = G, h = h)
			cl[, h]  = as.factor(obsLevels[cl.id])		
		}
	}	

	# ----    Sr distance -----------------
	
	if (method == "Sr.dist") {
	
		Sr.fun = function(x, G, Yprim, h) {
			q = nrow(G)
			Xe = Yprim %*% G[, 1:h]
			Xr = object$variates$X[, 1:h] - Xe
			Sr = t(Xr) %*% Xr / nrow(Y)
			Sr.inv = solve(Sr)
			x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
			if (h > 1) {
				mat = (x - G[, 1:h]) %*% Sr.inv %*% t(x - G[, 1:h])
				d = apply(mat^2, 1, sum)
			}
			else {
				d = drop(Sr.inv) * (x - G[, 1])^2
			}
			cl.id = which.min(d)
		}
		
		for (h in 1:ncomp) {
			cl.id = apply(as.matrix(t.pred[, 1:h]), 1, Sr.fun, G = G, Yprim = Yprim, h = h)
			cl[, h]  = as.factor(obsLevels[cl.id])		
		}
	}
	
    #-- valeurs sortantes --#
    rownames(t.pred) = rownames(newdata)
    colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
    rownames(Y.hat) = rownames(newdata)
    colnames(Y.hat) = colnames(Y)
    colnames(G) = paste("dim", c(1:ncomp), sep = " ")
    colnames(cl) = paste("dim", c(1:ncomp), sep = " ")
     
    return(invisible(list(predict = Y.hat, variates = t.pred, B.hat = B.hat, 
		                  centroids = G, method = method, class = cl)))
}




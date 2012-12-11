# Copyright (C) 2009 
# Seébastien Deéjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzàlez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research, Toulouse France and 
# The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD 
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


pls <-
function(X, 
         Y, 
         ncomp = 2, 
         mode = c("regression", "canonical", "invariant", "classic"),
         max.iter = 500, 
         tol = 1e-06,
         near.zero.var = TRUE,
         ...)
{

    #-- validation des arguments --#
    if (length(dim(X)) != 2) 
        stop("'X' must be a numeric matrix.")
     
    X = as.matrix(X)
    Y = as.matrix(Y)
     
    if (!is.numeric(X) || !is.numeric(Y)) 
        stop("'X' and/or 'Y' must be a numeric matrix.")
     
    n = nrow(X)
    q = ncol(Y)
     
    if ((n != nrow(Y))) 
        stop("unequal number of rows in 'X' and 'Y'.")
     
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
        stop("invalid number of variates, 'ncomp'.")
		
if(near.zero.var == TRUE){ 
    nzv = nearZeroVar(X, ...)
    if (length(nzv$Position > 0)) {
        warning("Zero- or near-zero variance predictors. 
  Reset predictors matrix to not near-zero variance predictors.
  See $nzv for problematic predictors.")
        X = X[, -nzv$Position]
    }
    }
	p = ncol(X)
	
    ncomp = round(ncomp)
    if(ncomp > p) {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", p, ".")
        ncomp = p
    }
     
    mode = match.arg(mode)
     
    #-- initialisation des matrices --#
    X.names = dimnames(X)[[2]]
    if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")
     
    if (dim(Y)[2] == 1) Y.names = "Y"
    else {
        Y.names = dimnames(Y)[[2]]
        if (is.null(Y.names)) Y.names = paste("Y", 1:q, sep = "")
    }
     
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names)) {
        ind.names = dimnames(Y)[[1]]
        rownames(X) = ind.names
    }
     	
    if (is.null(ind.names)) {
        ind.names = 1:n
        rownames(X) = rownames(Y) = ind.names
    }		
    	
    #-- centrer et r?duire les donn?es --#
    X = scale(X, center = TRUE, scale = TRUE)
    Y = scale(Y, center = TRUE, scale = TRUE) 
     
    X.temp = X
    Y.temp = Y
    mat.t = matrix(nrow = n, ncol = ncomp)
    mat.u = matrix(nrow = n, ncol = ncomp)
    mat.a = matrix(nrow = p, ncol = ncomp)
    mat.b = matrix(nrow = q, ncol = ncomp)
    mat.c = matrix(nrow = p, ncol = ncomp)
    mat.d = matrix(nrow = q, ncol = ncomp)
    mat.e = matrix(nrow = q, ncol = ncomp)
	n.ones = rep(1, n)
	p.ones = rep(1, p)
	q.ones = rep(1, q)
	na.X = FALSE
    na.Y = FALSE
    is.na.X = is.na(X)
    is.na.Y = is.na(Y)
	if (any(is.na.X)) na.X = TRUE
    if (any(is.na.Y)) na.Y = TRUE
     
    #-- boucle sur h --#
    for (h in 1:ncomp) {
        
        #-- initialisation --#
        u = Y.temp[, 1] 
		if (any(is.na(u))) u[is.na(u)] = 0
        a.old = 0
        b.old = 0
        iter = 1
         	
        if (na.X) {
            X.aux = X.temp
            X.aux[is.na.X] = 0
        }
         
        if (na.Y) {
            Y.aux = Y.temp
            Y.aux[is.na.Y] = 0
        }
         
        repeat {
            #--compute loading vectors and variates associated to X
            if (na.X) {
                a = crossprod(X.aux, u)				
                U = drop(u) %o% p.ones
                U[is.na.X] = 0
                u.norm = crossprod(U)				
                a = a / diag(u.norm)
                a = a / drop(sqrt(crossprod(a)))
                t = X.aux %*% a
                A = drop(a) %o% n.ones
                A[t(is.na.X)] = 0
                a.norm = crossprod(A)
                t = t / diag(a.norm)				
            }
            else {			
                a = crossprod(X.temp, u) #/ drop(crossprod(u)), not useful as a is scaled below
                a = a / drop(sqrt(crossprod(a)))
                t = X.temp %*% a / drop(crossprod(a))
            }
             
            #--compute loading vectors and variates associated to Y
            if (na.Y) {
                b = crossprod(Y.aux, t)
                T = drop(t) %o% q.ones
                T[is.na.Y] = 0
                t.norm = crossprod(T)				
                b = b / diag(t.norm)
                u = Y.aux %*% b
                B = drop(b) %o% n.ones
                B[t(is.na.Y)] = 0
                b.norm = crossprod(B)
                u = u / diag(b.norm)					
            }
            else {			
                b = crossprod(Y.temp, t) #/ drop(crossprod(t)), not useful as b is scaled below
                b=b / drop(sqrt(crossprod(b)))
                u = Y.temp %*% b / drop(crossprod(b))
            }
				
            if (crossprod(a - a.old) < tol) break
             
            if (iter == max.iter) {
                warning(paste("Maximum number of iterations reached for dimension", h),
                        call. = FALSE)
                break
            }
             
            a.old = a
            b.old = b
            iter = iter + 1
        }
         
        #-- deflation des matrices --#
        if (na.X) {
            X.aux = X.temp
            X.aux[is.na.X] = 0
            c = crossprod(X.aux, t)				
            T = drop(t) %o% p.ones
            T[is.na.X] = 0
            t.norm = crossprod(T)				
            c = c / diag(t.norm)
        }
        else {
            c = crossprod(X.temp, t) / drop(crossprod(t))
        }	
		
        X.temp = X.temp - t %*% t(c)   
         
        #-- mode canonique --#
        if (mode == "canonical") {
            if (na.Y) {
                Y.aux = Y.temp
                Y.aux[is.na.Y] = 0
                e = crossprod(Y.aux, u)
                U = drop(u) %o% q.ones
                U[is.na.Y] = 0
                u.norm = crossprod(U)				
                e = e / diag(u.norm)					
            }
            else {
                e = crossprod(Y.temp, u) / drop(crossprod(u))
            }
			
            Y.temp = Y.temp - u %*% t(e)
        }
         
        #-- mode classic --#
        if(mode == "classic") Y.temp = Y.temp - t %*% t(b)                 
         
        #-- mode regression --#
        if(mode == "regression") {
            if (na.Y) {
                Y.aux = Y.temp
                Y.aux[is.na.Y] = 0
                d = crossprod(Y.aux, t)
                T = drop(t) %o% q.ones
                T[is.na.Y] = 0
                t.norm = crossprod(T)				
                d = d / diag(t.norm)
            }
            else {				
                d = crossprod(Y.temp, t) / drop(crossprod(t))
            }
             
            Y.temp = Y.temp - t %*% t(d)
        }
		
        #-- mode invariant --#
        if (mode == "invariant") Y.temp = Y 
         
        mat.t[, h] = t
        mat.u[, h] = u
        mat.a[, h] = a
        mat.b[, h] = b
        mat.c[, h] = c
        if (mode == "regression") mat.d[, h] = d
	if (mode == "canonical") mat.e[, h] = e
     
    } #-- fin boucle sur h --#

    #-- valeurs sortantes --#
    rownames(mat.a) = rownames(mat.c) = X.names
    rownames(mat.b) = Y.names
    rownames(mat.t) = rownames(mat.u) = ind.names

    comp = paste("comp", 1:ncomp)
    colnames(mat.t) = colnames(mat.u) = comp
    colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = comp 

    cl = match.call()
    cl[[1]] = as.name('pls')
     
    result = list(call = cl,
	              X = X, 
	              Y = Y, 
	              ncomp = ncomp, 
	              mode = mode, 
	              mat.c = mat.c,
		      mat.d = mat.d,
		      mat.e = mat.e, 
	              variates = list(X = mat.t, Y = mat.u),
	              loadings = list(X = mat.a, Y = mat.b), 
	              names = list(X = X.names, Y = Y.names, indiv = ind.names),
          tol = tol,
          max.iter = max.iter,
          iter=iter
                )
    if (near.zero.var == TRUE & length(nzv$Position > 0)) result$nzv = nzv  
	
    class(result) = "pls"
    return(invisible(result))
}


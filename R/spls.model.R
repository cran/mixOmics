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


spls.model <-
function(X, Y, 
         ncomp = 2, 
         mode = c("regression", "canonical"),
         max.iter = 500, 
         tol = 1e-06,
         keepX, keepY)
{
     
    X = as.matrix(X)
    Y = as.matrix(Y)
     
    n = nrow(X)
    q = ncol(Y)
    p = ncol(X)
     
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
         
        #-- svd de M = t(X)*Y --#
        X.aux = X.temp        
        if (na.X) X.aux[is.na.X] = 0
         
        Y.aux = Y.temp         
        if (na.Y) Y.aux[is.na.Y] = 0
         
        M = crossprod(X.aux, Y.aux)
        svd.M = svd(M, nu = 1, nv = 1)
        a.old = svd.M$u
        b.old = svd.M$v
         
        #-- latent variables --#
        if (na.X) {
            t = X.aux %*% a.old
            A = drop(a.old) %o% n.ones
            A[is.na.X] = 0
            a.norm = crossprod(A)
            t = t / diag(a.norm)
            t = t / drop(sqrt(crossprod(t)))			
        }
        else {
            t = X.temp %*% a.old / drop(crossprod(a.old))
            t = t / drop(sqrt(crossprod(t)))
        }
         
        if (na.Y) {
            u = Y.aux %*% b.old
            B = drop(b.old) %o% n.ones
            B[is.na.Y] = 0
            b.norm = crossprod(B)
            u = u / diag(b.norm)
            u = u / drop(sqrt(crossprod(u)))			
        }
        else {
            u = Y.temp %*% b.old / drop(crossprod(b.old))
            u = u / drop(sqrt(crossprod(u)))
        }
         
        iter = 1
         
        #-- boucle jusqu'à convergence de a et de b --#
        repeat {
            if (na.X) a = t(X.aux) %*% u
            else a = t(X.temp) %*% u
			
            if (na.Y) b = t(Y.aux) %*% t
            else b = t(Y.temp) %*% t
             
            a[!keepX[, h]] = 0
            a = a / drop(crossprod(u))
            a = a / drop(sqrt(crossprod(a)))
		     
            b[!keepY[, h]] = 0
            b = b / drop(crossprod(t))
			 
            if (na.X) {
                t = X.aux %*% a
                A = drop(a) %o% n.ones
                A[is.na.X] = 0
                a.norm = crossprod(A)
                t = t / diag(a.norm)
                t = t / drop(sqrt(crossprod(t)))			
            }
            else {
                t = X.temp %*% a / drop(crossprod(a))
                t = t / drop(sqrt(crossprod(t)))
            }
             
            if (na.Y) {
                u = Y.aux %*% b
                B = drop(b) %o% n.ones
                B[is.na.Y] = 0
                b.norm = crossprod(B)
                u = u / diag(b.norm)
                u = u / drop(sqrt(crossprod(u)))			
            }
            else {
                u = Y.temp %*% b / drop(crossprod(b))
                u = u / drop(sqrt(crossprod(u)))
            }
           
            if (crossprod(a - a.old) < tol) break
             
            if (iter == max.iter) {
                warning(paste("Maximum number of iterations reached for the component", h),
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
         
        mat.t[, h] = t
        mat.u[, h] = u         
        mat.a[, h] = a
        mat.b[, h] = b
        mat.c[, h] = c
        if (mode == "regression") mat.d[, h] = d
         
    } #-- fin boucle sur h --#
     
    #-- valeurs sortantes --#     
    result = list(X = X, Y = Y, ncomp = ncomp, mode = mode, 
                  mat.c = mat.c, 
                  mat.t = mat.t, 
                  variates = list(X = mat.t, Y = mat.u),
                  loadings = list(X = mat.a, Y = mat.b))

    class(result) = c("spls", "pls") 
    return(invisible(result))
}

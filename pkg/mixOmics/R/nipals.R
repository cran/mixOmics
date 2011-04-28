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


nipals <-
function (X, ncomp = 1, reconst = FALSE, max.iter = 500, tol = 1e-09) 
{

    #-- validation des arguments --#
    if (length(dim(X)) != 2) 
        stop("'X' must be a numeric matrix.")
     
    X = as.matrix(X)
     
    if (!is.numeric(X)) 
        stop("'X' must be a numeric matrix.")
     
    nc = ncol(X)
    nr = nrow(X)
     
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || ncomp > nc)
        stop("invalid number of variates, 'ncomp'.")    
    ncomp = round(ncomp)
     
    #-- initialisation des matrices --#
    p = matrix(nrow = nc, ncol = ncomp)
    t = matrix(nrow = nr, ncol = ncomp)
    eig = vector("numeric", length = ncomp)
    nc.ones = rep(1, nc)
    nr.ones = rep(1, nr)
    is.na.X = is.na(X)
    na.X = FALSE
    if (any(is.na.X)) na.X = TRUE
     
    #-- boucle sur h --#	
    for (h in 1:ncomp) {
        th = X[, 1]
        if (any(is.na(th))) th[is.na(th)] = 0
        ph.old = rep(1 / sqrt(nc), nc)
        ph.new = vector("numeric", length = nc)
        iter = 1
        diff = 1
		 
        if (na.X) {
            X.aux = X
            X.aux[is.na.X] = 0
        }
         
        while (diff > tol & iter <= max.iter) {
            if (na.X) {
                ph.new = crossprod(X.aux, th)				
                T = drop(th) %o% nc.ones
                T[is.na.X] = 0
                th.cross = crossprod(T)				
                ph.new = ph.new / diag(th.cross)				
            }
            else {			
                ph.new = crossprod(X, th) / drop(crossprod(th))
            }
             
            ph.new = ph.new / drop(sqrt(crossprod(ph.new)))
    		 
            if (na.X) {
                th = X.aux %*% ph.new
                P = drop(ph.new) %o% nr.ones
                P[t(is.na.X)] = 0
                ph.cross = crossprod(P)
                th = th / diag(ph.cross)				
            }
            else {			
                th = X %*% ph.new / drop(crossprod(ph.new))
            }
             
            diff = drop(sum((ph.new - ph.old)^2, na.rm = TRUE))
            ph.old = ph.new
            iter = iter + 1
        }
         
        if (iter > max.iter) 
            warning(paste("Maximum number of iterations reached for comp.", h))
         
        X = X - th %*% t(ph.new)
        p[, h] = ph.new
        t[, h] = th
        eig[h] = sum(th * th, na.rm = TRUE)
    }
     
    eig = sqrt(eig)
    t = scale(t, center = FALSE, scale = eig)
    attr(t, "scaled:scale") = NULL
    result = list(eig = eig, p = p, t = t)
     
    if (reconst) {
        X.hat = matrix(0, nrow = nr, ncol = nc)
         
        for (h in 1:ncomp) {
            X.hat = X.hat + eig[h] * t[, h] %*% t(p[, h])
        }
         
        colnames(X.hat) = colnames(X)
        rownames(X.hat) = rownames(X)
        result$rec = X.hat
    }
     
    return(invisible(result))
}

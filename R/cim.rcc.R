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


cim.rcc <-
function(object, 
         comp = 1, 
         X.names = NULL, 
         Y.names = NULL, 
         ...) 
{

    p = ncol(object$X)
    q = ncol(object$Y)
    dim = object$ncomp	
	
    if (length(comp) == 1) {
	    if (is.null(comp) || !is.numeric(comp) || comp <= 0 || comp > dim)
            stop("invalid value for 'comp'.")
    }
	
    if (length(comp) > 1) {
        if(length(comp) > dim) 
            stop("the length of 'comp' must be smaller or equal than ", dim, ".")
        if (!is.numeric(comp) || any(comp < 1))
            stop("invalid vector for 'comp'.")
        if (any(comp > dim)) 
            stop("the elements of 'comp' must be smaller or equal than ", dim, ".")
    }
	
    if (length(X.names) != p & !is.null(X.names))
        stop("'X.names' must be a character vector of length ", p, ".")
		
    if (length(Y.names) != q & !is.null(Y.names))
        stop("'Y.names' must be a character vector of length ", q, ".")

    comp = round(comp)

    if (is.null(X.names)) X.names = object$names$X
    if (is.null(Y.names)) Y.names = object$names$Y

    bisect = object$variates$X[, comp] + object$variates$Y[, comp]
    cord.X = cor(object$X, bisect, use = "pairwise")
    cord.Y = cor(object$Y, bisect, use = "pairwise")
    simMat = as.matrix(cord.X %*% t(cord.Y))

    if (ncol(simMat) < nrow(simMat)) {
        simMat = t(simMat)
        aux.names = X.names
        X.names = Y.names
        Y.names = aux.names
    }  
	
    result = cim(simMat, labRow = X.names, labCol = Y.names, ...)
    return(invisible(result))
}


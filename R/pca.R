# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# Leigh Coonan, Student, University of Quuensland, Australia
# Fangzhou Yao, Student, University of Queensland, Australia
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


pca <-
function(X, 
         ncomp = 3,
#         retx = TRUE, 
         center = TRUE, 
         scale. = FALSE, 
         comp.tol = NULL,
         max.iter = 500, 
         tol = 1e-09) 
{

    retx  =TRUE  # to return pc's
    X = as.matrix(X)
	
    X.names = dimnames(X)[[2]]
    if (is.null(X.names)) X.names = paste("X", 1:ncol(X), sep = "")

    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names)) ind.names = 1:nrow(X)
	
    X = scale(X, center = center, scale = scale.)
    cen = attr(X, "scaled:center")
    sc = attr(X, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance.")

## added warning
    if (ncomp > min(ncol(X),nrow(X))) {
        stop("Use smaller 'ncomp'")
    }

    is.na.X = is.na(X)
    na.X = FALSE
    if (any(is.na.X)) na.X = TRUE
    NA.X = any(is.na.X)       

    if (is.null(ncomp)) {
        ncomp = min(nrow(X),ncol(X))
    }
# REMOVED 16-2-11 the estimation of the rank (too long or crashes when p >> n)
#        comp = 1
#        rank.old = 0
#        repeat {
#            s = nipals(x, ncomp = comp, reconst = TRUE, max.iter = max.iter, tol = tol)
#            rank = sum(s$eig > (s$eig[1L] * .Machine$double.eps))
#
#            if ((rank - rank.old) == 0) break
#            comp = comp + 1
#
#            if (comp > (ncol(x) - 1)) break	
#            rank.old = rank			
#        }
#		ncomp = comp
#    }
 

# If there are missing values use NIPALS agorithm
    if(na.X){
        result = nipals(X, ncomp = ncomp, reconst = retx, max.iter = max.iter, tol = tol)
        result$eig = (result$eig/sqrt(max(1, nrow(X) - 1)))^2
        result$rotation = result$p
 
  
    dimnames(result$rotation) = list(X.names, paste("PC", 1:ncol(result$rotation), sep = ""))
    dimnames(result$p) = list(X.names, paste("PC", 1:ncol(result$p), sep = ""))
    r = list(X = X, ncomp = ncomp, NA.X = NA.X, sdev = result$eig, rotation = result$rotation, center = if (is.null(cen)) FALSE else cen, 
             scale = if (is.null(sc)) FALSE else sc)

    if (retx) {
        r$x = result$rec %*% result$p
        dimnames(r$x) = list(ind.names, paste("X", 1:ncol(result$p), sep = ""))
    }


    }
# If data is complete use PCASVD, singular value decomposition
    if(!na.X){
        result = pcasvd(X, ncomp=ncomp, center=center, scale.=scale.)
        result$eig = (result$sdev^2)

# REMOVED 16-2-11 rank estimation (too long or crashes when p>>n)
#    if (!is.null(comp.tol)) {
#        rank = sum(s$eig > (s$eig[1L] * comp.tol))
#        if (rank < ncol(x)) {
#            s$p = s$p[, 1L:rank, drop = FALSE]
#            s$eig = s$eig[1L:rank]

    dimnames(result$rotation) = list(X.names, paste("PC", 1:ncol(result$rotation), sep = ""))
    r = list(X = X, ncomp = ncomp, NA.X = NA.X, sdev = (result$eig), rotation = (result$rotation), center = if (is.null(cen)) FALSE else cen, 
             scale = if (is.null(sc)) FALSE else sc)
    if (retx){
        r$x = X %*% result$rotation
        dimnames(r$x) = list(ind.names, paste("X", 1:ncol(result$rotation), sep = ""))
    }
    }

    class(r) = c("pca", "prcomp")
    return(invisible(r))
}

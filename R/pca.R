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


pca <-
function(x, 
         ncomp = NULL,
         retx = TRUE, 
         center = TRUE, 
         scale. = FALSE, 
         comp.tol = NULL,
         max.iter = 500, 
         tol = 1e-09) 
{
    x = as.matrix(x)
	
    x.names = dimnames(x)[[2]]
    if (is.null(x.names)) x.names = paste("X", 1:ncol(x), sep = "")

    ind.names = dimnames(x)[[1]]
    if (is.null(ind.names)) ind.names = 1:nrow(x)
	
    x = scale(x, center = center, scale = scale.)
    cen = attr(x, "scaled:center")
    sc = attr(x, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance.")

    if (is.null(ncomp)) {
        comp = 1
        rank.old = 0
        repeat {
            s = nipals(x, ncomp = comp, reconst = TRUE, max.iter = max.iter, tol = tol)
            rank = sum(s$eig > (s$eig[1L] * .Machine$double.eps))

            if ((rank - rank.old) == 0) break
            comp = comp + 1

            if (comp > (ncol(x) - 1)) break	
            rank.old = rank			
        }
		ncomp = comp
    }
    else {
        s = nipals(x, ncomp = ncomp, reconst = retx, max.iter = max.iter, tol = tol)
    }
	
    s$eig = s$eig/sqrt(max(1, nrow(x) - 1))
    if (!is.null(comp.tol)) {
        rank = sum(s$eig > (s$eig[1L] * comp.tol))
        if (rank < ncol(x)) {
            s$p = s$p[, 1L:rank, drop = FALSE]
            s$eig = s$eig[1L:rank]
        }
    }
	
    dimnames(s$p) = list(x.names, paste("PC", 1:ncol(s$p), sep = ""))
    r = list(ncomp = ncomp, sdev = s$eig, rotation = s$p, center = if (is.null(cen)) FALSE else cen, 
             scale = if (is.null(sc)) FALSE else sc)
			 
    if (retx) {
        r$x = s$rec %*% s$p
        dimnames(r$x) = list(ind.names, paste("X", 1:ncol(s$p), sep = ""))
    }
	
    class(r) = c("pca", "prcomp")
    return(invisible(r))
}

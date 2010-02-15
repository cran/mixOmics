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


plsda <-
function(X, 
         Y, 
         ncomp = 2, 
         max.iter = 500, 
         tol = 1e-06)
{

    #-- validation des arguments --#
    if (length(dim(X)) != 2 || !is.numeric(X)) 
        stop("'X' must be a numeric matrix.")
     
    X = as.matrix(X)
    n = nrow(X)
     
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
        stop("invalid number of variates, 'ncomp'.")

    if (is.factor(Y)) {
        obsLevels = levels(Y)
        nrow.Y = length(Y)
            if ((n != nrow.Y)) 
                stop("unequal number of rows in 'X' and 'Y'.")
        Y.old = Y
        Y = matrix(0, n, length(levels(Y.old)))
        Y[(1:n) + n * (as.vector(unclass(Y.old)) - 1)] = 1
        dimnames(Y) = list(Y.old, levels(Y.old))
    } 
    else {
        if (is.matrix(Y)) {
            test = apply(abs(Y), 1, sum)
                if (any(test != 1)) 
                    stop("the rows of 'Y' must be 0/1 and sum to 1.")
            obsLevels = colnames(Y)
            if(is.null(obsLevels)) obsLevels = paste("g", 1:ncol(Y), sep = "")
            colnames(Y) = obsLevels
            if(is.null(rownames(Y))) rownames(Y) = obsLevels[apply(Y, 1, which.max)]
        } else stop("'Y' must be a matrix or a factor.")
    }

    result = pls(X, Y, ncomp = ncomp, mode = "regression", max.iter = max.iter, 
                 tol = tol, scaleY = FALSE)
	
    class(result) = c("pls") 
    return(invisible(result))
}

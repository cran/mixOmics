# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
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
         center = TRUE, 
         scale = FALSE, 
         comp.tol = NULL,
         max.iter = 500, 
         tol = 1e-09) 
{

    retx = TRUE  # to return pc's
    X = as.matrix(X)
	
    X.names = dimnames(X)[[2]]
    if (is.null(X.names)) X.names = paste("X", 1:ncol(X), sep = "")

    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names)) ind.names = 1:nrow(X)
	
    X = scale(X, center = center, scale = scale)
    cen = attr(X, "scaled:center")
    sc = attr(X, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance.")

## added warning
    if (ncomp > min(ncol(X),nrow(X))) {
        stop("Use smaller 'ncomp'")
    }
    
    # check that the user did not enter extra arguments #
    # --------------------------------------------------#
    # what the user has entered
    match.user =names(match.call())
    # what the function is expecting
    match.function = c('X', 'ncomp', 'center', 'scale', 'comp.tol', 'max.iter', 'tol')
    
    #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
    if(length(setdiff(match.user[-1], match.function)) != 0) warning('Some of the input arguments do not match the function arguments, see ?plotVar')
    
    

    is.na.X = is.na(X)
    na.X = FALSE
    if (any(is.na.X)) na.X = TRUE
    NA.X = any(is.na.X)       

    if (is.null(ncomp)) {
        ncomp = min(nrow(X),ncol(X))
    }
    
    cl = match.call()
		cl[[1]] = as.name('pca')

# If there are missing values use NIPALS agorithm
    if(na.X){
        result = nipals(X, ncomp = ncomp, reconst = retx, max.iter = max.iter, tol = tol)
        result$eig = (result$eig/sqrt(max(1, nrow(X) - 1)))^2
        result$rotation = result$p
        dimnames(result$rotation) = list(X.names, paste("PC", 1:ncol(result$rotation), sep = ""))
        dimnames(result$p) = list(X.names, paste("PC", 1:ncol(result$p), sep = ""))
        r = list(call = cl, X = X, ncomp = ncomp, NA.X = NA.X, sdev = result$eig, rotation = result$rotation, 
		         center = if (is.null(cen)) FALSE else cen, scale = if (is.null(sc)) FALSE else sc)

        if (retx) {
            r$x = result$rec %*% result$p
            dimnames(r$x) = list(ind.names, paste("X", 1:ncol(result$p), sep = ""))
        }
    }
	
# If data is complete use PCASVD, singular value decomposition
    if(!na.X){
        result = pcasvd(X, ncomp = ncomp, center = center, scale = scale)
        result$eig = (result$sdev^2)
        dimnames(result$rotation) = list(X.names, paste("PC", 1:ncol(result$rotation), sep = ""))
        r = list(call = cl, X = X, ncomp = ncomp, NA.X = NA.X, sdev = (result$eig), 
                 rotation = (result$rotation), center = if (is.null(cen)) FALSE else cen, 
                 scale = if (is.null(sc)) FALSE else sc)
				 
        if (retx){
            r$x = X %*% result$rotation
            dimnames(r$x) = list(ind.names, paste("X", 1:ncol(result$rotation), sep = ""))
        }
    }

    class(r) = c("pca", "prcomp")
    return(invisible(r))
}

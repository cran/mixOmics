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





"rcc" <-
function(X, Y, lambda1 = 0, lambda2 = 0, ...) 
{

# validation des arguments #
#--------------------------#
if (length(dim(X)) != 2 || length(dim(Y)) != 2) 
        stop("'X' and/or 'Y' must be a numeric matrix.")

X = as.matrix(X)
Y = as.matrix(Y)

if (!is.numeric(X) || !is.numeric(Y)) 
        stop("'X' and/or 'Y' must be a numeric matrix.")

    if ((n = nrow(X)) != nrow(Y)) 
        stop("unequal number of rows in 'X' and 'Y'.")

    p = ncol(X)
    q = ncol(Y)

if(lambda1 < 0 || lambda2 < 0 || !is.numeric(lambda1) || !is.numeric(lambda2)) 
stop("invalid value for the regularization parameters.")

    X.names = dimnames(X)[[2]]
if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")

    Y.names = dimnames(Y)[[2]]
if (is.null(X.names)) Y.names = paste("Y", 1:q, sep = "")
    ind.names = dimnames(X)[[1]]

    if (is.null(ind.names)) ind.names = dimnames(Y)[[1]]
if (is.null(ind.names)) ind.names = 1:n

# calcul des matrices des covariances #
#-------------------------------------#
    Cxx = var(X, na.rm = TRUE, use = "pairwise") + 
diag(lambda1, ncol(X))
    Cyy = var(Y, na.rm = TRUE, use = "pairwise") + 
diag(lambda2, ncol(Y))
    Cxy = cov(X, Y, use = "pairwise")

# calcul des corrélations canoniques, des variables #
# canoniques et des vecteurs bisecteurs             #
#---------------------------------------------------#
Cxx.fac = chol(Cxx)
    Cyy.fac = chol(Cyy)
    Cxx.fac.inv = solve(Cxx.fac)
    Cyy.fac.inv = solve(Cyy.fac)
    mat = t(Cxx.fac.inv) %*% Cxy %*% Cyy.fac.inv

    if (p >= q) {
        result = svd(mat)
        cor = result$d
        xcoef = Cxx.fac.inv %*% result$u
        ycoef = Cyy.fac.inv %*% result$v
    }
    else {
        result = svd(t(mat))
        cor = result$d
        xcoef = Cxx.fac.inv %*% result$v
        ycoef = Cyy.fac.inv %*% result$u
    }

names(cor) = 1:length(cor)
    X.aux = scale(X, center = TRUE, scale = FALSE)
    Y.aux = scale(Y, center = TRUE, scale = FALSE)
    X.aux[is.na(X.aux)] = 0
    Y.aux[is.na(Y.aux)] = 0

    U = X.aux %*% xcoef
    V = Y.aux %*% ycoef

# valeurs sortantes #
#-------------------#
cl = match.call()
cl[[1]] = as.name('rcc')

    result = list(call = cl, X = X, Y = Y, lambda = c(lambda1, lambda2), cor = cor,
loadings = list(X = xcoef, Y = ycoef),
variates = list(X = U, Y = V), 
names = list(X = X.names, Y = Y.names, 
indiv = ind.names))

class(result) = "rcc"
return(invisible(result))
}


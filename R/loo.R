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



`loo` <-
function (X, Y, lambda1, lambda2) 
{

if (length(dim(X)) != 2 || length(dim(Y)) != 2) 
        stop("'X' and/or 'Y' must be a numeric matrix.")

X = as.matrix(X)
Y = as.matrix(Y)

if (!is.numeric(X) || !is.numeric(Y)) 
        stop("'X' and/or 'Y' must be a numeric matrix.")

p = ncol(X)
q = ncol(Y)

if ((nr = nrow(X)) != nrow(Y)) 
        stop("unequal number of rows in 'X' and 'Y'.")

if(lambda1 < 0 || lambda2 < 0 || !is.numeric(lambda1) || !is.numeric(lambda2)) 
stop("invalid value for the regularization parameters.")

    xscore = vector(mode = "numeric")
    yscore = vector(mode = "numeric")

X[is.na(X)] = 0
Y[is.na(Y)] = 0

    for (i in 1:nr) {
        rcc.v = rcc(X[-i, ], Y[-i, ], lambda1, lambda2)
        xscore[i] = crossprod(X[i, ], rcc.v$loadings$X[, 1])
        yscore[i] = crossprod(Y[i, ], rcc.v$loadings$Y[, 1])
    }

    cv.score = cor(xscore, yscore)
    return(invisible(cv.score))
}


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



`nipals` <-
function (X, ncomp, reconst = FALSE, max.iter = 500, tol = 1e-09) 
{
    if (length(dim(X)) != 2) 
        stop("'X' must be a numeric matrix.")

    X = as.matrix(X)

    if (!is.numeric(X)) 
        stop("'X' must be a numeric matrix.")

    nc = ncol(X)
    nr = nrow(X)
    if (missing(ncomp)) ncomp = nc

    p = matrix(nrow = nc, ncol = ncomp)
    t = matrix(nrow = nr, ncol = ncomp)
    eig = vector("numeric", length = ncomp)

    for (h in 1:ncomp) {
        th = X[, 1]
        ph.old = rep(1 / sqrt(nc), nc)
        ph.new = vector("numeric", length = nc)
        iter = 1
        diff = 1

        while (diff > tol & iter <= max.iter) {
            for (j in 1:nc) {
                ph.new[j] = sum(X[, j] * th, na.rm = TRUE) / 
                            sum(th * th, na.rm = TRUE)
            }

            ph.new = ph.new / sqrt(sum(ph.new * ph.new, na.rm = TRUE))

            for (i in 1:nr) {
                th[i] = sum(X[i, ] * ph.new, na.rm = TRUE) / 
                        sum(ph.new * ph.new, na.rm = TRUE)
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


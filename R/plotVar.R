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


#----------------------------------------------#
#-- Includes plotIndiv for PLS, sPLS and rCC --#
#----------------------------------------------#

plotVar <-
function(object, ...) UseMethod("plotVar")


#--------------------- PLS and sPLS ---------------------#
plotVar.pls <- plotVar.spls <- 
function(object, 
         comp = 1:2, 
         rad.in = 0.5, 
         X.label = FALSE, 
         Y.label = FALSE, 
         keep.var = FALSE, 
         pch = NULL, 
         cex = NULL, 
         col = NULL, 
         font = NULL,
	 ...) 
{

    # validation des arguments #
    #--------------------------#
    if (length(comp) != 2)
        stop("'comp' must be a numeric vector of length 2.")

    if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.")
		
    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])

    # calcul des coordonnées #
    #------------------------#
    if (isTRUE(keep.var)) {
        keep.X = apply(abs(object$loadings$X), 1, sum) > 0
        keep.Y = apply(abs(object$loadings$Y), 1, sum) > 0

        if (object$mode == "canonical") {
            cord.X = cor(object$X[, keep.X], object$variates$X[, c(comp1, comp2)], 
                     use = "pairwise")
            cord.Y = cor(object$Y[, keep.Y], object$variates$Y[, c(comp1, comp2)], 
                     use = "pairwise")
        }
        else {
            cord.X = cor(object$X[, keep.X], object$variates$X[, c(comp1, comp2)], 
                     use = "pairwise")
            cord.Y = cor(object$Y[, keep.Y], object$variates$X[, c(comp1, comp2)], 
                     use = "pairwise")
        }
    }
    else {
        if (object$mode == "canonical") {
            cord.X = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")
            cord.Y = cor(object$Y, object$variates$Y[, c(comp1, comp2)], use = "pairwise")
        }
        else {
            cord.X = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")
            cord.Y = cor(object$Y, object$variates$X[, c(comp1, comp2)], use = "pairwise")
        }
    }

    p = ncol(object$X)
    q = ncol(object$Y)

    # le plot des variables #
    #-----------------------#
    if (length(X.label) > 1 & length(X.label) != p)
        stop("'X.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

    if (length(Y.label) > 1 & length(Y.label) != q)
        stop("'Y.label' must be a vector of length 'ncol(Y)' or a boolean atomic vector.")

    if (is.null(pch)) {
        pch = list(rep(16, p), rep(17, q))
    }
    else {
        if (is.list(pch)) {
            if (length(pch[[1]]) != p || length(pch[[2]]) != q) 
                stop("'pch' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
        }
        else { 
            if (length(pch) == 2) { 
                pch = list(rep(pch[1], p), rep(pch[2], q))
            }
            else {
                stop("'pch' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
            }
        }
    }

    if (is.null(cex)) {
        cex = list(rep(1, p), rep(1, q))
    }
    else {
        if (is.list(cex)) {
            if (length(cex[[1]]) != p || length(cex[[2]]) != q) 
                stop("'cex' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
        }
        else { 
            if (length(cex) == 2) { 
                cex = list(rep(cex[1], p), rep(cex[2], q))
                }
            else {
                    stop("'cex' must be a vector of length 2 or a list of two
                         vector components of length ", p, " and ", q, " respectively.")
            }
        }
    }

    if (is.null(col)) {
        col = list(rep("red", p), rep("blue", q))
    }
    else {
        if (is.list(col)) {
            if (length(col[[1]]) != p || length(col[[2]]) != q) 
                stop("'col' must be a vector of length 2 or a list of two
                vector components of length ", p, " and ", q, " respectively.")
        }
        else { 
            if (length(col) == 2) { 
                col = list(rep(col[1], p), rep(col[2], q))
            }
            else {
                stop("'col' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
            }
        }
    }

    if (is.null(font)) {
        font = list(rep(2, p), rep(3, q))
    }
    else {
        if (is.list(font)) {
            if (length(font[[1]]) != p || length(font[[2]]) != q) 
                stop("'font' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
        }
        else { 
            if (length(font) == 2) { 
                font = list(rep(font[1], p), rep(font[2], q))
            }
            else {
                stop("'font' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
            }
        }
    }

    if (isTRUE(keep.var)) {
        pch[[1]] = pch[[1]][keep.X]
        pch[[2]] = pch[[2]][keep.Y]
        col[[1]] = col[[1]][keep.X]
        col[[2]] = col[[2]][keep.Y]
        cex[[1]] = cex[[1]][keep.X]
        cex[[2]] = cex[[2]][keep.Y]
        font[[1]] = font[[1]][keep.X]
        font[[2]] = font[[2]][keep.Y]
    }

    def.par = par(no.readonly = TRUE)

    if (isTRUE(X.label)) X.label = object$names$X
    if (isTRUE(Y.label)) Y.label = object$names$Y

    if (isTRUE(keep.var)) {
        if (length(X.label) == p) X.label = X.label[keep.X]
        if (length(Y.label) == q) Y.label = Y.label[keep.Y]
    }

    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste("Comp ", comp1), ylab = paste("Comp ", comp2))

    if (length(X.label) > 1) {
        text(cord.X[, 1], cord.X[, 2], X.label, col = col[[1]], 
             font = font[[1]], cex = cex[[1]])
    }
    else {
        points(cord.X[, 1], cord.X[, 2], pch = pch[[1]], 
               cex = cex[[1]], col = col[[1]])
    }

    if (length(Y.label) > 1) {
        text(cord.Y[, 1], cord.Y[, 2], Y.label, col = col[[2]], 
             font = font[[2]], cex = cex[[2]])
    }
    else {
        points(cord.Y[, 1], cord.Y[, 2], pch = pch[[2]], 
               cex = cex[[2]], col = col[[2]]) 
    }

    abline(v = 0, h = 0)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
  
    par(def.par)  
}
	

#-------------------------- rCC -------------------------#
plotVar.rcc <-
function(object, 
         comp = 1:2, 
         rad.in = 0.5, 
         cutoff = NULL, 
         X.label = FALSE, 
         Y.label = FALSE, 
         pch = NULL, 
         cex = NULL, 
         col = NULL, 
         font = NULL,
	 ...) 
{

    # validation des arguments #
    #--------------------------#
    if (length(comp) != 2)
        stop("'comp' must be a numeric vector of length 2.")

    if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.")

    p = ncol(object$X)
    q = ncol(object$Y)
    dim = min(p, q)

    if (any(comp > dim)) 
        stop("the elements of 'comp' must be smaller or equal than ", dim, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])

    # coordonnées des variables #
	#---------------------------#
    bisect = object$variates$X[, c(comp1, comp2)] + object$variates$Y[, c(comp1, comp2)]
    cord.X = cor(object$X, bisect, use = "pairwise")
    cord.Y = cor(object$Y, bisect, use = "pairwise")

    # choix des variables avec au moins une coordonnée supérieur au cutoff #
    #----------------------------------------------------------------------#
    if (!is.null(cutoff)) {
        gp.X = vector(mode = "numeric")
        gp.Y = vector(mode = "numeric")

        k = 1
        for(i in 1:p){  
            if(any(abs(cord.X[i, ]) > cutoff)) { 
                gp.X[k] = i
                k = k + 1
            }
        }

        k = 1
        for(i in 1:q){  
            if(any(abs(cord.Y[i, ]) > cutoff)) { 
                gp.Y[k] = i
                k = k + 1
            }
        }

        if(length(gp.X) == 0 || length(gp.Y) == 0) 
            stop("Cutoff value very high for the components ", comp1,
                 " and ", comp2, ".No variable was selected.")

        cord.X = matrix(cord.X[gp.X, ], length(gp.X), 2)
        cord.Y = matrix(cord.Y[gp.Y, ], length(gp.Y), 2)
        rad.in = cutoff
    }

    # le plot des variables #
    #-----------------------#
    if (length(X.label) > 1 & length(X.label) != p)
        stop("'X.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

    if (length(Y.label) > 1 & length(Y.label) != q)
        stop("'Y.label' must be a vector of length 'ncol(Y)' or a boolean atomic vector.")

    if (is.null(pch)) {
        pch = list(rep(16, p), rep(17, q))
    }
    else {
        if (is.list(pch)) {
            if (length(pch[[1]]) != p || length(pch[[2]]) != q) 
                stop("'pch' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
        }
        else { 
            if (length(pch) == 2) { 
                pch = list(rep(pch[1], p), rep(pch[2], q))
            }
            else {
                stop("'pch' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
            }
        }
    }

    if (is.null(cex)) {
        cex = list(rep(1, p), rep(1, q))
    }
    else {
        if (is.list(cex)) {
            if (length(cex[[1]]) != p || length(cex[[2]]) != q) 
                stop("'cex' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
        }
        else { 
            if (length(cex) == 2) { 
                cex = list(rep(cex[1], p), rep(cex[2], q))
            }
            else {
                stop("'cex' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
            }
        }
    }

    if (is.null(col)) {
        col = list(rep("red", p), rep("blue", q))
    }
    else {
        if (is.list(col)) {
            if (length(col[[1]]) != p || length(col[[2]]) != q) 
                stop("'col' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
        }
        else { 
            if (length(col) == 2) { 
                col = list(rep(col[1], p), rep(col[2], q))
            }
            else {
                stop("'col' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
            }
        }
    }

    if (is.null(font)) {
        font = list(rep(2, p), rep(3, q))
    }
    else {
        if (is.list(font)) {
            if (length(font[[1]]) != p || length(font[[2]]) != q) 
                stop("'font' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
        }
        else { 
            if (length(font) == 2) { 
                font = list(rep(font[1], p), rep(font[2], q))
            }
            else {
                stop("'font' must be a vector of length 2 or a list of two
                     vector components of length ", p, " and ", q, " respectively.")
            }
        }
    }

    if (!is.null(cutoff)) {
        pch[[1]] = pch[[1]][gp.X]
        pch[[2]] = pch[[2]][gp.Y]
        col[[1]] = col[[1]][gp.X]
        col[[2]] = col[[2]][gp.Y]
        cex[[1]] = cex[[1]][gp.X]
        cex[[2]] = cex[[2]][gp.Y]
        font[[1]] = font[[1]][gp.X]
        font[[2]] = font[[2]][gp.Y]
    }

    def.par = par(no.readonly = TRUE)

    if (isTRUE(X.label)) X.label = object$names$X
    if (isTRUE(Y.label)) Y.label = object$names$Y

    if (!is.null(cutoff)) {
        if (length(X.label) == p) X.label = X.label[gp.X]
        if (length(Y.label) == q) Y.label = Y.label[gp.Y]
    }

    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste("Comp ", comp1), ylab = paste("Comp ", comp2))

    if (length(X.label) > 1) {
        text(cord.X[, 1], cord.X[, 2], X.label, col = col[[1]], 
             font = font[[1]], cex = cex[[1]])
    }
    else {
        points(cord.X[, 1], cord.X[, 2], pch = pch[[1]], 
               cex = cex[[1]], col = col[[1]])
    }

    if (length(Y.label) > 1) {
        text(cord.Y[, 1], cord.Y[, 2], Y.label, col = col[[2]], 
             font = font[[2]], cex = cex[[2]])
    }
    else {
        points(cord.Y[, 1], cord.Y[, 2], pch = pch[[2]], 
               cex = cex[[2]], col = col[[2]]) 
    }

    abline(v = 0, h = 0)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
  
par(def.par)  
}

# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
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


#------------------------------------------------------------------#
#-- Includes plotVar for PLS, sPLS, PLS-DA, SPLS-DA, rCC and PCA --#
#------------------------------------------------------------------#

plotVar <-
function(object, ...) UseMethod("plotVar")


#--------------------- PLS, (plsda, sPLS, sPLSDA below)---------------------#
plotVar.pls <-  
function(object, 
         comp = 1:2, 
         rad.in = 0.5, 
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

    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")

		
    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#     # what the function is expecting
#     match.function = c('object', 'comp', 'rad.in', 'X.label', 'Y.label', 'pch', 'cex', 'col', 'font')
#     
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
#     if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
#     

    # calcul des coordonnees #
    #------------------------#
      if (object$mode == "canonical") {
          cord.X = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")
          cord.Y = cor(object$Y, object$variates$Y[, c(comp1, comp2)], use = "pairwise")
      }
      else {
          cord.X = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")
          cord.Y = cor(object$Y, object$variates$X[, c(comp1, comp2)], use = "pairwise")
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

    if (isTRUE(X.label)) X.label = object$names$X
    if (isTRUE(Y.label)) Y.label = object$names$Y

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
    
	return(invisible(list(coord.X = cord.X, coord.Y = cord.Y)))
}

# ----------------------plsda ---------------------------#
plotVar.plsda <-  
function(object, 
         comp = 1:2, 
         rad.in = 0.5, 
         var.label = FALSE, 
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
    
    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")

		
    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])

    # calcul des coordonnees #
    #------------------------#
    cord.X = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")

    p = ncol(object$X)

    # le plot des variables #
    #-----------------------#
    if (length(var.label) > 1 & length(var.label) != p)
        stop("'var.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

    if (is.null(pch)) {
        pch = list(rep(16, p))

    }
    else {
        if (is.list(pch)) {
            if (length(pch[[1]]) != p) 
                stop("'pch'  must be a vector of length 1 or a vector of length ", p)

        }
        else { 
            if (length(pch) == 1) { 
                pch = list(rep(pch[1], p))
            }

            else {
                stop("'pch'  must be a vector of length 1 or a vector of length ", p)

            }
        }
    }

    if (is.null(cex)) {
        cex = list(rep(1, p))

    }
    else {
        if (is.list(cex)) {
            if (length(cex[[1]]) != p) 
                stop("'cex'  must be a vector of length 1 or a vector of length ", p)

        }
        else { 
            if (length(cex) == 1) { 
                cex = list(rep(cex[1], p))

                }
            else {
                    stop("'cex'  must be a vector of length 1 or a vector of length ", p)

            }
        }
    }

    if (is.null(col)) {
        col = list(rep("red", p))

    }
    else {
        if (is.list(col)) {
            if (length(col[[1]]) != p) 
                stop("'col' must be a vector of length 1 or a vector of length ", p)
        }
        else { 
            if (length(col) == 1) { 
                col = list(rep(col[1], p))

            }
            if (length(col) == p) {
                col = list(col)
              }
            else {
                stop("'col' must be a vector of length 1 or a vector of length ", p)
            }
        }
    }

    if (is.null(font)) {
        font = list(rep(2, p))

    }
    else {
        if (is.list(font)) {
            if (length(font[[1]]) != p) 
                stop("'font' must be a vector of length 1 or a vector of length ", p)

        }
        else { 
            if (length(font) == 2) { 
                font = list(rep(font[1], p), rep(font[2], q))

            }
            else {
                stop("'font'  must be a vector of length 1 or a vector of length ", p)

            }
        }
    }

    if (isTRUE(var.label)) var.label = object$names$X

    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste("Comp ", comp1), ylab = paste("Comp ", comp2))

    if (length(var.label) > 1) {
        text(cord.X[, 1], cord.X[, 2], var.label, col = col[[1]], 
             font = font[[1]], cex = cex[[1]])
    }
    else {
        points(cord.X[, 1], cord.X[, 2], pch = pch[[1]], 
               cex = cex[[1]], col = col[[1]])
    }

    abline(v = 0, h = 0)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
  
    return(invisible(list(coord.X = cord.X)))
}


#--------------------- sPLS ---------------------#
plotVar.spls <- 
function(object, 
         comp = 1:2, 
         rad.in = 0.5, 
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

    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")

		
    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#      # what the function is expecting
#     match.function = c('object', 'comp', 'rad.in', 'X.label', 'Y.label', 'pch', 'cex', 'col', 'font')
#     
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
#     if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
#     

    # calcul des coordonnees #
    #------------------------#
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

      pch[[1]] = pch[[1]][keep.X]
      pch[[2]] = pch[[2]][keep.Y]
      col[[1]] = col[[1]][keep.X]
      col[[2]] = col[[2]][keep.Y]
      cex[[1]] = cex[[1]][keep.X]
      cex[[2]] = cex[[2]][keep.Y]
      font[[1]] = font[[1]][keep.X]
      font[[2]] = font[[2]][keep.Y]

    if (isTRUE(X.label)) X.label = object$names$X
    if (isTRUE(Y.label)) Y.label = object$names$Y

      if (length(X.label) == p) X.label = X.label[keep.X]
      if (length(Y.label) == q) Y.label = Y.label[keep.Y]

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
    
    return(invisible(list(coord.X = cord.X, coord.Y = cord.Y)))
}

	

# ---------------------- sPLSDA -------------------------#
plotVar.splsda <- 
function(object, 
         comp = 1:2, 
         rad.in = 0.5, 
         var.label = FALSE, 
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

    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")

		
    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#     # what the function is expecting
#     match.function = c('object', 'comp', 'rad.in', 'var.label', 'pch', 'cex', 'col', 'font')
#     
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
#     if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
#     

    # calcul des coordonnees #
    #------------------------#
    keep.X = apply(abs(object$loadings$X), 1, sum) > 0
    cord.X = cor(object$X[, keep.X], object$variates$X[, c(comp1, comp2)], 
                 use = "pairwise")

    p = ncol(object$X)

    # le plot des variables #
    #-----------------------#
    if (length(var.label) > 1 & length(var.label) != p)
        stop("'var.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

    if (is.null(pch)) {
        pch = list(rep(16, p))

    }
    else {
        if (is.list(pch)) {
            if (length(pch[[1]]) != p) 
                stop("'pch' must be a vector of length 1 or a vector of length ", p)
        }

        else { 
            if (length(pch) == 2) { 
                pch = list(rep(pch[1], p))

            }
            else {
                stop("'pch' must be a vector of length 1 or a vector of length ", p)

            }
        }
    }

    if (is.null(cex)) {
        cex = list(rep(1, p))
    }
    else {
        if (is.list(cex)) {
            if (length(cex[[1]]) != p) 
                stop("'cex' must be a vector of length 1 or a vector of length ", p)
        }
        else { 
            if (length(cex) == 1) { 
                cex = list(rep(cex[1], p))

                }
            else {
                    stop("'cex' must be a vector of length 1 or a vector of length ", p)

            }
        }
    }

    if (is.null(col)) {
        col = list(rep("red", p))
    }
    else {
        if (is.list(col)) {
            if (length(col[[1]]) != p) 
                stop("'col' must be a vector of length 1 or a vector of length ", p)
        }

        else { 
            if (length(col) == 1) { 
                col = list(rep(col[1], p))
            }
            
              if (length(col) == p) {
                col = list(col)
              }
            else {
                stop("'col' must be a vector of length 1 or a vector of length ", p)
            }
        }
    }

    if (is.null(font)) {
        font = list(rep(2, p))
    }
    else {
        if (is.list(font)) {
            if (length(font[[1]]) != p) 
                stop("'font' must be a vector of length 1 or a vector of length ", p)
        }
        else { 
            if (length(font) == 1) { 
                font = list(rep(font[1], p))
            }
            else {
                stop("'font' must be a vector of length 1 or a vector of length ", p)
            }
        }
    }

        pch[[1]] = pch[[1]][keep.X]
        col[[1]] = col[[1]][keep.X]
        cex[[1]] = cex[[1]][keep.X]
        font[[1]] = font[[1]][keep.X]

    if (isTRUE(var.label)) var.label = object$names$X
    if (length(var.label) == p) var.label = var.label[keep.X]

    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste("Comp ", comp1), ylab = paste("Comp ", comp2))

    if (length(var.label) > 1) {
        text(cord.X[, 1], cord.X[, 2], var.label, col = col[[1]], 
             font = font[[1]], cex = cex[[1]])
    }
    else {
        points(cord.X[, 1], cord.X[, 2], pch = pch[[1]], 
               cex = cex[[1]], col = col[[1]])
    }

    abline(v = 0, h = 0)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
  
    return(invisible(list(coord.X = cord.X)))
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

    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")
    

    p = ncol(object$X)
    q = ncol(object$Y)
    dim = min(p, q)

    if (any(comp > dim)) 
        stop("the elements of 'comp' must be smaller or equal than ", dim, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#     # what the function is expecting
#     match.function = c('object', 'comp', 'rad.in', 'X.label', 'Y.label', 'pch', 'cex', 'col', 'font')
#     
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
#     if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
#     

    # coordonnees des variables #
	#---------------------------#
    bisect = object$variates$X[, c(comp1, comp2)] + object$variates$Y[, c(comp1, comp2)]
    cord.X = cor(object$X, bisect, use = "pairwise")
    cord.Y = cor(object$Y, bisect, use = "pairwise")

    # choix des variables avec au moins une coordonnee superieur au cutoff #
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
   
    return(invisible(list(coord.X = cord.X, coord.Y = cord.Y)))
}


# -------------------------------- sPCA ----------------------------------------
plotVar.spca <- 
function(object, 
         comp = 1:2, 
         rad.in = 0.5, 
         var.label = FALSE, 
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

    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")
    
		
    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#     # what the function is expecting
#     match.function = c('object', 'comp', 'rad.in', 'var.label', 'pch', 'cex', 'col', 'font')
#     
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
#     if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
#     


    # calcul des coordonnees #
    #------------------------#
    keep.X = apply(abs(object$rotation), 1, sum) > 0
    cord.X = cor(object$X[, keep.X], object$x[, c(comp1, comp2)], 
                 use = "pairwise")

    p = ncol(object$X)

    # le plot des variables #
    #-----------------------#
    if (length(var.label) > 1 & length(var.label) != p)
        stop("'var.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

    if (is.null(pch)) {
        pch = list(rep(16, p))

    }
    else {
        if (is.list(pch)) {
            if (length(pch[[1]]) != p) 
                stop("'pch' must be a vector of length 1 or a vector of length ", p)
        }

        else { 
            if (length(pch) == 2) { 
                pch = list(rep(pch[1], p))

            }
            else {
                stop("'pch' must be a vector of length 1 or a vector of length ", p)

            }
        }
    }


    if (is.null(cex)) {
        cex = list(rep(1, p))
    }
    else {
        if (is.list(cex)) {
            if (length(cex[[1]]) != p) 
                stop("'cex' must be a vector of length 1 or a vector of length ", p)
        }
        else { 
            if (length(cex) == 1) { 
                cex = list(rep(cex[1], p))

                }
            else {
                    stop("'cex' must be a vector of length 1 or a vector of length ", p)

            }
        }
    }

    if (is.null(col)) {
        col = list(rep("red", p))
    }
    else {
        if (is.list(col)) {
            if (length(col[[1]]) != p) 
                stop("'col' must be a vector of length 1 or a vector of length ", p)
        }

        else { 
            if (length(col) == 1) { 
                col = list(rep(col[1], p))
            }
              if (length(col) == p) {
                col = list(col)
              }
            else {
                stop("'col' must be a vector of length 1 or a vector of length ", p)
            }
        }
    }

    if (is.null(font)) {
        font = list(rep(2, p))
    }
    else {
        if (is.list(font)) {
            if (length(font[[1]]) != p) 
                stop("'font' must be a vector of length 1 or a vector of length ", p)
        }
        else { 
            if (length(font) == 1) { 
                font = list(rep(font[1], p))
            }
            else {
                stop("'font' must be a vector of length 1 or a vector of length ", p)
            }
        }
    }

        pch[[1]] = pch[[1]][keep.X]
        col[[1]] = col[[1]][keep.X]
        cex[[1]] = cex[[1]][keep.X]
        font[[1]] = font[[1]][keep.X]

    if (isTRUE(var.label)) var.label = object$names$X
    if (length(var.label) == p) var.label = var.label[keep.X]

    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste("Comp ", comp1), ylab = paste("Comp ", comp2))

    if (length(var.label) > 1) {
        text(cord.X[, 1], cord.X[, 2], var.label, col = col[[1]], 
             font = font[[1]], cex = cex[[1]])
    }
    else {
        points(cord.X[, 1], cord.X[, 2], pch = pch[[1]], 
               cex = cex[[1]], col = col[[1]])
    }

    abline(v = 0, h = 0)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
    
    return(invisible(list(coord.X = cord.X)))
}


# ------------------------------ PCA object ------------------------------------
plotVar.pca <-
function(object, 
         comp = 1:2,
         rad.in = 0.5, 		 
         var.label = FALSE,		 
         ...) 
{

    # validation des arguments #
    #--------------------------#
    if (length(comp) != 2)
        stop("'comp' must be a numeric vector of length 2.")

    if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.")

	q = nrow(object$rotation)
	
    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
	
    if (is.logical(var.label)) {
        if (isTRUE(var.label)) var.label = rownames(object$rotation)
    }
	
	if (length(var.label) > 1) {
        if (length(var.label) != q)
            stop("'var.label' must be a character vector of length ", q, " or a boolean atomic vector.")
    }
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#     # what the function is expecting
#     match.function = c('object', 'comp', 'rad.in', 'var.label')
#     
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
# 	if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
# 	
	
    # calcul des coordonnees #
    #------------------------#
    cord.X = cor(object$X, object$x[, c(comp1, comp2)], use = "pairwise")

    # le plot des variables #
    #-----------------------#	
    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste("Comp ", comp[1]), ylab = paste("Comp ", comp[2]))

    if (length(var.label) > 1) {
        text(cord.X[, 1], cord.X[, 2], var.label, ...)
    }
    else {
        points(cord.X[, 1], cord.X[, 2], ...)
    }

    abline(v = 0, h = 0, lty = 2)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
    
    return(invisible(list(coord.X = cord.X)))
}

# ------------------------------ IPCA object ------------------------------------
plotVar.ipca <-
function(object, 
         comp = 1:2,
         rad.in = 0.5, 		 
         var.label = FALSE,		 
         ...) 
{

    # validation des arguments #
    #--------------------------#
    if (length(comp) != 2)
        stop("'comp' must be a numeric vector of length 2.")

    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")
    

	q = nrow(object$loadings)
	
    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
	
    if (is.logical(var.label)) {
        if (isTRUE(var.label)) var.label = rownames(object$loadings)
    }
	
	if (length(var.label) > 1) {
        if (length(var.label) != q)
            stop("'var.label' must be a character vector of length ", q, " or a boolean atomic vector.")
    }
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#     # what the function is expecting
#     match.function = c('object', 'comp', 'rad.in', 'var.label')
#     
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
# 	if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
# 	
	
    # calcul des coordonnees #
    #------------------------#
    cord.X = cor(object$X, object$x[, c(comp1, comp2)], use = "pairwise")

    # le plot des variables #
    #-----------------------#	
    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste("Comp ", comp[1]), ylab = paste("Comp ", comp[2]))

    if (length(var.label) > 1) {
        text(cord.X[, 1], cord.X[, 2], var.label, ...)
    }
    else {
        points(cord.X[, 1], cord.X[, 2], ...)
    }

    abline(v = 0, h = 0, lty = 2)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
    
    return(invisible(list(coord.X = cord.X)))
}


# -------------------------------- sIPCA ----------------------------------------
plotVar.sipca <- 
function(object, 
         comp = 1:2, 
         rad.in = 0.5, 
         var.label = FALSE, 
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

    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")
    
		
    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#     # what the function is expecting
#     match.function = c('object', 'comp', 'rad.in', 'var.label', 'pch', 'cex', 'col', 'font')
#     
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
#     if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
#     


    # calcul des coordonnees #
    #------------------------#
    keep.X = apply(abs(object$loadings), 1, sum) > 0
    cord.X = cor(object$X[, keep.X], object$x[, c(comp1, comp2)], 
                 use = "pairwise")

    p = ncol(object$X)

    # le plot des variables #
    #-----------------------#
    if (length(var.label) > 1 & length(var.label) != p)
        stop("'var.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

    if (is.null(pch)) {
        pch = list(rep(16, p))

    }
    else {
        if (is.list(pch)) {
            if (length(pch[[1]]) != p) 
                stop("'pch' must be a vector of length 1 or a vector of length ", p)
        }

        else { 
            if (length(pch) == 2) { 
                pch = list(rep(pch[1], p))

            }
            else {
                stop("'pch' must be a vector of length 1 or a vector of length ", p)

            }
        }
    }


    if (is.null(cex)) {
        cex = list(rep(1, p))
    }
    else {
        if (is.list(cex)) {
            if (length(cex[[1]]) != p) 
                stop("'cex' must be a vector of length 1 or a vector of length ", p)
        }
        else { 
            if (length(cex) == 1) { 
                cex = list(rep(cex[1], p))

                }
            else {
                    stop("'cex' must be a vector of length 1 or a vector of length ", p)

            }
        }
    }

    if (is.null(col)) {
        col = list(rep("red", p))
    }
    else {
        if (is.list(col)) {
            if (length(col[[1]]) != p) 
                stop("'col' must be a vector of length 1 or a vector of length ", p)
        }

        else { 
            if (length(col) == 1) { 
                col = list(rep(col[1], p))
            }
              if (length(col) == p) {
                col = list(col)
              }
            else {
                stop("'col' must be a vector of length 1 or a vector of length ", p)
            }
        }
    }

    if (is.null(font)) {
        font = list(rep(2, p))
    }
    else {
        if (is.list(font)) {
            if (length(font[[1]]) != p) 
                stop("'font' must be a vector of length 1 or a vector of length ", p)
        }
        else { 
            if (length(font) == 1) { 
                font = list(rep(font[1], p))
            }
            else {
                stop("'font' must be a vector of length 1 or a vector of length ", p)
            }
        }
    }

        pch[[1]] = pch[[1]][keep.X]
        col[[1]] = col[[1]][keep.X]
        cex[[1]] = cex[[1]][keep.X]
        font[[1]] = font[[1]][keep.X]

    if (isTRUE(var.label)) var.label = object$names$X
    if (length(var.label) == p) var.label = var.label[keep.X]

    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste("Comp ", comp1), ylab = paste("Comp ", comp2))

    if (length(var.label) > 1) {
        text(cord.X[, 1], cord.X[, 2], var.label, col = col[[1]], 
             font = font[[1]], cex = cex[[1]])
    }
    else {
        points(cord.X[, 1], cord.X[, 2], pch = pch[[1]], 
               cex = cex[[1]], col = col[[1]])
    }

    abline(v = 0, h = 0)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
    
    return(invisible(list(coord.X = cord.X)))
}


# ---------------------- RGCCA/SGCCA -------------------------#

plotVar.sgcca <- 
  function(object, 
           comp = c(1,2), 
           #sgcca specific
           block = c(1,2),
           # the comp where the variables are selected
           ncomp.select = c(1,2),
           labels = FALSE,
           pch = c(16,17), 
           cex =  c(0.5, 0.5), 
           col =  color.mixo(2),   #c('green', 'blue'),
           font = c(2,3),
           rad.in = 0.5, 
           ...) 
{
    
    ## validation des arguments
    if (length(comp) != 2)
      stop("'comp' must be a numeric vector of length 2.")
    
    if (!is.numeric(comp) || any(comp < 1))
      stop("invalid vector for 'comp'.")
    
    if (any(comp > object$ncomp[block])) 
      stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")
    
    if(any(c(length(pch),length(cex),length(font),length(cex)) > length(block)) )
      warning("Will only take into account the first ", length(block), " arguments (pch, cex, font or cex) for the plot")
    
    if(all(ncomp.select != comp))
      stop("All argument from 'ncomp.select' differ from 'comp'")
    
    #KA changed condition
    #if(any(ncomp.select > comp))
    if(any(ncomp.select > object$ncomp[comp]))
      stop("At least one argument from 'ncomp.select' is greater than the actual number of components in the sgcca model")
    
    cat('PlotVar will only display variables selected on components', ncomp.select[1], 'and', ncomp.select[2], '\n')
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#     # what the function is expecting
#     match.function = c('object', 'comp', 'block', 'ncomp.select', 'labels', 'pch', 'cex', 'col', 'font', 'rad.in')    
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
#     if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
#     
    
    
    # define if single block is true or false
    single.block = ifelse(length(block) == 1, TRUE, FALSE)
    
    #extract data
    data = object$data
    
    # extraire la matrice de design
    design = object$design
    
    # choose the mode (canonical or regression) and if regression, decide which component to compute
    # the correlation on
    reg = NULL
    if(length(block) >1){
      if(design[block[1], block[2]] + design[block[1], block[2]] == 0) warning("There is no relationship designed between these two blocks")
      
      if(design[block[1], block[2]] == design[block[2], block[1]]){ 
        mode = 'canonical'
      } else {
        mode = 'regression'
        #define where to regress on:
        # data 1 on data 2 -> reg = 1
        # data 2 on data 1 -> reg = 0  # double check with Artur
        reg = ifelse(design[block[1], block[2]] == 1, 1, 0)
      }
    } else{ # single block case
      #check that there is a relationship between this block??
      mode = 'regression'
    }
    
    
    keep = list()
    j=1
    # identify the selected variables selected on a given component ncomp.select
    if(length(ncomp.select) > 1){
      for(k in block){
        keep[[j]] = apply(abs(object$loadings[[k]][,ncomp.select]), 1, sum) > 0
        j=j+1
      }
    }else{
      for(k in block){
        keep[[j]] = abs(object$loadings[[k]][,ncomp.select])> 0
        j=j+1        
      }   
    }
    
    # -------------------------
    # compute coordinates
    # -------------------------
    # ----- canonical mode type
    # correlation between original variables and latent variables

    coord = list()
    j=1
    # canonical mode or when representing only one block
    if((mode == 'canonical') | is.null(reg)){
      for(k in block){
        coord[[j]] = cor(data[[k]][, keep[[j]]], object$variates[[k]][,comp], use = "pairwise")
        j=j+1
      }
    }else{
      # ----- regression mode type
      # correlation between original variables and frst latent variables (because of the deflation
      # made for regression mode)
      # !!!! remove??!
      if(reg == 1){  #data 1 on data 2 
        for(j in 1:length(block)){
          coord[[j]] = cor(data[[block[2]]][, keep[[j]]], object$variates[[block[1]]][,comp], use = "pairwise")
        }
      }
      if(reg == 0){        #data 2 on data 1 
        for(j in 1:length(block)){
          coord[[j]] = cor(data[[block[1]]][, keep[[j]]], object$variates[[block[2]]][,comp], use = "pairwise")
        }
      }
    }
    
    
    # -------------------
    # input parameters
    # -------------------
    # labels
    name.labels = list()
    #if(IS.TRUE(labels)){
    for(j in 1:length(block)){
      name.labels[[j]] = colnames(data[[j]][, keep[[j]]])
    }
    
    #determine number of variables in each block
    num.var = unlist(lapply(data, ncol))
    
    pch.plot = col.plot = cex.plot = font.plot = list()
    # set arguments for plot
    for(j in 1:length(block)){
      pch.plot[[j]] = rep(pch[j], num.var[j])
      col.plot[[j]] = rep(col[j], num.var[j])
      cex.plot[[j]] = rep(cex[j], num.var[j])
      font.plot[[j]] = rep(font[j], num.var[j])
    }

    # this is to display only the selected variables
    for(j in 1:length(block)){
      pch.plot[[j]] = pch.plot[[j]][keep[[j]]]
      col.plot[[j]] = col.plot[[j]][keep[[j]]]
      cex.plot[[j]] = cex.plot[[j]][keep[[j]]]
      font.plot[[j]] = font.plot[[j]][keep[[j]]]
    }
    
    
    # ---------------------------
    # plot the correlation circles
    # --------------------------
    
    
    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste('sGCCA, component 1'), ylab = paste('sGCCA, component 2'))
    
    for(j in 1:length(block)){
      if(labels){
        text(coord[[j]][, 1], coord[[j]][, 2], name.labels[[j]], col = col.plot[[j]], 
             font = font.plot[[j]], cex = cex.plot[[j]])
      }else{
        points(coord[[j]][, 1], coord[[j]][, 2], pch = pch.plot[[j]], 
               cex = cex.plot[[j]], col = col.plot[[j]])
      }
    }
    
    
    abline(v = 0, h = 0)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
    
    # output coordinates
    return(invisible(list(coord = coord)))
    
  }


# --------------
# RGCCA
# -------------



plotVar.rgcca <- 
  function(object, 
           comp = c(1,2), 
           #sgcca specific
           block = c(1,2),
           labels = FALSE,
           pch = c(16,17), 
           cex =  c(0.5, 0.5), 
           col =  color.mixo(2),       #c('green', 'blue'),
           font = c(2,3),
           rad.in = 0.5, 
           ...) 
{
    
    ## validation des arguments
    if (length(comp) != 2)
      stop("'comp' must be a numeric vector of length 2.")
    
    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")
    
    
    if (any(comp > object$ncomp[block])) 
      stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")
    
    if(any(c(length(pch),length(cex),length(font),length(cex)) > length(block)) )
      warning("Will only take into account the first ", length(block), " arguments (pch, cex, font or cex) for the plot")
    
#     # check that the user did not enter extra arguments #
#     # --------------------------------------------------#
#     # what the user has entered
#     match.user =names(match.call())
#     # what the function is expecting
#     match.function = c('object', 'comp', 'block', 'ncomp.select', 'labels', 'pch', 'cex', 'col', 'font', 'rad.in')    
#     #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
#     if(length(setdiff(match.user[-1], match.function)) > 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
#     
    
    
    # define if single block is true or false
    single.block = ifelse(length(block) == 1, TRUE, FALSE)
    
    #extract data
    data = object$data
    
    # extraire la matrice de design
    design = object$design
    
    # choose the mode (canonical or regression) and if regression, decide which component to compute
    # the correlation on
    reg = NULL
    if(single.block == FALSE){
      if(design[block[1], block[2]] + design[block[1], block[2]] == 0) stop("There is no relationship designed between these two blocks")
      
      if(design[block[1], block[2]] == design[block[2], block[1]]){ 
        mode = 'canonical'
      } else {
        mode = 'regression'
        #define where to regress on:
        # data 1 on data 2 -> reg = 1
        # data 2 on data 1 -> reg = 0  # double check with Artur
        reg = ifelse(design[block[1], block[2]] == 1, 1, 0)
      }
    } else{ # single block case
      #check that there is a relationship between this block??
      mode = 'regression'
    }
    
    
    # -------------------------
    # compute coordinates
    # -------------------------
    # ----- canonical mode type
    # correlation between original variables and latent variables
    
    coord = list()
    j=1
    # canonical mode or when representing only one block
    if((mode == 'canonical') | is.null(reg)){
      for(k in block){
        coord[[j]] = cor(data[[k]], object$variates[[k]][,comp], use = "pairwise")
        j=j+1
      }
    }else{
      # ----- regression mode type
      # correlation between original variables and frst latent variables (because of the deflation
      # made for regression mode)
      # !!!! remove??!
      if(reg == 1){  #data 1 on data 2 
        for(j in 1:length(block)){
          coord[[j]] = cor(data[[block[2]]], object$variates[[block[1]]][,comp], use = "pairwise")
        }
      }
      if(reg == 0){        #data 2 on data 1 
        for(j in 1:length(block)){
          coord[[j]] = cor(data[[block[1]]], object$variates[[block[2]]][,comp], use = "pairwise")
        }
      }
    }
    
    
    # -------------------
    # input parameters
    # -------------------
    # labels
    name.labels = list()
    #if(IS.TRUE(labels)){
    for(j in 1:length(block)){
      name.labels[[j]] = colnames(data[[j]])
    }
    
    
    #determine number of variables in each block
    num.var = unlist(lapply(data, ncol))    
    pch.plot = col.plot = cex.plot = font.plot = list()
    # set arguments for plot
    for(j in 1:length(block)){
      pch.plot[[j]] = rep(pch[j], num.var[j])
      col.plot[[j]] = rep(col[j], num.var[j])
      cex.plot[[j]] = rep(cex[j], num.var[j])
      font.plot[[j]] = rep(font[j], num.var[j])
    }
    
    # ---------------------------
    # plot the correlation circles
    # --------------------------
    
    
    par(pty = "s")
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         xlab = paste('sGCCA, component 1'), ylab = paste('sGCCA, component 2'))
    
    for(j in 1:length(block)){
      if(labels){
        text(coord[[j]][, 1], coord[[j]][, 2], name.labels[[j]], col = col.plot[[j]], 
             font = font.plot[[j]], cex = cex.plot[[j]])
      }else{
        points(coord[[j]][, 1], coord[[j]][, 2], pch = pch.plot[[j]], 
               cex = cex.plot[[j]], col = col.plot[[j]])
      }
    }
    
    
    abline(v = 0, h = 0)
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)))
    lines(rad.in * cos(seq(0, 2 * pi, l = 100)), 
          rad.in * sin(seq(0, 2 * pi, l = 100)))
    
    # output coordinates
    return(invisible(list(coord = coord)))
    
  }
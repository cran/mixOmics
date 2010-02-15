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

plotIndiv <-
function(object, ...) UseMethod("plotIndiv")


#--------------------- PLS and sPLS ---------------------#
plotIndiv.pls <- plotIndiv.spls <- plotIndiv.plsda <- plotIndiv.splsda <- 
function(object, 
         comp = 1:2, 
         ind.names = TRUE,
         rep.space = "X-variate",
         x.label = NULL, 
         y.label = NULL,  
         col = "black", 
         cex = 1, 
         pch = 1, 
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
    rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate"))

    # l'espace de représentation #
    #----------------------------#
    if (rep.space == "X-variate"){
        x = object$variates$X[, comp1]
        y = object$variates$X[, comp2]
        if (is.null(x.label)) x.label = paste("X-variate", comp1)
        if (is.null(y.label)) y.label = paste("X-variate", comp2)
    }
	
    if (rep.space == "Y-variate"){
        x = object$variates$Y[, comp1]
        y = object$variates$Y[, comp2]
        if (is.null(x.label)) x.label = paste("Y-variate", comp1)
        if (is.null(y.label)) y.label = paste("Y-variate", comp2)
    }
	
    if (rep.space == "XY-variate"){
        x = (object$variates$X[, comp1] + object$variates$Y[, comp1]) / 2
        y = (object$variates$X[, comp2] + object$variates$Y[, comp2]) / 2
        if (is.null(x.label)) x.label = paste("X-variate", comp1)
        if (is.null(y.label)) y.label = paste("Y-variate", comp2)
    }

    # le plot des individus #
    #-----------------------#
    if (length(ind.names) > 1) {
        plot(x, y, type = "n", xlab = x.label, ylab = y.label)
        text(x, y, ind.names, col = col, cex = cex, ...)
        abline(v = 0, h = 0, lty = 2)
    }
    else {
        if (isTRUE(ind.names)) {
            ind.names = object$names$indiv
            plot(x, y, type = "n", xlab = x.label, ylab = y.label)
            text(x, y, ind.names, col = col, cex = cex, ...)
            abline(v = 0, h = 0, lty = 2)
        }
        else {
            plot(x, y, xlab = x.label, ylab = y.label, 
            col = col, cex = cex, pch = pch)
            abline(v = 0, h = 0, lty = 2)
        }
    }
}


#-------------------------- rCC -------------------------#
plotIndiv.rcc <-
function (object, 
          comp = 1:2, 
          ind.names = TRUE,
          rep.space = "XY-variate",
          x.label = NULL, 
          y.label = NULL,
          col = "black", 
          cex = 1, 
          pch = 1, 
          ...) 
{

    # validation des arguments #
    #--------------------------#
    if (length(comp) != 2)
        stop("'comp' must be a numeric vector of length 2.")

    if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.")

    dim = min(ncol(object$X), ncol(object$Y))		
    if (any(comp > dim)) 
        stop("the elements of 'comp' must be smaller or equal than ", dim, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
    rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate"))

    # l'espace de représentation #
    #----------------------------#
    if (rep.space == "X-variate") {
        x = object$variates$X[, comp1]
        y = object$variates$X[, comp2]
    }

    if (rep.space == "Y-variate") {
        x = object$variates$Y[, comp1]
        y = object$variates$Y[, comp2]
    }

    if (rep.space == "XY-variate"){
        x = (object$variates$X[, comp1] + object$variates$Y[, comp1]) / 2
        y = (object$variates$X[, comp2] + object$variates$Y[, comp2]) / 2
    }

    if (is.null(x.label)) x.label = paste("Dimension ", comp1)
    if (is.null(y.label)) y.label = paste("Dimension ", comp2)

    # le plot des individus #
    #-----------------------#
    if (length(ind.names) > 1) {
        plot(x, y, type = "n", xlab = x.label, ylab = y.label)
        text(x, y, ind.names, col = col, cex = cex, ...)
        abline(v = 0, h = 0, lty = 2)
    }
    else {
        if (isTRUE(ind.names)) {
            ind.names = object$names$indiv
            plot(x, y, type = "n", xlab = x.label, ylab = y.label)
            text(x, y, ind.names, col = col, cex = cex, ...)
            abline(v = 0, h = 0, lty = 2)
        }
        else {
            plot(x, y, xlab = x.label, ylab = y.label, 
            col = col, cex = cex, pch = pch)
            abline(v = 0, h = 0, lty = 2)
        }
    }
}



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


#--------------------------------------------------------------------#
#-- Includes plotIndiv for PLS, sPLS, PLS-DA, SPLS-DA, rCC and PCA --#
#--------------------------------------------------------------------#

plotIndiv <-
function(object, ...) UseMethod("plotIndiv")


#--------------------- PLS and sPLS ---------------------#
plotIndiv.pls <- plotIndiv.spls <- plotIndiv.plsda <- plotIndiv.splsda <- 
function(object, 
         comp = 1:2, 
         ind.names = TRUE,
         rep.space = "X-variate",
         X.label = NULL, 
         Y.label = NULL,  
         col = "black", 
         cex = 1, 
         pch = 1, 
         abline.line = TRUE,
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
		
    if (is.logical(ind.names)) {
        if (isTRUE(ind.names)) ind.names = object$names$indiv
    }
	
	if (length(ind.names) > 1) {
        if (length(ind.names) != nrow(object$X))
            stop("'ind.names' must be a character vector of length ", nrow(object$X), " or a boolean atomic vector.")
    }

    comp1 = round(comp[1])
    comp2 = round(comp[2])
    rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate"))

    
    # check that the user did not enter extra arguments #
    # --------------------------------------------------#
    # what the user has entered
    match.user =names(match.call())
    # what the function is expecting
    match.function = c('object', 'comp', 'rep.space', 'X.label', 'Y.label', 'col', 'cex', 'pch')
    
    #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
    if(length(setdiff(match.user[-1], match.function)) != 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')

    
    
    # l'espace de représentation #
    #----------------------------#
    if (rep.space == "X-variate"){
        x = object$variates$X[, comp1]
        y = object$variates$X[, comp2]
        if (is.null(X.label)) X.label = paste("X-variate", comp1)
        if (is.null(Y.label)) Y.label = paste("X-variate", comp2)
    }
	
    if (rep.space == "Y-variate"){
        x = object$variates$Y[, comp1]
        y = object$variates$Y[, comp2]
        if (is.null(X.label)) X.label = paste("Y-variate", comp1)
        if (is.null(Y.label)) Y.label = paste("Y-variate", comp2)
    }
	
    if (rep.space == "XY-variate"){
        x = (object$variates$X[, comp1] + object$variates$Y[, comp1]) / 2
        y = (object$variates$X[, comp2] + object$variates$Y[, comp2]) / 2
        if (is.null(X.label)) X.label = paste("X-variate", comp1)
        if (is.null(Y.label)) Y.label = paste("Y-variate", comp2)
    }

    # le plot des individus #
    #-----------------------#
    if (length(ind.names) > 1) {
        plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
        text(x, y, ind.names, col = col, cex = cex, ...)
        if(abline.line) abline(v = 0, h = 0, lty = 2)
    }
    else {
        if (isTRUE(ind.names)) {
            plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
            text(x, y, ind.names, col = col, cex = cex, ...)
            if(abline.line) abline(v = 0, h = 0, lty = 2)
        }
        else {
            plot(x, y, xlab = X.label, ylab = Y.label, 
            col = col, cex = cex, pch = pch)
            if(abline.line) abline(v = 0, h = 0, lty = 2)
        }
    }
}


#-------------------------- rCC -------------------------#
plotIndiv.rcc <-
function (object, 
          comp = 1:2, 
          ind.names = TRUE,
          rep.space = "XY-variate",
          X.label = NULL, 
          Y.label = NULL,
          col = "black", 
          cex = 1, 
          pch = 1, 
          abline.line = TRUE,
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


    dim = min(ncol(object$X), ncol(object$Y))		
    if (any(comp > dim)) 
        stop("the elements of 'comp' must be smaller or equal than ", dim, ".")
		
    if (is.logical(ind.names)) {
        if (isTRUE(ind.names)) ind.names = object$names$indiv
    }
	
	if (length(ind.names) > 1) {
        if (length(ind.names) != nrow(object$X))
            stop("'ind.names' must be a character vector of length ", nrow(object$X), " or a boolean atomic vector.")
    }

    comp1 = round(comp[1])
    comp2 = round(comp[2])
    rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate"))

    
    # check that the user did not enter extra arguments #
    # --------------------------------------------------#
    # what the user has entered
    match.user =names(match.call())
    # what the function is expecting
    match.function = c('object', 'comp', 'rep.space', 'X.label', 'Y.label', 'col', 'cex', 'pch')
    
    #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
    if(length(setdiff(match.user[-1], match.function)) != 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
    
    
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

    if (is.null(X.label)) X.label = paste("Dimension ", comp1)
    if (is.null(Y.label)) Y.label = paste("Dimension ", comp2)

    # le plot des individus #
    #-----------------------#
    if (length(ind.names) > 1) {
        plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
        text(x, y, ind.names, col = col, cex = cex, ...)
        if(abline.line) abline(v = 0, h = 0, lty = 2)
    }
    else {
        if (isTRUE(ind.names)) {
            plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
            text(x, y, ind.names, col = col, cex = cex, ...)
            if(abline.line) abline(v = 0, h = 0, lty = 2)
        }
        else {
            plot(x, y, xlab = X.label, ylab = Y.label, 
            col = col, cex = cex, pch = pch)
            if(abline.line) abline(v = 0, h = 0, lty = 2)
        }
    }
}

#-------------------------- PCA -------------------------#
plotIndiv.pca <- plotIndiv.spca <-
function (object, 
          comp = 1:2, 
          ind.names = TRUE,
          X.label = NULL, 
          Y.label = NULL,
          abline.line = TRUE,
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
	
    if (is.logical(ind.names)) {
        if (isTRUE(ind.names)) ind.names = rownames(object$x)
    }
	
	if (length(ind.names) > 1) {
        if (length(ind.names) != nrow(object$x))
            stop("'ind.names' must be a character vector of length ", nrow(object$x), " or a boolean atomic vector.")
    }
    
    # check that the user did not enter extra arguments #
    # --------------------------------------------------#
    # what the user has entered
    match.user =names(match.call())
    # what the function is expecting
    match.function = c('object', 'comp', 'ind.names', 'X.label', 'Y.label')
    
    #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
    if(length(setdiff(match.user[-1], match.function)) != 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
	
    # l'espace de représentation #
    #----------------------------#
    x = object$x[, comp[1]]
    y = object$x[, comp[2]]

    if (is.null(X.label)) X.label = paste("Dimension ", comp1)
    if (is.null(Y.label)) Y.label = paste("Dimension ", comp2)

    # le plot des individus #
    #-----------------------#
    if (length(ind.names) > 1) {
        plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
        text(x, y, ind.names, ...)
        if(abline.line) abline(v = 0, h = 0, lty = 2)
    }
    else {
        if (isTRUE(ind.names)) {
            plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
            text(x, y, ind.names, ...)
            if(abline.line) abline(v = 0, h = 0, lty = 2)
        }
        else {
            plot(x, y, xlab = X.label, ylab = Y.label, ...)
            if(abline.line) abline(v = 0, h = 0, lty = 2)
        }
    }
}


#-------------------------- IPCA -------------------------#
plotIndiv.ipca <- plotIndiv.sipca <-
function (object, 
          comp = 1:2, 
          ind.names = TRUE,
          X.label = NULL, 
          Y.label = NULL,
          abline.line = TRUE,
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
        stop("the elements of 'comp' must be smaller than or equal to ", object$ncomp, ".")

    comp1 = round(comp[1])
    comp2 = round(comp[2])
	
    if (is.logical(ind.names)) {
        if (isTRUE(ind.names)) ind.names = rownames(object$x)
    }
	
	if (length(ind.names) > 1) {
        if (length(ind.names) != nrow(object$x))
            stop("'ind.names' must be a character vector of length ", nrow(object$x), " or a boolean atomic vector.")
    }
    
    # check that the user did not enter extra arguments #
    # --------------------------------------------------#
    # what the user has entered
    match.user =names(match.call())
    # what the function is expecting
    match.function = c('object', 'comp', 'ind.names', 'X.label', 'Y.label')
    
    #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
    if(length(setdiff(match.user[-1], match.function)) != 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
	
    # l'espace de représentation #
    #----------------------------#
    x = object$x[, comp[1]]
    y = object$x[, comp[2]]

    if (is.null(X.label)) X.label = paste("Dimension ", comp1)
    if (is.null(Y.label)) Y.label = paste("Dimension ", comp2)

    # le plot des individus #
    #-----------------------#
    if (length(ind.names) > 1) {
        plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
        text(x, y, ind.names, ...)
        if(abline.line) abline(v = 0, h = 0, lty = 2)
    }
    else {
        if (isTRUE(ind.names)) {
            plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
            text(x, y, ind.names, ...)
            if(abline.line) abline(v = 0, h = 0, lty = 2)
        }
        else {
            plot(x, y, xlab = X.label, ylab = Y.label, ...)
            if(abline.line) abline(v = 0, h = 0, lty = 2)
        }
    }
}

# ------------------ RGCCA / SGCCA --------------------------
plotIndiv.sgcca <- plotIndiv.rgcca <-
  function(
    object, 
    comp = 1:2, 
    ind.names = TRUE,
    rep.space = 1,
    X.label = NULL, 
    Y.label = NULL,  
    col = "black", 
    cex = 1, 
    pch = 1, 
    abline.line = TRUE,
    ...){
    
     # validation des arguments #
     #--------------------------#
    if (length(comp) != 2)
      stop("'comp' must be a numeric vector of length 2.")
    
        if (!is.numeric(comp))
        stop("invalid vector for 'comp'.")
        
        if (length(comp) == 1)
        stop("Need at least 2 components to plot the graph")

    
    if (any(comp > object$ncomp[rep.space])) 
       stop("the elements of 'comp' must be smaller or equal than ", object$ncomp[rep.space], ".")
     
    
    if (is.logical(ind.names)) {
      if (isTRUE(ind.names)) ind.names = object$names$indiv
    }
     
    if (length(ind.names) > 1) {
      if (length(ind.names) !=  nrow(object$variates[[rep.space]]))
        stop("'ind.names' must be a character vector of length ", nrow(object$data[[comp[1]]]), " or a boolean atomic vector.")
    }
 
    comp1 = round(comp[1])
    comp2 = round(comp[2])
    ##rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate"))
        
        # check that the user did not enter extra arguments #
        # --------------------------------------------------#
        # what the user has entered
        match.user =names(match.call())
        # what the function is expecting
        match.function = c('object', 'comp', 'ind.names', 'rep.space', 'X.label', 'Y.label', 'col', 'cex', 'pch')
        
        #if arguments are not matching, put a warning (put a [-1] for match.user as we have a first blank argument)
        if(length(setdiff(match.user[-1], match.function)) != 0) warning('Some of the input arguments do not match the function arguments, see ?plotIndiv')
        
    
    # l'espace de représentation #
    #----------------------------#
#    if (rep.space == "X-variate"){
      x = object$variates[[rep.space]][, comp1]
      y = object$variates[[rep.space]][, comp2]
      if (is.null(X.label)) X.label = paste("Block ", rep.space,': component ', comp1, sep = '')
      if (is.null(Y.label)) Y.label = paste("Block ", rep.space, ': component ', comp2, sep = '')

    
    # le plot des individus #
    #-----------------------#
    if (length(ind.names) > 1) {
      plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
      text(x, y, ind.names, col = col, cex = cex, ...)
      if(abline.line) abline(v = 0, h = 0, lty = 2)
    }
    else {
      if (isTRUE(ind.names)) {
        plot(x, y, type = "n", xlab = X.label, ylab = Y.label)
        text(x, y, ind.names, col = col, cex = cex, ...)
        if(abline.line) abline(v = 0, h = 0, lty = 2)
      }
      else {
        plot(x, y, xlab = X.label, ylab = Y.label, 
             col = col, cex = cex, pch = pch)
        if(abline.line) abline(v = 0, h = 0, lty = 2)
      }
    }

}  # end function







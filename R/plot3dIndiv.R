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


plot3dIndiv <-
function(object, ...) UseMethod("plot3dIndiv")

# ------------------------- PLS, sPLS, PLSDA, sPLDA objects -------------------------------
plot3dIndiv.pls <- plot3dIndiv.spls <- plot3dIndiv.plsda <-  plot3dIndiv.splsda <-
function(object, 
         comp = 1:3, 
         ind.names = FALSE,
         rep.space = "X-variate",
         xlab = NULL, 
         ylab = NULL, 
         zlab = NULL,
         col = "blue", 
         cex = 1, 
         pch = "s", 
         font = 1,
         axes.box = "box", 
         ...) 
{

    # validation des arguments #
    #--------------------------#
    if (length(comp) != 3)
        stop("'comp' must be a numeric vector of length 3.")

    if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.")

    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    if (is.logical(ind.names)) {
        if (isTRUE(ind.names)) ind.names = object$names$indiv
    }
	
	if (length(ind.names) > 1) {
        if (length(ind.names) != nrow(object$X))
            stop("'ind.names' must be a character vector of length ", nrow(object$X), " or a boolean atomic vector.")
    }

    # l'espace de représentation #
    #----------------------------#
    rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate"))
		
    if (rep.space == "X-variate") {
        x = object$variates$X[, comp[1]]
        y = object$variates$X[, comp[2]]
		z = object$variates$X[, comp[3]]
    }

    if (rep.space == "Y-variate") {
        x = object$variates$Y[, comp[1]]
        y = object$variates$Y[, comp[2]]
		z = object$variates$Y[, comp[3]]
    }
	
    if (rep.space == "XY-variate"){
        x = (object$variates$X[, comp[1]] + object$variates$Y[, comp[1]]) / 2
        y = (object$variates$X[, comp[2]] + object$variates$Y[, comp[2]]) / 2
        z = (object$variates$X[, comp[3]] + object$variates$Y[, comp[3]]) / 2
    }

    # le plot des individus #
    #-----------------------#
	pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
	if (length(pch) == 1) pch = rep(pch, length(x))
    pchlevels = levels(as.factor(pch))	
	
    if (length(col) == 1) col = rep(col, length(x))
	
    if ((length(cex) != length(pchlevels)) & (length(cex) > 1) & !isTRUE(ind.names)) 
        stop("'cex' must be a numeric vector of length equal than number of 'pch' levels.")
    else cex = rep(cex, length(pchlevels))
	
	opw = open3d(windowRect = c(500, 30, 1100, 630))
    par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    if (length(ind.names) > 1) {
        text3d(x, y, z, text = ind.names,
               color = col, font = font, cex = cex)
    }
    else {
        if (isTRUE(ind.names)) {
            text3d(x, y, z, text = ind.names,
               color = col, font = font, cex = cex[1])    
        }
        else {
            k = 0
            for (level in pchlevels) {
			    id = (pch == level)
				k = k + 1
                switch(pch.name[pch.type[, level]], 
                    sphere = plot3d(x = x[id], y = y[id], z = z[id], type = "s", 
                            col = col[id], radius = cex[k]/20, add = TRUE),
                    tetra = shapelist3d(tetrahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/25),
                    cube = shapelist3d(cube3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/30),
                    octa = shapelist3d(octahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/17),
                    icosa = shapelist3d(icosahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/20),
                    dodeca = shapelist3d(dodecahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/20))
            }
        }
    }
	
	if (axes.box == "box") {
        axes3d(marklen = 25)
        box3d()
    }
    if (axes.box == "bbox") {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
            specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)
    }
	if (axes.box == "both")	{
        axes3d(marklen = 25); box3d()
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)			   
    }

	if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])
	
    mtext3d(xlab, "x-+", line = 1)
    mtext3d(ylab, "y-+", line = 1.5)
    mtext3d(zlab, "z+-", line = 1)

    title3d(...)	

}

# ----------------------- RCC objects --------------------------------------------
plot3dIndiv.rcc <-
function (object, 
          comp = 1:3, 
          ind.names = FALSE,
          rep.space = "XY-variate",
          xlab = NULL, 
          ylab = NULL, 
          zlab = NULL,
          col = "blue", 
          cex = 1, 
          pch = "s", 
          font = 1,
          axes.box = "box", 
          ...) 
{

    # validation des arguments #
    #--------------------------#
    if (length(comp) != 3)
        stop("'comp' must be a numeric vector of length 3.")

    if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.")

    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")

    if (is.logical(ind.names)) {
        if (isTRUE(ind.names)) ind.names = object$names$indiv
    }
	
	if (length(ind.names) > 1) {
        if (length(ind.names) != nrow(object$X))
            stop("'ind.names' must be a character vector of length ", nrow(object$X), " or a boolean atomic vector.")
    }		

    # l'espace de représentation #
    #----------------------------#
    rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate"))
		
    if (rep.space == "X-variate") {
        x = object$variates$X[, comp[1]]
        y = object$variates$X[, comp[2]]
		z = object$variates$X[, comp[3]]
    }

    if (rep.space == "Y-variate") {
        x = object$variates$Y[, comp[1]]
        y = object$variates$Y[, comp[2]]
		z = object$variates$Y[, comp[3]]
    }
	
    if (rep.space == "XY-variate"){
        x = (object$variates$X[, comp[1]] + object$variates$Y[, comp[1]]) / 2
        y = (object$variates$X[, comp[2]] + object$variates$Y[, comp[2]]) / 2
        z = (object$variates$X[, comp[3]] + object$variates$Y[, comp[3]]) / 2
    }

    # le plot des individus #
    #-----------------------#
	pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
	if (length(pch) == 1) pch = rep(pch, length(x))
    pchlevels = levels(as.factor(pch))	
	
    if (length(col) == 1) col = rep(col, length(x))
	
    if ((length(cex) != length(pchlevels)) & (length(cex) > 1) & !isTRUE(ind.names)) 
        stop("'cex' must be a numeric vector of length equal than number of 'pch' levels.")
    else cex = rep(cex, length(pchlevels))
	
	opw = open3d(windowRect = c(500, 30, 1100, 630))	
	par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    if (length(ind.names) > 1) {
        text3d(x, y, z, text = ind.names,
               color = col, font = font, cex = cex)
    }
    else {
        if (isTRUE(ind.names)) {
            text3d(x, y, z, text = ind.names,
               color = col, font = font, cex = cex[1])    
        }
        else {
            k = 0
            for (level in pchlevels) {
			    id = (pch == level)
				k = k + 1
                switch(pch.name[pch.type[, level]], 
                    sphere = plot3d(x = x[id], y = y[id], z = z[id], type = "s", 
                            col = col[id], radius = cex[k]/20, add = TRUE),
                    tetra = shapelist3d(tetrahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/25),
                    cube = shapelist3d(cube3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/30),
                    octa = shapelist3d(octahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/17),
                    icosa = shapelist3d(icosahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/20),
                    dodeca = shapelist3d(dodecahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/20))
            }
        }
    }
	
	if (axes.box == "box") {
        axes3d(marklen = 25)
        box3d()
    }
    if (axes.box == "bbox") {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
            specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)
    }
	if (axes.box == "both")	{
        axes3d(marklen = 25); box3d()
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)			   
    }

	if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])
	
    mtext3d(xlab, "x-+", line = 1)
    mtext3d(ylab, "y-+", line = 1.5)
    mtext3d(zlab, "z+-", line = 1)

    title3d(...)	

}

# ------------------------------- PCA objects ----------------------------------------
plot3dIndiv.pca <-
function (object, 
          comp = 1:3, 
          ind.names = FALSE,
          xlab = NULL, 
          ylab = NULL, 
          zlab = NULL,
          col = "blue", 
          cex = 1, 
          pch = "s", 
          font = 1,
          axes.box = "box", 
          ...) 
{

    # validation des arguments #
    #--------------------------#
    if (length(comp) != 3)
        stop("'comp' must be a numeric vector of length 3.")

    if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.")

    if (any(comp > object$ncomp)) 
        stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")
		
    if (is.logical(ind.names)) {
        if (isTRUE(ind.names)) ind.names = rownames(object$x)
    }
	
	if (length(ind.names) > 1) {
        if (length(ind.names) != nrow(object$x))
            stop("'ind.names' must be a character vector of length ", nrow(object$x), " or a boolean atomic vector.")
    }

    # l'espace de représentation #
    #----------------------------#
    x = object$x[, comp[1]]
    y = object$x[, comp[2]]
    z = object$x[, comp[3]]

    # le plot des individus #
    #-----------------------#
	pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
	if (length(pch) == 1) pch = rep(pch, length(x))
    pchlevels = levels(as.factor(pch))	
	
    if (length(col) == 1) col = rep(col, length(x))
	
    if ((length(cex) != length(pchlevels)) & (length(cex) > 1) & !isTRUE(ind.names)) 
        stop("'cex' must be a numeric vector of length equal than number of 'pch' levels.")
    else cex = rep(cex, length(pchlevels))
	
	opw = open3d(windowRect = c(500, 30, 1100, 630))
	par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    if (length(ind.names) > 1) {
        text3d(x, y, z, text = ind.names,
               color = col, font = font, cex = cex)
    }
    else {
        if (isTRUE(ind.names)) {
            text3d(x, y, z, text = ind.names,
               color = col, font = font, cex = cex[1])    
        }
        else {
            cex = 20 * cex
            k = 0
            for (level in pchlevels) {
			    id = (pch == level)
				k = k + 1
                switch(pch.name[pch.type[, level]], 
                    sphere = plot3d(x = x[id], y = y[id], z = z[id], type = "s", 
                            col = col[id], radius = cex[k]/20, add = TRUE),
                    tetra = shapelist3d(tetrahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/25),
                    cube = shapelist3d(cube3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/30),
                    octa = shapelist3d(octahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/17),
                    icosa = shapelist3d(icosahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/20),
                    dodeca = shapelist3d(dodecahedron3d(), x = x[id], y = y[id], z = z[id], 
                            col = col[id], size = cex[k]/20))
            }
        }
    }
	
	if (axes.box == "box") {
        axes3d(marklen = 25)
        box3d()
    }
    if (axes.box == "bbox") {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
            specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)
    }
	if (axes.box == "both")	{
        axes3d(marklen = 25); box3d()
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)			   
    }

	if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])
	
    mtext3d(xlab, "x-+", line = 1)
    mtext3d(ylab, "y-+", line = 1.5)
    mtext3d(zlab, "z+-", line = 1)

    title3d(...)	

}

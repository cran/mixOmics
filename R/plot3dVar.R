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


plot3dVar <-
function(object, ...) UseMethod("plot3dVar")



# -------------------- PLS  -------------
plot3dVar.pls <- 
function(object, 
         comp = 1:3, 
         rad.in = 1, 
         X.label = FALSE, 
         Y.label = FALSE,
         pch = NULL, 
         cex = NULL, 
         col = NULL, 
         font = NULL, 
         axes.box = "all", 
         label.axes.box = "both",		
         xlab = NULL, 
         ylab = NULL, 
         zlab = NULL, 
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

    comp = round(comp)	

    # calcul des coordonnées #
    #------------------------#

      if (object$mode == "canonical") {
          cord.X = cor(object$X, object$variates$X[, comp], use = "pairwise")
          cord.Y = cor(object$Y, object$variates$Y[, comp], use = "pairwise")
      }
      else {
          cord.X = cor(object$X, object$variates$X[, comp], use = "pairwise")
          cord.Y = cor(object$Y, object$variates$X[, comp], use = "pairwise")
      }

    p = ncol(object$X)
    q = ncol(object$Y)

    # le plot des variables #
    #-----------------------#
    pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
    if (length(X.label) > 1 & length(X.label) != p)
        stop("'X.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

    if (length(Y.label) > 1 & length(Y.label) != q)
        stop("'Y.label' must be a vector of length 'ncol(Y)' or a boolean atomic vector.")
		
    if (is.null(pch)) {
        pch = c("s", "s")
    }

    if (is.null(cex)) {
        cex = c(1, 1)
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
        font = c(2, 3)
    }


    if (X.label) X.label = object$names$X
    if (Y.label) Y.label = object$names$Y

    if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])

    # le plot 3d #
    #------------#
    opw = open3d(windowRect = c(500, 30, 1100, 630))
     
	par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    if (any(axes.box == "axes") || any(axes.box == "all"))
        axes3d(c('x','y','z'), pos = c(0, 0, 0), nticks = 2, at = c(-1.2, 1.2), 
		    tick = FALSE, labels = "")
	
    if (length(X.label) > 1) {
        text3d(cord.X[, 1] + 0.05, cord.X[, 2], cord.X[, 3] + 0.05, text = X.label,
               color = col[[1]], font = font[1], cex = cex[1])
    }
	
    if (length(Y.label) > 1) {
        text3d(cord.Y[, 1] + 0.05, cord.Y[, 2], cord.Y[, 3] + 0.05, text = Y.label,
               color = col[[2]], font = font[2], cex = cex[2])
    }
	par3d(cex = 0.8)	
	
    if (any(axes.box == "axes") || any(axes.box == "all")) { 
        if (any(label.axes.box == "axes") || any(label.axes.box == "both")) {	
            text3d(1.2, -0.05, 0, text = xlab, cex = 0.8, color = "black")
            text3d(0, 1.27, 0, text = ylab, cex = 0.8, color = "black")
            text3d(0, -0.05, 1.2, text = zlab, cex = 0.8, color = "black")
		}
		
        x =	c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09,  1.09,
	          0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4))
	  
        y = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09,
	          0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4))
	  
        z = c(0.0, 0.035, -0.035, 0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4),
	          1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09)	  

        triangles3d(x = x, y = y, z = z, col = "black")
    }
	
    points3d(1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, 1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, 1.2, size = 0.1, alpha = 0)
    points3d(-1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, -1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, -1.2, size = 0.1, alpha = 0)
	
	cex = 1.5*cex	
	
    if (length(X.label) == 1) {
        switch(pch.name[pch.type[, pch[1]]], 
            sphere = plot3d(x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], type = "s", 
                            col = col[[1]], radius = cex[1]/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/65),
            cube = shapelist3d(cube3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/70),
            octa = shapelist3d(octahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50))
    }
   
    if (length(Y.label) == 1) {
        switch(pch.name[pch.type[, pch[2]]], 
            sphere = plot3d(x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], type = "s", 
                            col = col[[2]], radius = cex[2]/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/65),
            cube = shapelist3d(cube3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/70),
            octa = shapelist3d(octahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/50))
    }
	
    spheres3d(0, 0, 0, radius = rad.in, front = "fill", back = "fill", emission = gray(0.9), alpha = 0.6)
	spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", emission = gray(0.9))
	
    if (any(axes.box == "box") || any(axes.box == "all")) {
        axes3d(marklen = 25)
        box3d()
		if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }
	
    if (any(axes.box == "bbox") || any(axes.box == "all")) {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25) 
        if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }

    title3d(...)	
}


# ------------ plsda ------------------------------

plot3dVar.plsda <- 
function(object, 
         comp = 1:3, 
         rad.in = 1, 
         X.label = FALSE, 
         pch = NULL, 
         cex = NULL, 
         col = NULL, 
         font = NULL, 
         axes.box = "all", 
         label.axes.box = "both",		
         xlab = NULL, 
         ylab = NULL, 
         zlab = NULL, 
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

    comp = round(comp)	

    # calcul des coordonnées #
    #------------------------#
      if (object$mode == "canonical") {
          cord.X = cor(object$X, object$variates$X[, comp], use = "pairwise")
#          cord.Y = cor(object$Y, object$variates$Y[, comp], use = "pairwise")
      }
      else {
          cord.X = cor(object$X, object$variates$X[, comp], use = "pairwise")
      }

    p = ncol(object$X)

    # le plot des variables #
    #-----------------------#
    pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
    if (length(X.label) > 1 & length(X.label) != p)
        stop("'X.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")


    if (is.null(pch)) {
        pch = c("s")
    }

    if (is.null(cex)) {
        cex = c(1)
    }

    if (is.null(col)) {
        col = list(rep("red", p))
    }
    else {
        if (is.list(col)) {
            if (length(col[[1]]) != p) 
                stop("'col' must be a vector of length 1 or of length ", p)
        }
        else { 
            if (length(col) == 1) { 
                col = list(rep(col[1], p))
            }
            else {
                stop("'col' must be a vector of length 1 or of length ", p)
            }
        }
    }

    if (is.null(font)) {
        font = c(2)
    }


    if (X.label) X.label = object$names$X
    if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])

    # le plot 3d #
    #------------#
    opw = open3d(windowRect = c(500, 30, 1100, 630))
     
	par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    if (any(axes.box == "axes") || any(axes.box == "all"))
        axes3d(c('x','y','z'), pos = c(0, 0, 0), nticks = 2, at = c(-1.2, 1.2), 
		    tick = FALSE, labels = "")
	
    if (length(X.label) > 1) {
        text3d(cord.X[, 1] + 0.05, cord.X[, 2], cord.X[, 3] + 0.05, text = X.label,
               color = col[[1]], font = font[1], cex = cex[1])
    }
	
	par3d(cex = 0.8)	
	
    if (any(axes.box == "axes") || any(axes.box == "all")) { 
        if (any(label.axes.box == "axes") || any(label.axes.box == "both")) {	
            text3d(1.2, -0.05, 0, text = xlab, cex = 0.8, color = "black")
            text3d(0, 1.27, 0, text = ylab, cex = 0.8, color = "black")
            text3d(0, -0.05, 1.2, text = zlab, cex = 0.8, color = "black")
		}
		
        x =	c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09,  1.09,
	          0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4))
	  
        y = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09,
	          0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4))
	  
        z = c(0.0, 0.035, -0.035, 0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4),
	          1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09)	  

        triangles3d(x = x, y = y, z = z, col = "black")
    }
	
    points3d(1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, 1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, 1.2, size = 0.1, alpha = 0)
    points3d(-1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, -1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, -1.2, size = 0.1, alpha = 0)
	
	cex = 1.5*cex	
	
    if (length(X.label) == 1) {
        switch(pch.name[pch.type[, pch[1]]], 
            sphere = plot3d(x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], type = "s", 
                            col = col[[1]], radius = cex[1]/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/65),
            cube = shapelist3d(cube3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/70),
            octa = shapelist3d(octahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50))
    }
   
    spheres3d(0, 0, 0, radius = rad.in, front = "fill", back = "fill", emission = gray(0.9), alpha = 0.6)
	spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", emission = gray(0.9))
	
    if (any(axes.box == "box") || any(axes.box == "all")) {
        axes3d(marklen = 25)
        box3d()
		if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }
	
    if (any(axes.box == "bbox") || any(axes.box == "all")) {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25) 
        if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }

    title3d(...)	
}



# -------------------- sPLS -------------
plot3dVar.spls <-  
function(object, 
         comp = 1:3, 
         rad.in = 1, 
         X.label = FALSE, 
         Y.label = FALSE,
         pch = NULL, 
         cex = NULL, 
         col = NULL, 
         font = NULL, 
         axes.box = "all", 
         label.axes.box = "both",		
         xlab = NULL, 
         ylab = NULL, 
         zlab = NULL, 
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

    comp = round(comp)	

    # calcul des coordonnées #
    #------------------------#
#    if (isTRUE(keep.var)) {
        keep.X = apply(abs(object$loadings$X), 1, sum) > 0
        keep.Y = apply(abs(object$loadings$Y), 1, sum) > 0

        if (object$mode == "canonical") {
            cord.X = cor(object$X[, keep.X], object$variates$X[, comp], 
                     use = "pairwise")
            cord.Y = cor(object$Y[, keep.Y], object$variates$Y[, comp], 
                     use = "pairwise")
        }
        else {
            cord.X = cor(object$X[, keep.X], object$variates$X[, comp], 
                     use = "pairwise")
            cord.Y = cor(object$Y[, keep.Y], object$variates$X[, comp], 
                     use = "pairwise")
        }

    p = ncol(object$X)
    q = ncol(object$Y)

    # le plot des variables #
    #-----------------------#
    pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
    if (length(X.label) > 1 & length(X.label) != p)
        stop("'X.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

    if (length(Y.label) > 1 & length(Y.label) != q)
        stop("'Y.label' must be a vector of length 'ncol(Y)' or a boolean atomic vector.")
		
    if (is.null(pch)) {
        pch = c("s", "s")
    }

    if (is.null(cex)) {
        cex = c(1, 1)
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
        font = c(2, 3)
    }

      col[[1]] = col[[1]][keep.X]
      col[[2]] = col[[2]][keep.Y]

    if (X.label) X.label = object$names$X
    if (Y.label) Y.label = object$names$Y

      if (length(X.label) == p) X.label = X.label[keep.X]
      if (length(Y.label) == q) Y.label = Y.label[keep.Y]

    if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])

    # le plot 3d #
    #------------#
    opw = open3d(windowRect = c(500, 30, 1100, 630))
     
	par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    if (any(axes.box == "axes") || any(axes.box == "all"))
        axes3d(c('x','y','z'), pos = c(0, 0, 0), nticks = 2, at = c(-1.2, 1.2), 
		    tick = FALSE, labels = "")
	
    if (length(X.label) > 1) {
        text3d(cord.X[, 1] + 0.05, cord.X[, 2], cord.X[, 3] + 0.05, text = X.label,
               color = col[[1]], font = font[1], cex = cex[1])
    }
	
    if (length(Y.label) > 1) {
        text3d(cord.Y[, 1] + 0.05, cord.Y[, 2], cord.Y[, 3] + 0.05, text = Y.label,
               color = col[[2]], font = font[2], cex = cex[2])
    }
	par3d(cex = 0.8)	
	
    if (any(axes.box == "axes") || any(axes.box == "all")) { 
        if (any(label.axes.box == "axes") || any(label.axes.box == "both")) {	
            text3d(1.2, -0.05, 0, text = xlab, cex = 0.8, color = "black")
            text3d(0, 1.27, 0, text = ylab, cex = 0.8, color = "black")
            text3d(0, -0.05, 1.2, text = zlab, cex = 0.8, color = "black")
		}
		
        x =	c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09,  1.09,
	          0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4))
	  
        y = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09,
	          0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4))
	  
        z = c(0.0, 0.035, -0.035, 0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4),
	          1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09)	  

        triangles3d(x = x, y = y, z = z, col = "black")
    }
	
    points3d(1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, 1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, 1.2, size = 0.1, alpha = 0)
    points3d(-1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, -1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, -1.2, size = 0.1, alpha = 0)
	
	cex = 1.5*cex	
	
    if (length(X.label) == 1) {
        switch(pch.name[pch.type[, pch[1]]], 
            sphere = plot3d(x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], type = "s", 
                            col = col[[1]], radius = cex[1]/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/65),
            cube = shapelist3d(cube3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/70),
            octa = shapelist3d(octahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50))
    }
   
    if (length(Y.label) == 1) {
        switch(pch.name[pch.type[, pch[2]]], 
            sphere = plot3d(x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], type = "s", 
                            col = col[[2]], radius = cex[2]/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/65),
            cube = shapelist3d(cube3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/70),
            octa = shapelist3d(octahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/50))
    }
	
    spheres3d(0, 0, 0, radius = rad.in, front = "fill", back = "fill", emission = gray(0.9), alpha = 0.6)
	spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", emission = gray(0.9))
	
    if (any(axes.box == "box") || any(axes.box == "all")) {
        axes3d(marklen = 25)
        box3d()
		if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }
	
    if (any(axes.box == "bbox") || any(axes.box == "all")) {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25) 
        if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }

    title3d(...)	
}

# ------------------------------ sPLS-DA ---------------------------------------------
plot3dVar.splsda <- 
function(object, 
         comp = 1:3, 
         rad.in = 1, 
         X.label = FALSE, 
         pch = NULL, 
         cex = NULL, 
         col = NULL, 
         font = NULL, 
         axes.box = "all", 
         label.axes.box = "both",		
         xlab = NULL, 
         ylab = NULL, 
         zlab = NULL, 
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

    comp = round(comp)	

    # calcul des coordonnées #
    #------------------------#
      keep.X = apply(abs(object$loadings$X), 1, sum) > 0

##      if (object$mode == "canonical") {
##          cord.X = cor(object$X[, keep.X], object$variates$X[, comp], 
##                   use = "pairwise")
##      }
##      else {
          cord.X = cor(object$X[, keep.X], object$variates$X[, comp], 
                   use = "pairwise")
##      }

    p = ncol(object$X)

    # le plot des variables #
    #-----------------------#
    pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
    if (length(X.label) > 1 & length(X.label) != p)
        stop("'X.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

   if (is.null(pch)) {
        pch = c("s")
    }

    if (is.null(cex)) {
        cex = c(1)
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
            else {
                stop("'col' must be a vector of length 1 or a vector of length ", p)
            }
        }
    }

    if (is.null(font)) {
        font = c(2)
    }

      col[[1]] = col[[1]][keep.X]

    if (X.label) X.label = object$names$X

      if (length(X.label) == p) X.label = X.label[keep.X]

    if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])

    # le plot 3d #
    #------------#
    opw = open3d(windowRect = c(500, 30, 1100, 630))
     
	par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    if (any(axes.box == "axes") || any(axes.box == "all"))
        axes3d(c('x','y','z'), pos = c(0, 0, 0), nticks = 2, at = c(-1.2, 1.2), 
		    tick = FALSE, labels = "")
	
    if (length(X.label) > 1) {
        text3d(cord.X[, 1] + 0.05, cord.X[, 2], cord.X[, 3] + 0.05, text = X.label,
               color = col[[1]], font = font[1], cex = cex[1])
    }
	
	par3d(cex = 0.8)	
	
    if (any(axes.box == "axes") || any(axes.box == "all")) { 
        if (any(label.axes.box == "axes") || any(label.axes.box == "both")) {	
            text3d(1.2, -0.05, 0, text = xlab, cex = 0.8, color = "black")
            text3d(0, 1.27, 0, text = ylab, cex = 0.8, color = "black")
            text3d(0, -0.05, 1.2, text = zlab, cex = 0.8, color = "black")
		}
		
        x =	c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09,  1.09,
	          0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4))
	  
        y = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09,
	          0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4))
	  
        z = c(0.0, 0.035, -0.035, 0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4),
	          1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09)	  

        triangles3d(x = x, y = y, z = z, col = "black")
    }
	
    points3d(1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, 1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, 1.2, size = 0.1, alpha = 0)
    points3d(-1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, -1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, -1.2, size = 0.1, alpha = 0)
	
	cex = 1.5*cex	
	
    if (length(X.label) == 1) {
        switch(pch.name[pch.type[, pch[1]]], 
            sphere = plot3d(x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], type = "s", 
                            col = col[[1]], radius = cex[1]/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/65),
            cube = shapelist3d(cube3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/70),
            octa = shapelist3d(octahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50))
    }
   
	
    spheres3d(0, 0, 0, radius = rad.in, front = "fill", back = "fill", emission = gray(0.9), alpha = 0.6)
	spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", emission = gray(0.9))
	
    if (any(axes.box == "box") || any(axes.box == "all")) {
        axes3d(marklen = 25)
        box3d()
		if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }
	
    if (any(axes.box == "bbox") || any(axes.box == "all")) {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25) 
        if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }

    title3d(...)	
}


# ------------------------------- for RCC ----------------------------------------------
plot3dVar.rcc <-
function(object, 
         comp = 1:3, 
         rad.in = 1, 
         cutoff = NULL, 
         X.label = FALSE, 
         Y.label = FALSE, 
         pch = NULL, 
         cex = NULL, 
         col = NULL, 
         font = NULL,
         axes.box = "all", 
         label.axes.box = "both",		
         xlab = NULL, 
         ylab = NULL, 
         zlab = NULL, 
         ...) 
{

    # validation des arguments #
    #--------------------------#
    if (length(comp) != 3)
        stop("'comp' must be a numeric vector of length 3.")

    if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.")

    dim = min(ncol(object$X), ncol(object$Y))
    if (any(comp > dim)) 
        stop("the elements of 'comp' must be smaller or equal than ", dim, ".")

    comp = round(comp)	
	
    # calcul des coordonnées #
    #------------------------#

    p = ncol(object$X)
    q = ncol(object$Y)
    dim = min(p, q)	

    bisect = object$variates$X[, comp] + object$variates$Y[, comp]
    cord.X = cor(object$X, bisect, use = "pairwise")
    cord.Y = cor(object$Y, bisect, use = "pairwise")

    if (!is.null(cutoff)) {
        # choix des variables avec au moins une coordonnée # 
        # supérieur au cutoff                              #
        #--------------------------------------------------#
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
            stop("Cutoff value very high for the components of 'comp'. No variable was selected.")

        cord.X = matrix(cord.X[gp.X, ], length(gp.X), 3)
        cord.Y = matrix(cord.Y[gp.Y, ], length(gp.Y), 3)
        rad.in = cutoff
    }

    # le plot des variables #
    #-----------------------#
    pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
    if (length(X.label) > 1 & length(X.label) != p)
        stop("'X.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

    if (length(Y.label) > 1 & length(Y.label) != q)
        stop("'Y.label' must be a vector of length 'ncol(Y)' or a boolean atomic vector.")

    if (is.null(pch)) {
        pch = c("s", "s")
    }

    if (is.null(cex)) {
        cex = c(1, 1)
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
        font = c(2, 3)
    }

    if (!is.null(cutoff)) {
        col[[1]] = col[[1]][gp.X]
        col[[2]] = col[[2]][gp.Y]
    }

    if (X.label) X.label = object$names$X
    if (Y.label) Y.label = object$names$Y

    if (!is.null(cutoff)) {
        if (length(X.label) == p) X.label = X.label[gp.X]
        if (length(Y.label) == q) Y.label = Y.label[gp.Y]
    }

    if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])

    # le plot 3d #
    #------------#
    opw = open3d(windowRect = c(500, 30, 1100, 630))
	par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    if (any(axes.box == "axes") || any(axes.box == "all"))
        axes3d(c('x','y','z'), pos = c(0, 0, 0), nticks = 2, at = c(-1.2, 1.2), 
		    tick = FALSE, labels = "")
	
    if (length(X.label) > 1) {
        text3d(cord.X[, 1] + 0.05, cord.X[, 2], cord.X[, 3] + 0.05, text = X.label,
               color = col[[1]], font = font[1], cex = cex[1])
    }
	
    if (length(Y.label) > 1) {
        text3d(cord.Y[, 1] + 0.05, cord.Y[, 2], cord.Y[, 3] + 0.05, text = Y.label,
               color = col[[2]], font = font[2], cex = cex[2])
    }
	par3d(cex = 0.8)	
	
    if (any(axes.box == "axes") || any(axes.box == "all")) { 
        if (any(label.axes.box == "axes") || any(label.axes.box == "both")) {	
            text3d(1.2, -0.05, 0, text = xlab, cex = 0.8, color = "black")
            text3d(0, 1.27, 0, text = ylab, cex = 0.8, color = "black")
            text3d(0, -0.05, 1.2, text = zlab, cex = 0.8, color = "black")
		}
		
        x =	c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09,  1.09,
	          0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4))
	  
        y = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09,
	          0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4))
	  
        z = c(0.0, 0.035, -0.035, 0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4),
	          1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09)	  

        triangles3d(x = x, y = y, z = z, col = "black")
    }
	
    points3d(1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, 1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, 1.2, size = 0.1, alpha = 0)
    points3d(-1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, -1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, -1.2, size = 0.1, alpha = 0)
	
	cex = 1.5*cex	
	
    if (length(X.label) == 1) {
        switch(pch.name[pch.type[, pch[1]]], 
            sphere = plot3d(x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], type = "s", 
                            col = col[[1]], radius = cex[1]/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/65),
            cube = shapelist3d(cube3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/70),
            octa = shapelist3d(octahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50))
    }
   
    if (length(Y.label) == 1) {
        switch(pch.name[pch.type[, pch[2]]], 
            sphere = plot3d(x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], type = "s", 
                            col = col[[2]], radius = cex[2]/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/65),
            cube = shapelist3d(cube3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/70),
            octa = shapelist3d(octahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.Y[, 1], y = cord.Y[, 2], z = cord.Y[, 3], 
                            col = col[[2]], size = cex[2]/50))
    }
	
    spheres3d(0, 0, 0, radius = rad.in, front = "fill", back = "fill", emission = gray(0.9), alpha = 0.6)
	spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", emission = gray(0.9))
	
    if (any(axes.box == "box") || any(axes.box == "all")) {
        axes3d(marklen = 25)
        box3d()
		if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }
	
    if (any(axes.box == "bbox") || any(axes.box == "all")) {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25) 
        if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }

    title3d(...)	
}

# ------------------------------ PCA object ------------------------------------
plot3dVar.pca <-
function(object, 
         comp = 1:3, 
	 rad.in = 1, 
         var.label = FALSE, 
         pch = NULL, 
         cex = NULL, 
         col = NULL, 
         font = NULL,
         axes.box = "all", 
         label.axes.box = "both",		
         xlab = NULL, 
         ylab = NULL, 
         zlab = NULL,		 
         ...) 
{


    # validation des arguments #
    #--------------------------#
    if (length(comp) != 3)
        stop("'comp' must be a numeric vector of length 3.")

    if (!is.numeric(comp) || any(comp < 1))
        stop("invalid vector for 'comp'.")

    p = ncol(object$rotation)
	q = nrow(object$rotation)
	
    if (any(comp > p)) 
        stop("the elements of 'comp' must be smaller or equal than ", p, ".")

    comp = round(comp)
	
    # calcul des coordonnées #
    #------------------------#
    cord.X = object$rotation[, comp] 

    # le plot des variables #
    #-----------------------#
    pch.type = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")

    if (is.null(pch)) pch = "s"
    if (is.null(cex)) cex = 1
    if (is.null(col)) col = rep("blue", q)
    if (is.null(font)) font = 2
    if (is.logical(var.label)) {
        if (isTRUE(var.label)) var.label = rownames(object$rotation)
    }
	
	if (length(var.label) > 1) {
        if (length(var.label) != q)
            stop("'var.label' must be a character vector of length ", q, ".")
    }

    if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])

    # le plot 3d #
    #------------#
    opw = open3d(windowRect = c(500, 30, 1100, 630))
	par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    MaxAbs = max(abs(cord.X))
    tickmark = round(seq(-MaxAbs, MaxAbs, length = 5), 2)
    tick = seq(-1, 1, length = 5)	
    cord.X = cord.X / MaxAbs

    if (any(axes.box == "axes") || any(axes.box == "all"))
        axes3d(c('x','y','z'), pos = c(0, 0, 0), nticks = 2, at = c(-1.2, 1.2), 
		    tick = FALSE, labels = "")
	
    if (length(var.label) > 1) {
        text3d(cord.X[, 1], cord.X[, 2], cord.X[, 3], text = var.label,
               color = col, font = font, cex = cex)
    }
	
	par3d(cex = 0.8)	
	
    if (any(axes.box == "axes") || any(axes.box == "all")) { 
        if (any(label.axes.box == "axes") || any(label.axes.box == "both")) {	
            text3d(1.2, -0.05, 0, text = xlab, cex = 0.8, color = "black")
            text3d(0, 1.27, 0, text = ylab, cex = 0.8, color = "black")
            text3d(0, -0.05, 1.2, text = zlab, cex = 0.8, color = "black")
		}
		
        x =	c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09,  1.09,
	          0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4))
	  
        y = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09,
	          0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4))
	  
        z = c(0.0, 0.035, -0.035, 0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4),
	          1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09)	  

        triangles3d(x = x, y = y, z = z, col = "black")
    }
	
    points3d(1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, 1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, 1.2, size = 0.1, alpha = 0)
    points3d(-1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, -1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, -1.2, size = 0.1, alpha = 0)
	
	cex = 1.5 * cex	

    if (length(var.label) == 1) {
        switch(pch.name[which(pch.type == pch)], 
            sphere = plot3d(x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], type = "s", 
                            col = col, radius = cex/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col, size = cex/65),
            cube = shapelist3d(cube3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col, size = cex/70),
            octa = shapelist3d(octahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col, size = cex/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col, size = cex/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col, size = cex/50))
    }

    spheres3d(0, 0, 0, radius = rad.in, front = "fill", back = "fill", emission = gray(0.9), alpha = 0.6)
	spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", emission = gray(0.9))


    if (any(axes.box == "box") || any(axes.box == "all")) {
        axes3d(marklen = 25, xat = tick, yat = tick, zat = tick,
               xlab = tickmark, ylab = tickmark, zlab = tickmark)
        box3d()
		if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }
	
    if (any(axes.box == "bbox") || any(axes.box == "all")) {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5),
               xat = tick, yat = tick, zat = tick,		
               xlab = tickmark, ylab = tickmark, zlab = tickmark, 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25) 
        if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }

    title3d(...)	
}










# ----------------------------- sPCA---------------------------------------------
plot3dVar.spca <- 
function(object, 
         comp = 1:3, 
         rad.in = 1, 
         var.label = FALSE, 
         pch = NULL, 
         cex = NULL, 
         col = NULL, 
         font = NULL, 
         axes.box = "all", 
         label.axes.box = "both",		
         xlab = NULL, 
         ylab = NULL, 
         zlab = NULL, 
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

    comp = round(comp)	

    # calcul des coordonnées #
    #------------------------#
      keep.X = apply(abs(object$rotation), 1, sum) > 0

          cord.X = cor(object$X[, keep.X], object$x[, comp], 
                   use = "pairwise")

    p = ncol(object$X)

    # le plot des variables #
    #-----------------------#
    pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
    if (length(var.label) > 1 & length(var.label) != p)
        stop("'var.label' must be a vector of length 'ncol(X)' or a boolean atomic vector.")

   if (is.null(pch)) {
        pch = c("s")
    }

    if (is.null(cex)) {
        cex = c(1)
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
            else {
                stop("'col' must be a vector of length 1 or a vector of length ", p)
            }
        }
    }

    if (is.null(font)) {
        font = c(2)
    }

      col[[1]] = col[[1]][keep.X]

    if (var.label) var.label = object$names$X

      if (length(var.label) == p) var.label = var.label[keep.X]

    if (is.null(xlab)) xlab = paste("Comp ", comp[1])
	if (is.null(ylab)) ylab = paste("Comp ", comp[2])
	if (is.null(zlab)) zlab = paste("Comp ", comp[3])

    # le plot 3d #
    #------------#
    opw = open3d(windowRect = c(500, 30, 1100, 630))
     
	par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
	
    if (any(axes.box == "axes") || any(axes.box == "all"))
        axes3d(c('x','y','z'), pos = c(0, 0, 0), nticks = 2, at = c(-1.2, 1.2), 
		    tick = FALSE, labels = "")
	
    if (length(var.label) > 1) {
        text3d(cord.X[, 1] + 0.05, cord.X[, 2], cord.X[, 3] + 0.05, text = var.label,
               color = col[[1]], font = font[1], cex = cex[1])
    }
	
	par3d(cex = 0.8)	
	
    if (any(axes.box == "axes") || any(axes.box == "all")) { 
        if (any(label.axes.box == "axes") || any(label.axes.box == "both")) {	
            text3d(1.2, -0.05, 0, text = xlab, cex = 0.8, color = "black")
            text3d(0, 1.27, 0, text = ylab, cex = 0.8, color = "black")
            text3d(0, -0.05, 1.2, text = zlab, cex = 0.8, color = "black")
		}
		
        x =	c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09,  1.09,
	          0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4))
	  
        y = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09,
	          0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4))
	  
        z = c(0.0, 0.035, -0.035, 0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
              0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4),
	          1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09)	  

        triangles3d(x = x, y = y, z = z, col = "black")
    }
	
    points3d(1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, 1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, 1.2, size = 0.1, alpha = 0)
    points3d(-1.2, 0, 0, size = 0.1, alpha = 0)	
	points3d(0, -1.2, 0, size = 0.1, alpha = 0)
	points3d(0, 0, -1.2, size = 0.1, alpha = 0)
	
	cex = 1.5*cex	
	
    if (length(var.label) == 1) {
        switch(pch.name[pch.type[, pch[1]]], 
            sphere = plot3d(x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], type = "s", 
                            col = col[[1]], radius = cex[1]/50, add = TRUE),
            tetra = shapelist3d(tetrahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/65),
            cube = shapelist3d(cube3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/70),
            octa = shapelist3d(octahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/35),
            icosa = shapelist3d(icosahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50),
            dodeca = shapelist3d(dodecahedron3d(), x = cord.X[, 1], y = cord.X[, 2], z = cord.X[, 3], 
                            col = col[[1]], size = cex[1]/50))
    }
   
	
    spheres3d(0, 0, 0, radius = rad.in, front = "fill", back = "fill", emission = gray(0.9), alpha = 0.6)
	spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", emission = gray(0.9))
	
    if (any(axes.box == "box") || any(axes.box == "all")) {
        axes3d(marklen = 25)
        box3d()
		if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }
	
    if (any(axes.box == "bbox") || any(axes.box == "all")) {
        bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
               specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25) 
        if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
            mtext3d(xlab, "x-+", line = 1)
            mtext3d(ylab, "y-+", line = 1.5)
            mtext3d(zlab, "z+-", line = 1)
        }
    }

    title3d(...)	
}


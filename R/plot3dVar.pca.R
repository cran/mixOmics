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


plot3dVar.pca <-
function(object, 
         comp = 1:3,  
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

    p = ncol(object$x)
    if (any(comp > p)) 
        stop("the elements of 'comp' must be smaller or equal than ", p, ".")

    comp = round(comp)
	
    # calcul des coordonnées #
    #------------------------#
    cord.X = object$rotation[, comp] 

    # le plot des variables #
    #-----------------------#
    pch.type = matrix(1:6, ncol = 6)
    colnames(pch.type) = c("s", "t", "c", "o", "i", "d")
    pch.name = c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")
	
    if (length(X.label) > 1 & length(X.label) != p)
        stop("'X.label' must be a vector of length ", p, " or a boolean atomic vector.")

    if (is.null(pch)) pch = "s"

    if (is.null(cex)) cex = 1

    if (is.null(col)) col = rep("blue", p)
   
    if (is.null(font)) font = 2

    if (X.label) X.label = rownames(object$x)

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
        switch(pch.name[pch.type[pch[1]]], 
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

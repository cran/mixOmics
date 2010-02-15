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



#`imgCor` <- function(X, ...) UseMethod("imgCor")

`imgCor` <-
function(X, Y, type = c("combine", "separated"), col = jet.colors(64), ...) 
{

if (length(dim(X)) != 2 || length(dim(Y)) != 2) 
        stop("'X' and/or 'Y' must be a numeric matrix.")

X = as.matrix(X)
Y = as.matrix(Y)

if (!is.numeric(X) || !is.numeric(Y)) 
        stop("'X' and/or 'Y' must be a numeric matrix.")

    type = match.arg(type)
    p = ncol(X)
    q = ncol(Y)

matcor = cor(cbind(X, Y), use = "pairwise")
breaks = seq(-1, 1, length = length(col) + 1)

# représentation de la matrice de corrélation de #
# la concatenationdes variables X et Y, [X Y]    #
#------------------------------------------------#
    def.par = par(no.readonly = TRUE)

    if (type == "combine") {
        layout(matrix(c(1, 1, 1, 1, 2, 2), ncol = 2, nrow = 3, 
             byrow = TRUE), widths = 1, heights = c(0.8, 1, 0.35))

#-- layout 1 --# 
        par(pty = "s")
        image(1:(p + q), 1:(p + q), t(matcor[(p + q):1, ]), 
            zlim = c(-1, 1), main = "Combine [X Y] correlation", 
            col = col, axes = FALSE, xlab = "", ylab = "",
breaks = breaks)
        box()
        abline(h = q + 0.5, v = p + 0.5, lwd = 1, lty = 2)

#-- layout 2 --#
par(pty = "m", mai = c(0.6, 1.2, 0.1, 1))  
z = seq(-1, 1, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = breaks, 
            xaxt = "n", yaxt = "n")
box()
        par(usr = c(-1, 1, -1, 1))
        axis(1, at = c(-1, -0.5, 0, 0.5, 1))
mtext(side = 1, "Value", line = 2.5, cex = 0.8)
    }

# représentation des matrices de corrélation de #
# X, Y et entre X et Y                          #
#-----------------------------------------------#
    if (type == "separated") {
Xcor = cor(X, use = "pairwise")
Ycor = cor(Y, use = "pairwise")
XYcor = cor(X, Y, use = "pairwise")
        layout(matrix(c(1, 2, 3, 3, 4, 4), ncol = 2, nrow = 3, 
             byrow = TRUE), widths = 1, heights = c(0.8, 1, 0.35))

#-- layout 1 --#
        par(pty = "s", mar = c(2, 2, 2, 1))
        image(1:p, 1:p, t(Xcor[p:1, ]), zlim = c(-1, 1), col = col,
            main = "X correlation", axes = FALSE, xlab = "", ylab = "",
breaks = breaks)
        box()

#-- layout 2 --#
        image(1:q, 1:q, t(Ycor[q:1, ]), zlim = c(-1, 1), col = col,
main = "Y correlation", axes = FALSE, xlab = "", ylab = "",
breaks = breaks)
        box()

#-- layout 3 --#
        if (p > q) {
            XYcor = t(XYcor)
            p = ncol(Ycor)
            q = ncol(Xcor)
        }

        par(pty = "m", mai = c(0.25, 0.5, 0.3, 0.4)) 
        image(1:q, 1:p, t(XYcor), zlim = c(-1, 1), col = col, 
            main = "Cross-correlation", axes = FALSE, xlab = "", ylab = "",
breaks = breaks)
        box()

#-- layout 4 --#
par(mai = c(0.6, 1.2, 0.1, 1))  
z = seq(-1, 1, length = length(col))
breaks = seq(-1, 1, length = length(col) + 1)

        image(z = matrix(z, ncol = 1), col = col, breaks = breaks, 
            xaxt = "n", yaxt = "n")
box()

        par(usr = c(-1, 1, -1, 1))
        axis(1, at = c(-1, -0.5, 0, 0.5, 1))
mtext(side = 1, "Value", line = 2.5, cex = 0.8)
    }
    par(def.par)
}


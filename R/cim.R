# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
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


cim <-
function(...) UseMethod("cim")

# -------------------------------------
# CIM default
# -------------------------------------

cim.default <-function(
		mat, 
         breaks, 
         col = jet.colors, 
         distfun = dist, 
         hclustfun = hclust,
		 dendrogram = c("both", "row", "column", "none"),
         labRow = NULL, 
         labCol = NULL,
         ColSideColors = NULL,
         RowSideColors = NULL,		 
         symkey = TRUE, 
         keysize = 1,
         zoom = FALSE, 
         main = NULL, 
         xlab = NULL, 
         ylab = NULL,  
         cexRow = min(1, 0.2 + 1/log10(nr)), 
         cexCol = min(1, 0.2 + 1/log10(nc)), 
         margins = c(5, 5), 
         lhei = NULL, 
         lwid = NULL,
         ...) 
{

    # validation des arguments #
	#--------------------------#
    if (length(dim(mat)) != 2) 
        stop("'mat' must be a numeric matrix.")

    mat = simMat = as.matrix(mat)

    if (!is.numeric(mat)) 
        stop("'mat' must be a numeric matrix.")

    nr = nrow(mat)
    nc = ncol(mat)
	
    if (!is.null(labRow)) {
        if (length(labRow) != nr)
            stop("the length of 'labRow' must be equal to nrow(mat) = ", nr, ".")
    }

    if (!is.null(labCol)) {
        if (length(labCol) != nc)
            stop("the length of 'labCol' must be equal to ncol(mat) = ", nc, ".")
    }	

    if (isTRUE(symkey)) {
        max.mat = max(abs(mat))
        min.mat = -max.mat
    }
    else {
        max.mat = max(mat)
        min.mat = min(mat)
    }

    if (!is.numeric(margins) || length(margins) != 2) 
        stop("'margins' must be a numeric vector of length 2.")
		
	dendrogram = match.arg(dendrogram)
	
    if (is.null(labRow)) 
        labRow = if (is.null(rownames(mat))) (1:nr)[rowInd] else rownames(mat)

    if (is.null(labCol)) 
        labCol = if (is.null(colnames(mat))) (1:nc)[colInd] else colnames(mat) 

	if (any(dendrogram == "both") || any(dendrogram == "row")) { 
        Rowv = rowMeans(mat)
        hcr = hclustfun(distfun(mat))
        ddr = as.dendrogram(hcr)
        ddr = reorder(ddr, Rowv)
        rowInd = order.dendrogram(ddr)
        mat = mat[rowInd, ]
        labRow = labRow[rowInd]
    }    
     
	if (any(dendrogram == "both") || any(dendrogram == "column")) {
        Colv = colMeans(mat)
        hcc = hclustfun(distfun(t(mat)))
        ddc = as.dendrogram(hcc)
        ddc = reorder(ddc, Colv)
        colInd = order.dendrogram(ddc)
        mat = mat[, colInd]
        labCol = labCol[colInd]
    }
	
    rownames(mat) = labRow
    colnames(mat) = labCol

    if (missing(breaks) || is.null(breaks)) {
        if (class(col) == "function") breaks = 33
            else breaks = length(col) + 1
    }

    if (length(breaks) == 1) {
        breaks = seq(min.mat, max.mat, length = breaks)
        if (missing(col)) col = jet.colors
    }

    nbr = length(breaks)
    ncol = length(breaks) - 1

    if (class(col) == "function") 
        col = col(ncol)

    min.breaks = min(breaks)
    max.breaks = max(breaks)

    mat[mat < min.breaks] = min.breaks
    mat[mat > max.breaks] = max.breaks
    mat = t(mat)

    lmat = matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
	csc = rsc = FALSE
	
    if (!missing(ColSideColors) || !is.null(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != nc) 
            stop("'ColSideColors' must be a colors character vector of length ncol(mat)")
        lmat = rbind(lmat[1, ], c(NA, 3), lmat[2, ] + 1)
        if (missing(lhei) || is.null(lhei)) 
            lhei = c(keysize, 0.15, 4)
        if (length(lhei) != nrow(lmat)) 
            stop("lhei must have length = ", nrow(lmat))
        csc = TRUE
        if (any(dendrogram == "both") || any(dendrogram == "column"))
            ColSideColors = ColSideColors[colInd]
    }
     
    if (!missing(RowSideColors) || !is.null(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != nr) 
            stop("'RowSideColors' must be a colors character vector of length nrow(mat)")
        lmat = cbind(lmat[, 1], c(rep(NA, nrow(lmat) - 1), nrow(lmat) + 2), lmat[, 2] + 
               c(rep(0, nrow(lmat) - 1), 1))
        if (missing(lwid) || is.null(lwid)) 
            lwid = c(keysize, 0.15, 4)
        if (length(lwid) != ncol(lmat)) 
            stop("lwid must have length = ", ncol(lmat))
        rsc = TRUE
        if (any(dendrogram == "both") || any(dendrogram == "row")) 
            RowSideColors = RowSideColors[rowInd]
    }

    lmat[is.na(lmat)] = 0

    if (missing(lhei) || is.null(lhei)) 
        lhei = c(keysize, 4)

    if (missing(lwid) || is.null(lwid)) 
        lwid = c(keysize, 4)

    if (length(lhei) != nrow(lmat)) 
        stop("'lhei' must have length = ", nrow(lmat))

    if (length(lwid) != ncol(lmat)) 
        stop("'lwid' must have length = ", ncol(lmat))

    if (isTRUE(zoom)) {
        getOption("device")("xpos" = 65)
    }

    op = par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    #-- layout 1 --#
    par(mar = c(5, 2, 2, 1), cex = 0.75)

    z = seq(0, 1, length = length(col))
    z = matrix(z, ncol = 1)
    image(z, col = col, xaxt = "n", yaxt = "n")
    box()
    par(usr = c(0, 1, 0, 1))
    lv = c(min.breaks, (3*min.breaks + max.breaks)/4, (min.breaks + max.breaks)/2,
           (3*max.breaks + min.breaks)/4, max.breaks)
    xv = (as.numeric(lv) - min.mat) / (max.mat - min.mat)
    axis(1, at = xv, labels = round(lv, 2))
    title("Color key", font.main = 1)
	
    #-- layout 2 --#   
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (any(dendrogram == "both") || any(dendrogram == "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else {
        plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")	
    }

    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])

    #-- layout 3 --#
    if (isTRUE(csc)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors, axes = FALSE)	
    }
	
    #-- layout 4 --#
    par(mar = c(margins[1], 0, 0, 0))
    if (any(dendrogram == "both") || any(dendrogram == "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")    
    }
    else {
        plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")	
    }

    #-- layout 5 --#
    if (isTRUE(rsc)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors, axes = FALSE)	
    }
	
    #-- layout 6 --#
    par(mar = c(margins[1], 0, 0, margins[2]))
    image(1:nc, 1:nr, mat, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks)
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)

    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)

    axis(4, 1:nr, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)

    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)

    #-- ZOOM --#
    #----------#
    if (isTRUE(zoom)) {
        getOption("device")() 
        nD = dev.cur()
        dev.set(dev.prev())
        DEV.x = grDevices::dev.cur()

        zone = FALSE
        repeat {
            grDevices::dev.set(DEV.x)

##            loc = .Internal(locator(1, type = "n"))
            loc = locator(1, type = "n")
            x1 = round(loc[[1]] - 0.5, 0) + 0.5
            y1 = round(loc[[2]] - 0.5, 0) + 0.5
            flag = loc[[3]]

            if (flag == 0 & zone == TRUE) break

            if (flag != 0 & zone == TRUE) {
                rect(xleft.old, ybottom.old, xright.old, ytop.old, border = "white")
                points(x1.old, y1.old, type = "p", pch = 3, 
                cex = 2, col = "white")
            }

            if (flag == 0 & zone == FALSE) {
                break
            }
            else {
                if (x1 < 0) x1 = 0.5
                if (x1 > nc) x1 = nc + 0.5
                if (y1 < 0) y1 = 0.5
                if (y1 > nr) y1 = nr + 0.5
                x1.old = x1
                y1.old = y1

                points(x1, y1, type = "p", pch = 3, cex = 2)

##                loc = .Internal(locator(1, type = "n"))
				loc = locator(1, type = "n")
                x2 = round(loc[[1]] - 0.5, 0) + 0.5
                y2 = round(loc[[2]] - 0.5, 0) + 0.5

                if (x2 < 0) x2 = 0.5
                if (x2 > nc) x2 = nc + 0.5
                if (y2 < 0) y2 = 0.5
                if (y2 > nr) y2 = nr + 0.5

                xleft.old = min(x1, x2) 
                xright.old = max(x1, x2) 
                ybottom.old = min(y1, y2) 
                ytop.old = max(y1, y2) 

                rect(xleft.old, ybottom.old, xright.old, ytop.old)
                zone = TRUE
            }

            grDevices::dev.set(nD) 

            if (zone) {
                xleft = xleft.old + 0.5
                ybottom = ybottom.old + 0.5
                xright = xright.old - 0.5
                ytop = ytop.old - 0.5
                nr.zoom = length(xleft:xright)
                nc.zoom = length(ybottom:ytop)
                mat.zoom = matrix(mat[xleft:xright, ybottom:ytop], 
                nrow = nr.zoom, ncol = nc.zoom)
                rlab.zoom = rownames(mat)[xleft:xright]
                clab.zoom = colnames(mat)[ybottom:ytop]
                cexRow = min(1, 0.2 + 1/log10(nr.zoom))
                cexCol = min(1, 0.2 + 1/log10(nc.zoom))

                lmat = matrix(c(1, 0, 0, 2), 2, 2, byrow = TRUE)
                if (isTRUE(csc)) lmat = matrix(c(1, 0, 0, 2, 0, 3), 3, 2, byrow = TRUE)
                if (isTRUE(rsc)) lmat = matrix(c(1, 0, 0, 0, 2, 3), 2, 3, byrow = TRUE)
                if (isTRUE(csc) && isTRUE(rsc)) lmat = matrix(c(1, 0, 0, 0, 0, 2, 0, 3, 4), 3, 3, byrow = TRUE)
                layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

                # layout 1
                par(mar = c(5, 4, 2, 1), cex = 0.75)				
                image(z, col = col, xaxt = "n", yaxt = "n")
                box()
                par(usr = c(0, 1, 0, 1))
                lv = c(min.breaks, (3*min.breaks + max.breaks)/4, (min.breaks + max.breaks)/2,
                      (3*max.breaks + min.breaks)/4, max.breaks)
                xv = (as.numeric(lv) - min.mat) / (max.mat - min.mat)
                axis(1, at = xv, labels = round(lv, 2))
                title("Color key", font.main = 1)

                # layout 2
                if (isTRUE(csc)) {
                    par(mar = c(0.5, 0, 0, margins[2] + 1))
                    image(cbind(1:nr.zoom), col = ColSideColors[xleft:xright], axes = FALSE)	
                }
				
                # layout 3
                if (isTRUE(rsc)) {
                    par(mar = c(margins[1] + 1, 0, 0, 0.5))
                    image(rbind(1:nc.zoom), col = RowSideColors[ybottom:ytop], axes = FALSE)	
                }
				
                # layout 4
                par(mar = c(margins[1] + 1, 0, 0, margins[2] + 1))
                image(1:nr.zoom, 1:nc.zoom, mat.zoom, col = col, 
                breaks = breaks, axes = FALSE, xlab = "", ylab = "")

                axis(1, 1:nr.zoom, labels = rlab.zoom, las = 2, 
                line = -0.5, tick = 0, cex.axis = cexRow)

                if (!is.null(xlab)) 
                    mtext(xlab, side = 1, line = margins[1] - 1.25)

                axis(4, 1:nc.zoom, labels = clab.zoom, las = 2, 
                line = -0.5, tick = 0, cex.axis = cexCol)

                if (!is.null(ylab)) 
                    mtext(ylab, side = 4, line = margins[2] - 1.25)
            }

            grDevices::dev.set(DEV.x)
        }
    }

    par(op)
    res = list(simMat = simMat, labCol = labCol, labRow = labRow)
	
    if (any(dendrogram == "both") || any(dendrogram == "row")) {
        res$rowInd = rowInd
        res$ddr = ddr		
    }
	
    if (any(dendrogram == "both") || any(dendrogram == "column")) {
        res$colInd = colInd
        res$ddc = ddc
    }
	
    return(invisible(res))
}


# -------------------------------------
# cim for rcc object 
# -------------------------------------


cim.rcc <-
function(object, 
         comp = 1, 
         X.names = NULL, 
         Y.names = NULL, 
         ...) 
{

    p = ncol(object$X)
    q = ncol(object$Y)
    ncomp = object$ncomp	
	
    if (length(comp) == 1) {
	    if (is.null(comp) || !is.numeric(comp) || comp <= 0 || comp > ncomp)
            stop("invalid value for 'comp'.")
    }
	
    if (length(comp) > 1) {
        if(length(comp) > ncomp) 
            stop("the length of 'comp' must be smaller or equal than ", ncomp, ".")
        if (!is.numeric(comp) || any(comp < 1))
            stop("invalid vector for 'comp'.")
        if (any(comp > ncomp)) 
            stop("the elements of 'comp' must be smaller or equal than ", ncomp, ".")
    }
	
    if (length(X.names) != p & !is.null(X.names))
        stop("'X.names' must be a character vector of length ", p, ".")
		
    if (length(Y.names) != q & !is.null(Y.names))
        stop("'Y.names' must be a character vector of length ", q, ".")

    comp = round(comp)

    if (is.null(X.names)) X.names = object$names$X
    if (is.null(Y.names)) Y.names = object$names$Y

    bisect = object$variates$X[, comp] + object$variates$Y[, comp]
    cord.X = cor(object$X, bisect, use = "pairwise")
    cord.Y = cor(object$Y, bisect, use = "pairwise")
    simMat = as.matrix(cord.X %*% t(cord.Y))

    if (ncol(simMat) < nrow(simMat)) {
        simMat = t(simMat)
        aux.names = X.names
        X.names = Y.names
        Y.names = aux.names
    }  
	
    result = cim(simMat, labRow = X.names, labCol = Y.names, ...)
    return(invisible(result))
}

# --------------------------------
# CIM for sPLS object 
# ---------------------------------

cim.spls <-
function(object, 
         comp = 1, 
         X.names = NULL, 
         Y.names = NULL, 
         keep.var = TRUE, 
         ...) 
{

    # validation des arguments #
	#--------------------------#
    dim = object$ncomp	
    	
    if (length(comp) == 1) {
	    if (is.null(comp) || !is.numeric(comp) || comp <= 0)
            stop("invalid value for 'comp'.")
        if (comp > dim) 
            stop("'comp' must be smaller or equal than ", dim, ".")
    }
    	
    if (length(comp) > 1) {
        if(length(comp) > dim) 
            stop("the length of 'comp' must be smaller or equal than ", dim, ".")
        if (!is.numeric(comp) || any(comp < 1))
            stop("invalid vector for 'comp'.")
        if (any(comp > dim)) 
            stop("the elements of 'comp' must be smaller or equal than ", dim, ".")
    }
    	
    p = ncol(object$X)
    q = ncol(object$Y)
		
    if (length(X.names) != p & !is.null(X.names))
        stop("'X.names' must be a character vector of length ", p, ".")
		
    if (length(Y.names) != q & !is.null(Y.names))
        stop("'Y.names' must be a character vector of length ", q, ".")
		
    comp = round(comp)
     
    # Calcul de la matrice des associations entre les variables X et Y #
    #------------------------------------------------------------------#
    if (isTRUE(keep.var)) {
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
		
        if (is.null(X.names)) X.names = object$names$X[keep.X]
	    if (is.null(Y.names)) Y.names = object$names$Y[keep.Y]
    }
    else {
        if (object$mode == "canonical") {
            cord.X = cor(object$X, object$variates$X[, comp], use = "pairwise")
            cord.Y = cor(object$Y, object$variates$Y[, comp], use = "pairwise")
        }
        else {
            cord.X = cor(object$X, object$variates$X[, comp], use = "pairwise")
            cord.Y = cor(object$Y, object$variates$X[, comp], use = "pairwise")
        }
     		
        if (is.null(X.names)) X.names = object$names$X
        if (is.null(Y.names)) Y.names = object$names$Y
    }
     	
    simMat = cord.X %*% t(cord.Y)
    if (ncol(simMat) < nrow(simMat)) {
        simMat = t(simMat)
        aux.names = X.names
        X.names = Y.names
        Y.names = aux.names
    }  
     	
    result = cim(simMat, labRow = X.names, labCol = Y.names, ...)
    return(invisible(result))
}

# ---------------------------------
# CIM for PLS object 
# ---------------------------------
cim.pls <-
function(object, 
         comp = 1, 
         X.names = NULL, 
         Y.names = NULL, 
         ##keep.var = TRUE, 
         ...) 
{

    # validation des arguments #
	#--------------------------#
    dim = object$ncomp	
    	
    if (length(comp) == 1) {
	    if (is.null(comp) || !is.numeric(comp) || comp <= 0)
            stop("invalid value for 'comp'.")
        if (comp > dim) 
            stop("'comp' must be smaller or equal than ", dim, ".")
    }
    	
    if (length(comp) > 1) {
        if(length(comp) > dim) 
            stop("the length of 'comp' must be smaller or equal than ", dim, ".")
        if (!is.numeric(comp) || any(comp < 1))
            stop("invalid vector for 'comp'.")
        if (any(comp > dim)) 
            stop("the elements of 'comp' must be smaller or equal than ", dim, ".")
    }
    	
    p = ncol(object$X)
    q = ncol(object$Y)
		
    if (length(X.names) != p & !is.null(X.names))
        stop("'X.names' must be a character vector of length ", p, ".")
		
    if (length(Y.names) != q & !is.null(Y.names))
        stop("'Y.names' must be a character vector of length ", q, ".")
		
    comp = round(comp)
     
    # Calcul de la matrice des associations entre les variables X et Y #
    #------------------------------------------------------------------#
#    if (isTRUE(keep.var)) {
#        keep.X = apply(abs(object$loadings$X), 1, sum) > 0
#        keep.Y = apply(abs(object$loadings$Y), 1, sum) > 0

#        if (object$mode == "canonical") {
#            cord.X = cor(object$X[, keep.X], object$variates$X[, comp], 
#                     use = "pairwise")
#            cord.Y = cor(object$Y[, keep.Y], object$variates$Y[, comp], 
#                     use = "pairwise")
#        }
#        else {
#            cord.X = cor(object$X[, keep.X], object$variates$X[, comp], 
#                     use = "pairwise")
#            cord.Y = cor(object$Y[, keep.Y], object$variates$X[, comp], 
#                     use = "pairwise")
#        }		
		
#        if (is.null(X.names)) X.names = object$names$X[keep.X]
#	    if (is.null(Y.names)) Y.names = object$names$Y[keep.Y]
#    }
#    else {
        if (object$mode == 'canonical'){
            cord.X = cor(object$X, object$variates$X[, comp], use = 'pairwise')
            cord.Y = cor(object$Y, object$variates$Y[, comp], use = 'pairwise')
        }
        else{
            cord.X = cor(object$X, object$variates$X[, comp], use = 'pairwise')
            cord.Y = cor(object$Y, object$variates$X[, comp], use = 'pairwise')
        }
     		
        if (is.null(X.names)) X.names = object$names$X
        if (is.null(Y.names)) Y.names = object$names$Y
#    }
     	
    simMat = cord.X %*% t(cord.Y)
    if (ncol(simMat) < nrow(simMat)) {
        simMat = t(simMat)
        aux.names = X.names
        X.names = Y.names
        Y.names = aux.names
    }  
     	
    result = cim(simMat, labRow = X.names, labCol = Y.names, ...)
    return(invisible(result))
}


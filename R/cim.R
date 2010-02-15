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


cim <-
function(...) UseMethod("cim")

cim.default <-
function(mat, 
         breaks, 
         col = jet.colors, 
         distfun = dist, 
         hclustfun = hclust,
         labRow = NULL, 
         labCol = NULL, 
         symkey = TRUE, 
         zoom = FALSE, 
         main = NULL, 
         xlab = NULL, 
         ylab = NULL, 
         keysize = 1, 
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

    mat = as.matrix(mat)

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

    Rowv = rowMeans(mat)
    hcr = hclustfun(distfun(mat))
    ddr = as.dendrogram(hcr)
    ddr = reorder(ddr, Rowv)
    rowInd = order.dendrogram(ddr)
        
    Colv = colMeans(mat)
    hcc = hclustfun(distfun(t(mat)))
    ddc = as.dendrogram(hcc)
    ddc = reorder(ddc, Colv)
    colInd = order.dendrogram(ddc)
        
    mat = mat[rowInd, colInd]

    if (is.null(labRow)) 
        labRow = if (is.null(rownames(mat))) (1:nr)[rowInd] else rownames(mat)
    else labRow = labRow[rowInd]

    if (is.null(labCol)) 
        labCol = if (is.null(colnames(mat))) (1:nc)[colInd] else colnames(mat)
    else labCol = labCol[colInd]
	
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

    if (missing(lhei) || is.null(lhei)) 
        lhei = c(keysize, 4)

    if (missing(lwid) || is.null(lwid)) 
        lwid = c(keysize, 4)
 
    lmat = matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)

    if (length(lhei) != nrow(lmat)) 
        stop("'lhei' must have length = 2.")

    if (length(lwid) != ncol(lmat)) 
        stop("'lwid' must have length = 2.")

    if (isTRUE(zoom)) {
        getOption("device")("xpos" = 65)
    }

    op = par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    #-- layout 1 --#
    par(mar = c(5, 4, 2, 1), cex = 0.75)

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
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")

    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])

    #-- layout 3 --#
    par(mar = c(margins[1], 0, 0, 0))
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")

    #-- layout 4 --#
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

            loc = .Internal(locator(1, type = "n"))
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

                loc = .Internal(locator(1, type = "n"))
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
    invisible(list(rowInd = rowInd, colInd = colInd, ddc = ddc, ddr = ddr,
              labCol = labCol, labRow = labRow))
}

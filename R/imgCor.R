# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
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


imgCor <-
function(X, 
         Y,  
         type = c("combine", "separate"), 
         col = color.jet, 
         X.names = TRUE, 
         Y.names = TRUE,		 
         XsideColor = "blue",
         YsideColor = "red",
         symkey = TRUE, 
         keysize = 1, 
         interactive.dev = TRUE,		 
         cexRow = NULL, 
         cexCol = NULL, 
         margins = c(5, 5), 
         lhei = NULL, 
         lwid = NULL) 
{

    #-- validation des arguments --#
    if (length(dim(X)) != 2 || length(dim(Y)) != 2) 
        stop("'X' and/or 'Y' must be a numeric matrix.")
    
    X = as.matrix(X)
    Y = as.matrix(Y)
     
    if (!is.numeric(X) || !is.numeric(Y)) 
        stop("'X' and/or 'Y' must be a numeric matrix.")
		
    if (class(col) == "function") breaks = seq(-1, 1, length = 26)
    else breaks = seq(-1, 1, length = length(col) + 1)
     
    p = ncol(X)
    q = ncol(Y)	
	
    if (!is.logical(X.names)) {
        if (!is.vector(X.names) || (length(X.names) != p))
            stop("'X.names' must be a character vector of length ", p, ".")
    }
    else {
        if (isTRUE(X.names)) X.names = NULL else X.names = rep(" ", p)	
    }
	
    if (!is.logical(Y.names)) {
        if (!is.vector(Y.names) || (length(Y.names) != q))
            stop("'Y.names' must be a character vector of length ", q, ".")
    }
    else {
        if (isTRUE(Y.names)) Y.names = NULL else Y.names = rep(" ", q)	
    }
	
    if (!is.null(X.names)) colnames(X) = X.names
    if (!is.null(Y.names)) colnames(Y) = Y.names
    type = match.arg(type)
     
    # representation de la matrice de correlation de #
    # la concatenationdes variables X et Y, [X Y]    #
    #------------------------------------------------#

    if (type == "combine") {		
        matcor = cor(cbind(X, Y), use = "pairwise")
        matcor = t(matcor[(p + q):1, ])
		
        ColSideColors = c(rep(YsideColor, q), rep(XsideColor, p))
        RowSideColors = c(rep(XsideColor, p), rep(YsideColor, q))
		
        cim(matcor, col = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "[X,Y] correlation matrix",
            ColSideColors = ColSideColors,
            RowSideColors = RowSideColors,
            cexRow = if(is.null(cexRow)) min(1, 0.2 + 1/log10(p + q)) else cexRow, 
            cexCol = if(is.null(cexCol)) min(1, 0.2 + 1/log10(p + q)) else cexCol,			
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks)
    }
     
    # representation des matrices de correlation de #
    # X, Y et entre X et Y                          #
    #-----------------------------------------------#
    if (type == "separate") {
        Xcor = cor(X, use = "pairwise")
        Ycor = cor(Y, use = "pairwise")
        XYcor = cor(X, Y, use = "pairwise")

        Xcor = t(Xcor[p:1, ])
        Ycor = t(Ycor[q:1, ])
        XYcor = XYcor[, q:1]
		
        if (interactive.dev) {   
            cim(Xcor, col = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "X correlation matrix",			
            cexRow = if(is.null(cexRow)) min(1, 0.2 + 1/log10(p)) else cexRow, 
            cexCol = if(is.null(cexRow)) min(1, 0.2 + 1/log10(p)) else cexRow, 
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks)
			
            devAskNewPage(TRUE)
            cim(Ycor, col = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "Y correlation matrix",			
            cexRow = if(is.null(cexCol)) min(1, 0.2 + 1/log10(q)) else cexCol, 
            cexCol = if(is.null(cexCol)) min(1, 0.2 + 1/log10(q)) else cexCol, 
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks)

		    cim(XYcor, col = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "XY correlation matrix",			
            cexRow = if(is.null(cexRow)) min(1, 0.2 + 1/log10(p)) else cexRow, 
            cexCol = if(is.null(cexCol)) min(1, 0.2 + 1/log10(q)) else cexCol,
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks)
        }
		else {
            getOption("device")("xpos" = 0, "ypos" = 0)
            cim(XYcor, col = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "XY correlation matrix",			
            cexRow = if(is.null(cexRow)) min(1, 0.2 + 1/log10(p)) else cexRow, 
            cexCol = if(is.null(cexCol)) min(1, 0.2 + 1/log10(q)) else cexCol,
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks)
			
            getOption("device")("xpos" = 34, "ypos" = 34)
            cim(Ycor, col = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "Y correlation matrix",			
            cexRow = if(is.null(cexCol)) min(1, 0.2 + 1/log10(q)) else cexCol, 
            cexCol = if(is.null(cexCol)) min(1, 0.2 + 1/log10(q)) else cexCol,
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks)
			
            getOption("device")("xpos" = 64, "ypos" = 64)
            cim(Xcor, col = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "X correlation matrix",			
            cexRow = if(is.null(cexRow)) min(1, 0.2 + 1/log10(p)) else cexRow, 
            cexCol = if(is.null(cexRow)) min(1, 0.2 + 1/log10(p)) else cexRow,
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks)
		}
    }	
}

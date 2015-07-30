# Copyright (C) 2009 
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
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

# -----------------------------------------
# --- small example but see help file
# ----------------------------------------
#library(mixOmics)
#data(liver.toxicity)
#X <- liver.toxicity$gene
#Y <- liver.toxicity$clinic
#imgCor(X=X[,1:10],Y=Y[,1:4])


imgCor <-
function(X, 
         Y,  
         type = c("combine", "separate"), 
         col = color.jet(25),
         X.var.names = TRUE, 
         Y.var.names = TRUE,		 
         x.sideColors = "blue",
         y.sideColors = "red",
         symkey = TRUE, 
         keysize = c(1,1),
         interactive.dev = TRUE,		 
         row.cex = NULL, 
         col.cex = NULL, 
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
	
    if (!is.logical(X.var.names)) {
        if (!is.vector(X.var.names) || (length(X.var.names) != p))
            stop("'X.var.names' must be a character vector of length ", p, ".")
    }
    else {
        if (isTRUE(X.var.names)) X.var.names = NULL else X.var.names = rep(" ", p)	
    }
	
    if (!is.logical(Y.var.names)) {
        if (!is.vector(Y.var.names) || (length(Y.var.names) != q))
            stop("'Y.var.names' must be a character vector of length ", q, ".")
    }
    else {
        if (isTRUE(Y.var.names)) Y.var.names = NULL else Y.var.names = rep(" ", q)	
    }
	
    if (!is.null(X.var.names)) colnames(X) = X.var.names
    if (!is.null(Y.var.names)) colnames(Y) = Y.var.names
    type = match.arg(type)
     
    # representation de la matrice de correlation de #
    # la concatenationdes variables X et Y, [X Y]    #
    #------------------------------------------------#

    if (type == "combine") {		
        matcor = cor(cbind(X, Y), use = "pairwise")
        #matcor = t(matcor[(p + q):1, ])
        matcor=matcor[(p + q):1,]
		
        ColSideColors = c(rep(x.sideColors, p), rep(y.sideColors, q))
        RowSideColors = c(rep(y.sideColors, q), rep(x.sideColors, p))
		
        cim(matcor, color = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "[X,Y] correlation matrix",
            col.sideColors = ColSideColors,
            row.sideColors = RowSideColors,
            row.cex = if(is.null(row.cex)) min(1, 0.2 + 1/log10(p + q)) else row.cex,
            col.cex = if(is.null(col.cex)) min(1, 0.2 + 1/log10(p + q)) else col.cex,			
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks,cluster="none")
    }
     
    # representation des matrices de correlation de #
    # X, Y et entre X et Y                          #
    #-----------------------------------------------#
    if (type == "separate") {
        Xcor = cor(X, use = "pairwise")
        Ycor = cor(Y, use = "pairwise")
        XYcor = cor(X, Y, use = "pairwise")

        Xcor = Xcor[p:1,]#t(Xcor[p:1, ])
        Ycor = Ycor[q:1,]#t(Ycor[q:1, ])
        XYcor = XYcor[p:1,]
		
        if (interactive.dev) {   
            cim(Xcor, color = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "X correlation matrix",			
            row.cex = if(is.null(row.cex)) min(1, 0.2 + 1/log10(p)) else row.cex, 
            col.cex = if(is.null(row.cex)) min(1, 0.2 + 1/log10(p)) else row.cex, 
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks,cluster="none")
			
            devAskNewPage(TRUE)
            cim(Ycor, color = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "Y correlation matrix",			
            row.cex = if(is.null(col.cex)) min(1, 0.2 + 1/log10(q)) else col.cex, 
            col.cex = if(is.null(col.cex)) min(1, 0.2 + 1/log10(q)) else col.cex, 
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks,cluster="none")

		    cim(XYcor, color = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "XY correlation matrix",			
            row.cex = if(is.null(row.cex)) min(1, 0.2 + 1/log10(p)) else row.cex, 
            col.cex = if(is.null(col.cex)) min(1, 0.2 + 1/log10(q)) else col.cex,
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks,cluster="none")
        }
		else {
            getOption("device")("xpos" = 0, "ypos" = 0)
            cim(XYcor, color = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "XY correlation matrix",			
            row.cex = if(is.null(row.cex)) min(1, 0.2 + 1/log10(p)) else row.cex, 
            col.cex = if(is.null(col.cex)) min(1, 0.2 + 1/log10(q)) else col.cex,
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks,cluster="none")
			
            getOption("device")("xpos" = 34, "ypos" = 34)
            cim(Ycor, color = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "Y correlation matrix",			
            row.cex = if(is.null(col.cex)) min(1, 0.2 + 1/log10(q)) else col.cex, 
            col.cex = if(is.null(col.cex)) min(1, 0.2 + 1/log10(q)) else col.cex,
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks,cluster="none")
			
            getOption("device")("xpos" = 64, "ypos" = 64)
            cim(Xcor, color = col, dendrogram = "none",
            labRow = NULL, labCol = NULL,		 
            symkey = symkey, keysize = keysize, 
            main = "X correlation matrix",			
            row.cex = if(is.null(row.cex)) min(1, 0.2 + 1/log10(p)) else row.cex, 
            col.cex = if(is.null(row.cex)) min(1, 0.2 + 1/log10(p)) else row.cex,
            margins = margins, lhei = lhei, lwid = lwid,
            breaks = breaks,cluster="none")
		}
    }	
}

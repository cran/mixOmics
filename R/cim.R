# Copyright (C) 2015
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# This program is free software; you can redistribute it and/or
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


# --------------------------------------------
# CIM for objects "pca","spca","ipca","sipca","mlsplsda","splsda","plsda","rcc","pls","spls","mlspls"
# --------------------------------------------
cim <-
  function(mat,
           color = NULL,
           row.names = TRUE,
           col.names = TRUE,
           row.sideColors = NULL,
           col.sideColors = NULL,
           row.cex = NULL,
           col.cex = NULL,
           cluster = "both",
           dist.method = c("euclidean", "euclidean"),
           clust.method = c("complete", "complete"),
           cut.tree = c(0, 0),
           transpose = FALSE,
           comp = NULL,
           symkey = TRUE, 
           keysize = c(1, 1),            
           zoom = FALSE, 
           main = NULL,
           xlab = NULL,
           ylab = NULL,
           margins = c(5, 5),
           lhei = NULL,
           lwid = NULL,
           sample.names = TRUE,
           var.names = TRUE,
           sample.sideColors = NULL,
           var.sideColors = NULL,
           center = TRUE,
           scale = FALSE,
           X.var.names = TRUE, 
           Y.var.names = TRUE,
           x.sideColors = NULL,
           y.sideColors = NULL,
           mapping = "XY",
           legend=NULL,
           ...)
{
    class.object=class(mat)
    
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- check that the user did not enter extra arguments
    arg.call = match.call()
    user.arg = names(arg.call)[-1]
    
    err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())), 
                   error = function(e) e)
    
    if ("simpleError" %in% class(err))
      stop(err[[1]], ".", call. = FALSE)
    
    #     function.arg = c(names(mget(names(formals()), sys.frame(sys.nframe()))))
    #     not.arg = !(user.arg %in% function.arg)
    #     
    #     if (any(not.arg)) {
    #       unused.arg = user.arg[not.arg]
    #       not.arg = which(not.arg) + 1
    #       output = rep("", length(not.arg))
    #       
    #       for (i in 1:length(not.arg)) {
    #         output[i] = paste0(unused.arg[i], " = ", arg.call[[not.arg[i]]])
    #       }
    #       
    #       output = paste0("(", paste(output, collapse = ", "), ").")
    #       msg = "unused argument "
    #       if (length(not.arg) > 1) msg = "unused arguments "  
    #       stop(msg, output, call. = FALSE)
    #     }
    
    
    
    
    
    #-- color
    if (is.null(color)) 
      color = color.spectral(25)
    
    
    
    #-- cluster
    choices = c("both", "row", "column", "none")
    cluster = choices[pmatch(cluster, choices)]
    
    if (is.na(cluster)) 
      stop("'cluster' should be one of 'both', 'row', 'column' or 'none'.", 
           call. = FALSE)
    
    #-- cluster method
    if (!is.character(clust.method) | length(as.vector(clust.method)) != 2)
      stop("'clust.method' must be a character vector of length 2.", call. = FALSE)
    
    choices = c("ward.D", "single", "complete", "average", "mcquitty",
                "median", "centroid")
    clust.method = choices[c(pmatch(clust.method[1], choices),
                             pmatch(clust.method[2], choices))]
    
    if (any(is.na(clust.method))) 
      stop("invalid clustering method.", call. = FALSE)
    
    #-- distance method
    if (!is.character(dist.method) | length(as.vector(dist.method)) != 2)
      stop("'dist.method' must be a character vector of length 2.", call. = FALSE)
    
    choices = c("euclidean", "correlation", "maximum", "manhattan", 
                "canberra", "binary", "minkowski")
    dist.method = choices[c(pmatch(dist.method[1], choices),
                            pmatch(dist.method[2], choices))]
    
    if (any(is.na(dist.method))) 
      stop("invalid distance method.", call. = FALSE)
    
    
    #-- checking general input arguments ---------------------------------------#
    #---------------------------------------------------------------------------#  
    
    #-- color
    if (any(!sapply(color, function(color) { 
      tryCatch(is.matrix(col2rgb(color)), error = function(e) FALSE) }))) 
      stop("'color' must be a character vector of recognized colors.", 
           call. = FALSE)
    
    #-- row.sideColors
    if (any(!sapply(row.sideColors, function(row.sideColors) { 
      tryCatch(is.matrix(col2rgb(row.sideColors)), error = function(e) FALSE) }))) 
      stop("color names for vertical side bar must be a character vector of recognized colors.", 
           call. = FALSE)
    
    #-- col.sideColors
    if (any(!sapply(col.sideColors, function(col.sideColors) { 
      tryCatch(is.matrix(col2rgb(col.sideColors)), error = function(e) FALSE) }))) 
      stop("color names for horizontal side bar must be a character vector of recognized colors.", 
           call. = FALSE)
    
    #-- row.cex
    if (!is.null(row.cex)) {
      if (!is.numeric(row.cex) || length(row.cex) != 1)
        stop("'row.cex' must be a numerical value.", call. = FALSE)
    }
    
    #-- col.cex
    if (!is.null(col.cex)) {
      if (!is.numeric(col.cex) || length(col.cex) != 1)
        stop("'col.cex' must be a numerical value.", call. = FALSE)
    }
    
    #-- transpose
    if (!is.logical(transpose))
      stop("'transpose' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)
    
    #-- cut.tree
    if (!is.numeric(cut.tree) || length(cut.tree) != 2) 
      stop("'cut.tree' must be a numeric vector of length 2.",
           call. = FALSE)
    else {
      if (!(all(0 <= cut.tree & cut.tree <= 1)))
        stop("Components of 'cut.tree' must be between 0 and 1.",
             call. = FALSE) 
    }
    
    #-- keysize
    if (length(keysize) != 2 || any(!is.finite(keysize))) 
      stop("'keysize' must be a numeric vector of length 2.",
           call. = FALSE)
    
    #-- zoom
    if (!is.logical(zoom))
      stop("'zoom' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)
    
    #-- margins
    if (!is.numeric(margins) || length(margins) != 2) 
      stop("'margins' must be a numeric vector of length 2.",
           call. = FALSE)
    
    #-- symkey
    if (!is.logical(symkey))
      stop("'symkey' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)
    
    #-- lhei
    if (!is.null(lhei)) {
      if (is.null(col.sideColors)) {
        if (length(lhei) != 2 | !is.numeric(lhei) | any(is.na(lhei))) 
          stop("'lhei' must be a numeric vector of length 2.",
               call. = FALSE)
      }
      else {
        if (length(lhei) != 3 | !is.numeric(lhei) | any(is.na(lhei))) 
          stop("'lhei' must be a numeric vector of length 3.",
               call. = FALSE)
      }
    }
    
    #-- lwid
    if (!is.null(lwid)) {
      if (is.null(row.sideColors)) {
        if (length(lwid) != 2 | !is.numeric(lwid) | any(is.na(lwid))) 
          stop("'lwid' must be a numeric vector of length 2.",
               call. = FALSE)
      }
      else {
        if (length(lwid) != 3 | !is.numeric(lwid) | any(is.na(lwid))) 
          stop("'lwid' must be a numeric vector of length 3.",
               call. = FALSE)
      }
    }
    
    #-- xlab
    xlab = as.graphicsAnnot(xlab)
    
    #-- ylab
    ylab = as.graphicsAnnot(ylab)
    
    #-- main
    main = as.graphicsAnnot(main)
    
    #-- end checking --#
    #------------------#
    object.list1=c("pca","spca","ipca","sipca","mlsplsda","splsda","plsda")
    object.list2=c("rcc")
    object.list3=c("pls","spls","mlspls")
    object.list=c("pca","spca","ipca","sipca","mlsplsda","splsda","plsda","rcc","pls","spls","mlspls")
    if(class.object[1] %in%  object.list)
      
    {p = ncol(mat$X)
     q = ncol(mat$Y)
     n = nrow(mat$X)
     ncomp = mat$ncomp
     #-- comp
     if(is.null(comp))
     {comp=1:mat$ncomp}
     if (length(comp) > 1) {
       comp=unique(comp)
       if (!is.numeric(comp) || any(comp < 1))
         stop("invalid vector for 'comp'.", call. = FALSE)
       if (any(comp > ncomp)) 
         stop("the elements of 'comp' must be smaller or equal than ", ncomp, ".", 
              call. = FALSE)
     }
     
     if (length(comp) == 1) {
       if (is.null(comp) || !is.numeric(comp) || comp <= 0 || comp > ncomp)
         stop("invalid value for 'comp'.", call. = FALSE)
       comp=c(comp,comp)
     }
     
     
     
     comp = round(comp)
     
     #-- sample.names
     if (is.logical(sample.names)) {
       if (isTRUE(sample.names)) sample.names = mat$names$indiv else sample.names = rep("", n)
     }
     else {
       if (length(sample.names) != n)
         stop("'sample.names' must be a character vector of length ", n, ".", 
              call. = FALSE)
     }
     
     #-- var.names
     if (is.logical(var.names)) {
       if(isTRUE(var.names)) var.names = mat$names$X else var.names = rep("", p)
     }
     else {
       var.names = as.vector(var.names)
       if (length(var.names) != p)
         stop("'var.names' must be a character vector of length ", p, ".", 
              call. = FALSE)
     }
     
     #-- sample.sideColors
     if (!is.null(sample.sideColors)) {
       sample.sideColors = as.matrix(sample.sideColors)
       if (nrow(sample.sideColors) != n)
         stop("'sample.sideColors' must be a colors character vector (matrix) of length (nrow) ", n, ".", 
              call. = FALSE)
     }
     row.sideColors = sample.sideColors
     
     #-- var.sideColors
     if (!is.null(var.sideColors)) {
       var.sideColors = as.matrix(var.sideColors)
       if (nrow(var.sideColors) != p)
         stop("'var.sideColors' must be a colors character vector (matrix) of length (nrow) ", p, ".", 
              call. = FALSE)
     }
     col.sideColors = var.sideColors
     
     #-- X.var.names
     if (is.logical(X.var.names)) {
       if(isTRUE(X.var.names)) X.var.names = mat$names$X else X.var.names = rep("", p)
     }
     else {
       X.var.names = as.vector(X.var.names)
       if (length(X.var.names) != p)
         stop("'X.var.names' must be a character vector of length ", p, ".", 
              call. = FALSE)
     }
     
     #-- Y.var.names
     if (is.logical(Y.var.names)) {
       if(isTRUE(Y.var.names)) Y.var.names = mat$names$Y else Y.var.names = rep("", q)
     }
     else {
       Y.var.names = as.vector(Y.var.names)
       if (length(Y.var.names) != q)
         stop("'Y.var.names' must be a character vector of length ", q, ".", 
              call. = FALSE)
     }
     if( ! class.object[1] %in%  object.list1)
     {#-- x.sideColors
       if (!is.null(x.sideColors)) {
         x.sideColors = as.matrix(x.sideColors)
         if (nrow(x.sideColors) != p)
           stop("'x.sideColors' must be a colors character vector (matrix) of length (nrow) ", p, ".", 
                call. = FALSE)
       }
       
       #-- y.sideColors
       if (!is.null(y.sideColors)) {
         y.sideColors = as.matrix(y.sideColors)
         if (nrow(y.sideColors) != q)
           stop("'y.sideColors' must be a colors character vector (matrix) of length (nrow) ", q, ".", 
                call. = FALSE)
       }
       
       #-- mapping
       choices = c("XY", "X", "Y")
       mapping = choices[pmatch(mapping, choices)]
       
       if (is.na(mapping)) 
         stop("'mapping' should be one of 'XY', 'X' or 'Y'.", call. = FALSE)
       
       if (mapping == "XY") {
         row.sideColors = x.sideColors
         col.sideColors = y.sideColors
       }
       
       if (mapping == "X") {
         row.sideColors = sample.sideColors
         col.sideColors = x.sideColors
       }
       
       if (mapping == "Y") {
         row.sideColors = sample.sideColors
         col.sideColors = y.sideColors
       }
     }
     
     
     if(class.object[1] %in%  object.list1)
     {
       
       #-- clustering -------------------------------------------------------------#
       #---------------------------------------------------------------------------#
       if(class.object[1] %in%  c("splsda","plsda",'mlsplsda'))
       {
         if(class.object[1] %in%  c("splsda",'mlsplsda'))
           keep.X = apply(abs(mat$loadings$X[,comp]), 1, sum) > 0
         else
           keep.X = apply(abs(mat$loadings$X), 1, sum) > 0
         cord.X = cor(mat$X[, keep.X], mat$variates$X[, comp], use = "pairwise")
         X.mat = as.matrix(mat$variates$X[, comp])
       }
       else{
         if(class.object[1] %in%  c("spca","sipca"))
           keep.X = apply(abs(mat$rotation[,comp]), 1, sum) > 0
         else
           keep.X = apply(abs(mat$rotation), 1, sum) > 0
         cord.X = cor(mat$X[, keep.X], mat$x[, comp], use = "pairwise")
         X.mat = as.matrix(mat$x[, comp])
       }
       #-- cheking center and scale
       if (!is.logical(center)) {
         if (!is.numeric(center) || (length(center) != ncol(mat$X)))
           stop("'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.", 
                call. = FALSE)
       }
       if (!is.logical(scale)) {
         if (!is.numeric(scale) || (length(scale) != ncol(mat$X)))
           stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.", 
                call. = FALSE)
       }
       
       object = scale(mat$X[, keep.X], center = center, scale = scale)
       
       row.names = sample.names
       col.names = var.names[keep.X]
       
       if (!is.null(col.sideColors))
         col.sideColors = as.matrix(col.sideColors[keep.X, ])
       
       if ((cluster == "both") || (cluster == "row")) { 
         Rowv = rowMeans(X.mat)
         
         if (dist.method[1] == "correlation") 
           dist.mat = as.dist(1 - cor(t(as.matrix(X.mat)), method = "pearson"))
         else
           dist.mat = dist(X.mat, method = dist.method[1])
         
         hcr = hclust(dist.mat, method = clust.method[1])
         ddr = as.dendrogram(hcr)
         ddr = reorder(ddr, Rowv)
         rowInd = order.dendrogram(ddr)
         object = object[rowInd, ]
         row.names = row.names[rowInd]
         
         if (!is.null(row.sideColors)) 
           row.sideColors = as.matrix(row.sideColors[rowInd, ])
       }    
       
       if ((cluster == "both") || (cluster == "column")) {
         Colv = rowMeans(cord.X)
         
         if (dist.method[2] == "correlation") 
           dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
         else
           dist.mat = dist(cord.X, method = dist.method[2])
         
         hcc = hclust(dist.mat, method = clust.method[2])
         ddc = as.dendrogram(hcc)
         ddc = reorder(ddc, Colv)
         colInd = order.dendrogram(ddc)
         object = object[, colInd]
         col.names = col.names[colInd]
         
         if (!is.null(col.sideColors)) 
           col.sideColors = as.matrix(col.sideColors[colInd, ])
       }
       
       #-- calling the image.map function -----------------------------------------#
       #---------------------------------------------------------------------------#
       
       
       
       #-- output -----------------------------------------------------------------#
       #---------------------------------------------------------------------------#
       res = list(mat = object, row.names = row.names, col.names = col.names)
       
       if ((cluster == "both") || (cluster == "row")) {
         res$rowInd = rowInd
         res$ddr = ddr    
       }
       
       if ((cluster == "both") || (cluster == "column")) {
         res$colInd = colInd
         res$ddc = ddc
       }
       
       class(res) = paste("cim",class.object[1],sep="_")
       
     }
     else if(class.object[1] %in%  object.list2)
     {
       
       bisect = mat$variates$X[, comp] + mat$variates$Y[, comp]
       cord.X = cor(mat$X, bisect, use = "pairwise")
       cord.Y = cor(mat$Y, bisect, use = "pairwise")
       XY.mat = as.matrix(cord.X %*% t(cord.Y))
       
       #-- if mapping = "XY"
       if (mapping == "XY") {
         object = XY.mat
         row.names = X.var.names
         col.names = Y.var.names
         
         if ((cluster == "both") || (cluster == "row")) {
           #Rowv = rowMeans(XY.mat)
           Rowv = rowMeans(cord.X)
           
           if (dist.method[1] == "correlation") 
             dist.mat = as.dist(1 - cor(t(as.matrix(object)), method = "pearson"))
           else
             #dist.mat = dist(mat, method = dist.method[1])
             dist.mat = dist(cord.X, method = dist.method[1])
           
           hcr = hclust(dist.mat, method = clust.method[1])
           ddr = as.dendrogram(hcr)
           ddr = reorder(ddr, Rowv)
           rowInd = order.dendrogram(ddr)
           object = object[rowInd, ]
           row.names = row.names[rowInd]
           
           if (!is.null(row.sideColors)) 
             row.sideColors = as.matrix(row.sideColors[rowInd, ])
         }    
         
         if ((cluster == "both") || (cluster == "column")) {
           #Colv = colMeans(mat)
           Colv = rowMeans(cord.Y)
           
           if (dist.method[2] == "correlation") 
             dist.mat = as.dist(1 - cor(as.matrix(object), method = "pearson"))
           else
             #dist.mat = dist(t(mat), method = dist.method[2])
             dist.mat = dist(cord.Y, method = dist.method[2])
           
           hcc = hclust(dist.mat, method = clust.method[2])
           ddc = as.dendrogram(hcc)
           ddc = reorder(ddc, Colv)
           colInd = order.dendrogram(ddc)
           object = object[, colInd]
           col.names = col.names[colInd]
           
           if (!is.null(col.sideColors)) 
             col.sideColors = as.matrix(col.sideColors[colInd, ])
         }
       }
       
       #-- if mapping = "X"
       if (mapping == "X") {
         
         #-- cheking center and scale
         if (!is.logical(center)) {
           if (!is.numeric(center) || (length(center) != ncol(mat$X)))
             stop("'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.", 
                  call. = FALSE)
         }
         if (!is.logical(scale)) {
           if (!is.numeric(scale) || (length(scale) != ncol(mat$X)))
             stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.", 
                  call. = FALSE)
         }
         
         object = scale(mat$X, center = center, scale = scale)
         X.mat = as.matrix(mat$variates$X[, comp])
         row.names = sample.names
         col.names = X.var.names
         
         if ((cluster == "both") || (cluster == "row")) { 
           Rowv = rowMeans(X.mat)
           
           if (dist.method[1] == "correlation") 
             dist.mat = as.dist(1 - cor(t(as.matrix(X.mat)), method = "pearson"))
           else
             dist.mat = dist(X.mat, method = dist.method[1])
           
           hcr = hclust(dist.mat, method = clust.method[1])
           ddr = as.dendrogram(hcr)
           ddr = reorder(ddr, Rowv)
           rowInd = order.dendrogram(ddr)
           object = object[rowInd, ]
           row.names = row.names[rowInd]
           
           if (!is.null(row.sideColors)) 
             row.sideColors = as.matrix(row.sideColors[rowInd, ])
         }    
         
         if ((cluster == "both") || (cluster == "column")) {
           Colv = rowMeans(cord.X)
           
           if (dist.method[2] == "correlation") 
             dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
           else
             dist.mat = dist(cord.X, method = dist.method[2])
           
           hcc = hclust(dist.mat, method = clust.method[2])
           ddc = as.dendrogram(hcc)
           ddc = reorder(ddc, Colv)
           colInd = order.dendrogram(ddc)
           object = object[, colInd]
           col.names = col.names[colInd]
           
           if (!is.null(col.sideColors)) 
             col.sideColors = as.matrix(col.sideColors[colInd, ])
         }
       }
       
       #-- if mapping = "Y"
       if (mapping == "Y") {
         
         #-- cheking center and scale
         if (!is.logical(center)) {
           if (!is.numeric(center) || (length(center) != ncol(mat$Y)))
             stop("'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'Y'.", 
                  call. = FALSE)
         }
         if (!is.logical(scale)) {
           if (!is.numeric(scale) || (length(scale) != ncol(mat$Y)))
             stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'Y'.", 
                  call. = FALSE)
         }
         
         object = scale(mat$Y, center = center, scale = scale)
         Y.mat = as.matrix(mat$variates$Y[, comp])
         row.names = sample.names
         col.names = Y.var.names
         
         if ((cluster == "both") || (cluster == "row")) { 
           Rowv = rowMeans(Y.mat)
           
           if (dist.method[1] == "correlation") 
             dist.mat = as.dist(1 - cor(t(as.matrix(Y.mat)), method = "pearson"))
           else
             dist.mat = dist(Y.mat, method = dist.method[1])
           
           hcr = hclust(dist.mat, method = clust.method[1])
           ddr = as.dendrogram(hcr)
           ddr = reorder(ddr, Rowv)
           rowInd = order.dendrogram(ddr)
           object = object[rowInd, ]
           row.names = row.names[rowInd]
           
           if (!is.null(row.sideColors)) 
             row.sideColors = as.matrix(row.sideColors[rowInd, ])
         }    
         
         if ((cluster == "both") || (cluster == "column")) {
           Colv = rowMeans(cord.Y)
           
           if (dist.method[2] == "correlation") 
             dist.mat = as.dist(1 - cor(t(cord.Y), method = "pearson"))
           else
             dist.mat = dist(cord.Y, method = dist.method[2])
           
           hcc = hclust(dist.mat, method = clust.method[2])
           ddc = as.dendrogram(hcc)
           ddc = reorder(ddc, Colv)
           colInd = order.dendrogram(ddc)
           object = object[, colInd]
           col.names = col.names[colInd]
           
           if (!is.null(col.sideColors)) 
             col.sideColors = as.matrix(col.sideColors[colInd, ])
         }
       }
       
       #-- calling the image.map function -----------------------------------------#
       #---------------------------------------------------------------------------#
       
       
       
       #-- output -----------------------------------------------------------------#
       #---------------------------------------------------------------------------#
       res = list(mat = object, row.names = row.names, col.names = col.names)
       
       if ((cluster == "both") || (cluster == "row")) {
         res$rowInd = rowInd
         res$ddr = ddr    
       }
       
       if ((cluster == "both") || (cluster == "column")) {
         res$colInd = colInd
         res$ddc = ddc
       }
       
       class(res) = "cim_rcc"
       
     }
     else if(class.object[1] %in%  object.list3)
     {
       if(class.object[1] %in% c("spls","mlspls"))
       {
         keep.X = apply(abs(mat$loadings$X[,comp]), 1, sum) > 0
         keep.Y = apply(abs(mat$loadings$Y[,comp]), 1, sum) > 0}
       else
       {
         keep.X = apply(abs(mat$loadings$X), 1, sum) > 0
         keep.Y = apply(abs(mat$loadings$Y), 1, sum) > 0}
         
       
       
       if (mat$mode == "canonical") {
         cord.X = cor(mat$X[, keep.X], mat$variates$X[, comp], use = "pairwise")
         cord.Y = cor(mat$Y[, keep.Y], mat$variates$Y[, comp], use = "pairwise")
       }
       else {
         cord.X = cor(mat$X[, keep.X], mat$variates$X[, comp], use = "pairwise")
         cord.Y = cor(mat$Y[, keep.Y], mat$variates$X[, comp], use = "pairwise")
       }
       
       XY.mat = as.matrix(cord.X %*% t(cord.Y))
       
       #-- if mapping = "XY"
       if (mapping == "XY") {
         object = XY.mat
         row.names = X.var.names[keep.X]
         col.names = Y.var.names[keep.Y]
         
         if (!is.null(row.sideColors))
           row.sideColors = as.matrix(row.sideColors[keep.X, ])
         if (!is.null(col.sideColors))
           col.sideColors = as.matrix(col.sideColors[keep.Y, ])
         
         if ((cluster == "both") || (cluster == "row")) { 
           #Rowv = rowMeans(XY.mat)
           Rowv = rowMeans(cord.X)
           
           if (dist.method[1] == "correlation") 
             dist.mat = as.dist(1 - cor(t(as.matrix(object)), method = "pearson"))
           else
             #dist.mat = dist(mat, method = dist.method[1])
             dist.mat = dist(cord.X, method = dist.method[1])
           
           hcr = hclust(dist.mat, method = clust.method[1])
           ddr = as.dendrogram(hcr)
           ddr = reorder(ddr, Rowv)
           rowInd = order.dendrogram(ddr)
           object = object[rowInd, ]
           row.names = row.names[rowInd]
           
           if (!is.null(row.sideColors)) 
             row.sideColors = as.matrix(row.sideColors[rowInd, ])
         }    
         
         if ((cluster == "both") || (cluster == "column")) {
           #Colv = colMeans(mat)
           Colv = rowMeans(cord.Y)
           
           if (dist.method[2] == "correlation") 
             dist.mat = as.dist(1 - cor(as.matrix(object), method = "pearson"))
           else
             #dist.mat = dist(t(mat), method = dist.method[2])
             dist.mat = dist(cord.Y, method = dist.method[2])
           
           hcc = hclust(dist.mat, method = clust.method[2])
           ddc = as.dendrogram(hcc)
           ddc = reorder(ddc, Colv)
           colInd = order.dendrogram(ddc)
           object = object[, colInd]
           col.names = col.names[colInd]
           
           if (!is.null(col.sideColors)) 
             col.sideColors = as.matrix(col.sideColors[colInd, ])
         }
       }
       
       #-- if mapping = "X"
       if (mapping == "X") {
         
         #-- cheking center and scale
         if (!is.logical(center)) {
           if (!is.numeric(center) || (length(center) != ncol(mat$X)))
             stop("'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.", 
                  call. = FALSE)
         }
         if (!is.logical(scale)) {
           if (!is.numeric(scale) || (length(scale) != ncol(mat$X)))
             stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.", 
                  call. = FALSE)
         }
         
         object = scale(mat$X[, keep.X], center = center, scale = scale)
         X.mat = as.matrix(mat$variates$X[, comp])
         row.names = sample.names
         col.names = X.var.names[keep.X]
         
         if (!is.null(col.sideColors))
           col.sideColors = as.matrix(col.sideColors[keep.X, ])
         
         if ((cluster == "both") || (cluster == "row")) { 
           Rowv = rowMeans(X.mat)
           
           if (dist.method[1] == "correlation") 
             dist.mat = as.dist(1 - cor(t(as.matrix(X.mat)), method = "pearson"))
           else
             dist.mat = dist(X.mat, method = dist.method[1])
           
           hcr = hclust(dist.mat, method = clust.method[1])
           ddr = as.dendrogram(hcr)
           ddr = reorder(ddr, Rowv)
           rowInd = order.dendrogram(ddr)
           object = object[rowInd, ]
           row.names = row.names[rowInd]
           
           if (!is.null(row.sideColors)) 
             row.sideColors = as.matrix(row.sideColors[rowInd, ])
         }    
         
         if ((cluster == "both") || (cluster == "column")) {
           Colv = rowMeans(cord.X)
           
           if (dist.method[2] == "correlation") 
             dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
           else
             dist.mat = dist(cord.X, method = dist.method[2])
           
           hcc = hclust(dist.mat, method = clust.method[2])
           ddc = as.dendrogram(hcc)
           ddc = reorder(ddc, Colv)
           colInd = order.dendrogram(ddc)
           object = object[, colInd]
           col.names = col.names[colInd]
           
           if (!is.null(col.sideColors)) 
             col.sideColors = as.matrix(col.sideColors[colInd, ])
         }
       }
       
       #-- if mapping = "Y"
       if (mapping == "Y") {
         
         #-- cheking center and scale
         if (!is.logical(center)) {
           if (!is.numeric(center) || (length(center) != ncol(mat$Y)))
             stop("'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'Y'.", 
                  call. = FALSE)
         }
         if (!is.logical(scale)) {
           if (!is.numeric(scale) || (length(scale) != ncol(mat$Y)))
             stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'Y'.", 
                  call. = FALSE)
         }
         
         object = scale(mat$Y[, keep.Y], center = center, scale = scale)
         Y.mat = as.matrix(mat$variates$Y[, comp])
         row.names = sample.names
         col.names = Y.var.names[keep.Y]
         
         if (!is.null(col.sideColors))
           col.sideColors = as.matrix(col.sideColors[keep.Y, ])
         
         if ((cluster == "both") || (cluster == "row")) { 
           Rowv = rowMeans(Y.mat)
           
           if (dist.method[1] == "correlation") 
             dist.mat = as.dist(1 - cor(t(as.matrix(Y.mat)), method = "pearson"))
           else
             dist.mat = dist(Y.mat, method = dist.method[1])
           
           hcr = hclust(dist.mat, method = clust.method[1])
           ddr = as.dendrogram(hcr)
           ddr = reorder(ddr, Rowv)
           rowInd = order.dendrogram(ddr)
           object = object[rowInd, ]
           row.names = row.names[rowInd]
           
           if (!is.null(row.sideColors)) 
             row.sideColors = as.matrix(row.sideColors[rowInd, ])
         }    
         
         if ((cluster == "both") || (cluster == "column")) {
           Colv = rowMeans(cord.Y)
           
           if (dist.method[2] == "correlation") 
             dist.mat = as.dist(1 - cor(t(cord.Y), method = "pearson"))
           else
             dist.mat = dist(cord.Y, method = dist.method[2])
           
           hcc = hclust(dist.mat, method = clust.method[2])
           ddc = as.dendrogram(hcc)
           ddc = reorder(ddc, Colv)
           colInd = order.dendrogram(ddc)
           object = object[, colInd]
           col.names = col.names[colInd]
           
           if (!is.null(col.sideColors)) 
             col.sideColors = as.matrix(col.sideColors[colInd, ])
         }
       }
       
       #-- calling the image.map function -----------------------------------------#
       #---------------------------------------------------------------------------#
       
       
       #-- output -----------------------------------------------------------------#
       #---------------------------------------------------------------------------#
       res = list(mat = object, row.names = row.names, col.names = col.names)
       
       if (!is.null(sample.sideColors)) {
         res$sample.sideColors = sample.sideColors
       }
       
       if ((cluster == "both") || (cluster == "row")) {
         res$rowInd = rowInd
         res$ddr = ddr    
       }
       
       if ((cluster == "both") || (cluster == "column")) {
         res$colInd = colInd
         res$ddc = ddc
       }
       class(res) = paste("cim",class.object[1],sep="_")
     }}
    else
    {
      #-- mat
      isMat = tryCatch(is.matrix(mat), error = function(e) e)
      
      if ("simpleError" %in% class(isMat))
        stop(isMat[[1]], ".", call. = FALSE)
      
      if (!is.matrix(mat) || !is.numeric(mat))
        stop("'mat' must be a numeric matrix.", call. = FALSE)
      
      p = nrow(mat)
      q = ncol(mat)
      #-- row.names
      if (is.logical(row.names)) {
        if(isTRUE(row.names)) row.names = rownames(mat) else row.names = rep("", p)
      }
      else {
        row.names = as.vector(row.names)
        if (length(row.names) != p)
          stop("'row.names' must be a character vector of length ", p, ".", call. = FALSE)
      }
      
      #-- col.names
      if (is.logical(col.names)) {
        if(isTRUE(col.names)) col.names = colnames(mat) else col.names = rep("", q)
      }
      else {
        col.names = as.vector(col.names)
        if (length(col.names) != q)
          stop("'col.names' must be a character vector of length ", q, ".", call. = FALSE)
      }
      
      #-- row.sideColors
      if (!is.null(row.sideColors)) {
        row.sideColors = as.matrix(row.sideColors)
        if (nrow(row.sideColors) != p)
          stop("'row.sideColors' must be a colors character vector (matrix) of length (nrow) ", p, ".", 
               call. = FALSE)
      }
      
      #-- col.sideColors
      if (!is.null(col.sideColors)) {
        col.sideColors = as.matrix(col.sideColors)
        if (nrow(col.sideColors) != q)
          stop("'col.sideColors' must be a colors character vector (matrix) of length (nrow) ", q, ".", 
               call. = FALSE)
      }
      
      #-- clustering -------------------------------------------------------------#
      #---------------------------------------------------------------------------#
      object=mat
      if ((cluster == "both") || (cluster == "row")) {
        Rowv = rowMeans(mat)
        
        if (dist.method[1] == "correlation") 
          dist.mat = as.dist(1 - cor(t(as.matrix(mat)), method = "pearson"))
        else
          dist.mat = dist(mat, method = dist.method[1])
        
        hcr = hclust(dist.mat, method = clust.method[1])
        ddr = as.dendrogram(hcr)
        ddr = reorder(ddr, Rowv)
        rowInd = order.dendrogram(ddr)
        object = mat[rowInd, ]
        row.names = row.names[rowInd]
        
        if (!is.null(row.sideColors)) 
          row.sideColors = as.matrix(row.sideColors[rowInd, ])
      }      
      
      if ((cluster == "both") || (cluster == "column")) {
        Colv = colMeans(mat)
        
        if (dist.method[2] == "correlation") 
          dist.mat = as.dist(1 - cor(as.matrix(mat), method = "pearson"))
        else
          dist.mat = dist(t(mat), method = dist.method[2])
        
        hcc = hclust(dist.mat, method = clust.method[2])
        ddc = as.dendrogram(hcc)
        ddc = reorder(ddc, Colv)
        colInd = order.dendrogram(ddc)
        object = mat[, colInd]
        col.names = col.names[colInd]
        
        if (!is.null(col.sideColors)) 
          col.sideColors = as.matrix(col.sideColors[colInd, ])
      }
      
      #-- calling the image.map function -----------------------------------------#
      
      
      
      #-- output -----------------------------------------------------------------#
      #---------------------------------------------------------------------------#
      res = list(mat = object, row.names = row.names, col.names = col.names, 
                 row.sideColors = row.sideColors, col.sideColors = col.sideColors)
      
      if ((cluster == "both") || (cluster == "row")) {
        res$rowInd = rowInd
        res$ddr = ddr    
      }
      
      if ((cluster == "both") || (cluster == "column")) {
        res$colInd = colInd
        res$ddc = ddc
      }
      
      class(res) = "cim_default"
      
    }
    #---------------------------------------------------------------------------#
    imageMap(object,
             color = color,
             row.names = row.names,
             col.names = col.names,
             row.sideColors = row.sideColors,
             col.sideColors = col.sideColors,             
             row.cex = row.cex,
             col.cex = col.cex,
             cluster = cluster,
             ddr = ddr,
             ddc = ddc,
             cut.tree = cut.tree,
             transpose = transpose,
             symkey = symkey, 
             keysize = keysize,            
             zoom = zoom, 
             main = main,
             xlab = xlab,
             ylab = ylab,
             margins = margins,
             lhei = lhei,
             lwid = lwid)
    if (!is.null(legend))
    {if(is.null(legend$x)) legend$x = "topright"
     if(is.null(legend$bty)) legend$bty = "n"
     if (is.null(legend$cex)) legend$cex = 0.8
     if(class.object[1] %in%  c("splsda","plsda"))
     {
       if (is.null(legend$legend)) legend$legend = mat$names$Y
       
       #-- col
       if (is.null(legend$col)) {
         if (!is.null(row.sideColors) & !is.null(sample.sideColors)) 
           legend$col = unique(as.matrix(sample.sideColors[order(map(mat$ind.mat)), 1]))
       }
       
     }
     else if(class.object[1] %in%  c("mlsplsda"))
     {
       if (is.null(legend$legend) && is.null(legend$col)) {
         if (ncol(mat$design) >= 2) {
           df = data.frame(mat$design[, 2], sample.sideColors[, 1])
           df = unique(df)
           legend$legend = as.character(df[, 1])
           legend$col = as.character(df[, 2])
         }
         if (ncol(mat$design) == 3) {
           df = data.frame(mat$design[, 3], sample.sideColors[, 2])
           df = unique(df)
           legend$legend = c(legend$legend, as.character(df[, 1]))
           legend$col = c(legend$col, as.character(df[, 2]))
         }
       }
     }
     else if(class.object[1] %in%  c("mlspls"))
     {
       if (mapping != "XY") {
         if (is.null(legend$legend) && is.null(legend$col)) {
           if (ncol(mat$design) >= 2) {
             df = data.frame(mat$design[, 2], sample.sideColors[, 1])
             df = unique(df)
             legend$legend = as.character(df[, 1])
             legend$col = as.character(df[, 2])
           }
           if (ncol(mat$design) == 3) {
             df = data.frame(mat$design[, 3], sample.sideColors[, 2])
             df = unique(df)
             legend$legend = c(legend$legend, as.character(df[, 1]))
             legend$col = c(legend$col, as.character(df[, 2]))
           }
         }
       }
     }
     if (is.null(legend$legend))
       stop("argument \"legend$legend\" is missing, with no default")
     
     #-- fill
     if (is.null(legend$fill)) legend$fill = legend$col
     
     par(mar = c(0, 0, 0, 0), new = TRUE)
     plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")
     
     if (!is.null(legend$title))
       legend(x = legend$x, y = legend$y, legend = legend$legend, 
              col = legend$col, fill = legend$fill, bty = legend$bty,title=legend$title,cex=legend$cex, ...)
     else
       legend(x = legend$x, y = legend$y, legend = legend$legend, 
              col = legend$col, fill = legend$fill, bty = legend$bty,cex=legend$cex, ...)
     
    }
    return(invisible(res))
  }



imageMap <- 
  function (mat,
            color,
            row.names,
            col.names,
            row.sideColors,
            col.sideColors,
            row.cex = NULL,
            col.cex = NULL,
            cluster,
            ddr,
            ddc,
            cut.tree = c(0, 0),
            transpose = FALSE,
            symkey = TRUE, 
            keysize = c(1, 1),            
            zoom = FALSE, 
            main = NULL,
            xlab = NULL,
            ylab = NULL,
            margins = c(5, 5),
            lhei = NULL,
            lwid = NULL)
{
    #-- image map --------------------------------------------------------------#
    #----------
    
    if (isTRUE(symkey)) {
      max.mat = max(abs(mat), na.rm = TRUE)
      min.mat = -max.mat
    }
    else {
      max.mat = max(mat, na.rm = TRUE)
      min.mat = min(mat, na.rm = TRUE)
    }
    
    if (isTRUE(transpose)) {
      mat = t(mat)
      
      temp = col.sideColors
      col.sideColors = row.sideColors
      row.sideColors = temp
      
      temp = col.names
      col.names = row.names
      row.names = temp
      
      if (cluster == "both")
      {temp = ddc
       ddc = ddr
       ddr = temp}
      else if (cluster == "column") 
        ddr=ddc
      else if (cluster == "row") 
        ddc=ddr
      
      lhei = NULL
      lwid = NULL
    }
    
    nr = nrow(mat)
    nc = ncol(mat)
    
    #-- row.cex and col.cex
    if (is.null(row.cex)) row.cex = min(1, 0.2 + 1/log10(nr))
    if (is.null(col.cex)) col.cex = min(1, 0.2 + 1/log10(nc))
    
    #-- breaks
    breaks = length(color) + 1
    breaks = seq(min.mat, max.mat, length = breaks)
    
    nbr = length(breaks)
    ncol = nbr - 1
    
    min.breaks = min(breaks)
    max.breaks = max(breaks)
    
    mat[mat < min.breaks] = min.breaks
    mat[mat > max.breaks] = max.breaks
    mat = t(mat)
    
    #-- layout matrix
    lmat = matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
    csc = rsc = FALSE
    
    if (!is.null(col.sideColors)) {
      lmat = rbind(lmat[1, ], c(NA, 3), lmat[2, ] + 1)
      if (is.null(lhei)) {
        n.csc = ncol(col.sideColors)
        lhei = c(keysize[2], 0.15 + 0.1 * (n.csc - 1), 4)
      }
      csc = TRUE
    }
    
    if (!is.null(row.sideColors)) {
      lmat = cbind(lmat[, 1], c(rep(NA, nrow(lmat) - 1), nrow(lmat) + 2), 
                   lmat[, 2] + c(rep(0, nrow(lmat) - 1), 1))
      if (is.null(lwid)) {
        n.rsc = ncol(row.sideColors)
        lwid = c(keysize[2], 0.15 + 0.1 * (n.rsc - 1), 4)
      }
      rsc = TRUE
    }
    
    lmat[is.na(lmat)] = 0
    
    if (is.null(lhei)) 
      lhei = c(keysize[2], 4)
    
    if (is.null(lwid)) 
      lwid = c(keysize[1], 4)
    
    if (isTRUE(zoom)) {
      graphics.off()
      dev.new(pos=-1)#x11(xpos = -1)
    }
    
    op = par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    #-- layout 1 --#
    par(mar = c(5, 2, 2, 1), cex = 0.75)
    
    z = seq(0, 1, length = length(color))
    z = matrix(z, ncol = 1)
    image(z, col = color, xaxt = "n", yaxt = "n")
    box()
    par(usr = c(0, 1, 0, 1))
    lv = c(min.breaks, (3*min.breaks + max.breaks)/4, (min.breaks + max.breaks)/2,
           (3*max.breaks + min.breaks)/4, max.breaks)
    xv = (as.numeric(lv) - min.mat) / (max.mat - min.mat)
    axis(1, at = xv, labels = round(lv, 2))
    title("Color key", font.main = 1)
    
    #-- layout 2 --#   
    par(mar = c(ifelse(cut.tree[2] != 0, 0.5, 0), 0, ifelse(!is.null(main), 5, 0), margins[2]))
    if ((cluster == "both") || (!transpose && cluster == "column") || (transpose && cluster == "row" )) {
      h = attr(ddc, "height")
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none", 
           ylim = c(cut.tree[2] * h, h))
    }
    else {
      plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")  
    }
    
    if (!is.null(main)) 
      title(main, cex.main = 1.5 * op[["cex.main"]])
    
    #-- layout 3 --#
    if (isTRUE(csc)) {
      par(mar = c(0.5, 0, 0, margins[2]))
      sideColors = as.vector(col.sideColors)
      img = matrix(c(1:(n.csc * nc)), ncol = n.csc, byrow = FALSE)
      
      image(1:nc, 1:n.csc, img, col = sideColors, axes = FALSE, xlab = "", ylab = "")
      abline(h = 1:(n.csc - 1) + 0.5, lwd = 2,
             col = ifelse(par("bg") == "transparent", "white", par("bg")))  
    }
    
    #-- layout 4 --#
    par(mar = c(margins[1], 0, 0, ifelse(cut.tree[1] != 0, 0.5, 0)))
    if ((cluster == "both") || (cluster == "row" & !transpose) || (cluster == "column" & transpose)) {
      h = attr(ddr, "height")
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none", 
           xlim = c(h, cut.tree[1] * h))    
    }
    else {
      plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")  
    }
    
    #-- layout 5 --#
    if (isTRUE(rsc)) {
      par(mar = c(margins[1], 0, 0, 0.5))
      n.rsc = ncol(row.sideColors)
      r.sideColors = row.sideColors[, n.rsc:1]
      sideColors = as.vector(r.sideColors)
      img = matrix(1:(n.rsc * nr), nrow = n.rsc, byrow = TRUE)
      
      image(1:n.rsc, 1:nr, img, col = sideColors, axes = FALSE, xlab = "", ylab = "")
      abline(v = 1:(n.rsc - 1) + 0.5, lwd = 2,
             col = ifelse(par("bg") == "transparent", "white", par("bg")))  
    }
    
    #-- layout 6 --#
    par(mar = c(margins[1], 0, 0, margins[2]))
    image(1:nc, 1:nr, mat, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
            c(0, nr), axes = FALSE, xlab = "", ylab = "", col = color, 
          breaks = breaks)
    axis(1, 1:nc, labels = col.names, las = 2, line = -0.5, tick = 0, 
         cex.axis = col.cex)
    
    if (!is.null(xlab)) 
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    
    axis(4, 1:nr, labels = row.names, las = 2, line = -0.5, tick = 0, 
         cex.axis = row.cex)
    
    if (!is.null(ylab)) 
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    
    #-- ZOOM --#
    #----------#
    flag1 = flag2 = FALSE
    
    if (isTRUE(zoom)) {
      nD = dev.cur()
      zone = FALSE
      
      repeat {
        dev.set(nD)
        
        repeat {
          loc = locator(1, type = "n")
          
          if (is.null(loc) & zone == TRUE) break
          if (is.null(loc) & !isTRUE(flag1)) break
          flag1 = TRUE
          
          x1 = round(loc[[1]] - 0.5) + 0.5
          y1 = round(loc[[2]] - 0.5) + 0.5
          
          if (!(x1 < 0 | x1 > nc + 0.5 | y1 < 0 | y1 > nr + 0.5)) break
        }
        
        if (is.null(loc) & zone == TRUE) break
        
        if (!is.null(loc) & zone == TRUE) {
          rect(xleft.old, ybottom.old, xright.old, ytop.old, border = "white")
          points(x1.old, y1.old, type = "p", pch = 3, 
                 cex = 2, col = "white")
        }
        
        if (is.null(loc) & zone == FALSE) {
          break
        }
        else {
          x1.old = x1
          y1.old = y1
          
          points(x1, y1, type = "p", pch = 3, cex = 2)
          
          repeat {
            loc = locator(1, type = "n")
            
            if (is.null(loc) & zone == TRUE) break
            if (is.null(loc) & !isTRUE(flag2)) break
            flag2 = TRUE
            
            x2 = round(loc[[1]] - 0.5) + 0.5
            y2 = round(loc[[2]] - 0.5) + 0.5
            
            if (!(x2 < 0 | x2 > nc + 0.5 | y2 < 0 | y2 > nr + 0.5)) { zone = TRUE; break }         
          }
          
          if (is.null(loc) & zone == TRUE) break
          if (is.null(loc) & !isTRUE(flag2)) break
          
          xleft.old = min(x1, x2) 
          xright.old = max(x1, x2) 
          ybottom.old = min(y1, y2) 
          ytop.old = max(y1, y2) 
          
          rect(xleft.old, ybottom.old, xright.old, ytop.old)
        }
        
        dev.new()#x11()
        plot.par = par(no.readonly = TRUE)
        
        if (isTRUE(zone)) {
          xleft = xleft.old + 0.5
          ybottom = ybottom.old + 0.5
          xright = xright.old - 0.5
          ytop = ytop.old - 0.5
          nr.zoom = length(xleft:xright)
          nc.zoom = length(ybottom:ytop)
          mat.zoom = matrix(mat[xleft:xright, ybottom:ytop], 
                            nrow = nr.zoom, ncol = nc.zoom)
          rlab.zoom = col.names[xleft:xright]
          clab.zoom = row.names[ybottom:ytop]
          r.cex = min(1.2, 0.2 + 1/log10(nr.zoom))
          c.cex = min(1.2, 0.2 + 1/log10(nc.zoom))
          
          layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
          
          # layout 1
          par(mar = c(5, 2, 2, 1), cex = 0.75)    		
          image(z, col = color, xaxt = "n", yaxt = "n")
          box()
          par(usr = c(0, 1, 0, 1))
          lv = c(min.breaks, (3*min.breaks + max.breaks)/4, (min.breaks + max.breaks)/2,
                 (3*max.breaks + min.breaks)/4, max.breaks)
          xv = (as.numeric(lv) - min.mat) / (max.mat - min.mat)
          axis(1, at = xv, labels = round(lv, 2))
          title("Color key", font.main = 1)
          
          #-- layout 2 --#   
          par(mar = c(ifelse(cut.tree[2] != 0, 0.5, 0), 0, ifelse(!is.null(main), 5, 0), margins[2]))
          if ((cluster == "both") || (cluster == "column" & !transpose) || (cluster == "row" & transpose)) {
            plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none", xlim = c(xleft - 0.5, xright + 0.5))
          }
          else {
            plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")  
          }
          
          if (!is.null(main)) 
            title(main, cex.main = 1.5 * op[["cex.main"]])
          
          # layout 3
          if (isTRUE(csc)) {
            par(mar = c(0.5, 0, 0, margins[2]))
            sideColors = as.vector(col.sideColors[xleft:xright, ])
            img = matrix(c(1:(n.csc * nr.zoom)), ncol = n.csc, byrow = FALSE)
            
            image(1:nr.zoom, 1:n.csc, img, col = sideColors, axes = FALSE, xlab = "", ylab = "")
            abline(h = 1:(n.csc - 1) + 0.5, lwd = 2,
                   col = ifelse(par("bg") == "transparent", "white", par("bg")))
          }
          
          #-- layout 4 --#
          par(mar = c(margins[1], 0, 0, ifelse(cut.tree[1] != 0, 0.5, 0)))
          if ((cluster == "both") || (cluster == "row" & !transpose) || (cluster == "column" & transpose)) {
            plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none", ylim = c(ybottom - 0.5, ytop + 0.5))    
          }
          else {
            plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")  
          }
          
          # layout 5
          if (isTRUE(rsc)) {
            par(mar = c(margins[1], 0, 0, 0.5))
            r.sideColors = row.sideColors[ybottom:ytop, n.rsc:1]
            sideColors = as.vector(r.sideColors)
            img = matrix(1:(n.rsc * nc.zoom), nrow = n.rsc, byrow = TRUE)
            
            image(1:n.rsc, 1:nc.zoom, img, col = sideColors, axes = FALSE, xlab = "", ylab = "")
            abline(v = 1:(n.rsc - 1) + 0.5, lwd = 2,
                   col = ifelse(par("bg") == "transparent", "white", par("bg")))
          }
          
          # layout 6
          par(mar = c(margins[1], 0, 0, margins[2]))
          image(1:nr.zoom, 1:nc.zoom, mat.zoom, col = color, 
                breaks = breaks, axes = FALSE, xlab = "", ylab = "")
          
          axis(1, 1:nr.zoom, labels = rlab.zoom, las = 2, 
               line = -0.5, tick = 0, cex.axis = r.cex)
          
          if (!is.null(xlab)) 
            mtext(xlab, side = 1, line = margins[1] - 1.25)
          
          axis(4, 1:nc.zoom, labels = clab.zoom, las = 2, 
               line = -0.5, tick = 0, cex.axis = c.cex)
          
          if (!is.null(ylab)) 
            mtext(ylab, side = 4, line = margins[2] - 1.25)
        }
        
        par(plot.par)
      }
    }
    
    par(op)
  }
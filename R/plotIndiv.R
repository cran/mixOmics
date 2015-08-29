# Copyright (C) 2015
# Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD


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


#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv for PLS, sPLS, PLS-DA, SPLS-DA, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#

plotIndiv <-
  function(object,
           comp = c(1, 2),
           ind.names = TRUE,
           rep.space = "X-variate",
           blocks = NULL, # to choose which block data to plot, when using GCCA module
           X.label = NULL,
           Y.label = NULL,
           abline.line = FALSE,
           col.per.group,
           col,
           cex,
           pch,
           plot.ellipse = FALSE,
           ellipse.level = 0.95,
           group,  # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
           main="plotIndiv",
           add.legend=FALSE,
           style="ggplot2", # can choose between graphics, lattice or ggplot2
           ...)
{
    
    class.object = class(object)
    
    object.pls=c("pls","spls","splsda","plsda","mlspls","mlsplsda","rcc")
    object.pca=c("ipca","sipca","pca","spca")
    object.blocks=c("sgcca","rgcca", "sgccda")
    
    ### Start: Validation of arguments
    ncomp = object$ncomp
    if (class.object[1] %in% object.blocks) {
        
        if (is.null(blocks)){
            blocks = object$names$blocks
        } else if (is.numeric(blocks) & min(blocks) > 0 &  max(blocks) <= length(object$names$blocks)) {
            blocks = object$names$blocks[blocks]
        } else if (is.character(blocks)) {
            if (!any(blocks %in% object$names$blocks))
            stop("One element of 'blocks' does not match with the names of the blocks")
        } else {
            stop("Incorrect value for 'blocks", call. = FALSE)
        }
        object$variates = object$variates[names(object$variates) %in% blocks]
        
        if (any(object$ncomp[blocks] == 1)) {
            stop(paste("The number of components for one selected block '", paste(blocks, collapse = " - "),"' is 1. The number of components must be superior or equal to 2."), call. = FALSE)
        }
        
        if (is.null(object$indY) & rep.space %in% c("Y-variate", "XY-variate")) {
            stop("For an object of class 'blocks', 'rep.space' must be 'X-variate", call. = FALSE)
        }
        ncomp = object$ncomp[blocks]
    }
    
    #-- ellipse.level
    if ((ellipse.level > 1) | (ellipse.level < 0))
    stop("The value taken by 'ellipse.level' must be between 0 and 1")
    
    #-- comp
    if (length(comp) != 2)
    stop("'comp' must be a numeric vector of length 2.")
    
    if (!is.numeric(comp))
    stop("Invalid vector for 'comp'.")
    
    if (any(ncomp < max(comp)))
    stop("Each element of 'comp' must be smaller or equal than ", max(object$ncomp), ".", call. = FALSE)
    
    comp1 = round(comp[1]); comp2 = round(comp[2])
    
    #-- Specific to pls object (choice betwee X, Y and XY)
    if (class.object[1] %in% object.pls){
        if (rep.space == "X-variate")
        {object$variates = object$variates["X"]; blocks = "X"}
        if (rep.space == "Y-variate")
        {object$variates = object$variates["Y"]; blocks = "Y"}
        if (rep.space == "XY-variate"){
            object$variates$XYvariates = do.call("cbind", sapply(c(comp1, comp2), function(k) {lapply(object$variates[names(object$variates) != "Y"],
                function(x){(x[, k, drop = FALSE] + object$variates$Y[, k, drop = FALSE]) / 2})}))
            object$variates = object$variates["XYvariates"]; blocks = "XY combined"
        }
    }
    
    if (class.object[1] %in% object.pca)
      blocks = "X"
    
    #-- rep.space
    rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate"))
    
    #-- Start: Retrieve variates from object
    x = y = list()
    if (class.object[1] %in%  c(object.pls, object.blocks)) {
        
        if (is.logical(ind.names) & isTRUE(ind.names)) {
            ind.names = object$names$indiv
        }
        if (length(ind.names) > 1) {
            if (length(ind.names) != length(object$names$indiv))
              stop("'ind.names' must be a character vector of length ", nrow(object$X), " or a boolean atomic vector.")
        }
        
        x = lapply(object$variates, function(x){x[, comp1, drop = FALSE]})
        y = lapply(object$variates, function(x){x[, comp2, drop = FALSE]})
        if (is.null(X.label)) {
          if (rep.space == "X-variate") {X.label = paste("X-variate", comp1)}
          if (rep.space == "Y-variate") {X.label = paste("Y-variate", comp1)}
          if (rep.space == "XY-variate") {X.label = paste("XY-variate", comp1)}
        }
        if (is.null(Y.label)) {
          if (rep.space == "X-variate") {Y.label = paste("X-variate", comp2)}
          if (rep.space == "Y-variate") {Y.label = paste("Y-variate", comp2)}
          if (rep.space == "XY-variate") {Y.label = paste("XY-variate", comp2)}
        }
        
    } else if (class.object[1] %in%  object.pca) {
        
        if (is.logical(ind.names)) {
            if (isTRUE(ind.names))
            ind.names = rownames(object$x)
        }
        if (length(ind.names) > 1) {
            if (length(ind.names) != nrow(object$x))
            stop("'ind.names' must be a character vector of length ", nrow(object$x), " or a boolean atomic vector.")
        }
        if (rep.space != "X-variate")
        stop("'rep.space' must be equal to 'X-variate' for an object belonging to 'pca' or 'ipca'")
        
        x[[1]] = object$x[, comp[1]]
        y[[1]] = object$x[, comp[2]]
        
        if (is.null(X.label)) X.label = paste("Dimension ", comp1)
        if (is.null(Y.label)) Y.label = paste("Dimension ", comp2)
    }
    #-- End: Retrieve variates from object
    
    #-- ind.names
    display.names = FALSE
    if (length(ind.names) == length(x[[1]])){
        display.names = TRUE
    }
    
    #-- Define group
    missing.group = FALSE
    if (missing(group) & any(class.object %in% c("plsda","splsda"))){
      group = factor(map(object$ind.mat), labels = object$names$Y)
    } else if (missing(group) & any(class.object %in% c("sgccda"))){
      group = factor(map(object$ind.mat), labels = object$names$colnames$Y)
    } else if (!missing(group)) {
      missing.group = TRUE
      if (!is.factor(group)){
        group = as.factor(group)
      }
      object$ind.mat = unmap(group)
      
      if (length(group) != length(x[[1]]))
        stop("Length of 'group' should be of length ", length(x[[1]]), ", the sample size of your data")
    } else {
      group = factor(rep("No group", length(x[[1]])))
      object$ind.mat = unmap(group)
    }
    
    #-- col.per.group argument
    if (missing(col.per.group)){
      if (nlevels(group) < 10) {
        #only 10 colors in color.mixo
        col.per.group = color.mixo(1:nlevels(group))
      } else {
        #use color.jet
        col.per.group = color.jet(nlevels(group))
      }
    } else {
      if (length(col.per.group) == 1) {
        col.per.group = rep(col.per.group, nlevels(group))
      } else if (length(col.per.group) != length(x[[1]]) & length(col.per.group) != nlevels(group)) {
        stop("Length of 'col.per.group' should be either of length 1 or of length ", nlevels(group), " (the number of groups) or of length ", length(x[[1]]), " (the sample size or your data).
          Alternatively, use the argument 'col' to give one color per sample")
      }
      missing.group = TRUE
    }
    
    levels.color = vector(, length(x[[1]]))
    if (length(col.per.group) != length(x[[1]])) {
      for (i in 1 : nlevels(group)){
        levels.color[group == levels(group)[i]] = col.per.group[i]
      }
    } else {
      levels.color = col.per.group
    }
    
    #-- col argument   
    missing.col = FALSE
    if (!missing(col)){
      if (length(col) > length(x[[1]]))
        stop("Length of 'col' should be of length inferior or equal to ", length(x[[1]]),".")
      
      col = factor(rep(col, ceiling(length(x[[1]])/length(col)))[1 : length(x[[1]])])
      if (!missing.group) {
        group = col
        levels.color = col
        col.per.group = levels(col)
        object$ind.mat = unmap(group)
      }
      missing.col = TRUE
    } else {
      col = levels.color
    }
      
    #-- cex argument
    if (missing(cex)){
      if (style == "ggplot2"){
        cex = rep(5, length(x[[1]]))
      } else {
        cex = rep(1, length(x[[1]]))
      }
    } else {
      if (length(cex) == 1){
        cex = rep(cex, length(x[[1]]))
      } else if (length(cex) > length(x[[1]])) {
        stop("Length of 'cex' should be of length inferior or equal to ", length(x[[1]]),".")
      } else {
        cex = rep(cex, ceiling(length(x[[1]])/length(cex)))[1 : length(x[[1]])]
      }
    }
      
    #-- pch argument
    if (missing(pch)){
      if (missing.col){
        pch = as.numeric(col)
      } else {
        pch = as.numeric(group)
      }
    } else {
      if (length(pch) == 1){
        pch = rep(pch, length(x[[1]]))
      } else if (length(pch) > length(x[[1]])){
        stop("Length of 'pch' should be of length inferior or equal to ", length(group),".")
      } else {
        pch = rep(pch, ceiling(length(x[[1]])/length(pch)))[1 : length(x[[1]])]
      }
    }
        
    if (plot.ellipse) {
      #-- Start: Computation ellipse
      min.ellipse = max.ellipse = xlim.min = xlim.max = ylim.min = ylim.max = list()
      ind.gp = matrice = cdg = variance = list()
      ind.gp = lapply(1 : ncol(object$ind.mat), function(x){which(object$ind.mat[, x]==1)})
      matrice = lapply(1 : length(x), function(z1) {lapply(ind.gp, function(z2){matrix(c(x[[z1]][z2], y[[z1]][z2]), ncol = 2)})})
      cdg = lapply(1 : length(x), function(z){ lapply(matrice[[z]], colMeans)})
      variance = lapply(1 : length(x), function(z){lapply(matrice[[z]], var)})
      coord.ellipse = lapply(1 : length(x), function(z1){ lapply(1 : ncol(object$ind.mat), function(z2){ellipse(variance[[z1]][[z2]],
          centre = cdg[[z1]][[z2]],
          level = ellipse.level)})})
      max.ellipse = lapply(1 : length(x), function(z1) {sapply(coord.ellipse[[z1]], function(z2){apply(z2, 2, max)})})
      min.ellipse = lapply(1 : length(x), function(z1) {sapply(coord.ellipse[[z1]], function(z2){apply(z2, 2, min)})})
      #-- End: Computation ellipse
        
      xlim = lapply(1 : length(x), function(z) {c(min(x[[z]], min.ellipse[[z]][1, ]), max(x[[z]], max.ellipse[[z]][1, ]))})
      ylim = lapply(1 : length(x), function(z) {c(min(y[[z]], min.ellipse[[z]][2, ]), max(y[[z]], max.ellipse[[z]][2, ]))})
    } else {
      xlim = lapply(1 : length(x), function(z) {c(min(x[[z]]), max(x[[z]]))})
      ylim = lapply(1 : length(x), function(z) {c(min(y[[z]]), max(y[[z]]))})
    }
    
    #-- Start: data set
    df = list()
    for (i in 1 : length(x)) {
      df[[i]] = data.frame(x = x[[i]], y = y[[i]], group = group)
    }
    
    df = data.frame(do.call(rbind, df), "Block" = paste0("Block: ", unlist(lapply(1 : length(df), function(z){rep(blocks[z], nrow(df[[z]]))}))))
    names(df)[1:2] = c("x", "y")
    
    if (display.names)
      df$names = rep(ind.names, length(x))
    
    if (plot.ellipse == TRUE){
      df.ellipse = data.frame(do.call("rbind", lapply(1 : length(x), function(k){do.call("cbind", coord.ellipse[[k]])})), "Block" = paste0("Block: ", rep(blocks, each = 100)))
      names(df.ellipse)[1 : (2*nlevels(group))] = paste0("Col", 1 : (2*nlevels(group)))
    }
    
    df$pch = pch; df$cex = cex; df$col.per.group = levels.color; df$col = as.character(col)
    #-- End: data set
    
    #-- Start: ggplot2
    if (style == "ggplot2"){
        #-- Initialise ggplot2
        p = ggplot(df, aes(x = x, y = y, color = group),
                  main = main, xlab = X.label, ylab = Y.label) + theme_bw()
        
        #-- Display sample or row.names
        for (i in levels(group)){
          if (display.names) {
            p = p + geom_text(data = subset(df, group == i), aes(label = names), size = 0)
          } else {
            p = p + geom_point(data = subset(df, group == i), size = 0, shape = 0)
          }
        }
        
        #-- Modify scale colour - Change X/Ylabel - split plots into Blocks  
        p = p + scale_colour_manual(values = col.per.group[match(levels(factor(as.character(group))), levels(group))], name = "Legend", breaks = levels(group))
        p = p + labs(list(title = main, x = X.label, y = Y.label)) + facet_wrap(~ Block, ncol = 2, scales = "free", as.table = FALSE)
        
        #-- color samples according to col
        for (i in unique(col)){
          if (display.names) {
            p = p + geom_text(data = subset(df, col == i), 
                              aes(label = names), 
                              color = df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$col,
                              size = df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$cex)
          } else {
            p = p + geom_point(data = subset(df, col == i), 
                               color = df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$col,
                               size = df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$cex, 
                               shape = df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$pch)
          }
        }      
        
        #-- Legend
        if (!add.legend) {
          p = p + theme(legend.position="none")
        } else {
          p = p + guides(colour = guide_legend(override.aes = list(shape = "-", size = 10)))
        }
        
        #-- abline
        if (abline.line)
          p = p + geom_vline(aes(xintercept = 0), linetype = 2, colour = "darkgrey") + geom_hline(aes(yintercept = 0),linetype = 2,colour = "darkgrey")
        
        #-- ellipse
        if (plot.ellipse == TRUE) {
            for (i in 1 : nlevels(group)){
                p = p + geom_path(data = df.ellipse,
                aes_string(x = paste0("Col", 2*(i - 1) + 1), y = paste0("Col", 2 * i),
                           label = "Block", group = NULL), color = col.per.group[i])
            }
        }
        return(p)
    }
    #-- End: ggplot2
    
    #-- Start: Lattice
    if(style=="lattice") {
        p = xyplot(y ~ x | Block, data = df, xlab = X.label, ylab = Y.label, main = main,
        group = if (display.names) {names} else {group},
        scales= list(x = list(relation = "free", limits = xlim),
        y = list(relation = "free", limits = ylim)),
        
        #-- Legend
        key = if(add.legend == TRUE) {list(space = "right", title = "Legend", cex.title = 1.5,
            text = list(levels(group)), point = list(col = col.per.group), cex = 2, pch = "-")}
        else {NULL},
        
        panel = function(x, y, subscripts, groups, display = display.names,...) {
            #-- Abline
            if (abline.line) { panel.abline(v = 0, lty = 2, col = "darkgrey")
                panel.abline(h = 0, lty = 2, col = "darkgrey")}
            
            #-- Display sample or row.names
            for (i in 1 : nlevels(group)){
                if (display){
                    ltext(x = x[group == levels(group)[i]], y = y[group == levels(group)[i]],
                          labels = groups[subscripts & group == levels(group)[i]], col = "white", cex = 0) 
                } else {
                    lpoints(x = x[group == levels(group)[i]], y = y[group == levels(group)[i]], col = "white", cex = 0, pch = 0)
                }
            }
            
            #-- color samples according to col
            for (i in unique(col)){
              if (display) {
                ltext(x = x[col == i], y = y[col == i], labels =  groups[subscripts & col == i],
                      col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex)
              } else {
                lpoints(x = x[col == i],  y = y[col == i],
                        col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex, pch = df[df$col == i, ]$pch)
              }
            }
        })
        print(p) #-- the lattice plot needs to be printed in order to display the ellipse(s)
        
        #-- ellipse
        if (plot.ellipse) {
            panels = trellis.currentLayout(which = "panel")
            for (k in 1 : length(x)){
                ind = which(panels == k, arr.ind = TRUE)
                trellis.focus("panel",ind[2], ind[1])
                
                for (i in 1 : nlevels(group)) {
                    panel.lines(x = df.ellipse[df.ellipse$Block %in% paste0("Block: ", blocks[k]), paste0("Col", 2*(i - 1) + 1)],
                    y = df.ellipse[df.ellipse$Block %in% paste0("Block: ", blocks[k]), paste0("Col", 2 * i)],
                    col = col.per.group[i])
                }
            }
            trellis.unfocus()
        }
    }
    #-- End: Lattice
    
    #-- Start: graphics
    if(style=="graphics") {
        
        opar <- par()[! names(par()) %in% c("cin", "cra", "csi", "cxy", "din", "page")]
        #-- Define layout
        if (add.legend) {
            layout(cbind(matrix(1 : (ceiling(length(x)/2) * 2), ceiling(length(x)/2), min(length(x), 2), byrow = TRUE) + 1, 1),
                  if (ceiling(length(x)/2) == 1) {widths=c(0.6, 0.4)} else {widths=c(0.35, 0.35, 0.3)})
            plot(1,1, type = "n", axes = FALSE, ann = FALSE)
            legend(0.6, 1, col = col.per.group, legend = levels(group), pch = "-", title = 'Legend', cex = 1.5)
        } else {
            layout(matrix(1 : (ceiling(length(x)/2) * 2), ceiling(length(x)/2), min(length(x), 2), byrow = TRUE))
        }
        
        for (k in 1 : length(x)){
            #-- initialise plot
            plot(df[df$Block %in% paste0("Block: ", blocks[k]), "x" ],
                 df[df$Block %in% paste0("Block: ", blocks[k]), "y" ],
                 type = "n", xlab = X.label, ylab = Y.label, main = paste0("Block: ", blocks[k]),
                 xlim = c(xlim[[k]][1], xlim[[k]][2]), ylim = c(ylim[[k]][1], ylim[[k]][2]))
            
            #-- Display sample or row.names
            for (i in 1 : nlevels(group)){
              if (display.names) {
                text(x = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "x"],
                       y = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "y"],
                       labels = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "names"],
                       col = "white", cex = 0)
              } else {
                points(x = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "x"],
                       y = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "y"],
                       col = "white", cex = 0, pch = 0)
              }
            }  
            
            #-- color samples according to col
            for (i in unique(col)){
              if (display.names) {
                text(x = df[df$col == i & df$Block %in% paste0("Block: ", blocks[k]), "x"],
                     y = df[df$col == i & df$Block %in% paste0("Block: ", blocks[k]), "y"],
                     labels = df[df$col == i & df$Block %in% paste0("Block: ", blocks[k]), "names"],
                     col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex)
              } else {
                points(x = df[df$col == i & df$Block %in% paste0("Block: ", blocks[k]), "x"],
                       y = df[df$col == i & df$Block %in% paste0("Block: ", blocks[k]), "y"],
                       col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex, pch = df[df$col == i, ]$pch)
              }
            }
            
            #-- Abline
            if (abline.line)
            abline(v = 0, h = 0, lty = 2)
              
            #-- Ellipse
            if (plot.ellipse == TRUE) {
                for (i in 1 : nlevels(group)){
                  lines(x = df.ellipse[df.ellipse$Block %in% paste0("Block: ", blocks[k]), paste0("Col", 2*(i - 1) + 1)],
                        y = df.ellipse[df.ellipse$Block %in% paste0("Block: ", blocks[k]), paste0("Col", 2 * i)],
                        col = col.per.group[i])
                }
            }
        }
        
        title(main, outer = TRUE, line = -1)
        if (length(x) != (round(length(x)/2) * 2) & length(x) != 1)
        plot(1,1, type = "n", axes = FALSE, ann = FALSE)
        par(opar)
    }
    #-- End: graphics
}
# Copyright (C) 2015
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.


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
cex,
pch,
plot.ellipse = FALSE,
ellipse.level = 0.95,
group = NULL,  # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
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
        if (is.null(X.label)) X.label = paste("X-variate", comp1)
        if (is.null(Y.label)) Y.label = paste("X-variate", comp2)
        
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
    if (is.null(group) & any(class.object %in% c("plsda","splsda"))){
        group = factor(map(object$ind.mat), labels = object$names$Y)
    }
    if (is.null(group) & any(class.object %in% c("sgccda"))){
        group = factor(map(object$ind.mat), labels = object$names$colnames$Y)
    }
    if (!is.null(group)) {
        if (!is.factor(group)){
            group = as.factor(group)
        }
        object$ind.mat = unmap(group)
    } else {
        group = factor(rep("No group", length(x[[1]])))
        object$ind.mat = unmap(group)
    }
    
    #-- col argument
    if (missing(col.per.group)){
        if (nlevels(group) < 10) {
            #only 10 colors in color.mixo
            levels.color = color.mixo(1:nlevels(group))
        } else{
            #use color.jet
            levels.color = color.jet(nlevels(group))
        }
    } else if (length(col.per.group) == 1) {
        levels.color = rep(col.per.group, nlevels(group))
    } else if (length(col.per.group) == nlevels(group)){
        levels.color = col.per.group
    } else if (length(col.per.group) == length(x[[1]])){
        stop("Length of 'col.per.group' should be of length = ", nlevels(group), " the number of groups.
        Alternatively, use the argument 'group' to give one color per sample")
    } else {
        stop("Length of 'col.per.group' should be of length = ", nlevels(group), " the number of groups.
        Alternatively, use the argument 'group' to give one color per sample")
    }
    
    #-- cex argument
    if (missing(cex)){
        if (style == "ggplot2"){
            cex = rep(3, nlevels(group))
        } else if(style == "graphics"){
            cex = rep(1, nlevels(group))
        } else if (style == "lattice") {
            cex = rep(1, nlevels(group))
        }
    } else if (length(cex) == 1){
        cex = rep(cex, nlevels(group))
    } else if (length(cex) != nlevels(group)){
        stop("'cex' must be a character vector of length ", nlevels(group) ," or one size")
    }
    
    #-- pch argument
    if (missing(pch)) {
        pch = 1 : nlevels(group)
    } else if (length(pch) == 1) {
        pch = rep(pch, nlevels(group))
    } else if (length(pch) != nlevels(group)) {
        stop("'pch' must be a character vector of length ", nlevels(group) ," or one size")
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
        names(df.ellipse)[1 : (2*nlevels(group))] = gsub(" ", "", paste(c("ellipse.x", "ellipse.y"), rep(levels(group), time = 1, each = 2),  sep = "."))
    }
    #-- End: data set
    
    #-- Start: ggplot2
    if (style == "ggplot2"){
        #-- Initialise ggplot2
        p = ggplot(df, aes(x = x, y = y, color = group),
        main = main,
        xlab = X.label,
        ylab = Y.label) + theme_bw()
        
        #-- Display sample or row.names
        for (i in 1 : nlevels(group)){
            if (display.names) {
                p = p + geom_text(data = subset(df, group == levels(group)[i]),
                aes(label = names), size = cex[i])
            } else {
                p = p + geom_point(data = subset(df, group == levels(group)[i]),
                size = cex[i], shape = pch[i])
            }
        }
        
        #-- Modify scale colour - Change X/Ylabel - split plots into Blocks
        p = p + scale_colour_manual(values = levels.color[match(levels(factor(as.character(group))), levels(group))],
        name = "Legend", breaks = levels(group))
        
        p = p + labs(list(title = main, x = X.label, y = Y.label)) + facet_wrap(~ Block, ncol = 2, scales = "free")
        
        #-- Legend
        if (!add.legend) {
            p = p + theme(legend.position="none")
        } else if (!(display.names)) {
            p = p + guides(colour = guide_legend(override.aes = list(shape = pch, size = cex)))
        }
        
        #-- abline
        if (abline.line)
        p = p + geom_vline(aes(xintercept = 0), linetype = 2, colour = "darkgrey") + geom_hline(aes(yintercept = 0),linetype = 2,colour = "darkgrey")
        
        #-- ellipse
        if (plot.ellipse == TRUE) {
            for (i in 1 : nlevels(group)){
                p = p + geom_path(data = df.ellipse,
                aes_string(x = gsub(" ", "", paste("ellipse.x", levels(group)[i], sep = ".")),
                y = gsub(" ", "", paste("ellipse.y", levels(group)[i], sep = ".")),
                label = "Block", group = NULL), color = levels.color[i])
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
        key = if(add.legend == TRUE) {list(space = "right", title = "Legend", cex.title = 1.25, cex = cex,
            text = list(levels(group)),
            point = list(col = levels.color),
            pch = if (display.names){15} else {pch})}
        else {NULL},
        
        panel = function(x, y, subscripts, groups, display = display.names, ...) {
            #-- Abline
            if (abline.line) { panel.abline(v = 0, lty = 2, col = "darkgrey")
                panel.abline(h = 0, lty = 2, col = "darkgrey")}
            
            #-- Display sample or row.names
            for (i in 1 : nlevels(group)){
                if (display){
                    ltext(x = x[group == levels(group)[i]],
                    y = y[group == levels(group)[i]],
                    cex = cex[i], col = levels.color[i],
                    labels = groups[subscripts & group == levels(group)[i]])
                } else {
                    lpoints(x = x[group == levels(group)[i]],
                    y = y[group == levels(group)[i]],
                    cex = cex[i], col = levels.color[i], pch = pch[i])
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
                    panel.lines(x = df.ellipse[df.ellipse$Block %in% paste0("Block: ", blocks[k]), gsub(" ", "", paste("ellipse.x", levels(group)[i], sep = "."))],
                    y = df.ellipse[df.ellipse$Block %in% paste0("Block: ", blocks[k]), gsub(" ", "", paste("ellipse.y", levels(group)[i], sep = "."))],
                    col = levels.color[i])
                }
            }
            trellis.unfocus()
        }
    }
    #-- End: Lattice
    
    #-- Start: graphics
    if(style=="graphics") {
        
        #-- Define layout
        if (add.legend) {
            layout(cbind(matrix(1 : (ceiling(length(x)/2) * 2), ceiling(length(x)/2), min(length(x), 2), byrow = TRUE) + 1, 1),
            if (ceiling(length(x)/2) == 1) {widths=c(0.6, 0.4)} else {widths=c(0.35, 0.35, 0.3)})
            plot(1,1, type = "n", axes = FALSE, ann = FALSE)
            if (length(ind.names) == length(x[[1]])){
                legend(0.6, 1, col = levels.color, legend = levels(group), pch = 15, title = 'Legend', cex = 1)
            } else {
                legend(0.6, 1, col = levels.color, legend = levels(group), pch = pch, title = 'Legend', cex = 1)
            }
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
                if (length(ind.names) == length(x[[1]])) {
                    text(x = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "x"],
                    y = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "y"],
                    labels = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "names"],
                    cex = cex[i], col = levels.color[i])
                } else {
                    points(x = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "x"],
                    y = df[group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "y"],
                    cex = cex[i], col = levels.color[i], pch = pch[i])
                }
            }
            
            #-- Abline
            if (abline.line)
            abline(v = 0, h = 0, lty = 2)
            
            #-- Ellipse
            if (plot.ellipse == TRUE) {
                for (i in 1 : nlevels(group)){
                    lines(x = df.ellipse[df.ellipse$Block %in% paste0("Block: ", blocks[k]), gsub(" ", "", paste("ellipse.x", levels(group)[i], sep = "."))],
                    y = df.ellipse[df.ellipse$Block %in% paste0("Block: ", blocks[k]), gsub(" ", "", paste("ellipse.y", levels(group)[i], sep = "."))],
                    col = levels.color[i])
                }
            }
        }
        
        title(main, outer = TRUE, line = -1)
        if (length(x) != (round(length(x)/2) * 2) & length(x) != 1)
        plot(1,1, type = "n", axes = FALSE, ann = FALSE)
    }
    #-- End: graphics
}
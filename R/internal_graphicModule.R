#############################################################################################################
# Author:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 16-03-2016
# last modified: 12-04-2016
#
# Copyright (C) 2016
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
#############################################################################################################

# --------------------------------------------------------------------------------------
# Internal helpers functions to run "plotIndiv" functions
# --------------------------------------------------------------------------------------

# df: data frame with all the information needed: coordinates (x,y,z), grouping factor 'group', 'Block' that indicates the block, names (ind.names), 'pch', 'cex' or each point, 'col.per.group' that gives the color of each point, 'pch.legend' that gives the pch of each point for the legend (same as pch?)
# as well as: x0 and y0 if plot centroid==TRUE
# centroid
# star
# ellipse
# df.ellipse
# xlim
# ylim
# title
# X.label
# Y.label
# legend
# display.names

internal_graphicModule=function(df,
centroid,
col.per.group,
title,
X.label,
Y.label,
Z.label,
xlim,
ylim,
zlim,
class.object,
display.names,
legend,
abline,
star,
ellipse,
df.ellipse,
style,
layout=NULL,
missing.col,
axes.box,
study.levels,
plot_parameters,
alpha)
{
    object.pls = c("pls", "spls", "mlspls", "rcc")
    object.pca = c("ipca", "sipca", "pca", "spca", "prcomp")
    object.blocks = c("sgcca", "rgcca")
    object.mint = c("mint.pls", "mint.spls", "mint.plsda", "mint.splsda")
    #class.object=class(object)
    
    # to satisfy R CMD check that doesn't recognise x, y and group as variable (in aes)
    x = y = group = pch = studyname = NULL
    
    size.title = plot_parameters$size.title
    size.subtitle = plot_parameters$size.subtitle
    size.xlabel = plot_parameters$size.xlabel
    size.ylabel = plot_parameters$size.ylabel
    size.axis = plot_parameters$size.axis
    size.legend = plot_parameters$size.legend
    size.legend.title = plot_parameters$size.legend.title
    legend.position = plot_parameters$legend.position
    point.lwd = plot_parameters$point.lwd
    
    #-- Start: ggplot2
    if (style == "ggplot2")
    {
        nResp = nlevels(df$Block)
        if (is.null(layout))
        {
            nRows = min(c(3, ceiling(nResp/2)))
            nCols = min(c(3, ceiling(nResp/nRows)))
            layout = c(nRows, nCols)
        } else {
            if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
            stop("'layout' must be a numeric vector of length 2.")
            nRows = layout[1]
            nCols = layout[2]
        }

        #note: at this present time, ggplot2 does not allow xlim to be changed per subplot, so cannot use xlim properly
        
        #-- Initialise ggplot2
        p = ggplot(df, aes(x = x, y = y, color = group),
        main = title, xlab = X.label, ylab = Y.label) +
        theme_bw() + theme(strip.text = element_text(size = size.subtitle, face = "bold"))
        

        #}
        
        #-- Display sample or row.names
        for (i in levels(df$group))
        {
            if (display.names)
            {
                p = p +geom_point(data = subset(df, df$group == i),size = 0, shape = 0)+ geom_text(data = subset(df, df$group == i), aes(label = names), size = 0,show.legend  = F)
            } else {
                p = p + geom_point(data = subset(df, df$group == i), size = 0, shape = 0)
            }
            if (centroid == TRUE)
            {
                p = p + geom_point(data = subset(df[, c("col", "x0", "y0", "Block", "cex", "pch", "group")], df$group == i), aes(x=x0,y=y0), size = 0, shape = 0)
            }
        }
        
        
        #-- Modify scale colour - Change X/Ylabel - split plots into Blocks
        p = p + scale_colour_manual(values = unique(col.per.group)[match(levels(factor(as.character(df$group))), levels(df$group))], name = "Legend", breaks = levels(df$group))
        
        p = p + labs(list(title = title, x = X.label, y = Y.label)) + facet_wrap(~ Block, ncol = nCols, scales = "free", as.table = TRUE) #as.table to plot in the same order as the factor
        p = p + theme(plot.title=element_text(size=size.title),axis.title.x=element_text(size=size.xlabel),axis.title.y=element_text(size=size.ylabel),axis.text=element_text(size=size.axis))# bigger title
        
        #-- xlim, ylim
        p = p + coord_cartesian(xlim=xlim,ylim=ylim)
        
        #-- color samples according to col
        
        for (i in unique(df$col))
        {
            for(j in 1:nlevels(df$Block))
            {
                if (display.names)
                {
                    p = p +geom_point(data = subset(df, col == i & df$Block == levels(df$Block)[j]),size = 0, shape = 0,
                    color = unique(df[df$col == i & df$Block == levels(df$Block)[j], ]$col))+
                    geom_text(data = subset(df, col == i & df$Block == levels(df$Block)[j]),
                    aes(label = names),
                    color = df[df$col == i & df$Block == levels(df$Block)[j], ]$col,
                    size = df[df$col == i & df$Block == levels(df$Block)[j], ]$cex,show.legend  = F)
                } else {
                    p = p + geom_point(data = subset(df, col == i & df$Block == levels(df$Block)[j]),
                    color = unique(df[df$col == i & df$Block == levels(df$Block)[j], ]$col),
                    size = df[df$col == i & df$Block == levels(df$Block)[j], ]$cex,
                    shape = df[df$col == i & df$Block == levels(df$Block)[j],  ]$pch, stroke = point.lwd)# unique(df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$pch))
                }
                if (centroid == TRUE)
                {
                    p = p + geom_point(data = subset(df[, c("col", "x0", "y0", "Block", "cex", "pch", "group")], col == i), aes(x = x0, y = y0),
                    color = unique(df[df$col == i & df$Block == levels(df$Block)[1], ]$col),
                    size = unique(df[df$col == i & df$Block == levels(df$Block)[1], ]$cex),
                    shape = 8, stroke = point.lwd)
                }
            
            }
        }
        
        
        
        #-- Legend
        if (!legend)
        {
            p = p + theme(legend.position="none")
        } else {
            p = p + guides(colour = guide_legend(override.aes = list(shape = if(display.names | any(class.object%in%object.mint) ) {19} else unique(df$pch.legend), size = unique(df$cex),stroke=point.lwd)))    +
            theme(legend.title=element_text(size=size.legend.title),legend.text=element_text(size=size.legend)) +
            theme(legend.position=legend.position)
        }
        
        #-- abline
        if (abline)
        p = p + geom_vline(aes(xintercept = 0), linetype = 2, colour = "darkgrey") + geom_hline(aes(yintercept = 0), linetype = 2, colour = "darkgrey")
        
        #-- star
        if (star == TRUE)
        {
            for (i in 1 : nlevels(df$group))
            {
                p = p + geom_segment(data = subset(df, group == levels(df$group)[i]),
                aes(x = x0, y = y0, xend = x, yend = y,
                label = "Block"), color = unique(col.per.group)[i],size = point.lwd)
            }
        }
        
        #-- ellipse
        if (ellipse == TRUE)
        {
            for (i in 1 : nlevels(df$group))
            {
                p = p + geom_path(data = df.ellipse,
                aes_string(x = paste0("Col", 2*(i - 1) + 1), y = paste0("Col", 2 * i),
                label = "Block", group = NULL), color = unique(col.per.group)[i], size = point.lwd)
            }
        }
        
        plot(p)
    }
    #-- End: ggplot2
    
    #-- Start: ggplot2
    if (style=="ggplot2-MINT")
    {
        
        nResp = nlevels(df$Block)
        if (is.null(layout))
        {
            nRows = min(c(3, ceiling(nResp/2)))
            nCols = min(c(3, ceiling(nResp/nRows)))
            layout = c(nRows, nCols)
        } else {
            if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
            stop("'layout' must be a numeric vector of length 2.")
            nRows = layout[1]
            nCols = layout[2]
        }
        
        #note: at this present time, ggplot2 does not allow xlim to be changed per subplot, so cannot use xlim properly
        df$studyname = factor(df$pch, labels = study.levels)
        
        #-- Initialise ggplot2
        p = ggplot(df, aes(x = x, y = y, color = group, shape = studyname),
        main = title, xlab = X.label, ylab = Y.label) +
        theme_bw() + theme(strip.text = element_text(size = size.subtitle, face = "bold"))
        
        #-- Display sample or row.names
        for (i in levels(df$group))
        {
                p = p + geom_point(data = subset(df, df$group == i), size = subset(df, df$group == i)$cex[1], stroke = point.lwd)
        }
        
        #-- Modify scale colour - Change X/Ylabel - split plots into Blocks
        p = p + scale_colour_manual(values = unique(col.per.group)[match(levels(factor(as.character(df$group))), levels(df$group))], name = "Outcome", breaks = levels(df$group)) +
            labs(shape = "Study")#levels(object$study)[study.ind])

        p = p + scale_shape_manual(values = as.numeric(levels(factor(df$pch)))) # replace the shape/pch by the input, it's converted by default to 1,2,3.. by ggplots

        p = p + labs(list(title = title, x = X.label, y = Y.label)) + facet_wrap(~ Block, ncol = nCols, scales = "free", as.table = TRUE) #as.table to plot in the same order as the factor
        p = p + theme(plot.title = element_text(size=size.title), axis.title.x = element_text(size=size.xlabel), axis.title.y = element_text(size = size.ylabel), axis.text = element_text(size = size.axis))# bigger title
        
        #-- xlim, ylim
        p = p + coord_cartesian(xlim = xlim, ylim = ylim)
        
        #-- Legend
        if (!legend)
        {
            p = p + theme(legend.position="none")
        }else{
            p = p + guides(colour = guide_legend(override.aes = list(size = unique(df$cex)))) +
            theme(legend.title = element_text(size = size.legend.title),legend.text = element_text(size = size.legend))+
            theme(legend.position = legend.position)#,legend.direction="vertical")
        }
        
        #-- abline
        if (abline)
        p = p + geom_vline(aes(xintercept = 0), linetype = 2, colour = "darkgrey") + geom_hline(aes(yintercept = 0),linetype = 2,colour = "darkgrey")


        #-- centroid
        if (centroid == TRUE) #only when one block
        {
            for (i in levels(df$group))
            {
                p = p + geom_point(data = subset(df,  df$group == i), aes(x = x0, y = y0),
                color = subset(df, df$group == i)$col[1],
                size = subset(df, df$group == i)$cex[1],
                shape = 8, stroke = point.lwd)
            }
        }


        #-- star
        if (star == TRUE) #only when one block
        {
            for (i in 1 : nlevels(df$group))
            {
                p = p + geom_segment(data = subset(df, group == levels(df$group)[i]),
                aes(x = x0, y = y0, xend = x, yend = y,
                label = "Block"), color = unique(col.per.group)[i],size = point.lwd)
            }
        }

        #-- ellipse
        if (ellipse == TRUE) #only when one block
        {
            for (i in 1 : nlevels(df$group))
            {
                p = p + geom_path(data = df.ellipse,
                aes_string(x = paste0("Col", 2*(i - 1) + 1), y = paste0("Col", 2 * i),
                label = "Block", group = NULL, shape = NULL), color = unique(col.per.group)[i], size = point.lwd)
            }
        }

        plot(p)
    }
    #-- End: ggplot2
    
    #internal_lattice=function(df,group,blocks,names,centroid,x0,y0,col.per.group,title,X.label,Y.label,lim.X,xlim,lim.Y,ylim,class.object,
    #col,display.names,legend,abline,pch.legend,cex,star,x,y,ellipse,df.ellipse)
    if (style=="lattice")
    {
        #-- Start: Lattice
        p = xyplot(y ~ x | Block, data = df, xlab = list(label=X.label,cex=size.xlabel), ylab = list(label=X.label,cex=size.ylabel),
        main = list(label = title, cex = size.title), as.table = TRUE, #as.table plot in order
        groups = if (display.names) {names} else {group},
        scales= list(x = list(relation = "free", limits = xlim,cex=size.axis),
        y = list(relation = "free", limits = ylim,cex=size.axis)),
        
        #-- Legend
        key = if(legend == TRUE)
        {
            if (!any(class.object%in%object.mint))
            {
                list(space = legend.position, title = "Legend", cex.title = size.legend.title,
                point = list(col =  col.per.group),cex=size.legend, pch = if(display.names | any(class.object%in%object.mint)) {16} else unique(df$pch.legend),text = list(levels(df$group)))
            } else {#we add the shape legend
                list(space = legend.position, cex.title = size.legend.title,
                point = list(
                col =  c(NA, col.per.group, NA, NA, rep("black", length(study.levels)))),
                cex = c(size.legend.title, rep(size.legend, length(col.per.group)), size.legend, size.legend.title, rep(size.legend,nlevels(factor(df$pch)))),
                pch = c(NA, rep(16, length(col.per.group)), NA, NA, as.numeric(levels(factor(df$pch)))),
                text = list(outcome = c("Outcome", levels(df$group), "", "Study", study.levels))
                )
            }
        } else {
            NULL
        },
        
        panel = function(x, y, subscripts, groups, display = display.names,...)
        {
            #-- Abline
            if (abline)
            {
                panel.abline(v = 0, lty = 2, col = "darkgrey")
                panel.abline(h = 0, lty = 2, col = "darkgrey")
            }
            
            #-- Display sample or row.names
            for (i in 1 : nlevels(df$group))
            {
                
                if (display)
                {
                    ltext(x = df$x[df$group == levels(df$group)[i]], y = df$y[df$group == levels(df$group)[i]],
                    labels = groups[subscripts & df$group == levels(df$group)[i]], col = "white", cex = 0)
                } else {
                    lpoints(x = df$x[df$group == levels(df$group)[i]], y = df$y[df$group == levels(df$group)[i]], col = "white", cex = 0, pch = 0)
                }
            }
            
            
            #-- color samples according to col
            for (i in unique(df$col))
            {
                if (display)
                {
                    ltext(x = df$x[subscripts] [df$col[subscripts] == i], y = df$y[subscripts] [df$col[subscripts] == i],
                    labels =  groups[subscripts & df$col == i],
                    col = df[df$col == i, ]$col, cex = df$cex[subscripts][df$col[subscripts] == i])
                } else {
                    lpoints(x = df$x[subscripts] [df$col[subscripts] == i],  y = df$y[subscripts] [df$col[subscripts] == i],
                    col = df[df$col == i, ]$col, cex = df$cex[subscripts][df$col[subscripts] == i], pch = df$pch[subscripts][df$col[subscripts] == i])
                }
            }
        }
        )
        print(p) #-- the lattice plot needs to be printed in order to display the ellipse(s)
        
        #-- centroid
        if (centroid)
        {
            panels = trellis.currentLayout(which = "panel")
            for (k in 1 : nlevels(df$Block))
            {

                other = df$Block %in% levels(df$Block)[k] #paste0("Block: ", blocks[k])
                ind = which(panels == k, arr.ind = TRUE)
                trellis.focus("panel", ind[2], ind[1], highlight = FALSE)
                
                for (i in 1 : nlevels(df$group))
                {
                    x0 = mean(df[other & df$group == levels(df$group)[i], ]$x)
                    y0 = mean(df[other & df$group == levels(df$group)[i], ]$y)
                    panel.points(x = x0, y = y0, col = unique(col.per.group)[i],
                        pch = 8, cex = df[other & df$group == levels(df$group)[i], ]$cex)
                }
            }
            trellis.unfocus()
        }
        
        
        
        #-- star
        if (star)
        {
            panels = trellis.currentLayout(which = "panel")
            for (k in 1 : nlevels(df$Block))
            {
     
                other = df$Block %in% levels(df$Block)[k]#paste0("Block: ", blocks[k])
                ind = which(panels == k, arr.ind = TRUE)
                trellis.focus("panel", ind[2], ind[1], highlight = FALSE)
                
                for (i in 1 : nlevels(df$group))
                {
                    for (q in 1: length(df[other & df$group == levels(df$group)[i]  , "x"]))
                    {
                        x0 = mean(df[other & df$group == levels(df$group)[i] , ]$x)
                        y0 = mean(df[other & df$group == levels(df$group)[i] , ]$y)
                        panel.segments(x0, y0, df[other & df$group == levels(df$group)[i], ]$x[q],
                        df[other & df$group == levels(df$group)[i], ]$y[q],
                        col = unique(col.per.group)[i], cex = df[other & df$group == levels(df$group)[i], ]$cex,
                        pch = df[other & df$group == levels(df$group)[i], ]$pch)
                    }
                }
            }
            trellis.unfocus()
        }
        
        
        
        
        
        #-- ellipse
        if (ellipse)
        {
            panels = trellis.currentLayout(which = "panel")
            for (k in 1 : nlevels(df$Block))
            {
  
                other.ellipse = df.ellipse$Block %in% levels(df$Block)[k]#paste0("Block: ", blocks[k])
                ind = which(panels == k, arr.ind = TRUE)
                trellis.focus("panel",ind[2], ind[1], highlight = FALSE)
                
                for (i in 1 : nlevels(df$group))
                {
                    panel.lines(x = df.ellipse[other.ellipse, paste0("Col", 2*(i - 1) + 1)],
                    y = df.ellipse[other.ellipse, paste0("Col", 2 * i)],
                    col = unique(col.per.group)[i])
                }
            }
            trellis.unfocus()
        }
    }
    #-- End: Lattice
    
    #internal_graphics=function(df,group,blocks,names,centroid,x0,y0,col.per.group,title,X.label,Y.label,lim.X,xlim,lim.Y,ylim,class.object,
    #col,display.names,legend,abline,pch.legend,cex,star,x,y,ellipse,df.ellipse,layout,rep.space,missing.col,...)
    if (style=="graphics")
    {
        #-- Start: graphics
        
        opar = par(c("mai","mar","usr","cxy","xaxp","yaxp"))
        
        reset.mfrow = FALSE # if set to TRUE, the algorithm ends up with  par(mfrow=reset.mfrow)
        
        #-- Define layout
        nResp = nlevels(df$Block)
        if (is.null(layout))
        {
            # check if there are enough plots in mfrow
            omfrow = par("mfrow")
            available.plots = prod(omfrow)
            if (available.plots<nResp) # if not enough plots available, we create our new plot
            {
                nRows = min(c(3, ceiling(nResp/2)))
                nCols = min(c(3, ceiling(nResp/nRows)))
                layout = c(nRows, nCols)
                par(mfrow = layout)
                
                if (nRows * nCols < nResp)
                devAskNewPage(TRUE)
                
                reset.mfrow=TRUE # we changed mfrow to suits our needs, so we reset it at the end
            }
        } else {
            if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
            stop("'layout' must be a numeric vector of length 2.")
            
            nRows = layout[1]
            nCols = layout[2]
            par(mfrow = layout)
            
            if (nRows * nCols < nResp)
            devAskNewPage(TRUE)
        }


        for (k in 1 : nlevels(df$Block))
        {
            if (legend)
            par(mai=c(1.360000, 1.093333, 1.093333, (max(strwidth(levels(df$group),"inches"))) + 0.6), xpd=TRUE)

            other = df$Block %in% levels(df$Block)[k]

            plot(df[other, "x" ], df[other, "y" ],
            type = "n", xlab = X.label, ylab = Y.label,
            xlim = c(xlim[[k]][1], xlim[[k]][2]), ylim = c(ylim[[k]][1], ylim[[k]][2]),
            cex.axis = size.axis, cex.lab = size.xlabel, lwd = point.lwd)#,...)
            
            #-- initialise plot
            if (any(class.object %in% c("ipca", "sipca", "pca", "spca", "prcomp", "splsda", "plsda")) & nlevels(df$Block) == 1 & !any(class.object%in%object.mint)) # avoid double title when only one block is plotted
            {
                titlemain = NULL
                if (ellipse)
                other.ellipse = TRUE
                
            }else{
                titlemain = levels(df$Block)[k]
                if (ellipse)
                other.ellipse = df.ellipse$Block %in% levels(df$Block)[k]
            }
            
            #add title of the 'blocks'
            title(main = titlemain, line = 1, cex.main = size.subtitle)

            #-- Display sample or row.names
            for (i in 1 : nlevels(df$group))
            {
                if (display.names)
                {
                    text(x = df[df$group == levels(df$group)[i] & other, "x"],
                    y = df[df$group == levels(df$group)[i] & other, "y"],
                    labels = df[df$group == levels(df$group)[i] & other, "names"],
                    col = "white", cex = 0,lwd=point.lwd)#,...)
                } else {
                    points(x = df[df$group == levels(df$group)[i] & other, "x"],
                    y = df[df$group == levels(df$group)[i] & other, "y"],
                    col = "white", cex = 0, pch = 0,lwd=point.lwd)#,...)
                }
            }
            
            #-- color samples according to col
            for (i in unique(df$col))
            {
                if (display.names)
                {
                    text(x = df[df$col == i & other, "x"],
                    y = df[df$col == i & other, "y"],
                    labels = df[df$col == i & other, "names"],
                    col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex,lwd=point.lwd)#,...)
                } else {
                    points(x = df[df$col == i & other, "x"],
                    y = df[df$col == i & other, "y"],
                    col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex, pch = df[df$col == i, ]$pch,lwd=point.lwd)#,...)
                }
            }
            
            pch.legend = NULL
            if (missing.col)
            {
                
                for (i in 1:nlevels(factor(df$col)))
                pch.legend = c(pch.legend, df[df$col == levels(factor(df$col))[i], ]$pch)
            } else {
                for (i in 1:nlevels(df$group))
                pch.legend = c(pch.legend, df[df$group == levels(df$group)[i], ]$pch)
            }
            
            if (legend)
            {
                legend(par()$usr[2]+0.1, par()$usr[4] - (par()$usr[4]-par()$usr[3])/2, col = col.per.group, legend = levels(df$group), pch = if(display.names) {16} else unique(df$pch.legend), title = 'Legend', cex = size.legend, lty = 0,lwd = point.lwd)
                
            }
            if (legend)
            par(xpd=FALSE) # so the abline does not go outside the plot
            
            #-- Abline
            if (abline)
            abline(v = 0, h = 0, lty = 2,lwd=point.lwd)#,...)
            
            #-- Star
            if (star == TRUE)
            {
                for (i in 1 : nlevels(df$group))
                {
                    x0 = mean(df[df$group == levels(df$group)[i] & other, "x"])
                    y0 = mean(df[df$group == levels(df$group)[i] & other, "y"])
                    for (q in 1: length(df[df$group == levels(df$group)[i] & other, "x"]))
                    {
                        segments(x0, y0, df[df$group == levels(df$group)[i] & other, "x"][q], df[df$group == levels(df$group)[i] & other, "y"][q],
                        cex = df$df[df$group == levels(df$group)[i] & other, "cex"], col = df[df$group == levels(df$group)[i] & other, "col"], lwd = point.lwd)#,...)
                    }
                }
            }
            
            #-- Centroid
            if (centroid == TRUE)
            {
                for (i in 1 : nlevels(df$group))
                {
                    x0 = mean(df[df$group == levels(df$group)[i] & other, "x"])
                    y0 = mean(df[df$group == levels(df$group)[i] & other, "y"])
                    points(cbind(x0,y0), pch = 8, cex = df$df[df$group == levels(df$group)[i] & other, "cex"],
                    col = unique(col.per.group)[i], lwd = point.lwd)#,...)
                }
            }
            
            
            #-- Ellipse
            if (ellipse == TRUE)
            {
                for (i in 1 : nlevels(df$group))
                {
                    lines(x = df.ellipse[other.ellipse, paste0("Col", 2*(i - 1) + 1)],
                    y = df.ellipse[other.ellipse, paste0("Col", 2 * i)],
                    col = unique(col.per.group)[i],lwd=point.lwd)#,...)
                }
            }
            
            if (any(class.object %in% c("ipca", "sipca", "pca", "spca", "prcomp", "splsda", "plsda")) & nlevels(df$Block)==1 & !any(class.object %in% object.mint) )
            {
                title(title, line = 1, cex.main = size.title)#,...)
            } else {
                title(title, outer=TRUE, line = -2,cex.main = size.title)#,...)
            }
        }
        
        
        #par(opar)
        #par(usr=opar["usr"])
        #par(xaxp=opar["xaxp"])
        #par(yaxp=opar["yaxp"])

        #par(mai=opar["mai"])
        if (reset.mfrow)
        par(mfrow = omfrow)
        
        #par(mar=opar["mar"])
        
    }
    #-- End: graphics
    
    
    #internal_3d=function(df,group,blocks,names,centroid,x0,y0,col.per.group,title,X.label,Y.label,lim.X,xlim,lim.Y,ylim,class.object,
    #col,display.names,legend,abline,pch.legend,cex,star,x,y,ellipse,df.ellipse,axes.box,Z.label,z)
    if (style=="3d")
    {
        
        #-- Start: 3d
        
        for (k in 1 : nlevels(df$Block))
        {
            #if (nlevels(df$Block)>1) # removing the popping up window when there's only one block (for shiny)
            open3d()
            
            par3d(windowRect = c(500, 30, 1100, 630))
            Sys.sleep(0.1)
            
            if (!is.null(title))
            {
                mat = matrix(1:2, 2)
                layout3d(mat, heights = c(1, 10), model = "inherit")
                next3d()
                text3d(0, 0, 0, title)
                next3d()
            }
            
            if(any(class.object %in% c("ipca", "sipca", "pca", "spca", "prcomp", "splsda", "plsda", "mlsplsda")))
            {
                other = TRUE
                if (ellipse)
                other.ellipse = TRUE
            } else {
                other = df$Block %in% levels(df$Block)[k]
                if (ellipse)
                other.ellipse = df.ellipse$Block %in% levels(df$Block)[k]
            }
            
            par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
            
            if (legend)
            {
                legend3d(x = "right",
                legend = levels(df$group),
                col = col.per.group,
                pch = rep(16,length(unique(df$pch))),
                pt.cex = unique(df$cex),
                bty = "n")
            }
            
            #-- Display sample or row.names
            
            for (i in unique(df$col))
            {
                if (display.names)
                {
                    text3d(x = df[df$col == i &  other, "x"],
                    y = df[df$col == i &  other, "y"],
                    z = df[df$col == i &  other, "z"],
                    texts = df[df$col == i &  other, "names"],
                    color = df[df$col == i, ]$col, size = unique(df[df$col == i, ]$cex))
                }else{
                    cex=unique(df[df$col == i, ]$cex)*20
                    switch(unique(df[df$col == i, ]$pch),
                    sphere = plot3d(x = df[df$col == i & other, "x"],
                    y = df[df$col == i & other, "y"],
                    z = df[df$col == i & other, "z"],
                    col = df[df$col == i, ]$col, size = cex, radius = cex/20, add = TRUE),
                    tetra = shapelist3d(tetrahedron3d(), x = df[df$col == i &other, "x"],
                    y = df[df$col == i & other, "y"],
                    z = df[df$col == i & other, "z"],
                    col = df[df$col == i, ]$col, size = cex/25),
                    cube = shapelist3d(cube3d(),x = df[df$col == i & other, "x"],
                    y = df[df$col == i & other, "y"],
                    z = df[df$col == i & other, "z"],
                    col = df[df$col == i, ]$col, size = cex/30),
                    octa = shapelist3d(octahedron3d(), x = df[df$col == i & other, "x"],
                    y = df[df$col == i & other, "y"],
                    z = df[df$col == i & other, "z"],
                    col = df[df$col == i, ]$col, size = cex/17),
                    icosa = shapelist3d(icosahedron3d(), x = df[df$col == i & other, "x"],
                    y = df[df$col == i & other, "y"],
                    z = df[df$col == i &other, "z"],
                    col = df[df$col == i, ]$col, size = cex/20),
                    dodeca = shapelist3d(dodecahedron3d(), x = df[df$col == i &other, "x"],
                    y = df[df$col == i & other, "y"],
                    z = df[df$col == i & other, "z"],
                    col = df[df$col == i, ]$col, size = cex/20))
                }
            }
            
            #-- Ellipse
            if (ellipse)
            {
                coords = matrix(cbind(df[other, "x"],
                df[other, "y"],
                df[other,"z"]),ncol = 3)
                centr.coords = apply(coords, 2, function(x) tapply(x, df$group, mean))
                if (length(unique(df$group)) == 1)
                centr.coords = matrix(centr.coords, nrow=1)
                
                rownames(centr.coords) = levels(df$group)
                lg = levels(df$group)
                for(i in 1:length(lg))
                {
                    g   = lg[i]
                    sel = df$group == g
                    s   = cov(coords[sel, , drop = FALSE])
                    cc  = centr.coords[i,]
                    # lines(ellipse(s, centre=cc), col=unique(col.per.group)[i])
                    shade3d(ellipse3d(s, centre = cc, level = df.ellipse$ellipse.level[1]), col = unique(col.per.group)[i], alpha = alpha)
                    
                }
            }
            
            #-- Centroid
            if (centroid == TRUE)
            {
                for (i in 1 : nlevels(df$group))
                {
                    x0 = mean(df[df$group == levels(df$group)[i] & other, "x"])
                    y0 = mean(df[df$group == levels(df$group)[i] & other, "y"])
                    z0 = mean(df[df$group == levels(df$group)[i] & other, "z"])
                    points3d(x=x0, y=y0,z=z0, cex=df$df[df$group == levels(df$group)[i] & other, "cex"], col = unique(col.per.group)[i])
                }
            }
            
            #-- Star
            if (star == TRUE)
            {
                for (i in 1 : nlevels(df$group))
                {
                    x0 = mean(df[df$group == levels(df$group)[i] & other, "x"])
                    y0 = mean(df[df$group == levels(df$group)[i] & other, "y"])
                    z0 = mean(df[df$group == levels(df$group)[i] & other, "z"])
                    for (q in 1: length(df[df$group == levels(df$group)[i] & other, "x"]))
                    {
                        segments3d(x=c(x0, df[df$group == levels(df$group)[i] & other, "x"][q]), y=c(y0,df[df$group == levels(df$group)[i] & other, "y"][q]),
                        z=c(z0, df[df$group == levels(df$group)[i] & other, "z"][q]),
                        cex=df$df[df$group == levels(df$group)[i] & other, "cex"], col=df[df$group == levels(df$group)[i] & other, "col"])
                    }
                }
            }
            
            
            #-- draws axes/box --#
            if (axes.box == "box")
            {
                axes3d(marklen = 25)
                box3d()
            }
            if (axes.box == "bbox")
            {
                bbox3d(color = c("#333377", "black"), emission = gray(0.5),
                specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)
            }
            if (axes.box == "both")
            {
                axes3d(marklen = 25); box3d()
                bbox3d(color = c("#333377", "black"), emission = gray(0.5),
                specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)
            }
            #-- add axes labels --#
            mtext3d(X.label, "x-+", line = 1)
            mtext3d(Y.label, "y-+", line = 1.5)
            mtext3d(Z.label, "z+-", line = 1)
            if (! any(class.object%in% c("ipca", "sipca", "pca", "spca", "prcomp", "splsda", "plsda", "mlsplsda")))
            title3d(main = levels(df$Block)[k])
        }
        #-- output --#
        return(invisible(cbind(df$x, df$y, df$z)))
    }
    
    if (style%in%c("graphics","3d"))
    p = NULL
    
    return(p)
}



# Copyright (C) 2015
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
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


# ---------------------------------------------------------

# to do later
# create a S3 function?

# would need to expaned for other types of objects
# for sPCA, then a Y 'outcome' should be provided as input
# for sPLS, then block = X or Y should be provided as input
# for sGCCA, then block = k should be provided as input


plotContrib = function(object,
                       contrib=c("max"),  # choose between 'max" or "min"
                       method = c("mean", "median"), 
                       comp = 1,
                       ties = TRUE, 
                       ndisplay = NULL, 
                       cex.name = 0.7, 
                       cex.legend = 0.8,
                       name.var = NULL,
                       legend = TRUE,
                       legend.color = NULL, 
                       main = NULL,
                       legend.title = 'Outcome',
                       plot = TRUE 
                       ) {
  
  # -------------------
  # input checks 
  # ------------------
  
  if(!any(class(object)=="plsda")) stop("plotContrib is only available for 'plsda' and 'splsda' objects for now")

  
  # method
  # ----
  if(length(method) == 2){
    method = 'median'
    warning('Argument method: we expect either mean or median, set to median by default')
  }
  # cex
  # --
  if(cex.name <= 0){cex.name = 0.7}
  if(cex.legend <= 0){cex.name = 0.8}
  
  #comp
  # ----
  if(length(comp)>1){
    stop('Comp should have a length of 1')
  }
  if(comp > object$ncomp){
    stop(paste('Argument comp should be less or equal to', object$ncomp))
  }
  
  # ndisplay
  # ------
  # if null set by default to all variables
  if(is.null(ndisplay)){
    ndisplay = length(selectVar(object, comp = comp)$name)
  }else if(ndisplay > length(selectVar(object, comp = comp)$name)){
    stop('ndisply value is larger than the number of variables that you have selected!')
  }  
  # check on name.var
  # -----------------
  if(!is.null(name.var)){
    #check we have a vector or factor
    if(!is.null(ncol(name.var))){
      stop('name.var should be a vector or a factor')
    }
    # we match the names of name.var to the X colnames
    if(is.null(names(name.var))){
      if(length(name.var) == ncol(object$X)){
        warning('The names of the vector name.var will be set as the colnames of the data frame X')
        names(name.var) = colnames(object$X)
      }else{
        stop('The vector name.var has the wrong length and has no names, name your vector with your variable names.')
      }
    }
  } # end checked when name.var provided
  
  #title
  #-----
  if(!is.null(main) & !is.character(main)){
    warning('main needs to be of type character')
  }
  
  #legend.color
  #-----
  if(!is.null(legend.color) & (length(legend.color) != nlevels(factor(map(object$ind.mat), labels= object$names$Y)))){
    warning('legend.color must be the same length than the number of group, by default set to default colors')
    legend.color = color.mixo(1:10)  # by default set to the colors in color.mixo (10 colors)
  }
  if(is.null(legend.color)) legend.color = color.mixo(1:10) # by default set to the colors in color.mixo (10 colors)
  
  # --------------------
  # end check inputs
  # ---------------------
  
  # ==================================================
  # First, calculate the contribution of the loading
  # =================================================
  # Retrieve selected variables
  list.select = selectVar(object, comp = comp)$name
  X = object$X
  Y = factor(map(object$ind.mat), labels= object$names$Y)
  
  # Start: Initialisation
  which.comp = method.group = list()
  which.contrib = data.frame(matrix(FALSE, 
                                ncol = nlevels(Y) + 2, nrow = length(list.select), 
                                dimnames = list(list.select, c(paste0("Contrib.", levels(Y)), "Contrib", "GroupContrib"))))
  # End: Initialisation
  
  # calculate the max.method per group for each variable, and identifies which group has the max max.method
  for(k in list.select){
    if (method == 'mean'){
      method.group[[k]] = tapply(X[, k], Y, mean,na.rm=TRUE)
    } else if (method == 'median'){
      method.group[[k]] = tapply(X[, k], Y, median,na.rm=TRUE)
    }
    # determine which group has the highest mean/median
    if(contrib=="max")
    {
        which.contrib[k, 1:nlevels(Y)] = abs(method.group[[k]]) == max(abs(method.group[[k]]))
    }else if(contrib=="min")
    {
        which.contrib[k, 1:nlevels(Y)] = abs(method.group[[k]]) == min(abs(method.group[[k]]))
    }
  }
  
  # if ties, we set the color to white
  which.contrib$Contrib = apply(which.contrib, 1, function(x){if (length(which(x)) > 1){
    return("white")
  } else { # otherwise we use legend color provided
    ##return(color.mixo(1:nlevels(Y))[which(x)])
    return(legend.color[1:nlevels(Y)][which(x)])
  }})
 
  # we also add an output column indicating the group that is max
  which.contrib$GroupContrib = apply(which.contrib[, 1:(nlevels(Y))], 1, function(x){if (length(which(x)) > 1){
    return("tie")
  } else { 
    return(levels(Y)[which(x)])
  }})

  method.group = do.call(rbind, method.group)
  
  
  #   # run a kruskall wallis test and adjust the pvalues for multiple testing using Benjamini Hochberg
  # this has been removed as the type of stat test is particular to the data
  #   kruskall.pval = apply(X[, list.select, drop = FALSE], 2, function(x){kruskal.test(x ~ Y)$p.value})
  #   kruskall.pval.adj = p.adjust(kruskall.pval, method = 'BH')
  
  if (!is.null(name.var)){
    #removed output KW test
    contrib = data.frame(name = as.character(name.var[list.select]), method.group, which.contrib, importance = selectVar(object, comp = comp)$value)
  } else {
    contrib = data.frame(method.group, which.contrib, importance = selectVar(object, comp = comp)$value)
  }
  # End contribution calculation
  
  # =====================================
  # then determine the colors/groups matching max contribution
  # =======================================
  
  # when working with sparse counts in particular and using the median to measure contribution
  # ties to determine the contribution of a variable may happen, in that case remove them, otherwise they are showns as blank
  if (ties == TRUE){
    contrib = contrib[!contrib$Contrib %in% "white", ]
  }
  contrib = contrib[1 : min(nrow(contrib), ndisplay), , drop = FALSE]
  
  #display barplot with names of variables
  if(!is.null(name.var)){
    name.var.output = as.character(contrib$name) 
  }else{
    name.var.output = row.names(contrib)
  }
  
  # ==================================
  # represent the barplot
  # ==================================
  #added condition if all we need is the contribution stats
  if(plot){
  # in case the user does not want to show the legend, then margins can be reduced
  if(legend){
    layout(matrix(c(1, 2), 1, 2, byrow = TRUE), c(0.7, 0.3), TRUE)
  }else{
    layout(matrix(c(1, 2), 1, 2, byrow = TRUE), c(0.7, 0.1), TRUE)
  }
  # barplot with contributions
  par(mar = c(5, min(9, max(sapply(name.var.output, nchar))), 4, 0))
  mp = barplot(contrib$value.var, horiz = T, las = 1, col = contrib$Contrib, axisnames = TRUE, names.arg = name.var.output, #names.arg = row.names(contrib),
               cex.names = cex.name, cex.axis = 0.7, beside = TRUE,border=NA)
  if(is.null(main)){
    title(paste('Contribution on comp', comp))
  }else{
    title(paste(main))
  }
  
  # legend
  if(legend){
    par(mar = c(5, 0, 4, 3) + 0.1)
    plot(1,1, type = "n", axes = FALSE, ann = FALSE)
    legend(0.8, 1, col = legend.color[1:nlevels(Y)], legend = levels(Y), pch = 19, 
           title = paste(legend.title), 
           cex = cex.legend)
  }
  par(mfrow=c(1,1)) 
  } # end if plot
  
  # ===================================
  # return the contribution matrix
  # ===================================
  return(invisible(list(contrib = contrib)))
}



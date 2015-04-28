# Copyright (C) 2009 
# Kim-Anh Le Cao, 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
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


# object: a pls/spls object
# comp: to display the variables selected on dimension 'comp'
# names.X, names.Y: set to true means that the X and Y data frames have row names (see example below with srbct)

selectVar <-
function(...) UseMethod("selectVar")



# ------------------ for sPLS object --------------------
selectVar.spls = function(object, comp =1, ...){ 
  
  if(comp > object$ncomp) stop('The comp value you indicated is larger than the fitted model')
  
  # variables from data set X
  # name of selected variables
  name.var.X = names(sort(abs(object$loadings$X[,comp]), decreasing = T)[1:object$keepX[comp]])
  #value on the loading vector
  value.var.X = object$loadings$X[name.var.X,comp]
  
  # variables from data set Y
  # name of selected variables
  name.var.Y = names(sort(abs(object$loadings$Y[,comp]), decreasing = T)[1:object$keepY[comp]])
  #value on the loading vector
  value.var.Y = object$loadings$Y[name.var.Y,comp]
  
  return(
    list(name.X = name.var.X, value.X = data.frame(value.var.X), name.Y = name.var.Y, value.Y = data.frame(value.var.Y), comp = comp)
  )
  
}


# ------------------ for sPLS-DA object --------------------
selectVar.splsda =  function(object, comp=1, ...){ 
  
  if(comp > object$ncomp) stop('The comp value you indicated is larger than the fitted model')
  
  # variables from data set X
  name.var = names(sort(abs(object$loadings$X[,comp]), decreasing = T)[1:object$keepX[comp]])
  #value on the loading vector
  value.var = object$loadings$X[name.var,comp]
  
  return(
    list(name = name.var, value = data.frame(value.var), comp = comp)
  )
}



# ------------------ for sPCA object --------------------
selectVar.spca = function(object, comp=1, ...){ 
  
  if(comp > object$ncomp) stop('The comp value you indicated is larger than the fitted model')
  
  # variables from data set X
  name.var = names(sort(abs(object$rotation[,comp]), decreasing = T)[1:object$keepX[comp]])
  #value on the loading vector
  value.var = object$rotation[name.var,comp]
  
  return(
    list(name = name.var, value = data.frame(value.var), comp = comp)
  )
}

# ------------------ for siPCA object --------------------
selectVar.sipca = function(object, comp=1, ...){ 
  
  if(comp > object$ncomp) stop('The comp value you indicated is larger than the fitted model')
  
  # variables from data set X
  name.var = names(sort(abs(object$loadings[,comp]), decreasing = T)[1:object$keepX[comp]])
  #value on the loading vector
  value.var = object$loadings[name.var,comp]
  
  return(
    list(name = name.var, value = data.frame(value.var), comp = comp)
  )
}



# ------------------ for sgcca object --------------------

selectVar.sgcca = function(object, block = NULL, comp = 1, ...){ 
  
  # check arguments
  # -----------------
  if(length(comp) > 1)
    stop("Expecting one single value for 'comp'")
  if(is.null(block)){
    if(any(comp > object$ncomp))
      stop("'comp' is greater than the number of components in the fitted model,
         you need to specify the variable 'block' ")
  }else{
    if(any(comp > object$ncomp[block]))
      stop("'comp' is greater than the number of components in the fitted model for the block you specified")
    
  }
  
  
  # extract selected variables
  # --------------------------
  keep = list()
  name.var = value.var = list()
  if(is.null(block)){
    # identify the selected variables selected on a given component comp
    for(k in 1:length(object$data)){
      keep[[k]] = abs(object$loadings[[k]][,comp])> 0      
    }   
    #store name and value of the selected variables
    for(k in 1:length(object$data)){
      name.var[[k]] = names(which(keep[[k]] == TRUE))   #object$names$var[keep[[k]]]
      value.var[[k]] = object$loadings[[k]][keep[[k]], comp]
    }
  }else{ #end is.null(block)
    j=1
    # identify the selected variables selected on a given component ncomp.select
    for(k in block){
      keep[[j]] = abs(object$loadings[[k]][,comp])> 0
      j=j+1        
    }   
    l=1
    for(k in block){
      name.var[[l]] = names(which(keep[[l]] == TRUE))  
      value.var[[l]] = object$loadings[[k]][keep[[l]], comp]
      l = l+1
    }
  } # end is.null (block)
  
  return(
    list(name.var = name.var, value.var = value.var, comp = comp)
  )
  
}


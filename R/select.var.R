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

select.var <-
function(object, ...) UseMethod("select.var")




# # ------------------ for PLS object --------------------
# select.var.pls = function(object, comp){
#   
#   if(comp > object$ncomp) stop('The comp value you indicated is larger than the fitted model')
#   
#   # basically here we output all variables ranked by absolute value of the loading
#   # variables from data set X
#   # name of selected variables
#   name.var.X = names(sort(abs(object$loadings$X[,comp]), decreasing = T)[1:ncol(object$X)])
#   #value on the loading vector
#   value.var.X = object$loadings$X[name.var.X,comp]
#   
#   # variables from data set Y
#   # name of selected variables
#   name.var.Y = names(sort(abs(object$loadings$Y[,comp]), decreasing = T)[1:ncol(object$Y)])
#   #value on the loading vector
#   value.var.Y = object$loadings$Y[name.var.Y,comp]
#   
#   return(
#     list(name.X = name.var.X, value.X = data.frame(value.var.X), name.Y = name.var.Y, value.Y = data.frame(value.var.Y), comp = comp)
#   )
# }


# ------------------ for sPLS object --------------------
select.var.spls = function(object, comp =1, ...){ 
  
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

# # ------------------ for PLS-DA object --------------------
# select.var.plsda = function(object, comp, names = TRUE){
#   
#   if(comp > object$ncomp) stop('The comp value you indicated is larger than the fitted model')
#   
#   
#   # variables from data set X
#   # name of selected variables
#   name.var = names(sort(abs(object$loadings$X[,comp]), decreasing = T)[1:ncol(object$X)])
#   #value on the loading vector
#   value.var = object$loadings$X[name.var,comp]
#   
#   return(
#     list(name = name.var, value = data.frame(value.var), comp = comp)
#   )
# }


# ------------------ for sPLS-DA object --------------------
select.var.splsda =  function(object, comp=1, ...){ 
  
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
select.var.spca = function(object, comp=1, ...){ 
  
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
select.var.sipca = function(object, comp=1, ...){ 
  
  if(comp > object$ncomp) stop('The comp value you indicated is larger than the fitted model')
  
  # variables from data set X
  name.var = names(sort(abs(object$loadings[,comp]), decreasing = T)[1:object$keepX[comp]])
  #value on the loading vector
  value.var = object$loadings[name.var,comp]
  
  return(
    list(name = name.var, value = data.frame(value.var), comp = comp)
  )
}


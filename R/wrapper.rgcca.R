# Copyright (C) 2013 
# Kim-Anh Le Cao, University of Queensland, Brisbane, Australia
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

wrapper.rgcca = function(
  data,
  design = 1 - diag(length(data)),
  tau = rep(1, length(data)),
  ncomp = rep(1, length(data)), 
  scheme = "centroid",
  scale = TRUE, 
  init = "svd", 
  bias = TRUE,
  tol = .Machine$double.eps, 
  verbose = FALSE  
){
  
  # call function
  #rgcca <- function(A, C = 1-diag(length(A)), tau = rep(1, length(A)), ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE , init="svd", bias = TRUE, tol = .Machine$double.eps, verbose=TRUE)
    
  result.rgcca = rgcca(A = data, C = design, tau = tau, 
                       ncomp = ncomp,
                       scheme = scheme, scale = scale,
                       init = init, bias = bias, tol = tol, verbose = verbose)
  
  # outputs
#   out <- list(Y = shave.matlist(Y, ncomp),
#               a = shave.matlist(a, ncomp), 
#               astar = shave.matlist(astar, ncomp),
#               C = C, tau = tau_mat, scheme = scheme,
#               ncomp=ncomp, crit = crit,
#               mode = mode,
#               AVE=list(AVE_X=AVE_X,
#                        AVE_outer=AVE_outer,
#                        AVE_inner=AVE_inner),
#               #KA added names of rows and cols for plotIndiv and plotVar
#               names = list(indiv = rownames(A[[1]]))
#   )
#   class(out) <- "rgcca"
#   return(out)
  
  cl = match.call()
  cl[[1]] = as.name('rgcca')
  
  output = list(
    class = cl,
    data = data,
    variates = result.rgcca$Y,
    loadings = result.rgcca$a,
    loadings.star = result.rgcca$astar,
    design = design,
    tau = result.rgcca$tau,
    scheme = scheme,
    ncomp = ncomp, 
    crit = result.rgcca$crit,
    AVE = list(AVE.data = result.rgcca$AVE$AVE_X, result.rgcca$AVE$AVE_outer, result.rgcca$AVE$AVE_inner), #rename?
    names = list(indiv = rownames(data[[1]]), var = sapply(data, colnames))
    
  )

  class(output) = 'rgcca'
  return(invisible(output))
  
}

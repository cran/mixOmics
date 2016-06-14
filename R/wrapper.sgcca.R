#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2013
# last modified: 01-03-2016
#
# Copyright (C) 2013
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



wrapper.sgcca = function(
X,
design = 1 - diag(length(X)),
penalty = NULL,
ncomp = 1,
keepX.constraint,
keepX,
scheme = "centroid",
mode = "canonical",
scale = TRUE,
bias = TRUE,
init = "svd",
tol = .Machine$double.eps,
verbose = FALSE,
max.iter = 1000,
near.zero.var = FALSE
){
    
    
    check=Check.entry.sgcca(X = X, design = design ,ncomp = ncomp , scheme = scheme , scale = scale ,  bias = bias,
    init = init , tol = tol , verbose = verbose,mode = mode, max.iter = max.iter,near.zero.var = near.zero.var,keepX = keepX,keepX.constraint = keepX.constraint)
    
    
    A = check$A
    design = check$design
    ncomp = check$ncomp
    init = check$init
    scheme = check$scheme
    verbose = check$verbose
    bias = check$bias
    near.zero.var = check$near.zero.var
    keepA.constraint = check$keepA.constraint
    keepA = check$keepA
    nzv.A = check$nzv.A
    
    
    result.sgcca = internal_mint.block(A = A, design = design, tau = NULL,
    ncomp = ncomp,
    scheme = scheme, scale = scale,
    init = init, bias = bias, tol = tol, verbose = verbose,
    keepA.constraint = keepA.constraint,
    keepA = keepA,
    max.iter = max.iter,
    study = factor(rep(1,nrow(A[[1]]))),#mint.sgcca not coded yet
    mode = mode,penalty = penalty
    )
    
   
    out = list(
    call = match.call(),
    X = result.sgcca$X,
    variates = result.sgcca$variates,
    loadings = result.sgcca$loadings,
    loadings.star = result.sgcca$loadings.star,
    design = result.sgcca$design,
    penalty = penalty,
    scheme = result.sgcca$scheme,
    ncomp = result.sgcca$ncomp,
    crit = result.sgcca$crit,
    AVE = list(AVE.X = result.sgcca$AVE$AVE_X, result.sgcca$AVE$AVE_outer, result.sgcca$AVE$AVE_inner), #rename?
    names = result.sgcca$names,#names = list(indiv = rownames(X[[1]]), var = sapply(X, colnames)),
    init = result.sgcca$init,
    bias = result.sgcca$bias,
    tol = result.sgcca$tol,
    iter = result.sgcca$iter,
    max.iter = result.sgcca$max.iter,
    nzv = result.sgcca$nzv,
    scale = result.sgcca$scale,
    design = result.sgcca$design,
    scheme = result.sgcca$scheme,
    explained_variance = result.sgcca$explained_variance
    )
    
    class(out) = 'sgcca'
    return(invisible(out))
}


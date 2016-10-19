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


wrapper.rgcca = function(
X,
design = 1 - diag(length(X)),
tau = rep(1, length(X)),
ncomp = 1,
keepX.constraint,
keepX,
scheme = "horst",
scale = TRUE,
bias = FALSE,
init = "svd.single",
tol = .Machine$double.eps,
verbose = FALSE,
max.iter = 1000,
near.zero.var = FALSE)
{


    check = Check.entry.rgcca(X = X, design = design, tau = tau, ncomp = ncomp, scheme = scheme, scale = scale,
    init = init, bias = bias, tol = tol, verbose = verbose, max.iter = max.iter, near.zero.var = near.zero.var,keepX = keepX,
    keepX.constraint = keepX.constraint)
    X = check$A
    ncomp = check$ncomp
    design = check$design
    init = check$init
    scheme = check$scheme
    verbose = check$verbose
    bias = check$bias
    nzv.A = check$nzv.A
    keepA.constraint = check$keepA.constraint
    keepA = check$keepA



    result.rgcca = internal_mint.block(A = X, design = design, tau = tau,
    ncomp = ncomp,
    scheme = scheme, scale = scale,
    init = init, bias = bias, tol = tol, verbose = verbose,
    keepA.constraint = keepA.constraint, keepA = keepA,
    max.iter = max.iter,
    study = factor(rep(1,nrow(X[[1]]))),#mint.rgcca not coded yet
    mode = "canonical"
    )
    

    out = list(
    call = match.call(),
    X = result.rgcca$X,
    variates = result.rgcca$variates,
    loadings = result.rgcca$loadings,
    loadings.star = result.rgcca$loadings.star,
    keepX=result.rgcca$keepA,keepX.constraint = result.rgcca$keepA.constraint,
    design = result.rgcca$design,
    tau = result.rgcca$tau,
    scheme = result.rgcca$scheme,
    ncomp = result.rgcca$ncomp,
    crit = result.rgcca$crit,
    AVE = list(AVE.X = result.rgcca$AVE$AVE_X, result.rgcca$AVE$AVE_outer, result.rgcca$AVE$AVE_inner), #rename?
    names = result.rgcca$names,#names = list(indiv = rownames(X[[1]]), var = sapply(X, colnames)),
    init = result.rgcca$init,
    bias = result.rgcca$bias,
    tol = result.rgcca$tol,
    iter = result.rgcca$iter,
    max.iter = result.rgcca$max.iter,
    nzv = result.rgcca$nzv,
    scale = result.rgcca$scale,
    design = result.rgcca$design,
    scheme = result.rgcca$scheme,
    explained_variance = result.rgcca$explained_variance
    )
    
    class(out) = c("sparse.rgcca","rgcca")
    return(invisible(out))
}


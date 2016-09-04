#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2015
# last modified: 24-02-2016
#
# Copyright (C) 2015
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


# ========================================================================================================
# mint.splsda: perform a vertical sPLS-DA on a combination of experiments, input as a matrix in X
# this function is a particular setting of internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint, which then call 'internal_mint.block'
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 2.
# study: grouping factor indicating which samples are from the same study
# keepX.constraint: A list containing which variables of X are to be kept on each of the first PLS-components.
# keepX: number of \eqn{X} variables kept in the model on the last components (once all keepX.constraint[[i]] are used).
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations


mint.splsda = function(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
keepX.constraint=NULL,
keepX = rep(ncol(X), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE)
{
    
    #-- validation des arguments --#
    # most of the checks are done in 'internal_wrapper.mint'
    if (is.null(Y))
    stop("'Y' has to be something else than NULL.")
    
    if (is.null(dim(Y)))
    {
        Y = factor(Y)
    }  else {
        stop("'Y' should be a factor or a class vector.")
    }
    
    Y.mat = unmap(Y)
    colnames(Y.mat) = levels(Y)

    # call to 'internal_wrapper.mint'
    result = internal_wrapper.mint(X = X, Y = Y.mat, ncomp = ncomp, near.zero.var = near.zero.var, study = study, mode = mode,
    keepX = keepX, keepX.constraint = keepX.constraint, max.iter = max.iter, tol = tol, scale = scale)
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$X[-result$indY][[1]],
        Y = Y,
        ind.mat = result$X[result$indY][[1]],
        ncomp = result$ncomp,
        study = result$study,
        mode = result$mode,
        keepX = result$keepA[[1]],
        keepY = result$keepA[[2]],
        keepX.constraint = result$keepA.constraint[[1]],
        keepY.constraint = result$keepA.constraint[[2]],
        variates = result$variates,
        loadings = result$loadings,
        variates.partial = result$variates.partial,
        loadings.partial = result$loadings.partial,
        names  =  result$names,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = result$scale,
        explained_variance = result$explained_variance)

    class(out) = c("mint.splsda","splsda","spls","DA")
    return(invisible(out))
    
}

#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2009
# last modified: 24-05-2016
#
# Copyright (C) 2009
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
# spls: perform a sPLS
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

# X: numeric matrix of predictors
# Y: numeric vector or matrix of responses
# ncomp: the number of components to include in the model. Default to 2.
# keepX.constraint: A list containing which variables of X are to be kept on each of the first PLS-components.
# keepY.constraint: A list containing which variables of Y are to be kept on each of the first PLS-components
# keepX: number of \eqn{X} variables kept in the model on the last components (once all keepX.constraint[[i]] are used).
# keepY: number of \eqn{Y} variables kept in the model on the last components.
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations


spls = function(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
keepX.constraint=NULL,
keepY.constraint=NULL,
keepX,
keepY,
scale = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
logratio = "none",   # one of "none", "CLR"
multilevel = NULL)    # multilevel is passed to multilevel(design = ) in withinVariation. Y is ommited and shouldbe included in multilevel design
{

    # call to 'internal_wrapper.mint'
    result = internal_wrapper.mint(X = X, Y = Y, ncomp = ncomp, scale = scale, near.zero.var = near.zero.var, mode = mode,
    keepX = keepX, keepY = keepY, keepX.constraint = keepX.constraint, keepY.constraint = keepY.constraint, max.iter = max.iter,
    tol = tol, logratio = logratio,
    multilevel = multilevel, DA = FALSE)
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$X[-result$indY][[1]],
        Y = result$X[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        keepX = result$keepA[[1]],
        keepY = result$keepA[[2]],
        keepX.constraint = result$keepA.constraint[[1]],
        keepY.constraint = result$keepA.constraint[[2]],
        variates = result$variates,
        loadings = result$loadings,
        names = result$names,
        tol = result$tol,iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = scale,
        logratio = logratio,
        explained_variance = result$explained_variance,
        input.X = result$input.X,
        mat.c = result$mat.c,
        defl.matrix = result$defl.matrix
        )
    
   
    class(out) = c("spls")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        out$multilevel = multilevel
        class(out) = c("mlspls",class(out))
    }

    return(invisible(out))
}



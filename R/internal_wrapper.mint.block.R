#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2015
# last modified: 13-04-2016
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


# perform the mint.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects

# ========================================================================================================
# internal_wrapper.mint.block: this function is a particular setting of internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint.block
# ========================================================================================================
# used in (mint).block approaches

internal_wrapper.mint.block = function(X,
Y,
indY,
study,
ncomp,
keepX.constraint,
keepY.constraint,
keepX,
keepY,
design,
scheme,
mode,
scale = TRUE,
bias,
init ,
tol = 1e-06,
verbose,
max.iter = 100,
near.zero.var = FALSE)
{
    if (missing(scheme))
    scheme= "horst"
    
    if (missing(bias))
    bias= FALSE
    
    if (missing(verbose))
    verbose= FALSE
    
    if (missing(mode))
    mode="regression"
    
    # checks (near.zero.var is done there)
    check=Check.entry.wrapper.mint.block(X = X, Y = Y, indY = indY, ncomp = ncomp,
    keepX = keepX, keepX.constraint = keepX.constraint,
    keepY = keepY, keepY.constraint = keepY.constraint,
    study = study, design = design, init = init, scheme = scheme, scale = scale,
    bias = bias, near.zero.var = near.zero.var, mode = mode, tol = tol,
    max.iter = max.iter, verbose = verbose)

    # get some values after checks
    A = check$A
    indY = check$indY
    study = check$study
    design = check$design
    ncomp = check$ncomp
    keepA = check$keepA
    keepA.constraint = check$keepA.constraint
    init = check$init
    nzv.A = check$nzv.A
    
    # A: list of matrices
    # indY: integer, pointer to one of the matrices of A
    # design: design matrix, links between matrices. Diagonal must be 0
    # ncomp: vector of ncomp, per matrix
    # scheme: a function "g", refer to the article (thanks Benoit)
    # scale: do you want to scale ? mean is done by default and cannot be changed (so far)
    # bias: scale the data with n or n-1
    # init: one of "svd" or "random", initialisation of the algorithm
    # tol: nobody cares about this
    # verbose: show the progress of the algorithm
    # mode: canonical, classic, invariant, regression
    # max.iter: nobody cares about this
    # study: factor for each matrix of A, must be a vector
    # keepA: keepX of spls for each matrix of A. must be a list. Each entry must be of the same length (max ncomp)
    # keepA.constraint: keepX.constraint, which variables are kept on the first num.comp-1 components. It is a list of characters
    # near.zero.var: do you want to remove variables with very small variance
    
    result=internal_mint.block(A = A,
    indY = indY,
    design = design,
    ncomp = ncomp,
    scheme = scheme,
    scale = scale,
    bias = bias,
    init = init,
    tol = tol,
    verbose = verbose,
    tau = NULL,
    mode = mode,
    max.iter = max.iter,
    study = study,
    keepA = keepA,
    keepA.constraint = keepA.constraint)
       
    if(near.zero.var)
    result$nzv=nzv.A

    class(result) = c("sparse.mint.block")
    return(invisible(result))
    

}




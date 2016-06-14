#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
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


# ========================================================================================================
# mixOmics: perform one of the package's function depending on the input data (list or matrix, vector or categerical data, etc)
# ========================================================================================================


mixOmics = function(X,
Y,
indY, #only use if Y not provided
study, #mint
ncomp,
keepX.constraint, #hybrid
keepY.constraint, #hybrid
keepX, #sparse
keepY, #sparse
design, #block
tau = NULL,# rgcca, number between 0,1 or "optimal"
scheme, #block
mode,
scale,
bias,
init,
tol =  1e-06,
verbose = FALSE,
max.iter = 500,
near.zero.var = FALSE)

{
    if (is.list(X) & !is.data.frame(X))# either rgcca, sgcca,sgcca-DA, mint.block, mint.block-DA
    {
        
        #need to check if Y or indY is a factor, unmap it and then do the checks (no other factors etc)
        if ((missing(indY)& missing(Y)) & is.null(tau))
        stop("Either 'Y', 'indY' or 'tau' is needed")
        
        if (is.null(tau)) # SGCCA/mint
        {
            
            isfactorY = FALSE


            if (!missing(Y))
            {
                if (is.list(Y) & !is.data.frame(X)) stop("Y must be a matrix or a factor")
            
                if (is.factor(Y)) {
                    #Y = as.factor(Y)
                    isfactorY = TRUE
                }
                
            }else if (!missing(indY)) {
                temp = X[[indY]] #not called Y to not be an input of the wrappers
                if (is.factor(temp)) {
                    #temp = as.factor(temp)
                    isfactorY = TRUE
                }
            }else if (missing(indY)) {
                stop("Either 'Y' or 'indY' is needed")
                
            }
            


            
            if (isfactorY)# either block.plsda/block.splsda/mint.block.plsda/mint.block.splsda
            {
                
                if (missing(keepX) & missing(keepX.constraint))
                {
                    if (missing(study)) #block.plsda
                    {
                        if (missing(scale))
                        scale = FALSE
                        
                        message("a block Partial Least Squares - Discriminant Analysis is being performed (block.PLS-DA)")
                        res = block.plsda(X = X, Y = Y, indY = indY, ncomp = ncomp,design = design,scheme = scheme,
                        mode = mode,scale = scale,bias = bias,init = init,tol = tol,verbose = verbose,max.iter = max.iter,near.zero.var = near.zero.var)
                        
                    } else {# mint.block.plsda
                        if (missing(scale))
                        scale = FALSE
                        
                        message("a mint block Partial Least Squares - Discriminant Analysis is being performed (mint.block.PLS-DA)")
                        res = mint.block.plsda(X = X, Y = Y, indY = indY,study = study, ncomp = ncomp,design = design,scheme = scheme,
                        mode = mode,scale = scale,bias = bias,init = init,tol = tol,verbose = verbose,max.iter = max.iter,near.zero.var = near.zero.var)
                    }
                    
                    
                } else {
                    if (missing(study))# block.splsda
                    {
                        if (missing(scale))
                        scale = FALSE
                        
                        message("a block sparse Partial Least Squares - Discriminant Analysis is being performed (block.sPLS-DA)")
                        res = mint.block.splsda(X = X, Y = Y, indY = indY, ncomp = ncomp,keepX = keepX,keepX.constraint = keepX.constraint,
                        design = design,scheme = scheme,mode =  mode,scale = scale,bias = bias,init = init,tol = tol,verbose = verbose,
                        max.iter = max.iter,near.zero.var = near.zero.var)
                        
                        
                    } else {# mint.block.splsda
                        if (missing(scale))
                        scale = FALSE
                       
                        message("a mint block sparse Partial Least Squares - Discriminant Analysis is being performed (mint.block.sPLS-DA)")
                        res = mint.block.splsda(X = X, Y = Y, indY = indY, ncomp = ncomp,study = study,keepX = keepX,keepX.constraint = keepX.constraint,
                        design = design,scheme = scheme,mode =  mode,scale = scale,bias = bias,init = init,tol = tol,verbose = verbose,
                        max.iter = max.iter,near.zero.var = near.zero.var)
                    }
                    
                }
                
            } else { # either block.pls/block.spls/mint.block.pls/mint.block.spls
                
                
                if (missing(keepX) & missing(keepX.constraint))
                {
                    if (missing(study)) #block.pls
                    {
                        if (missing(scale))
                        scale = FALSE
                        
                        message("a block Partial Least Squares is being performed (block.PLS)")
                        res = block.pls(X = X, Y = Y, indY = indY, ncomp = ncomp,design = design,scheme = scheme,
                        mode = mode,scale = scale,bias = bias,init = init,tol = tol,verbose = verbose,max.iter = max.iter,near.zero.var = near.zero.var)
                        
                    } else {# mint.block.pls
                        if (missing(scale))
                        scale = FALSE
                        
                        message("a mint block Partial Least Squares is being performed (mint.block.PLS)")
                        res = mint.block.pls(X = X, Y = Y, indY = indY,study = study, ncomp = ncomp,design = design,scheme = scheme,
                        mode = mode,scale = scale,bias = bias,init = init,tol = tol,verbose = verbose,max.iter = max.iter,near.zero.var = near.zero.var)
                    }
                    
                    
                } else {
                    if (missing(study))# block.spls
                    {
                        if (missing(scale))
                        scale = FALSE
                        
                        message("a block sparse Partial Least Squares is being performed (block.sPLS)")
                        res = block.spls(X = X, Y = Y, indY = indY, ncomp = ncomp,keepX = keepX,keepX.constraint = keepX.constraint,
                        design = design,scheme = scheme,mode =  mode,scale = scale,bias = bias,init = init,tol = tol,verbose = verbose,
                        max.iter = max.iter,near.zero.var = near.zero.var)
                        
                        
                    } else {# mint.block.spls
                        if (missing(scale))
                        scale = FALSE
                        
                        message("a mint block sparse Partial Least Squares is being performed (mint.block.sPLS)")
                        res = mint.block.spls(X = X, Y = Y, indY = indY, ncomp = ncomp,study = study,keepX = keepX,keepX.constraint = keepX.constraint,
                        design = design,scheme = scheme,mode =  mode,scale = scale,bias = bias,init = init,tol = tol,verbose = verbose,
                        max.iter = max.iter,near.zero.var = near.zero.var)
                        
                    }
                    
                }
                
            }
            
        } else { # RGCCA
            
            
            if (!missing(study)) {message("'study' is not used")}
            
            if (missing(keepX) & missing(keepX.constraint)) #RGCCA
            {
                message("A RGCCA analysis is being performed")
                if (missing(scale))
                scale = FALSE
                
                res = wrapper.rgcca(X = X,design = design,tau = tau,ncomp = ncomp,
                max.iter = max.iter,scheme = scheme,scale = scale,init = init,bias = bias,tol = tol,verbose = verbose)
                
            } else { #sparse RGCCA
                if (missing(scale))
                scale = FALSE
                
                message("A sparse RGCCA analysis is being performed")
                res = wrapper.rgcca(X = X,design = design,tau = tau,ncomp = ncomp,keepX = keepX,keepX.constraint = keepX.constraint,
                max.iter = max.iter,scheme = scheme,scale = scale,init = init,bias = bias,tol = tol,verbose = verbose)
                
                
            }
        }
        
        
        
        #end if (is.list(X))
    } else {#either pls,spls, plsda, splsda or mint. pls/spls/plsda/splsda
        if (missing(Y))
        stop("Y is missing")
        if (is.list(Y) & !is.data.frame(X))
        stop("Y must be a matrix or a factor")
        
        if (missing(mode)) mode = "regression"
        #check for unused inputs (scheme, etc etc)
        if (!is.null(tau) | !missing(design) | !missing(init) | !missing(scheme) | !missing(bias) | !missing(verbose))
        {
            if (!is.null(tau))
            message("'tau' is not used")
            if (!missing(design))
            message("'design' is not used")
            if (!missing(init))
            message("'init' is not used")
            if (!missing(scheme))
            message("'scheme' is not used")
            if (!missing(bias))
            message("'bias' is not used")
            if (!missing(verbose))
            message("'verbose' is not used")
            
            stop("unused input parameters")
        }
        
        
        if (is.factor(Y))#either plsda, splsda
        {
            
            #Check.entry.pls.single(X, ncomp, keepX,keepX.constraint) # to have the warnings relative to X and Y, instead of blocks
            if (length(Y)!=nrow(X))
            stop("unequal number of rows in 'X' and 'Y'.")
            
            if (missing(keepX) & missing(keepX.constraint) & missing(keepY) & missing(keepY.constraint))  #plsda, mint.plsda
            {
                if (missing(study))
                {
                    if (missing(scale))
                    scale = TRUE
                    
                    message("a Partial Least Squares - Discriminant Analysis is being performed (PLS-DA)")
                    res = plsda(X = X, Y = Y, ncomp = ncomp, mode = mode,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                } else {# mint
                    if (missing(scale))
                    scale = FALSE
                    
                    message("a mint Partial Least Squares - Discriminant Analysis is being performed (mint.PLS-DA)")
                    res = mint.plsda(X = X, Y = Y, ncomp = ncomp, mode = mode, study = study,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
                
            } else {#splsda, mint.splsda
                if (missing(study))
                {
                    if (missing(scale))
                    scale = TRUE
                    
                    message("a sparse Partial Least Squares - Discriminant Analysis is being performed (sPLS-DA)")
                    res = splsda(X = X, Y = Y, ncomp = ncomp, mode = mode, keepX = keepX,keepX.constraint = keepX.constraint,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                } else {# mint
                    if (missing(scale))
                    scale = FALSE
                    
                    message("a mint sparse Partial Least Squares - Discriminant Analysis is being performed (mint.sPLS-DA)")
                    res = mint.splsda(X = X, Y = Y, ncomp = ncomp, mode = mode, study = study,keepX = keepX,keepX.constraint = keepX.constraint,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
            }
            
        } else { #pls or spls
            
            
            #Check.entry.pls(X, Y, ncomp, keepX, keepY,keepX.constraint,keepY.constraint) # to have the warnings relative to X and Y, instead of blocks
            
            if (missing(keepX) & missing(keepX.constraint) & missing(keepY) & missing(keepY.constraint))  #pls, mint.pls
            {
                if (missing(study))
                {
                    if (missing(scale))
                    scale = TRUE
                    
                    message("a Partial Least Squares is being performed (PLS)")
                    res = pls(X = X, Y = Y, ncomp = ncomp, mode = mode,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                } else { # mint
                    if (missing(scale))
                    scale = FALSE
                    
                    message("a mint Partial Least Squares is being performed (mint.PLS)")
                    res = mint.pls(X = X, Y = Y, ncomp = ncomp, mode = mode, study = study,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
            } else {
                if (missing(study))
                {
                    if (missing(scale))
                    scale = TRUE
                    
                    message("a sparse Partial Least Squares is being performed (sPLS)")
                    res = spls(X = X, Y = Y, ncomp = ncomp, mode = mode, keepX = keepX,keepY = keepY,
                    keepX.constraint = keepX.constraint,keepY.constraint = keepY.constraint,max.iter = max.iter, tol = tol,
                    near.zero.var = near.zero.var,scale = scale)
                } else {
                    if (missing(scale))
                    scale = FALSE
                    
                    message("a mint sparse Partial Least Squares is being performed (mint.sPLS)")
                    res = mint.spls(X = X, Y = Y, ncomp = ncomp, mode = mode, study = study,keepX = keepX,keepY = keepY,
                    keepX.constraint = keepX.constraint,keepY.constraint = keepY.constraint,max.iter = max.iter, tol = tol,
                    near.zero.var = near.zero.var,scale = scale)
                    
                }
            }
            
            
        }
    }
    cl = match.call()
    res$call = cl
    class(res) = c("mixOmics",class(res))
    return(invisible(res))
    
    
}

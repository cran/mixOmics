#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2015
# last modified: 12-04-2016
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

# --------------------------------------
# check and construct keepA and keepA.constraint
# --------------------------------------
get.keepA.and.keepA.constraint=function(X,keepX,keepX.constraint,ncomp)
{
    # X:data
    # keepA
    # keepA.constraint
    
    keepA = list()
    keepA.constraint = list()
    if(missing(keepX.constraint) || length(keepX.constraint) == 0)
    {
        if (missing(keepX) || length(keepX)==0)
        {
            #if both keepX.constraint and keepX are missing, pls-like: keepX=ncol(X)
            for (q in 1:length(X))
            keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)) #keepX
        
            names(keepA)=names(X)
        } else {
            
            if (!is.list(keepX))
            stop("'keepX' must be a list")
            
            if (length(keepX) > length(X))
            stop(paste0("length(keepX) is higher than the number of blocks in X, which is ", length(X), "."))
            
            # error if no names on keepX or not matching the names in X
            if(length(unique(names(keepX)))!=length(keepX) | sum(is.na(match(names(keepX),names(X)))) > 0)
            stop("Each entry of 'keepX' must have a unique name corresponding to a block of 'X'")

            # I want to match keepX to X by names
            ind.match = match(names(X), names(keepX))
            
            for (q in 1:length(X))
            {
                
                if (!is.na(ind.match[q])) # means there is a keepX with the same name as X[q] #(q <= length(keepX))
                {
                    #checking entries of keepX
                    if (is.list(keepX[[ind.match[q]]]))
                    stop(paste0("keepX[[",ind.match[q],"]]' must be a vector"))
                    
                    if (any(keepX[[ind.match[q]]] > ncol(X[[q]])))
                    stop(paste0("each component of 'keepX[[",ind.match[q],"]]' must be lower or equal to ncol(X[[",q,"]])=",ncol(X[[q]]),"."))
                    
                    if (any(keepX[[ind.match[q]]] < 0))
                    stop(paste0("each component of 'keepX[[",ind.match[q],"]]' must be non negative."))
                    
                    if (length(keepX[[ind.match[q]]]) > ncomp[q])
                    stop(paste0("length of 'keepX[[",ind.match[q],"]]' must be lower or equal to ncomp[",q,"]=",ncomp[q], "."))
                    
                    keepA[[q]] = keepX[[ind.match[q]]]
                    if (length(keepA[[q]]) < max(ncomp))
                    keepA[[q]] = c(keepA[[q]], rep(ncol(X[[q]]), max(ncomp) - length(keepA[[q]]))) #complete the keepX already provided
                
                }else{
                    keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))##
                }
            }
            
        }
        
        for(q in 1:length(X))
        keepA.constraint[[q]]=list() #keepX.constraint
        #print(keepA)
        
    }else{
        
        #check entries keepX.constraint
        if (length(keepX.constraint)>length(X))
        stop(paste0("length(keepX.constraint) is higher than the number of blocks in X, which is ", length(X), "."))
        
        # error if no names on keepX.constraint or not matching the names in X
        if (length(unique(names(keepX.constraint))) != length(keepX.constraint) | sum(is.na(match(names(keepX.constraint),names(X)))) > 0)
        stop("Each entry of 'keepX.constraint' must have a unique name corresponding to a block of 'X'")

        # I want to match keepX.constraint to X by names
        ind.match = match(names(X), names(keepX.constraint))
       
        # check that ncomp>=length(keepX.constraint)
        for(q in 1:length(X))
        {
            if (!is.na(ind.match[q])) # means there is a keepX with the same name as X[q] #(q <= length(keepX))
            {
                if (length(keepX.constraint[[ind.match[q]]]) > ncomp[q])
                stop(paste0("you should have length(keepX.constraint[[", ind.match[q], "]]) lower or equal to ncomp[", q, "]=",ncomp[q], "."))
                
            }
        }
        
        if (missing(keepX) || length(keepX)==0)
        {
            #if not missing keepX.constraint but missing keepX, we complete keepX pls-like to have length(keepX.constraint)+length(keepX)=ncomp
            for(q in 1:length(X))
            {
                if (!is.na(ind.match[q]))
                {
                    keepA[[q]] = rep(ncol(X[[q]]), max(ncomp) - length(keepX.constraint[[ind.match[q]]]))
                } else {
                    keepA[[q]] = rep(ncol(X[[q]]), max(ncomp))
                }
                
            }
            names(keepA) = names(X)
        } else { #complete keepA so that length(keepX.constraint[[q]])+length(keepX[[q]])=max(ncomp)
            
            if (!is.list(keepX))
            stop("'keepX' must be a list")
            
            if(length(keepX)>length(X))
            stop(paste0("length(keepX) is higher than the number of blocks in X, which is ",length(X),"."))
            
            ind.match.keepX = match(names(X), names(keepX))
            ind.match.keepX.constraint = match(names(X), names(keepX.constraint))

            for(q in 1:length(X))
            {
                
                #check the entries provided before completed by pls-like
                if (!is.na(ind.match.keepX[q]) & !is.na(ind.match.keepX.constraint[q])) # keepX and keepX.constraint are available for block X[q]
                {
                    if ((length(keepX.constraint[[ind.match.keepX.constraint[q]]]) + length(keepX[[ind.match.keepX[q]]])) > ncomp[q])
                    stop(paste0("length (keepX.constraint[[", ind.match.keepX.constraint[q], "]]) + length(keepX[[", ind.match.keepX[q], "]]) = ", (length(keepX.constraint[[ind.match.keepX.constraint[q]]]) + length(keepX[[ind.match.keepX[q]]])), "; it should be lower or equal to ncomp[", q, "]=", ncomp[q], "."))
                }

                if (!is.na(ind.match.keepX[q])) # # keepX is available for block X[q], maybe keepX.constraint
                {
                    #checking entries of keepX
                    if (is.list(keepX[[ind.match.keepX[q]]]))
                    stop(paste0("keepX[[", ind.match.keepX[q], "]]' must be a vector"))
                    
                    if (any(keepX[[ind.match.keepX[q]]] > ncol(X[[q]])))
                    stop(paste0("each component of 'keepX[[", ind.match.keepX[q], "]]' must be lower or equal to ncol(X[[", q, "]])=", ncol(X[[q]]), "."))
                    
                    if (any(keepX[[ind.match.keepX[q]]]<0))
                    stop(paste0("each component of 'keepX[[", ind.match.keepX[q], "]]' must be non negative."))
                    
                    if (length(keepX[[ind.match.keepX[q]]]) > ncomp[q])
                    stop(paste0("length of 'keepX[[", ind.match.keepX[q], "]]' must be lower or equal to ncomp[", q, "]=", ncomp[q], "."))
                    
                    keepA[[q]] = keepX[[ind.match.keepX[q]]]
                    
                    if (length(keepA[[q]]) < max(ncomp))
                    {
                        #complete the keepX already provided, depending if keepX.constraint is available for block X[q]
                        if (!is.na(ind.match.keepX.constraint[q]))
                        {
                            keepA[[q]] = c(keepA[[q]], rep(ncol(X[[q]]), max(ncomp) - length(keepA[[q]]) - length(keepX.constraint[[ind.match.keepX.constraint[q]]])))
                        } else {
                            keepA[[q]] = c(keepA[[q]], rep(ncol(X[[q]]), max(ncomp) - length(keepA[[q]])))
                        }
                    }
                    
                }else{ # keepX is not available for block X[q], maybe keepX.constraint
                    
                    #complete the keepX already provided, depending if keepX.constraint is available for block X[q]
                    if (!is.na(ind.match.keepX.constraint[q]))
                    {
                        keepA[[q]] = rep(ncol(X[[q]]),max(ncomp)-length(keepX.constraint[[ind.match.keepX.constraint[q]]]))
                    } else {
                        keepA[[q]] = rep(ncol(X[[q]]),max(ncomp))
                    }
                }
            }
            
        }
        
        # match keepX.constraint and the colnames of X in order for keepX.constraint to be a list of character
        # safety if keepX.constraint contains a mixed of character/numeric. It should one or the other, not a mix
        for(q in 1:length(X))
        {
            if (!is.na(ind.match[q])) # means there is a keepX.constraint with the same name as X[q]
            {
                if(length(keepX.constraint[[ind.match[q]]]) > 0)
                {
                    if (is.numeric(unlist(keepX.constraint[[ind.match[q]]])) && any(unlist(keepX.constraint[[ind.match[q]]]) > ncol(X[[q]])))
                    stop(paste0("each entry of 'keepX.constraint[[", ind.match[q], "]]' must be lower or equal than ", ncol(X[[q]]), "."))
                    
                    if (!is.numeric(unlist(keepX.constraint[[ind.match[q]]])))
                    {
                        ind=match(unlist(keepX.constraint[[ind.match[q]]]), colnames(X[[q]]))
                        if(sum(is.na(ind)) > 0)
                        stop("'keepX.constraint' must contain a subset of colnames(X) or the position of the X-variables you wish to keep.")
                    }
                    X.indice = X[[q]][, unlist(keepX.constraint[[ind.match[q]]]), drop=FALSE]
                    keepX.constraint[[ind.match[q]]]=relist(colnames(X.indice), skeleton = keepX.constraint[[ind.match[q]]])
                }
                
                # we need numbers in keepX.constraint from now on
                keepX.constraint[[ind.match[q]]] = sapply(keepX.constraint[[ind.match[q]]], function(x){match(x, colnames(X[[q]]))}, simplify = FALSE)
            
                keepA.constraint[[q]] = keepX.constraint[[ind.match[q]]] #match keepA.constraint with blocks in 'X'
            }else{
                keepA.constraint[[q]] = list()
            }
        }
    }
    names(keepA) = names(keepA.constraint) = names(X)
    return(list(keepA=keepA, keepA.constraint=keepA.constraint))
}


# --------------------------------------
# match.keepX.constraint after removing variables with near.zer.var for instance
# --------------------------------------
match.keepX.constraint=function(names.remove,keepX.constraint)
{
    #matching the new X (after removing some variables) and keepX.constraint
    if (length(names.remove) > 0)
    {
        
        ind.match = lapply(keepX.constraint, function(x){match(x, names.remove)})
        
        if (sum(!is.na(unlist(ind.match))) > 0)
        {
            warning("at least one variable was removed from keepX.constraint because of a null variance. Please check object$keepX.constraint to see which variables are used.")
            #remove from keepX.constraint
            keepX.constraint = lapply(keepX.constraint, function(x)
            {
                temp = match(x, names.remove)
                if (sum(!is.na(temp))>0)
                {
                    x = x[-which(!is.na(temp))]
                } else {
                    x = x
                }
            })
        }
        
        keepX = unlist(lapply(keepX.constraint, length))
        if (any(keepX == 0))
        {
            ind.min = min(which(keepX == 0))
            ncomp = ind.min - 1
            stop(paste("keepX.constraint was reduced from", length(keepX.constraint), " components to", ncomp, " components. Please change keepX.constraint or put near.zero.var=FALSE and restart the function"))
            #construction of the new keepX.constraint, using ncomp components
            keepX.constraint.temp = keepX.constraint
            for (i in (ncomp + 1) : length(keepX.constraint))
            keepX.constraint.temp[[i]] = NULL
            
            keepX.constraint = keepX.constraint.temp
        }
    }
    out=keepX.constraint
}


# --------------------------------------
# Check.entry.single: check input parameter for a single matrix that comes from a list of matrix
# --------------------------------------
# X: a single numeric matrix
# ncomp: the number of components to include in the model
# q: position of X in a list of matrix, as used in horizontal analysis (e.g. block.pls)

Check.entry.single = function(X,  ncomp, q)
{

    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop(paste0("'X[[", q, "]]' must be a numeric matrix."))
    
    X = as.matrix(X)
    
    if (!is.numeric(X))
    stop(paste0("'X[[", q, "]]'  must be a numeric matrix."))
    
    N = nrow(X)
    P = ncol(X)
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop(paste0("invalid number of variates 'ncomp' for matrix 'X[[", q, "]]'."))
    
    ncomp = round(ncomp)
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]] = X.names
    }
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X)  = ind.names
    }
    
    if (length(unique(rownames(X))) != nrow(X))
        stop("samples should have a unique identifier/rowname")
    if (length(unique(X.names)) != P)
        stop("Unique indentifier is needed for the columns of X")
    
    return(list(X=X, ncomp=ncomp, X.names=X.names, ind.names=ind.names))
}



# --------------------------------------
# Check.entry.pls: check the entry associated with (s)PLS. Usde in 'internal_wrapper.mint.R'
# --------------------------------------
# X: numeric matrix of predictors
# Y: numeric vector or matrix of responses
# ncomp: the number of components to include in the model. Default to 2.
# keepX: number of \eqn{X} variables kept in the model on the last components (once all keepX.constraint[[i]] are used).
# keepY: number of \eqn{Y} variables kept in the model on the last components.
# keepX.constraint: A list containing which variables of X are to be kept on each of the first PLS-components.
# keepY.constraint: A list containing which variables of Y are to be kept on each of the first PLS-components
# mode: deflation mode
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# max.iter: integer, the maximum number of iterations.
# tol: Convergence stopping value.
# logratio: whether a log transformation will be performed, and which one,
# DA: whether the input is used for a Discrimant Analysis
# multilevel: whether a multilevel analysis will be performed (repeated measurements)


Check.entry.pls = function(X, Y, ncomp, keepX, keepY, keepX.constraint, keepY.constraint, mode, scale,
     near.zero.var, max.iter, tol, logratio, DA, multilevel)
{

    if (missing(mode))
    mode = "regression"
    
    if (length(mode)>1)
    mode = mode[1]
    
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
    stop("Choose one of the four following modes: canonical, invariant, classic or regression")
    
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")
    
    X = as.matrix(X)
    
    if (!(logratio %in% c("none", "CLR")))
    stop("Choose one of the two following logratio transformation: none or CLR")
    
    if(!is.null(multilevel))
    {
        #multilevel analysis: withinVariation and then pls-like
        # if it's DA analysis, Y and 'multilevel' are combined
        if(DA)
        {
            Y = multilevel
        }else{
            if ((nrow(X) != nrow(multilevel)))
            stop("unequal number of rows in 'X' and 'multilevel'.")
            
            Y = as.matrix(Y)
            if (!is.numeric(X) || !is.numeric(Y))
            stop("'X' and/or 'Y' must be a numeric matrix.")
        }
    }else{
        Y = as.matrix(Y)
        if (!is.numeric(X) || !is.numeric(Y))
        stop("'X' and/or 'Y' must be a numeric matrix.")
    }
    N = nrow(X)
    Q = ncol(Y)
    P= ncol(X)
    
    if ((N != nrow(Y)))
    stop("Unequal number of rows in 'X' and 'Y'.")
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || length(ncomp)>1)
    stop("invalid number of variates, 'ncomp'.")
    
    ncomp = round(ncomp)
    if(ncomp > P)
    {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", P, ".")
        ncomp = P
    }
    
    if (!is.numeric(tol) | tol<=0)
    stop("tol must be non negative")
    
    if (!is.numeric(max.iter) | max.iter<=0)
    stop("max.iter must be non negative")
    
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]] = X.names
    }
    

    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = dimnames(Y)[[1]]
        rownames(X) = ind.names
    }
    
    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X) = rownames(Y) = ind.names
    }
    
    rownames(X) = rownames(Y) = ind.names
    
    
    #if (dim(Y)[2] == 1) Y.names = "Y"
    Y.names = dimnames(Y)[[2]]
    if (is.null(Y.names))
    {
        if (dim(Y)[2] == 1)
        {
            Y.names = "Y"
        } else {
            Y.names = paste("Y", 1:Q, sep = "")
        }
        
        dimnames(Y)[[2]]=Y.names
    }
    
    if (length(unique(X.names)) != P)
    stop("Unique indentifier is needed for the columns of X")
    
    if (length(unique(Y.names)) != Q)
    stop("Unique indentifier is needed for the columns of Y")
    
    
    # check on keepX and keepX.constraint
    if (missing(keepX.constraint))
    {
        if (missing(keepX))
        {
            keepX = rep(P, ncomp)
        } else {
            if (length(keepX)<ncomp)
            keepX = c(keepX, rep(P, ncomp - length(keepX))) #complete (with ncomp) the keepX already provided
        }
        keepX.constraint=list()
    } else {
        if (length(keepX.constraint)>ncomp)
        stop(paste0("You should have length(keepX.constraint) lower or equal to 'ncomp' = ", ncomp, "."))
        
        if (missing(keepX))
        {
            keepX = rep(P, ncomp - length(keepX.constraint))
        } else {
            
            if ((length(keepX.constraint) + length(keepX)) < ncomp)
            keepX = c(keepX, rep(P, ncomp - length(keepX) - length(keepX.constraint)))
            
            if ((length(keepX.constraint) + length(keepX)) > ncomp)
            stop(paste0("length (keepX.constraint) + length(keepX) should be lower than 'ncomp' = ", ncomp, "."))
            
        }
    }
    
    # check on keepY and keepY.constraint
    if (missing(keepY.constraint))
    {
        if (missing(keepY))
        {
            keepY = rep(Q, ncomp)
        } else {
            if (length(keepY) < ncomp)
            keepY = c(keepY, rep(Q, ncomp - length(keepY))) #complete the keepY already provided
        }
        keepY.constraint = list()
    } else {
        if (length(keepY.constraint)>ncomp)
        stop(paste0("you should have length(keepY.constraint) lower or equal to 'ncomp' = ", ncomp, "."))
        
        if (missing(keepY))
        {
            keepY = rep(Q, ncomp - length(keepY.constraint))
        } else {
            
            if ((length(keepY.constraint) + length(keepY)) < ncomp)
            keepY = c(keepY, rep(Q, ncomp - length(keepY) - length(keepY.constraint)))
            
            if ((length(keepY.constraint) + length(keepY)) > ncomp)
            stop(paste0("length (keepY.constraint) + length(keepY) should be lower than 'ncomp' = ", ncomp,"."))
        }
    }
    
    
    
    if (any(keepX<0))
    stop("each component of 'keepX' must be non negative ")
    if (any(keepY<0))
    stop("each component of 'keepY' must be non negative ")
    
    if (any(keepX > ncol(X)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    if (any(keepY > ncol(Y)))
    stop("each component of 'keepY' must be lower or equal than ", Q, ".")
    
    if (is.numeric(unlist(keepX.constraint)) && any(unlist(keepX.constraint) > ncol(X)))
    stop("each entry of 'keepX.constraint' must be lower or equal than ", P, ".")
    if ( is.numeric(unlist(keepY.constraint)) && any(unlist(keepY.constraint) > ncol(Y)))
    stop("each entry of 'keepY.constraint' must be lower or equal than ", Q, ".")
    
    if (!is.logical(scale))
    stop("'scale' must be either TRUE or FALSE")
    
    if (!is.logical(near.zero.var))
    stop("'near.zero.var' must be either TRUE or FALSE")
    
    # match keepX.constraint and the colnames of X in order for keepX.constraint to be a list of character
    # safety if keepX.constraint contains a mixed of character/numeric. It should one or the other, not a mix
    if (length(keepX.constraint) > 0)
    {
        if (!is.numeric(unlist(keepX.constraint)))
        {
            ind = match(unlist(keepX.constraint), colnames(X))
            if (sum(is.na(ind)) > 0)
            stop("'keepX.constraint' must contain a subset of colnames(X) or the position of the X-variables you wish to keep.")
        }
        X.indice = X[, unlist(keepX.constraint), drop=FALSE]
        keepX.constraint = relist(colnames(X.indice), skeleton=keepX.constraint)
    }
    
    # same for keepY.constraint
    if (length(keepY.constraint) > 0)
    {
        if (!is.numeric(unlist(keepY.constraint)))
        {
            ind = match(unlist(keepY.constraint),colnames(Y))
            if (sum(is.na(ind)) > 0)
            stop("'keepY.constraint' must contain a subset of colnames(Y) or the position of the Y-variables you wish to keep.")
        }
        Y.indice = Y[, unlist(keepY.constraint), drop=FALSE]
        keepY.constraint = relist(colnames(Y.indice), skeleton=keepY.constraint)
    }
    

    # at this stage keepA.constraint needs to be character, to easily remove variables with near zero variance
    ### near.zero.var, remove the variables with very small variances
    if (near.zero.var == TRUE)
    {
        nzv.A = nearZeroVar(X)
        
        if (length(nzv.A$Position) > 0)
        {
            names.remove.X = colnames(X)[nzv.A$Position]
            X = X[, -nzv.A$Position, drop=FALSE]
            warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
            if (ncol(X) == 0)
            stop("No more variables in X")
            
            # at this stage, keepA.constraint needs to be numbers
            if (length(keepX.constraint) > 0)
            {
                #remove the variables from keepA.constraint if removed by near.zero.var
                keepX.constraint = match.keepX.constraint(names.remove.X, keepX.constraint)
            }
            #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
            if (any(keepX > ncol(X)))
            {
                ind = which(keepX > ncol(X))
                keepX[ind] = ncol(X)
            }
        }
        
    }else{nzv.A=NULL}
    
    # we need numbers in keepX.constraint from now on
    keepX.constraint = lapply(keepX.constraint,function(x){match(x, colnames(X))})
    keepY.constraint = lapply(keepY.constraint,function(x){match(x, colnames(Y))})
    
    return(list(X=X, Y=Y, ncomp=ncomp, X.names=X.names, Y.names=Y.names, ind.names=ind.names, mode=mode, keepX.constraint=keepX.constraint,
    keepY.constraint=keepY.constraint, keepX=keepX, keepY=keepY, nzv.A=nzv.A))
}


# --------------------------------------
# Check.entry.wrapper.mint.block: used in 'internal_wrapper.mint.block.R'
# --------------------------------------
# X: a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets. \code{NA}s are not allowed.
# Y: a factor or a class vector for the discrete outcome.
# indY: to supply if Y is missing, indicate the position of the outcome in the list X.
# ncomp: numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.
# keepX: A vector of same length as X.  Each entry keepX[i] is the number of X[[i]]-variables kept in the model on the last components (once all keepX.constraint[[i]] are used).
# keepX.constraint: A list of same length as X. Each entry keepX.constraint[[i]] is a list containing which variables of X[[i]] are to be kept on each of the first PLS-components
# keepY: Only used if Y is provided. Each entry keepY[i] is the number of Y-variables kept in the model on the last components.
# keepY.constraint: Only used if Y is provided, otherwise extracted from keepX.constraint. A list containing which variables of Y are to be kept on each of the first PLS-components
# study: grouping factor indicating which samples are from the same study
# design: the input design.
# init: intialisation of the algorithm, one of "svd" or "svd.single". Default to "svd"
# scheme: the input scheme, one of "horst", "factorial" or ""centroid". Default to "centroid"
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# bias: boleean. A logical value for biaised or unbiaised estimator of the var/cov (defaults to FALSE).
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# mode: input mode, one of "canonical", "classic", "invariant" or "regression". Default to "regression"
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# verbose: if set to \code{TRUE}, reports progress on computing.

Check.entry.wrapper.mint.block = function(X,
Y,
indY, #only use if Y not provided
ncomp,
keepX,
keepX.constraint,
keepY,
keepY.constraint,
study, #mint
design, #block
init,
scheme, #block
scale,
bias,
near.zero.var,
mode,
tol,
max.iter,
verbose)
{
    #need to give the default values of mint.block.spls to mixOmics
    
    if (!is.list(X))
    stop("X must be a list")
    
    if (length(unique(unlist(lapply(X, nrow)))) != 1)
    stop("Unequal number of rows among the blocks of 'X'")
    
    if ((missing(indY) & missing(Y)))
    stop("Either 'Y' or 'indY' is needed")
    
    if (missing(ncomp))
    ncomp = 1
    
    #check length(ncomp)=1
    if (length(ncomp) != 1)
    stop("'ncomp' must be a single value")
    
    # transform ncomp to length(X)
    ncomp = rep(ncomp, length(X))

    # check names on X are unique
    if (length(unique(names(X))) != length(X))
    stop("Each block of 'X' must have a unique name.")
    #names(X)=paste0("block", 1:length(X)) #add names to the blocks if no names or not unique name for each block


    #check dimnames and ncomp per block of A
    for (q in 1:length(X))
    {
        check = Check.entry.single(X[[q]], ncomp[q], q=q)
        X[[q]] = check$X
        ncomp[q] = check$ncomp
    }
    
    #check ncomp[q] < ncol(X[[q]])
    for (q in 1:length(X))
    {
        ncomp[q] = round(ncomp[q])
        if(ncomp[q] > ncol(X[[q]]))
        {
            warning(paste0("Reset maximum number of variates 'ncomp[", q, "]' to ncol(X[[", q, "]])= ", ncol(X[[q]]), "."))
            ncomp[q] = ncol(X[[q]])
        }
    }
    
    # construction of keepA and keepA.constraint for X, we include Y later on if needed (i.e. if Y is provided, otherwise Y is in X[[indY]])
    check = get.keepA.and.keepA.constraint(X=X, keepX=keepX, keepX.constraint=keepX.constraint, ncomp=ncomp)
    keepA = check$keepA
    keepA.constraint = check$keepA.constraint
    
    
    if (!missing(Y))# Y is not missing, so we don't care about indY
    {        
        if (!missing(indY))
        warning("'Y' and 'indY' are provided, 'Y' is used.")
        
        if (is.list(Y))
        stop("Y must be a matrix")
        
        Y = as.matrix(Y)
        if (!is.numeric(Y) )
        stop("'Y' must be a numeric matrix.")
        
        #check dimnames and ncomp per block of A
        check = Check.entry.single(Y, max(ncomp), q=1)
        Y = check$X
        
        if (nrow(Y)!=nrow(X[[1]]))
        stop("Different number of samples between the blocks and Y")
        
        # if not missing, we transform keepY and keepY.constraint in list to input them in get.keepA.and.keepA.constraint
        if (!missing(keepY))
        keepY = list(Y = keepY)
        
        if (!missing(keepY.constraint))
        keepY.constraint = list(Y = keepY.constraint)
        
        check.temp.Y = get.keepA.and.keepA.constraint(X = list(Y = Y), keepX = keepY, keepX.constraint = keepY.constraint, ncomp = ncomp)
        keepY.temp = check.temp.Y$keepA
        keepY.constraint.temp = check.temp.Y$keepA.constraint

        keepA[[length(X)+1]] = keepY.temp[[1]] #add keepY in keepA
        keepA.constraint[[length(X)+1]] = keepY.constraint.temp[[1]] #add keepY.constraint in keepA


        # check design matrix before adding Y in
        A = X
        ### Start check design matrix
        if (missing(design))
        {
            design = 1 - diag(length(A)+1)
            rownames(design) = colnames(design) = c(names(A), "Y")
        } else if (ncol(design) != nrow(design) || ncol(design) < length(X) || ncol(design) > (length(X) + 1) || any(!design %in% c(0,1))) {
            stop(paste0("'design' must be a square matrix with ", length(X), "columns."))
        } else if (ncol(design) == length(X)) {
            message("Design matrix has changed to include Y; each block will be linked to Y.")
            design = rbind(cbind(design, 1), 1)
            diag(design) = 0
            rownames(design) = colnames(design) = c(names(A), "Y")
        }
        rownames(design) = colnames(design) = c(names(A), "Y")
        ### End check design matrix

        # build the list A by adding Y, and creating indY
        A[[length(A)+1]] = Y
        names(A)[length(A)] = names(keepA)[length(A)] = names(keepA.constraint)[length(A)] = "Y"
        indY = length(A)
        
        if (mode == "canonical")
        ncomp = c(ncomp, min(ncomp, ncol(Y) - 1)) #adjust ncomp for Y
        
        if (mode == "regression")
        ncomp = c(ncomp, max(ncomp)) #adjust ncomp for Y
        
    } else {        #missing(Y) but indY not missing
        
        if (!missing(keepY))
        warning("indY is provided: 'keepY' is ignored and only 'keepX' is used.")
        if (!missing(keepY.constraint))
        warning("indY is provided: 'keepY.constraint' is ignored and only 'keepX.constraint' is used.")

        A = X #list provided as input
        
        ### Start check design matrix
        if (missing(design))
        {
            design = 1 - diag(length(A))
        } else if (ncol(design) != nrow(design) || ncol(design) < length(A) || ncol(design) > (length(A) + 1) || any( !design %in% c(0,1)))
        {
            stop(paste0("'design' must be a square matrix with ", length(A), "columns."))
        }
        rownames(design) = colnames(design) = names(A)
        ### End check design matrix
        
        # check indY
        if (!is.numeric(indY) | indY > length(A) | length(indY) > 1)
        stop(paste0("'indY' must be a numeric value lower or equal to ", length(A), ", the number of blocks in A."))

    }
    
    # --------------------------------------------------------------------------------
    # at this stage, we have A, indY, keepA, keepA.constraint, ncomp verified
    # --------------------------------------------------------------------------------
    
    # force all rownames to be the same
    ind.names = rownames(A[[1]])
    lapply(A, function(x){rownames(x) = ind.names})
    
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
    stop("Choose one of the four following modes: canonical, invariant, classic or regression")
    
    #set the default study factor
    if (missing(study))
    {
        study = factor(rep(1, nrow(A[[1]])))
    } else {
        study = as.factor(study)
    }
    if (length(study) != nrow(A[[1]]))
    stop(paste0("'study' must be a factor of length ", nrow(A[[1]]), "."))
    
    if (any(table(study) <= 1))
    stop("At least one study has only one sample, please consider removing before calling the function again")
    if (any(table(study) < 5))
    warning("At least one study has less than 5 samples, mean centering might not do as expected")
    
    if (missing(init))
    init = "svd"
    
    if (!init%in%c("svd","svd.single"))
    stop("init should be one of 'svd' or 'svd.single'")
    
    if (!abs(indY - round(indY) < 1e-25))
    stop ("indY must be an integer")
    if (indY > length(A))
    stop ("indY must point to a block of A")
    
    # =====================================================
    # with or without tau (RGGCA or mint.block.spls algo)
    # =====================================================
    x = unlist(lapply(A,nrow))
    if (!isTRUE(all.equal(max(x), min(x))))
    stop("The samplesize must be the same for all blocks")
    
    #check scheme
    if (!(scheme %in% c("horst", "factorial", "centroid")))
    {
        stop("Choose one of the three following schemes: horst, centroid or factorial")
    } else {
        if (verbose)
        cat("Computation of the SGCCA block components based on the", scheme, "scheme \n")
    }
    
    if (!is.numeric(tol) | tol <= 0)
    stop("tol must be non negative")
    if (!is.numeric(max.iter) | max.iter <= 0)
    stop("max.iter must be non negative")
    
    if (!is.logical(verbose))
    stop("verbose must be either TRUE or FALSE")
    if (!is.logical(scale))
    stop("scale must be either TRUE or FALSE")
    if (!is.logical(bias))
    stop("bias must be either TRUE or FALSE")
    if (!is.logical(near.zero.var))
    stop("near.zero.var must be either TRUE or FALSE")
    
    
    # at this stage keepA.constraint need to be character, to remove easily variables with near zero variance
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(A, nearZeroVar)
        for(q in 1:length(A))
        {
            if (length(nzv.A[[q]]$Position) > 0)
            {
                names.remove.X = colnames(A[[q]])[nzv.A[[q]]$Position]
                A[[q]] = A[[q]][, -nzv.A[[q]]$Position, drop=FALSE]
                if (verbose)
                warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
                if (ncol(A[[q]]) == 0)
                stop(paste0("No more variables in",A[[q]]))
                
                # at this stage, keepA.constraint need to be numbers
                if (length(keepA.constraint[[q]]) > 0)
                {
                    #remove the variables from keepA.constraint if removed by near.zero.var
                    keepA.constraint[[q]] = match.keepX.constraint(names.remove.X, keepA.constraint[[q]])
                    # replace character by numbers
                    keepA.constraint[[q]] = lapply(keepA.constraint[[q]], function(x){match(x, colnames(A[[q]]))})
                }
                #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
                if (any(keepA[[q]] > ncol(A[[q]])))
                {
                    ind = which(keepA[[q]] > ncol(A[[q]]))
                    keepA[[q]][ind] = ncol(A[[q]])
                }
            }
            
        }
    } else {
        nzv.A=NULL
    }
    
    return(list(A=A, ncomp=ncomp, study=study, keepA=keepA, keepA.constraint=keepA.constraint,
    indY=indY, design=design, init=init, nzv.A=nzv.A))
}


# --------------------------------------
# Check.entry.sgcca: used in 'wrapper.sgcca.R'
# --------------------------------------
# X: a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets. \code{NA}s are not allowed.
# design: the input design.
# ncomp: numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid". Default to "centroid"
# mode: input mode, one of "canonical", "classic", "invariant" or "regression". Default to "regression"
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# init: intialisation of the algorithm, one of "svd" or "svd.single". Default to "svd"
# bias: boleean. A logical value for biaised or unbiaised estimator of the var/cov (defaults to FALSE).
# tol: Convergence stopping value.
# verbose: if set to \code{TRUE}, reports progress on computing.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# keepX: A vector of same length as X.  Each entry keepX[i] is the number of X[[i]]-variables kept in the model on the last components (once all keepX.constraint[[i]] are used).
# keepX.constraint: A list of same length as X. Each entry keepX.constraint[[i]] is a list containing which variables of X[[i]] are to be kept on each of the first PLS-components

Check.entry.sgcca = function(X,
design,
ncomp,
scheme,
mode,
scale,
init,
bias,
tol,
verbose,
max.iter,
near.zero.var,
keepX,
keepX.constraint)
{
    #need to give the default values of mint.block.spls to mixOmics
    
    if (!is.list(X))
    stop("X must be a list of at list two matrices")
    
    if (length(X)<2)
    stop("X must be a list of at list two matrices")
    
    if (length(unique(unlist(lapply(X, nrow)))) != 1)
    stop("Unequal number of rows among the blocks of 'X'")
    
    # check names on X are unique
    if (length(unique(names(X))) != length(X))
    stop("Each block of 'X' must have a unique name.")
    #names(X)=paste0("block", 1:length(X)) #add names to the blocks if no names or not unique name for each block

    if (missing(ncomp))
    ncomp = 1

    #check length(ncomp)=1
    if (length(ncomp) != 1)
    stop("'ncomp' must be a single value")

    # transform ncomp to length(X)
    ncomp = rep(ncomp, length(X))

    #check dimnames and ncomp per block of A
    for (q in 1:length(X))
    {
        check = Check.entry.single(X[[q]], ncomp[q], q=q)
        X[[q]] = check$X
        ncomp[q] = check$ncomp
    }
    
    #check ncomp[q]<ncol(X[[q]])
    for (q in 1:length(X))
    {
        ncomp[q] = round(ncomp[q])
        if(ncomp[q] > ncol(X[[q]]))
        {
            warning(paste0("Reset maximum number of variates 'ncomp[", q, "]' to ncol(X[[", q, "]])= ", ncol(X[[q]]), "."))
            ncomp[q] = ncol(X[[q]])
        }
    }
    
    A = X#input
    
    if (missing(init))
    init="svd"
    
    if (!init%in%c("svd", "svd.single"))
    stop("init should be one of 'svd' or 'svd.single'")
    
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
    stop("Choose one of the four following modes: canonical, invariant, classic or regression")
    
    
    #check scheme
    if (missing(scheme))
    scheme= "centroid"
    
    if (!(scheme %in% c("horst", "factorial","centroid")))
    {
        stop("Choose one of the three following schemes: horst, centroid or factorial")
    } else {
        if (verbose)
        cat("Computation of the SGCCA block components based on the", scheme, "scheme \n")
    }
    
    
    if(missing(design))
    design = 1 - diag(length(A))
    
    # check design matrix
    if (nrow(design) != ncol(design))
    stop(paste0("'design' must be a square matrix."))
    
    if (nrow(design) != length(A))
    stop(paste0("'design' must be a square matrix with", length(A), "columns."))
    
    
    if (missing(bias))
    bias = FALSE
    if (missing(verbose))
    verbose = FALSE
    
    if (tol <= 0)
    stop("tol must be non negative")
    
    if (max.iter <= 0)
    stop("max.iter must be non negative")
    
    if (!is.logical(verbose))
    stop("verbose must be either TRUE or FALSE")
    if (!is.logical(scale))
    stop("scale must be either TRUE or FALSE")
    if (!is.logical(bias))
    stop("bias must be either TRUE or FALSE")
    if (!is.logical(near.zero.var))
    stop("near.zero.var must be either TRUE or FALSE")
    
    # construction of keepA and keepA.constraint
    check = get.keepA.and.keepA.constraint(X=A, keepX=keepX, keepX.constraint=keepX.constraint, ncomp=ncomp)
    keepA = check$keepA
    keepA.constraint = check$keepA.constraint
    
    # at this stage keepA.constraint need to be character, to easily remove variables with near zero variance
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(A,nearZeroVar)
        for (q in 1:length(A))
        {
            if (length(nzv.A[[q]]$Position) > 0)
            {
                names.remove.X = colnames(A[[q]])[nzv.A[[q]]$Position]
                A[[q]] = A[[q]][, -nzv.A[[q]]$Position, drop=FALSE]
                if (verbose)
                warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
                if (ncol(A[[q]]) == 0)
                stop(paste0("No more variables in", A[[q]]))
                
                # at this stage, keepA.constraint need to be numbers
                if (length(keepA.constraint[[q]])>0)
                {
                    #remove the variables from keepA.constraint if removed by near.zero.var
                    keepA.constraint[[q]] = match.keepX.constraint(names.remove.X, keepA.constraint[[q]])
                    # replace character by numbers
                    keepA.constraint[[q]] = lapply(keepA.constraint[[q]], function(x){match(x, colnames(A[[q]]))})
                }
                #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
                if (any(keepA[[q]]>ncol(A[[q]])))
                {
                    ind = which(keepA[[q]] > ncol(A[[q]]))
                    keepA[[q]][ind] = ncol(A[[q]])
                }
            }
            
        }
    } else {
        nzv.A=NULL
    }
    
    
    return(list(A=A, ncomp=ncomp, design=design, init=init, scheme=scheme, verbose=verbose, bias=bias, nzv.A=nzv.A,
    keepA=keepA, keepA.constraint=keepA.constraint))
    
}



# --------------------------------------
# Check.entry.rgcca: used in 'wrapper.rgcca.R'
# --------------------------------------
# X: a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets. \code{NA}s are not allowed.
# design: the input design.
# tau:
# ncomp: numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid". Default to "centroid"
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# init: intialisation of the algorithm, one of "svd" or "svd.single". Default to "svd"
# bias: boleean. A logical value for biaised or unbiaised estimator of the var/cov (defaults to FALSE).
# tol: Convergence stopping value.
# verbose: if set to \code{TRUE}, reports progress on computing.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# keepX: A vector of same length as X.  Each entry keepX[i] is the number of X[[i]]-variables kept in the model on the last components (once all keepX.constraint[[i]] are used).
# keepX.constraint: A list of same length as X. Each entry keepX.constraint[[i]] is a list containing which variables of X[[i]] are to be kept on each of the first PLS-components
Check.entry.rgcca = function(X,
design,
tau,
ncomp,
scheme,
scale,
init,
bias,
tol,
verbose,
max.iter,
near.zero.var,
keepX,
keepX.constraint)
{
    #need to give the default values of mint.block.spls to mixOmics
    
    if (!is.list(X))
    stop("X must be a list of at list two matrices")
    
    if (length(X)<2)
    stop("X must be a list of at list two matrices")
    
    if (length(unique(unlist(lapply(X, nrow)))) != 1)
    stop("Unequal number of rows among the blocks of 'X'")

    # check names on X are unique
    if (length(unique(names(X))) != length(X))
    stop("Each block of 'X' must have a unique name.")
    #names(X)=paste0("block", 1:length(X)) #add names to the blocks if no names or not unique name for each block

    if (is.null(tau))
    stop("'tau' is needed")
    
    if (missing(ncomp))
    ncomp = 1
    
    #check length(ncomp)=1
    if (length(ncomp) != 1)
    stop("'ncomp' must be a single value")
    
    # transform ncomp to length(X)
    ncomp = rep(ncomp, length(X))
    
    #check dimnames and ncomp per block of A
    for (q in 1:length(X))
    {
        check = Check.entry.single(X[[q]], ncomp[q], q=q)
        X[[q]] = check$X
        ncomp[q] = check$ncomp
    }
    
    
    #check ncomp[q]<ncol(X[[q]])
    for(q in 1:length(X))
    {
        ncomp[q] = round(ncomp[q])
        if(ncomp[q] > ncol(X[[q]]))
        {
            warning(paste0("Reset maximum number of variates 'ncomp[", q, "]' to ncol(X[[", q, "]])= ", ncol(X[[q]]), "."))
            ncomp[q] = ncol(X[[q]])
        }
    }
    
    A = X#input
    
    if (is.numeric(tau))
    {
        if (any(tau < 0) | any(tau > 1))
        stop("'tau' contains either values between 0 and 1, or 'optimal'.")
        
        if (is.vector(tau))
        {
            if (length(tau) != length(A))
            stop(paste0("'tau' must be of length ", length(A), "."))
            
            tau = matrix(rep(tau, max(ncomp)), nrow = max(ncomp), ncol = length(tau), byrow = TRUE)
        }
    } else {
        if (tau != "optimal")
        stop("'tau' contains either values between 0 and 1, or 'optimal'.")
    }
    
    if (missing(init))
    init = "svd.single"
    
    if (init != "svd.single")
    stop("init should be 'svd.single'.")
    
    # check scheme
    if(missing(scheme)) scheme = "centroid"
    if (!(scheme %in% c("horst", "factorial", "centroid")))
    {
        stop("Choose one of the three following schemes: horst, centroid or factorial")
    } else {
        if (verbose)
        cat("Computation of the SGCCA block components based on the", scheme, "scheme \n")
    }
    
    
    if (missing(design))
    design = 1 - diag(length(A))
    
    #check design matrix
    if (nrow(design) != ncol(design))
    stop(paste0("'design' must be a square matrix."))
    if (nrow(design) != length(A))
    stop(paste0("'design' must be a square matrix with", length(A), "columns."))
    
    
    if (missing(bias))
    bias = FALSE
    if (missing(verbose))
    verbose = FALSE
    if (missing(near.zero.var))
    near.zero.var = FALSE
    
    if (tol <= 0)
    stop("tol must be non negative")
    
    if (max.iter <= 0)
    stop("max.iter must be non negative")
    
    if (!is.logical(verbose))
    stop("verbose must be either TRUE or FALSE")
    if (!is.logical(scale))
    stop("scale must be either TRUE or FALSE")
    if (!is.logical(bias))
    stop("bias must be either TRUE or FALSE")
    if (!is.logical(near.zero.var))
    stop("near.zero.var must be either TRUE or FALSE")
    
    # construction of keepA and keepA.constraint
    check = get.keepA.and.keepA.constraint(X=A, keepX=keepX, keepX.constraint=keepX.constraint, ncomp=ncomp)
    keepA = check$keepA
    keepA.constraint = check$keepA.constraint
    
    
    
    # at this stage keepA.constraint need to be character, to remove easily variables with near zero variance
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(A,nearZeroVar)
        for (q in 1:length(A))
        {
            if (length(nzv.A[[q]]$Position) > 0)
            {
                names.remove.X = colnames(A[[q]])[nzv.A[[q]]$Position]
                A[[q]] = A[[q]][, -nzv.A[[q]]$Position, drop=FALSE]
                if (verbose)
                warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
                if (ncol(A[[q]]) == 0)
                stop(paste0("No more variables in", A[[q]]))
                
                # at this stage, keepA.constraint need to be numbers
                if (length(keepA.constraint[[q]])>0)
                {
                    #remove the variables from keepA.constraint if removed by near.zero.var
                    keepA.constraint[[q]] = match.keepX.constraint(names.remove.X, keepA.constraint[[q]])
                    # replace character by numbers
                    keepA.constraint[[q]] = lapply(keepA.constraint[[q]], function(x){match(x, colnames(A[[q]]))})
                }
                #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
                if (any(keepA[[q]] > ncol(A[[q]])))
                {
                    ind = which(keepA[[q]] > ncol(A[[q]]))
                    keepA[[q]][ind] = ncol(A[[q]])
                }
            }
            
        }
    } else {
        nzv.A = NULL
    }
    
    
    return(list(A=A, ncomp=ncomp, design=design, init=init, scheme=scheme, verbose=verbose, bias=bias, nzv.A=nzv.A,
    keepA=keepA, keepA.constraint=keepA.constraint))
}





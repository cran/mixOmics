#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
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


# ========================================================================================================
# Internal helpers functions to run "mixOmics" and "internal_mint.block" functions
# ========================================================================================================

# Some of these functions have been borrowed from the RGCCA package, as indicated below

# --------------------------------------
# study_split: used in 'internal_mint.block.R' and 'predict.mint.block.pls.R'
# --------------------------------------
get.weights = function(variates, indY)
{
    ncomp = min(sapply(variates, ncol))
    x.xList <- list()
    compt = 1
    for(comp in 1:ncomp)
    {
        for(i in 1:length(variates)){
            corDat <- rep(0, length(variates))
            names(corDat) <- paste("cor", names(variates)[i], names(variates), sep = "_")
            for(j in 1:length(variates)){
                corDat[j] <- as.numeric(cor(variates[[i]][,comp], variates[[j]][,comp]))
            }
            x.xList[[compt]] <- corDat
            compt = compt +1
        }
    }
    corMat.diablo <- do.call(rbind, x.xList)
    rownames(corMat.diablo) <- paste(names(variates),".comp",rep(1:ncomp,each=length(variates)),sep="")
    colnames(corMat.diablo) <- names(variates)
    
    temp = matrix(corMat.diablo[,indY],ncol=ncomp)
    correlation = apply(temp, 1, function(x){mean(abs(x))})[1:length(variates)]
    names(correlation) = names(variates)
    
    correlation = correlation[-indY]
    return(correlation)
    
}


# --------------------------------------
# study_split: used in 'internal_mint.block.R' and 'predict.mint.block.pls.R'
# --------------------------------------
study_split = function(data, study)
{
    data = as.matrix(data)
    M = length(levels(study))
    P = ncol(data)
    
    #---------------------- split data
    data.list.study = split(data,study)
    if (!is.null(rownames(data)))
    study.name = split(rownames(data),study)
    
    for(m in 1:M)
    {
        data.list.study[[m]] = matrix(data.list.study[[m]], ncol=P)
        
        if (!is.null(colnames(data)))
        colnames(data.list.study[[m]]) = colnames(data)
        
        if (!is.null(rownames(data)))
        rownames(data.list.study[[m]]) = study.name[[m]]
    }
    result = data.list.study
    return(invisible(result))
}


# --------------------------------------
# soft_thresholding: used in sparsity (below)
# --------------------------------------
# x: vector
# nx: number of entries to put to zero

soft_thresholding_L1 = function(x,nx)
{
    #selection on a (loadings.X). modified on 19/02/15 to make sure that a!=0
    if (nx!=0)
    {
        absa = abs(x)
        if (any(rank(absa, ties.method = "max") <= nx))
        {
            x = ifelse(rank(absa, ties.method = "max") <= nx, 0,
                        sign(x) * (absa - max(absa[rank(absa, ties.method = "max") <= nx])))
        }
    }
    
    x
}

# ----------------------------------------------------------------------------------------------------------
# soft.threshold() - soft-thresholds a vector such that the L1-norm constraint is satisfied.
# ----------------------------------------------------------------------------------------------------------
soft.threshold = function (x, sumabs = 1)
return(soft(x, BinarySearch(x, sumabs)))

BinarySearch = function(argu,sumabs)
{
    if (norm2(argu)==0 || sum(abs(argu/norm2(argu)))<=sumabs)
    return(0)
    
    lam1 = 0
    lam2 = max(abs(argu))-1e-5
    iter = 1
    while (iter < 500)
    {
        su = soft(argu,(lam1+lam2)/2)
        if (sum(abs(su/norm2(su)))<sumabs)
        {
            lam2 = (lam1+lam2)/2
        } else {
            lam1 = (lam1+lam2)/2
        }
        if ((lam2-lam1)<1e-10)
        return((lam1+lam2)/2)
        
        iter = iter+1
    }
    warning("Didn't quite converge")
    return((lam1+lam2)/2)
}

soft = function(x,d) return(sign(x)*pmax(0, abs(x)-d))

norm2 = function(vec)
{
    a = sqrt(sum(vec^2))
    if (a == 0)
    a = .05
    
    return(a)
}


# --------------------------------------
# sparsity function: used in 'internal_mint.block.R'
# --------------------------------------
sparsity=function(loadings.A, keepA, keepA.constraint=NULL, penalty=NULL)
{
    
    if (!is.null(keepA.constraint))
    {
        loadings.A[-keepA.constraint] = 0
    } else if (!is.null(keepA)) {
        nx = length(loadings.A) - keepA
        loadings.A = soft_thresholding_L1(loadings.A, nx = nx)
    } else if (!is.null(penalty)) {
        loadings.A = soft.threshold(loadings.A, penalty)
    }

    return(loadings.A)
}



# --------------------------------------
# scaling with or without bias: used in mean_centering_per_study (below)
# --------------------------------------
scale.function=function(temp, scale = TRUE, bias = FALSE)
{
    meanX = colMeans(temp, na.rm = TRUE)
    data.list.study.scale_i = t(t(temp) - meanX)
    if (scale)
    {
        if (bias)
        {
            sqrt.sdX = sqrt(colSums(data.list.study.scale_i^2, na.rm = TRUE) / (nrow(temp)))
        } else {
            sqrt.sdX = sqrt(colSums(data.list.study.scale_i^2, na.rm = TRUE) / (nrow(temp) - 1))
        }
        data.list.study.scale_i = t(t(data.list.study.scale_i) / sqrt.sdX)
    } else {
        sqrt.sdX = NULL
    }
    
    is.na.data = is.na(data.list.study.scale_i)
    #if (sum(is.na.data) > 0)
    #data.list.study.scale_i[is.na.data] = 0
    
    out = list(data_scale=data.list.study.scale_i, meanX=meanX, sqrt.sdX=sqrt.sdX)
    return(out)
}

# --------------------------------------
# Mean centering/scaling per study: used in 'internal_mint.block.R'
# --------------------------------------
mean_centering_per_study=function(data, study, scale, bias=FALSE)
{
    
    M = length(levels(study))   # number of groups
    # split the data
    data.list.study = study_split(data, study)

    # center and scale data per group, and concatene the data
    res = lapply(data.list.study, scale.function, scale = scale, bias = bias)
    concat.data = do.call("rbind", lapply(res,function(x){x[[1]]}))
    meanX = lapply(res, function(x){x[[2]]})
    sqrt.sdX = lapply(res, function(x){x[[3]]})
    rownames.study = lapply(res, function(x){rownames(x[[1]])})

    #rename rows and cols of concatenated centered (and/or scaled) data
    colnames(concat.data) = colnames(data)
    
    #sort the samples as in the original X
    indice.match = match(rownames(data),rownames(concat.data))
    concat.data = concat.data[indice.match, ,drop=FALSE]
    
    if (M > 1)
    {
        for (m in 1:M)
        {
            attr(concat.data,paste0("means:", levels(study)[m])) = meanX[[m]]
            if(scale)
            {
                attr(concat.data,paste0("sigma:", levels(study)[m])) = sqrt.sdX[[m]]
            } else {
                attr(concat.data,paste0("sigma:", levels(study)[m])) = NULL
            }
        }
    } else {
        attr(concat.data,"scaled:center") = meanX[[1]]
        if (scale)
        {
            attr(concat.data,"scaled:scale") = sqrt.sdX[[1]]
        } else {
            attr(concat.data,"scaled:scale") = NULL
        }
    }
    
    return(list(concat.data=concat.data, rownames.study=rownames.study))
}


# --------------------------------------
# l2.norm: used in 'internal_mint.block.R'
# --------------------------------------
l2.norm=function(x)
{
    if (!is.vector(x))
    stop("x has to be a vector")
    
    out = x / drop(sqrt(crossprod(x)))
}

# ---------------------------------------------------
# tau.estimate() - Estimation of tau accoring to Strimmer formula
# ---------------------------------------------------
#used in 'internal_mint.block.R'
tau.estimate = function (x)
{
    if (is.matrix(x) == TRUE && is.numeric(x) == FALSE)
    stop("The data matrix must be numeric!")
    
    p = NCOL(x)
    n = NROW(x)
    #covm = cov(x)
    corm = cor(x)
    xs = scale(x, center = TRUE, scale = TRUE)
    xs2 = xs^2
    v = (n/((n - 1)^3)) * (crossprod(xs2) - 1/n * (crossprod(xs))^2)
    diag(v) = 0
    m = matrix(rep(apply(xs2, 2, mean), p), p, p)
    I = diag(NCOL(x))
    d = (corm - I)^2
    tau = (sum(v))/sum(d)
    tau = max(min(tau, 1), 0)
    return(tau)
}


#############################################################################################################
# Functions acquired from RGCCA R-library
#############################################################################################################
# ----------------------------------------------------------------------------------------------------------
# cov2() - Compute biased and unbiased covariance and variance estimates
# ----------------------------------------------------------------------------------------------------------
# used in 'internal_mint.block.R'
cov2 = function (x, y = NULL, bias = TRUE) {
    n = NROW(x)
    if (is.null(y)) {
        x = as.matrix(x)
        if (bias) {
            C = ((n - 1)/n) * cov(x, use = "pairwise.complete.obs")
        } else {
            C = cov(x, use = "pairwise.complete.obs")
        }
    } else {
        if (bias) {
            C = ((n - 1)/n) * cov(x, y, use = "pairwise.complete.obs")
        } else {
            C = cov(x, y, use = "pairwise.complete.obs")
        }
    }
    return(C)
}

# ----------------------------------------------------------------------------------------------------------
# initsvd() - performs SVD on matrix X
# ----------------------------------------------------------------------------------------------------------
# used in 'internal_mint.block.R'
initsvd = function (X) {
    n = NROW(X)
    p = NCOL(X)
    ifelse(n >= p, return(svd(X, nu = 0, nv = 1)$v), return(svd(X, nu = 1, nv = 0)$u))
}

# ----------------------------------------------------------------------------------------------------------
# miscrossprod() - Compute cross-product between vectors x and y
# ----------------------------------------------------------------------------------------------------------
# used in 'internal_mint.block.R'
miscrossprod = function (x, y) {
    d.p = sum(drop(x) * drop(y), na.rm = TRUE)
    #d.p = as.vector(d.p)/norm2(d.p)     ## change made
    return(d.p)
}


# ----------------------------------------------------------------------------------------------------------
# deflation()
# ----------------------------------------------------------------------------------------------------------
# used in defl.select (below)
deflation = function(X, y){
    # Computation of the residual matrix R
    # Computation of the vector p.
    is.na.tX = is.na(t(X))
    if (any(is.na.tX))
    {
        #p = apply(t(X),1,miscrossprod,y)/as.vector(crossprod(y))
        
        #variates.A[, q] =  apply(A[[q]], 1, miscrossprod, loadings.A[[q]])
        A.temp = replace(t(X), is.na.tX, 0) # replace NA in A[[q]] by 0
        variates.A.temp = A.temp %*% y
        temp = drop(y) %o% rep(1, nrow(A.temp))
        temp[(t(is.na.tX))] = 0
        loadings.A.norm = crossprod(temp)
        p = variates.A.temp / diag(loadings.A.norm)
        # we can have 0/0, so we put 0
        a = is.na(p)
        if (any(a))
        p[a] = 0
        
    } else {
        p = t(X)%*%y/as.vector(crossprod(y))
    }
    
    R = X - y%*%t(p)
    return(list(p=p,R=R))
}

# ----------------------------------------------------------------------------------------------------------
# defl.select() - computes residual matrices
# ----------------------------------------------------------------------------------------------------------
# used in 'internal_mint.block.R'
defl.select = function(yy, rr, nncomp, nn, nbloc, indY = NULL, mode = "canonical", aa = NULL) { ### Start: Add new parameter for estimation classic mode
    resdefl = NULL
    pdefl = NULL
    for (q in 1 : nbloc) {
        ### Start: insertion of new deflations (See La regression PLS Theorie et pratique (page 139))
        if ( nn <= nncomp[q] ) {
            if ((mode == "canonical") || (q != indY)) { #deflation of each block independently from the others, except indY
                defltmp = deflation(rr[[q]], yy[ , q])
                resdefl[[q]] = defltmp$R
                pdefl[[q]]   = defltmp$p
            } else if (mode == "classic") {
                resdefl[[q]] = Reduce("+", lapply(c(1:nbloc)[-q], function(x) {rr[[q]] - yy[ ,x]%*%t(aa[[q]])}))/(nbloc-1)
                pdefl[[q]]   =  rep(0,NCOL(rr[[q]]))
            } else if (mode == "invariant") { #no deflation
                resdefl[[q]] = rr[[q]]
                pdefl[[q]]   =  rep(0,NCOL(rr[[q]]))
            } else if (mode == "regression") {
                resdefl[[q]] = Reduce("+", lapply(c(1:nbloc)[-q], function(x) {deflation(rr[[q]],yy[, x])$R}))/(nbloc-1)
                pdefl[[q]]   =  rep(0,NCOL(rr[[q]]))
            }
            ### End: insertion of new deflations (See La regression PLS Theorie et pratique (page 139))
        } else {
            resdefl[[q]] = rr[[q]]
            pdefl[[q]]   =  rep(0,NCOL(rr[[q]]))
        }
    }
    return(list(resdefl=resdefl,pdefl=pdefl))
}



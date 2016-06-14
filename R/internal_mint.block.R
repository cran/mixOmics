#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Amrit Singh, the University of British Columbia, Vancouver, Canada
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 20-07-2014
# last modified: 12-04-2016
#
# Copyright (C) 2014
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


#############################################################################################################
# Functions modified from RGCCA R-library
#   sgcca(), sgccak()
#
# Functions acquired from RGCCA R-library, see 'internal_mint.block_helpers.R'
#   cov2(), initsvd(), crossprod(),
#   defl.select(),
#############################################################################################################

internal_mint.block = function (A, indY = NULL,  design = 1 - diag(length(A)), tau=NULL,#rep(1, length(A)),
ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,  bias = FALSE,
init = "svd.single", tol = 1e-06, verbose = FALSE,
mode = "canonical", max.iter = 500,study = NULL, keepA = NULL,
keepA.constraint = NULL, penalty = NULL)
{
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
    
    names(ncomp) = names(A)
    
    # keepA is updated to be of length length(A) now, the first entries correspond to the keepA.constraint if it was provided
    for (q in 1:length(A))
    keepA[[q]] = c(unlist(lapply(keepA.constraint[[q]], length)), keepA[[q]]) #of length ncomp, can contains 0

    # center the data per study, per matrix of A, scale if scale=TRUE, option bias
    mean_centered = lapply(A, function(x){mean_centering_per_study(x, study, scale, bias)})
    A = lapply(mean_centered, function(x){as.matrix(x$concat.data)})
    ni = table(study) #number of samples per study
    
    ### Start: Initialization parameters
    pjs = sapply(A, NCOL)
    nb_ind = NROW(A[[1]])
    J = length(A)
    R = A # R: residuals matrices, will be a list of length ncomp
    N = max(ncomp)
    
    AVE_inner = AVE_outer = rep(NA, max(ncomp))
    defl.matrix = AVE_X = crit = loadings.partial.A = variates.partial.A = tau.rgcca = list()
    P = loadings.A = loadings.Astar = c = t = b = variates.A = vector("list", J)
    
    for (k in 1:J)
    t[[k]] = variates.A[[k]] = matrix(NA, nb_ind, N)
    
    for (k in 1:J)
    P[[k]] = loadings.A[[k]] = loadings.Astar[[k]]= matrix(NA, pjs[[k]], N)
    
    for (k in 1:J)
    {
        loadings.partial.A[[k]] = variates.partial.A[[k]] = vector("list", length = nlevels(study))
        for(m in 1:nlevels(study))
        {
            loadings.partial.A[[k]][[m]] = matrix(nrow = NCOL(A[[k]]), ncol = N)
            variates.partial.A[[k]][[m]] = matrix(nrow = ni[m], ncol = N)
        }
    }
    
    defl.matrix[[1]] = A
    ndefl = ncomp - 1
    J2 = J-1
    
    if (is.vector(tau))
    tau = matrix(rep(tau, N), nrow = N, ncol = length(tau), byrow = TRUE)
    ### End: Initialization parameters
    
    
    iter=NULL
    for (n in 1 : N)
    {
        ### Start: Estimation ai
        if (is.null(tau))
        {
            mint.block.result = sparse.mint.block_iteration(R, design, study = study,
            keepA.constraint = if (!is.null(keepA.constraint)) {lapply(keepA.constraint, function(x){unlist(x[n])})} else {NULL} ,
            keepA = if (!is.null(keepA)) {lapply(keepA, function(x){x[n]})} else {NULL},
            scheme = scheme, init = init, max.iter = max.iter, tol = tol, verbose = verbose,penalty = penalty)
        } else {
            mint.block.result = sparse.rgcca_iteration(R, design, tau = if (is.matrix(tau)){tau[n, ]} else {"optimal"},
            scheme = scheme, init = init, tol = tol,
            verbose = verbose, max.iter = max.iter, penalty = penalty,
            keepA = if (!is.null(keepA)) {lapply(keepA, function(x){x[n]})} else {NULL})
        }
        ### End: Estimation ai
        
        if(is.null(tau))
        {
            #recording loadings.partials, $Ai$study[,ncomp]
            # recording variates.partials, $Ai[,ncomp]
            for(k in 1:J)
            {
                for(m in 1:nlevels(study))
                {
                    loadings.partial.A[[k]][[m]][, n] = matrix(mint.block.result$loadings.partial.A.comp[[k]][[m]], ncol=1)
                    variates.partial.A[[k]][[m]][, n] = matrix(mint.block.result$variates.partial.A.comp[[k]][[m]], ncol=1)
                }
                #variates.partial.A[[k]][,n]=matrix(unlist(mint.block.result$variates.partial.A.comp[[k]]),ncol=1)
            }
        }
        
        AVE_inner[n] = mint.block.result$AVE_inner
        crit[[n]] = mint.block.result$crit
        tau.rgcca[[n]] = mint.block.result$tau
        
        for (k in 1:J)
        variates.A[[k]][, n] = mint.block.result$variates.A[, k]
        
        # deflation if there are more than 1 component and if we haven't reach the max number of component (N)
        if (N != 1 & n != N)
        {
            defla.result = defl.select(mint.block.result$variates.A, R, ndefl, n, nbloc = J, indY = indY, mode = mode, aa = mint.block.result$loadings.A)
            R = defla.result$resdefl
            defl.matrix[[n + 1]] = R
        }
        
        for (k in 1 : J)
        {
            if (N != 1)
            {
                P[[k]][, n - 1] = defla.result$pdefl[[k]]
            }
            loadings.A[[k]][, n] = mint.block.result$loadings.A[[k]]
        }
        
        if (n == 1)
        {
            for (k in 1 : J)
            loadings.Astar[[k]][, n] = mint.block.result$loadings.A[[k]]
        } else {
            for (k in 1 : J)
            loadings.Astar[[k]][, n] = mint.block.result$loadings.A[[k]] - loadings.Astar[[k]][, (1 : n - 1), drop = F] %*% drop(t(loadings.A[[k]][, n]) %*% P[[k]][, 1 : (n - 1), drop = F])
        }
        iter = c(iter, mint.block.result$iter)
    }
    
    if (verbose)
    cat(paste0("Computation of the SGCCA block components #", N , " is under progress...\n"))
    
    shave.matlist = function(mat_list, nb_cols) mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY = FALSE)
    shave.veclist = function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
    
    for (k in 1:J)
    {
        rownames(loadings.A[[k]]) = rownames(loadings.Astar[[k]])=colnames(A[[k]])
        rownames(variates.A[[k]]) = rownames(A[[k]])
        colnames(variates.A[[k]]) = colnames(loadings.A[[k]]) = paste0("comp ", 1:max(ncomp))
        AVE_X[[k]] = apply(cor(A[[k]], variates.A[[k]])^2, 2, mean)
        if (is.null(tau))
        {
            names(loadings.partial.A[[k]]) = names(variates.partial.A[[k]]) = levels(study)
            for (m in 1:nlevels(study))
            {
                rownames(loadings.partial.A[[k]][[m]]) = colnames(A[[k]])
                colnames(loadings.partial.A[[k]][[m]]) = paste0("comp ", 1:max(ncomp))
                rownames(variates.partial.A[[k]][[m]]) = rownames(mean_centered[[1]]$rownames.study[[m]])
                colnames(variates.partial.A[[k]][[m]]) = paste0("comp ", 1:max(ncomp))
            }
        }
    }
    
    outer = matrix(unlist(AVE_X), nrow = max(ncomp))
    for (j in 1 : max(ncomp))
    AVE_outer[j] = sum(pjs * outer[j, ])/sum(pjs)
    
    variates.A = shave.matlist(variates.A, ncomp)
    AVE_X = shave.veclist(AVE_X, ncomp)
    AVE = list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)
    
    #calcul explained variance
    A_split=lapply(A, study_split, study) #split the data per study

    expl.A=lapply(1:length(A),function(x){
        if (nlevels(study) == 1)
        {
            temp = suppressWarnings(explained_variance(A[[x]], variates = variates.A[[x]], ncomp = ncomp[[x]]))
        }else{
            temp = lapply(1:nlevels(study), function(y){
                suppressWarnings(explained_variance(A_split[[x]][[y]], variates = variates.partial.A[[x]][[y]], ncomp = ncomp[[x]]))})
            temp[[length(temp)+1]] = explained_variance(A[[x]], variates = variates.A[[x]], ncomp = ncomp[[x]])
            names(temp) = c(levels(study), "all data")
        }
        temp
        })
    names(expl.A) = names(A)
    
    ### Start: output
    names(loadings.A) = names(variates.A) = names(A)
    
    if (is.null(tau))
    names(loadings.partial.A) = names(variates.partial.A) = names(A)
    
    names = lapply(1:J, function(x) {colnames(A[[x]])})
    names(names) = names(A)
    names[[length(names) + 1]] = row.names(A[[1]])
    names(names)[length(names)] = "indiv"
    
    out = list(X = A, indY = indY, ncomp = ncomp, mode = mode,
    keepA = keepA, keepA.constraint = keepA.constraint,
    variates = variates.A, loadings = shave.matlist(loadings.A, ncomp),
    variates.partial= if(is.null(tau)) {variates.partial.A} ,loadings.partial= if(is.null(tau)) {loadings.partial.A},
    loadings.star = shave.matlist(loadings.Astar, ncomp),
    names = list(sample = row.names(A[[1]]), colnames = lapply(A, colnames), blocks = names(A)),
    tol = tol, iter=iter, max.iter=max.iter,
    design = design,
    scheme = scheme,  crit = crit, AVE = AVE, defl.matrix = defl.matrix,
    init = init, bias = bias,
    scale = scale, tau = if(!is.null(tau)) tau.rgcca, study = study,
    explained_variance = expl.A)
    ### End: Output
    
    return(out)
}

# ----------------------------------------------------------------------------------------------------------
# sgccak - Runs sgccak() modified from RGCCA
#   inputs: A - list of datasets each with the same number of rows (samples)
#           design - design matrix
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ----------------------------------------------------------------------------------------------------------

sparse.mint.block_iteration = function (A, design, study = NULL, keepA.constraint = NULL, keepA = NULL,
scheme = "centroid", init = "svd", max.iter = 500, tol = 1e-06, verbose = TRUE, bias = FALSE,
penalty=NULL)
{
    
    # keepA.constraint is a list of positions in A of the variables to keep on the component
    # keepA is a vector of numbers
    # study is a vector
    # == == == == no check needed as this function is only used in internal_mint.block, in which the checks are conducted
    
    ### Start: Initialization parameters
    J = length(A)
    J2 = J-1
    pjs = sapply(A, NCOL)
    AVE_X = rep(0, J)
    if (!is.null(penalty))
    penalty = penalty * sqrt(pjs)
    
    
    iter = 1
    converg = crit = numeric()
    variates.A = Z = matrix(0, NROW(A[[1]]), J)
    misdata = any(sapply(A, function(x){any(is.na(x))})) # Detection of missing data
    
    if (misdata)
    is.na.A = lapply(A, is.na)
    
    g = function(x) switch(scheme, horst = x, factorial = x^2, centroid = abs(x))
    
    
    # study split
    A_split = lapply(A, study_split, study)
    if (misdata)
    is.na.A_split = lapply(is.na.A, study_split, study)
    
    nlevels_study = nlevels(study)
    ### End: Initialization parameters
    
    ### Start: Initialisation "loadings.A" vector
    if (init == "svd")
    {
        ### Start: Change initialization of loadings.A
        if (misdata)
        {
            M = lapply(c(1:(J-1)), function(x){crossprod(replace(A[[x]], is.na(A[[x]]), 0), replace(A[[J]], is.na(A[[J]]), 0))})
        } else {
            M = lapply(c(1:(J-1)), function(x){crossprod(A[[x]], A[[J]])})
        }
        
        svd.M = lapply(M, function(x){svd(x, nu = 1, nv = 1)})
        loadings.A = lapply(c(1:(J-1)), function(x){svd.M[[x]]$u})
        loadings.A[[J]] = svd.M[[1]]$v
    } else if (init=="svd.single") {
        alpha =  lapply(1 : J, function(y){initsvd(lapply(y, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])})
        
        loadings.A = list()
        for (j in 1:J)
        {
            if (nrow(A[[j]]) >= ncol(A[[j]]))
            {
                loadings.A[[j]] = alpha[[j]]
            } else {
                K = as.matrix(A[[j]]) %*% as.matrix(t(A[[j]]))
                N = ifelse(bias, nrow(A[[j]]), nrow(A[[j]]) - 1)
                alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% K %*% alpha[[j]])) * alpha[[j]]
                loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
            }
        }
        
        ### End: Change initialization of a
    } else {
        stop("init should be either 'svd' or 'svd.single'.")
    }
    ### End: Initialisation "a" vector
    variates.partial.A.comp = NULL
    loadings.partial.A.comp = list()
    for (q in 1:J)
    {
        if(misdata)
        {
            #variates.A[, q] =  apply(A[[q]], 1, miscrossprod, loadings.A[[q]])
            A.temp = replace(A[[q]], is.na.A[[q]], 0) # replace NA in A[[q]] by 0
            variates.A.temp = A.temp %*% loadings.A[[q]]
            temp = drop(loadings.A[[q]]) %o% rep(1, nrow(A[[q]]))
            temp[(t(is.na.A[[q]]))] = 0
            loadings.A.norm = crossprod(temp)
            variates.A[, q] = variates.A.temp / diag(loadings.A.norm)
            # we can have 0/0, so we put 0
            a = is.na(variates.A[, q])
            if (any(a))
            variates.A[a, q] = 0
        }else{
            variates.A[, q] = A[[q]]%*%loadings.A[[q]]
        }
        loadings.A[[q]] = l2.norm(as.vector(loadings.A[[q]]))
        loadings.partial.A.comp[[q]] = list()
    }
    loadings.A_old = loadings.A
    
    
    ### Start Algorithm 1 Sparse generalized canonical analysis (See Variable selection for generalized canonical correlation analysis (Tenenhaus))
    repeat {
        #variates.Aold can be used for the convergence of the algorithm, not used at the moment, so please disregard variates.Aold
        variates.Aold = variates.A
        for (q in 1:J)
        {
            ### Start : !!! Impact of the diag of the design matrix !!! ###
            if (scheme == "horst")
            CbyCovq = design[q, ]
            
            if (scheme == "factorial")
            CbyCovq = design[q, ] * cov2(variates.A, variates.A[, q], bias = bias)
            
            if (scheme == "centroid")
            CbyCovq = design[q, ] * sign(cov2(variates.A, variates.A[, q], bias = bias))
            ### End : !!! Impact of the diag of the design matrix !!! ###
            
            ### Step A start: Compute the inner components
            Z[, q] = rowSums(mapply("*", CbyCovq, as.data.frame(variates.A)))
            Z_split = study_split(Z[,q,drop=FALSE],study)  # split Z by the study factor
            ### Step A end: Compute the inner components
            
            
            ### Step B start: Computer the outer weight ###
            # possibility of removing NA (replacing by 0) and use crossprod, further development
            temp=0
            for (m in 1:nlevels_study)
            {
                if(misdata)
                {
                    loadings.partial.A.comp[[q]][[m]] = apply(t(A_split[[q]][[m]]), 1, miscrossprod, Z_split[[m]])
                }else{
                    loadings.partial.A.comp[[q]][[m]] = t(A_split[[q]][[m]])%*%Z_split[[m]]
                }
                temp=temp+loadings.partial.A.comp[[q]][[m]]
            }
            loadings.A[[q]] = temp
            
            
            # sparse using keepA / penalty
            if (!is.null(penalty))
            {
                loadings.A[[q]] = sparsity(loadings.A[[q]], keepA = NULL,keepA.constraint = NULL, penalty = penalty[q])
            }else{
                loadings.A[[q]] = sparsity(loadings.A[[q]], keepA[[q]], keepA.constraint[[q]], penalty = NULL)
            }
            
            loadings.A[[q]]=l2.norm(as.vector(loadings.A[[q]]))
            
            ### Step B end: Computer the outer weight ###
            if(misdata)
            {
                #variates.A[, q] =  apply(A[[q]], 1, miscrossprod, loadings.A[[q]])
                A.temp = replace(A[[q]], is.na.A[[q]], 0) # replace NA in A[[q]] by 0
                variates.A.temp = A.temp %*% loadings.A[[q]]
                temp = drop(loadings.A[[q]]) %o% rep(1, nrow(A[[q]]))
                temp[(t(is.na.A[[q]]))] = 0
                loadings.A.norm = crossprod(temp)
                variates.A[, q] = variates.A.temp / diag(loadings.A.norm)
                # we can have 0/0, so we put 0
                a = is.na(variates.A[, q])
                if (any(a))
                variates.A[a, q] = 0
                
            }else{
                variates.A[, q] =  A[[q]]%*%loadings.A[[q]]
            }
            
        }
        
        crit[iter] = sum(design * g(cov2(variates.A, bias = bias)))
        
        if (iter > max.iter)
        warning("The SGCCA algorithm did not converge", call. = FALSE)#cat("The SGCCA algorithm did not converge after", max.iter ,"iterations."))
        
        ### Start: Match algorithm with mixOmics algo (stopping point)
        ### if ((converg[iter] < tol & sum(stationnary_point) == J) | iter > max.iter)
        if (max(sapply(1:J, function(x){crossprod(loadings.A[[x]] - loadings.A_old[[x]])})) < tol | iter > max.iter)
        break
        ### End: Match algorithm with mixOmics algo (stopping point)
        
        loadings.A_old = loadings.A
        iter = iter + 1
    }
    ### End Algorithm 1 (See Variable selection for generalized canonical correlation analysis (Tenenhaus))
    
    
    #calculation variates.partial.A.comp
    variates.partial.A.comp = lapply(1 : J, function(x){
        lapply(1 : nlevels_study, function(y){
            if(misdata)
            {
                #apply(A_split[[x]][[y]], 1, miscrossprod, loadings.A[[x]])
                A.temp = replace(A_split[[x]][[y]], is.na.A_split[[x]][[y]], 0) # replace NA in A_split[[x]][[y]] by 0
                variates.part.A.temp = A.temp %*% loadings.A[[x]]
                temp = drop(loadings.A[[x]]) %o% rep(1, nrow(A_split[[x]][[y]]))
                temp[(t(is.na.A[[x]][[y]]))] = 0
                loadings.A.norm = crossprod(temp)
                variates.part.A = variates.part.A.temp / diag(loadings.A.norm)
                # we can have 0/0, so we put 0
                a = is.na(variates.part.A)
                if (any(a))
                variates.part.A[a] = 0
                
                return(variates.part.A)
            }else{
                A_split[[x]][[y]] %*% loadings.A[[x]]
            }
        })
    })
    
    if (verbose)
    plot(crit, xlab = "iteration", ylab = "criteria")
    
    AVE_inner = sum(design * cor(variates.A)^2/2)/(sum(design)/2)
    
    result = list(variates.A = variates.A, loadings.A = loadings.A, crit = crit[which(crit != 0)],
    AVE_inner = AVE_inner, loadings.partial.A.comp = loadings.partial.A.comp, variates.partial.A.comp = variates.partial.A.comp, iter = iter)
    return(result)
}

# ----------------------------------------------------------------------------------------------------------
# rgccak - Runs sgccak() modified from RGCCA
#   inputs: A - list of datasets each with the same number of rows (samples)
#           design - design matrix
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ----------------------------------------------------------------------------------------------------------


sparse.rgcca_iteration = function (A, design, tau = "optimal", scheme = "centroid", scale = FALSE, max.iter = 500,
verbose = FALSE, init = "svd.single", bias = FALSE, tol = .Machine$double.eps, keepA = NULL, penalty = NULL)
{
    ### Start: Initialisation parameters
    A = lapply(A, as.matrix)
    J = length(A)
    n = NROW(A[[1]])
    pjs = sapply(A, NCOL)
    variates.A = matrix(0, n, J)
    if (!is.null(penalty))
    penalty = penalty * sqrt(pjs)
    ### End: Initialisation parameters
    
    if (!is.numeric(tau))
    tau = sapply(A, tau.estimate)
    
    loadings.A = alpha = M = Minv = K = list()
    which.primal = which((n >= pjs) == 1)
    which.dual = which((n < pjs) == 1)
    
    if (init == "svd.single")
    {
        for (j in which.primal)
        loadings.A[[j]] = initsvd(lapply(j, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])

        for (j in which.dual)
        {
            alpha[[j]] = initsvd(lapply(j, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])
            K[[j]] = A[[j]] %*% t(A[[j]])
        }
    } else {
        stop("init should be 'svd.single'.")
    }
    
    N = ifelse(bias, n, n - 1)
    for (j in 1 : J)
    {
        if (j %in% which.primal)
        {
            M[[j]] = ginv(tau[j] * diag(pjs[j]) + (1 - tau[j]) * cov2(A[[j]], bias = bias))
            loadings.A[[j]] = drop(1/sqrt(t(loadings.A[[j]]) %*% M[[j]] %*% loadings.A[[j]])) * M[[j]] %*% loadings.A[[j]]
        }
        
        if (j %in% which.dual)
        {
            M[[j]] = tau[j] * diag(n) + (1 - tau[j])/(N) * K[[j]]
            Minv[[j]] = ginv(M[[j]])
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% M[[j]] %*% K[[j]] %*% alpha[[j]])) * alpha[[j]]
            loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
        }
        variates.A[, j] = A[[j]] %*% loadings.A[[j]]
    }
    
    iter = 1
    converg = crit = numeric()
    Z = matrix(0, NROW(A[[1]]), J)
    loadings.A_old = loadings.A
    g = function(x) switch(scheme, horst = x, factorial = x^2, centroid = abs(x))
    
    repeat {
        variates.Aold = variates.A
        
        for (j in c(which.primal, which.dual))
        {
            
            if (scheme == "horst")
            CbyCovq = design[j, ]
            
            if (scheme == "factorial")
            CbyCovq = design[j, ] * cov2(variates.A, variates.A[, j], bias = bias)
            
            if (scheme == "centroid")
            CbyCovq = design[j, ] * sign(cov2(variates.A, variates.A[, j], bias = bias))
            
            # Compute the inner components
            Z[, j] = rowSums(mapply("*", CbyCovq, as.data.frame(variates.A)))
            
            # Computer the outer weight
            if (j %in% which.primal)
            loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
            
            # Compute the outer weight
            if (j %in% which.dual)
            {
                alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
            }
            
            # sparse using keepA / penalty
            if (!is.null(keepA) || !is.null(penalty))
            {
                temp.norm = norm2(loadings.A[[j]])
                if (!is.null(keepA))
                {
                    loadings.A[[j]] = sparsity(loadings.A = loadings.A[[j]], keepA = keepA[[j]], penalty = NULL)
                } else if (!is.null(penalty)) {
                    loadings.A[[j]] = sparsity(loadings.A = loadings.A[[j]], keepA = NULL, penalty = penalty[j])
                }
                loadings.A[[j]] = (loadings.A[[j]]/norm2(loadings.A[[j]]))*temp.norm
            }
            
            # Update variate
            variates.A[, j] = A[[j]] %*% loadings.A[[j]]
        }
        
        crit[iter] = sum(design * g(cov2(variates.A, bias = bias)))
        
        if (iter > max.iter)
        warning("The RGCCA algorithm did not converge")#"cat("The RGCCA algorithm did not converge after", max.iter ,"iterations."))
        
        ### Start: Match algorithm with mixOmics algo (stopping point)
        if (max(sapply(1:J, function(x){crossprod(loadings.A[[x]] - loadings.A_old[[x]])})) < tol | iter > max.iter)
        break
        ### End: Match algorithm with mixOmics algo (stopping point)
        
        loadings.A_old = loadings.A
        iter = iter + 1
    }
    
    if (verbose) 
    plot(crit, xlab = "iteration", ylab = "criteria")
    
    AVE_inner = sum(design * cor(variates.A)^2/2)/(sum(design)/2)
    
    result = list(variates.A = variates.A, loadings.A = loadings.A, crit = crit[which(crit != 0)], 
    AVE_inner = AVE_inner, design = design, tau = tau, scheme = scheme,iter=iter, keepA = keepA)
    return(result)
}


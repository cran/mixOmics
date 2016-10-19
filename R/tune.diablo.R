#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2013
# last modified: 24-08-2016
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


# ========================================================================================================
# tune.splsda: chose the optimal number of parameters per component on a splsda method
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 1.
# test.keepX: grid of keepX among which to chose the optimal one
# already.tested.X: a vector giving keepX on the components that were already tuned
# validation: Mfold or loo cross validation
# folds: if validation=Mfold, how many folds?
# dist: distance to classify samples. see predict
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# progressBar: show progress,
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# nrepeat: number of replication of the Mfold process
# logratio = c('none','CLR'). see splsda
# verbose: if TRUE, shows component and nrepeat being tested.


tune.block.splsda = function (
X,
Y,
indY,
ncomp = 2,
test.keepX,
already.tested.X,
constraint = FALSE, #if TRUE, expect a list in already.tested.X, otherwise a number(keepX)
validation = "Mfold",
folds = 10,
dist = "max.dist",
measure = "BER", # one of c("overall","BER")
weighted = TRUE, # optimise the weighted or not-weighted prediction
progressBar = TRUE,
max.iter = 100,
near.zero.var = FALSE,
nrepeat = 1,
design,
scheme,
mode,
scale = TRUE,
bias,
init ,
tol = 1e-06,
verbose,
light.output = TRUE, # if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
cpus,
name.save = NULL)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #------------------#
    #-- check entries --#
    
    
    # check inpuy 'Y' and transformation in a dummy matrix
    if(!missing(Y))
    {
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
        stop("'Y' should be a factor with more than one level")
        
    } else if(!missing(indY)) {
        Y = X[[indY]]
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(temp) == 1)
        stop("'X[[indY]]' should be a factor with more than one level")
        
        X = X[-indY] #remove Y from X to pass the arguments simpler to block.splsda
        
    } else if(missing(indY)) {
        stop("Either 'Y' or 'indY' is needed")
        
    }
    
    #-- dist
    dist = match.arg(dist, choices = c("max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = FALSE)

    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    #-- ncomp
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    
    #-- validation
    choices = c("Mfold", "loo")
    validation = choices[pmatch(validation, choices)]
    if (is.na(validation))
    stop("'validation' must be either 'Mfold' or 'loo'")
    
    if (validation == 'loo')
    {
        if (nrepeat != 1)
        warning("Leave-One-Out validation does not need to be repeated: 'nrepeat' is set to '1'.")
        nrepeat = 1
    }
    
    
    #-- measure
    measure.input = measure
    if(measure == "BER")
    {
        measure = "Overall.BER"
    } else if (measure == "overall"){
        measure = "Overall.ER"
    } else {
        stop("'measure' must be 'overall' or 'BER'")
    }
    
    
    #-- already.tested.X
    if (missing(already.tested.X))
    {
        if(constraint == TRUE)
        {
            already.tested.X = list()
        } else {
            already.tested.X = NULL
        }
    } else {
        if(!is.list(already.tested.X))
        stop("'already.tested.X' must be a list, each entry corresponding to a block of X (Y excluded)")
        
        if(is.null(already.tested.X) | length(already.tested.X)==0)
        stop("'already.tested.X' must be a vector of keepX values (if 'constraint'= FALSE) or a list (if'constraint'= TRUE) ")
        
        # we require the same number of already tuned components on each block
        if(length(unique(sapply(already.tested.X, length))) > 1)
        stop("The same number of components must be already tuned for each block, in 'already.tested.X'")
        
        
        if(constraint == TRUE)
        {
            if(any(sapply(already.tested.X, function(x) is.list(x))) != TRUE)
            stop(" Each entry of 'already.tested.X' must be a list since 'constraint' is set to TRUE")
            
            #print(paste("A total of",lapply(already.tested.X, function(x){sapply(x,length)}),collapse=" "),"specific variables ('already.tested.X') were selected on the first ", length(already.tested.X[[1]]), "component(s)"))
        } else {
            if(any(sapply(already.tested.X, function(x) is.list(x))) == TRUE)
            stop(" Each entry of 'already.tested.X' must be a vector of keepX values since 'constraint' is set to FALSE")
            
            #print(paste("Number of variables selected on the first", length(already.tested.X), "component(s):", paste(already.tested.X,collapse = " ")))
        }
        if(length(already.tested.X[[1]]) >= ncomp)
        stop("'ncomp' needs to be higher than the number of components already tuned, which is length(already.tested.X)=",length(already.tested.X) , call. = FALSE)
    }
    
    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- test.keepX
    if(missing(test.keepX))
    {
        test.keepX = lapply(1:length(X),function(x){c(5,10,15)[which(c(5,10,15)<ncol(X[[x]]))]})
        names(test.keepX) = names(X)
        
    } else {
        if(length(test.keepX) != length(X))
        stop(paste("test.keepX should be a list of length ", length(X),", corresponding to the blocks: ", paste(names(X),collapse=", "), sep=""))
        
        aa = sapply(test.keepX, length)
        if (any(is.null(aa) | aa == 1 | !is.numeric(aa)))
        stop("Each entry of 'test.keepX' must be a numeric vector with more than two values", call. = FALSE)
        
    }
    
    l = sapply(test.keepX,length)
    n = names(test.keepX)
    temp = data.frame(l, n)
    
    
    message(paste("You have provided a sequence of keepX of length: ", paste(apply(temp, 1, function(x) paste(x,collapse=" for block ")), collapse= " and "), ".\nThis results in ",prod(sapply(test.keepX,length)), " models being fitted for each component and each nrepeat, this may take some time to run, be patient!",sep=""))
    
    if(missing(cpus))
    {
        parallel = FALSE
        message(paste("You can look into the 'cpus' argument to speed up computation time.",sep=""))

    } else {
        parallel = TRUE
        if(progressBar == TRUE)
        message(paste("As code is running in parallel, the progressBar will only show 100% upon completion of each component.",sep=""))

    }
    
    if(weighted == TRUE)
    {
        perfo = paste0("WeightedVote.error.rate.",dist)
    } else {
        perfo = paste0("MajorityVote.error.rate.",dist)
    }
    #-- end checking --#
    #------------------#
    
    test.keepX = lapply(test.keepX,sort) #sort test.keepX so as to be sure to chose the smallest in case of several minimum
    
    grid = expand.grid (test.keepX[length(test.keepX):1])[length(test.keepX):1] # each row is to be tested, the reordering is just a personal preference, works without it
    
    # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
    if ((!is.null(already.tested.X)) & length(already.tested.X) > 0)
    {
        comp.real = (length(already.tested.X[[1]]) + 1):ncomp
        #check and match already.tested.X to X
        if(constraint == TRUE & length(already.tested.X[[1]]) >0)
        {
            already.tested.X = get.keepA.and.keepA.constraint (X = X, keepX.constraint = already.tested.X, ncomp = rep(length(already.tested.X[[1]]), length(X)))$keepA.constraint
            
            # to get characters in already.tested.X
            already.tested.X = lapply(1:length(already.tested.X), function(y){temp = lapply(1:length(already.tested.X[[y]]), function(x){colnames(X[[y]])[already.tested.X[[y]][[x]]]}); names(temp) = paste("comp", 1:length(already.tested.X[[1]]), sep=""); temp})
            names(already.tested.X) = names(X)
            
        } else if(length(already.tested.X[[1]]) >0)
        {
            if(length(unique(names(already.tested.X)))!=length(already.tested.X) | sum(is.na(match(names(already.tested.X),names(X)))) > 0)
            stop("Each entry of 'already.tested.X' must have a unique name corresponding to a block of 'X'")
            
        }
        
    } else {
        comp.real = 1:ncomp
    }
    
    
    if (parallel == TRUE)
    cl <- makeCluster(cpus, type = "SOCK")
    
    
    fonction.indice.grid = function(indice.grid, mode, design, scheme, bias, init, verbose){
        test.keepX.comp = grid[indice.grid,]
        
        if(constraint == FALSE)
        {
            keepX.temp = lapply(1:length(X), function(x){c(already.tested.X[[x]],test.keepX.comp[[x]])})
            names(keepX.temp) = names(X)
        }
        
        # run block.splsda
        model = suppressMessages(block.splsda(X = X, Y = Y, ncomp=comp.real[comp],
        keepX.constraint = if(constraint){already.tested.X}else{NULL},
        keepX = if(constraint){test.keepX.comp}else{keepX.temp},
        design=design, scheme=scheme, mode=mode, scale=scale,
        bias=bias, init=init, tol=tol, verbose=verbose, max.iter=max.iter, near.zero.var=near.zero.var))
        
        
        # run perf on the model
        cvPerf = lapply(1 : nrepeat, function(u){out = suppressMessages(perf(model, validation = validation, folds = folds, dist = dist));
            if (progressBar ==  TRUE)
            setTxtProgressBar(pb, ((indice.grid-1)*nrepeat+u)/(nrow(grid)*nrepeat))
            out
        })
        names(cvPerf) = paste("nrepeat",1:nrepeat,sep="")
        
        # record results
        ## Majority Vote
        if(weighted == TRUE)
        {
            cvPerf2 = lapply(1 : nrepeat, function(x){unlist(cvPerf[[x]][names(cvPerf[[x]]) == "WeightedVote.error.rate"], recursive = FALSE)})
            names(cvPerf2) = paste("nrepeat",1:nrepeat,sep="")
        } else {
            cvPerf2 = lapply(1 : nrepeat, function(x){unlist(cvPerf[[x]][names(cvPerf[[x]]) == "MajorityVote.error.rate"], recursive = FALSE)})
            names(cvPerf2) = paste("nrepeat",1:nrepeat,sep="")
        }


    
        setTxtProgressBar(pb, indice.grid/nrow(grid))
        
        return(cvPerf2)
    }
    
    
    keepX = error.rate = mat.sd.error = NULL
    mat.error.rate = list()
    error.per.class.keepX.opt=list()
    # successively tune the components until ncomp: comp1, then comp2, ...
    for(comp in 1:length(comp.real))
    {
        if (progressBar == TRUE)
        cat("\ncomp",comp.real[comp], "\n")
        
        #-- set up a progress bar --#
        if (progressBar ==  TRUE & comp == 1)
        {
            pb = txtProgressBar(style = 3)
            nBar = 1
        }
        
        result.comp = matrix(nrow=nrow(grid),ncol=1)
        error.mean = error.sd = mat.error.rate.keepX = NULL
        error.per.class = array(0,c(nlevels(Y),nrepeat,nrow(grid)),dimnames = list(levels(Y), paste("nrep",1:nrepeat,sep=".")))
        

        if (parallel == TRUE)
        {
            clusterExport(cl, c("block.splsda","perf"))
            cvPerf3 = parLapply(cl, 1: nrow(grid), fonction.indice.grid)
        } else {
            cvPerf3 = lapply(1: nrow(grid), fonction.indice.grid)
            
        }
        names(cvPerf3) = paste("indice.grid",1:nrow(grid),sep="")
        

        
        
        for(indice.grid in 1 : nrow(grid))
        {
            # test.keepX.comp: keepX for each block on component "comp.real[comp]"
            # already.tested.X: either keepX (constraint=FALSE) or keepX.constraint.temp (constraint=TRUE) for all block on components 1:(comp.real[comp]-1)
            # keepX.temp: keepX for all block on all component 1:comp.real[comp], only used if constraint=FALSE
            mat.error.rate.temp = simplify2array( lapply(cvPerf3[[indice.grid]], function(x) x[[perfo]][measure,comp])) # error over the nrepeat
            mat.error.rate.keepX = rbind(mat.error.rate.keepX, mat.error.rate.temp)
            
            error.mean = c(error.mean, mean(mat.error.rate.temp))
            error.per.class[, , indice.grid] = simplify2array(lapply(cvPerf3[[indice.grid]], function(x){ x[[perfo]][1:nlevels(Y),comp]}))
            
            if(nrepeat > 1)
            error.sd = c(error.sd, apply(simplify2array(lapply(cvPerf3[[indice.grid]], function(x) x[[perfo]])), c(1,2), sd)[measure,comp])
            
            #if (progressBar ==  TRUE)
            #setTxtProgressBar(pb, (indice.grid)/nrow(grid))
            
        }
        
        names(error.mean) = apply(grid,1,function(x){paste(x, collapse = "_")})
        if(nrepeat > 1)
        names(error.sd) = names(error.mean)
        
        min.error = min(error.mean)
        min.keepX = names(which(error.mean == min.error)) # vector of all keepX combination that gives the minimum error
        
        a = strsplit(min.keepX,"_") # list of all keepX combination
        a = lapply(a, as.numeric) # transform characters of keepX into numbers
        
        #transform keepX in percentage of variable per dataset, so we choose the minimal overall
        p = sapply(X,ncol)
        percent = sapply(a, function(x) sum(x/p))
        ind.opt = which.min(percent) # we take only one
        a = a[[ind.opt]]# vector of each optimal keepX for all block on component comp.real[comp]
        
        # best keepX
        opt.keepX.comp = as.list(a)
        names(opt.keepX.comp) = names(X)
        
        error.per.class.keepX.opt[[comp]] = error.per.class[, , which.min(error.mean)] # error for the optimal keepXs
        error.rate = cbind(error.rate, error.mean)
        mat.sd.error = cbind(mat.sd.error, error.sd)
        mat.error.rate [[comp]] = mat.error.rate.keepX
        
        if(!constraint)
        {
            # add the optimal keepX to already.tested.X
            already.tested.X = lapply(1:length(X), function(x){c(already.tested.X[[x]],opt.keepX.comp[[x]])})
            
        } else {
            # get the variables selected by the optimal keepX, and add them in already.tested.X
            fit = suppressMessages(block.splsda(X = X, Y = Y, ncomp=comp.real[comp],
            keepX.constraint = already.tested.X,
            keepX = opt.keepX.comp,
            design=design, scheme=scheme, mode=mode, scale=scale,
            bias=bias, init=init, tol=tol, verbose=verbose, max.iter=max.iter, near.zero.var=near.zero.var))
            
            varselect = selectVar(fit, comp = comp.real[comp])
            varselect = varselect[which(names(varselect) %in% names(X))]
            
            if(length(already.tested.X) == 0)
            {
                already.tested.X = lapply(1:length(X), function(x){already.tested.X[[x]] = list()})
                already.tested.X = lapply(1:length(X),function(x){already.tested.X[[x]][[comp.real[comp]]] = list("comp1" = selectVar(fit, comp = 1)[[x]]$"name")})
                
            } else {
                already.tested.X = lapply(1:length(X),function(x){already.tested.X[[x]] = c(already.tested.X[[x]], list(selectVar(fit, comp = comp.real[comp])[[x]]$"name")); names(already.tested.X[[x]]) = paste("comp",1:comp.real[comp],sep=""); return(already.tested.X[[x]])})
            }
        }
        names(already.tested.X) = names(X)
        
        
        # prepping the results and save a file, if necessary
        if(!is.null(name.save))
        {
            colnames(error.rate) = paste("comp", comp.real[1:comp], sep='')
            names(mat.error.rate) = c(paste('comp', comp.real[1:comp], sep=''))
            mat.error.rate = lapply(mat.error.rate, function(x) {colnames(x) = paste("nrep",1:nrepeat,sep="."); rownames(x) = rownames(error.rate);x})
            names(error.per.class.keepX.opt) = c(paste('comp', comp.real[1:comp], sep=''))
            
            if(nrepeat > 1)
            colnames(mat.sd.error) = paste("comp", comp.real[1:comp], sep='')
            
            
            result = list(
            error.rate = error.rate,
            error.rate.sd = mat.sd.error,
            error.rate.all = mat.error.rate,
            choice.keepX = if(constraint){lapply(already.tested.X, function(x){sapply(x,length)})}else{already.tested.X},
            choice.keepX.constraint = if(constraint){already.tested.X}else{NULL},
            error.rate.class = error.per.class.keepX.opt)
            
            result$measure = measure.input
            result$call = match.call()
            
            class(result) = "tune.block.splsda"
        
            save(result, file = paste0(name.save,".comp",comp.real[1],"to",comp.real[comp],".Rdata"))
        }
        
        if (progressBar ==  TRUE)
        setTxtProgressBar(pb, 1)

    }
    #close the cluster after ncomp
    if (parallel == TRUE)
    stopCluster(cl)
    
    cat("\n")
    
    colnames(error.rate) = paste("comp", comp.real, sep='')
    names(mat.error.rate) = c(paste('comp', comp.real, sep=''))
    mat.error.rate = lapply(mat.error.rate, function(x) {colnames(x) = paste("nrep",1:nrepeat,sep="."); rownames(x) = rownames(error.rate);x})
    names(error.per.class.keepX.opt) = c(paste('comp', comp.real, sep=''))
    
    if(nrepeat > 1)
    colnames(mat.sd.error) = paste("comp", comp.real, sep='')
    
    
    result = list(
    error.rate = error.rate,
    error.rate.sd = mat.sd.error,
    error.rate.all = mat.error.rate,
    choice.keepX = if(constraint){lapply(already.tested.X, function(x){sapply(x,length)})}else{already.tested.X},
    choice.keepX.constraint = if(constraint){already.tested.X}else{NULL},
    error.rate.class = error.per.class.keepX.opt)
    
    result$measure = measure.input
    result$call = match.call()
    
    class(result) = "tune.block.splsda"
    
    return(result)
    
}


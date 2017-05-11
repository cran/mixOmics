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


tune.splsda = function (X, Y,
ncomp = 1,
test.keepX = c(5, 10, 15),
already.tested.X,
constraint = FALSE, #if TRUE, expect a list in already.tested.X, otherwise a number(keepX)
validation = "Mfold",
folds = 10,
dist = "max.dist",
measure = "BER", # one of c("overall","BER")
auc = FALSE,
progressBar = TRUE,
max.iter = 100,
near.zero.var = FALSE,
nrepeat = 1,
logratio = c('none','CLR'),
multilevel = NULL,
light.output = TRUE, # if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
cpus
)
{    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #------------------#
    #-- check entries --#
    X = as.matrix(X)
    
    if (length(dim(X)) != 2 || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
    
    
    # Testing the input Y
    if (is.null(multilevel))
    {
        if (is.null(Y))
        stop("'Y' has to be something else than NULL.")
        
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
        stop("'Y' should be a factor with more than one level")
        
    } else {
        # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
        multilevel = data.frame(multilevel)
        
        if ((nrow(X) != nrow(multilevel)))
        stop("unequal number of rows in 'X' and 'multilevel'.")
        
        if (ncol(multilevel) != 1)
        stop("'multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.")
        
        if (!is.null(ncol(Y)) && !ncol(Y) %in% c(0,1,2))# multilevel 1 or 2 factors
        stop("'Y' should either be a factor, a single column data.frame containing a factor, or a 2-columns data.frame containing 2 factors.")
        
        multilevel = data.frame(multilevel, Y)
        multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements
        
    }
    
    
    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    
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
    
    #-- logratio
    if (length(logratio) > 1)
    logratio = logratio[1]
    
    if (!logratio %in% c("none","CLR"))
    stop("'logratio must be one of 'none' or 'CLR'")
    
    #if ((!is.null(already.tested.X)) && (length(already.tested.X) != (ncomp - 1)) )
    #stop("The number of already tested parameters should be NULL or ", ncomp - 1, " since you set ncomp = ", ncomp)
    
    if (missing(already.tested.X))
    {
        if(constraint == TRUE)
        {
            already.tested.X = list()
        } else {
            already.tested.X = NULL
        }
    } else {
        if(is.null(already.tested.X) | length(already.tested.X)==0)
        stop("''already.tested.X' must be a vector of keepX values (if 'constraint'= FALSE) or a list (if'constraint'= TRUE) ")

        if(constraint == TRUE)
        {
            if(!is.list(already.tested.X))
            stop("''already.tested.X' must be a list since 'constraint' is set to TRUE")
            
            message(paste("A total of",paste(lapply(already.tested.X,length),collapse=" and "),"specific variables ('already.tested.X') were selected on the first ", length(already.tested.X), "component(s)"))
        } else {
            if(is.list(already.tested.X))
            stop("''already.tested.X' must be a vector of keepX values since 'constraint' is set to FALSE")

            message(paste("Number of variables selected on the first", length(already.tested.X), "component(s):", paste(already.tested.X,collapse = " ")))
        }
    }
    if(length(already.tested.X) >= ncomp)
    stop("'ncomp' needs to be higher than the number of components already tuned, which is length(already.tested.X)=",length(already.tested.X) , call. = FALSE)
    
    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- test.keepX
    if (is.null(test.keepX) | length(test.keepX) == 1 | !is.numeric(test.keepX))
    stop("'test.keepX' must be a numeric vector with more than two entries", call. = FALSE)
    
    if(!missing(cpus))
    {
        if(!is.numeric(cpus) | length(cpus)!=1)
        stop("'cpus' must be a numerical value")
        
        parallel = TRUE
        cl = makeCluster(cpus, type = "SOCK")
        clusterExport(cl, c("splsda","selectVar"))
        
        if(progressBar == TRUE)
        message(paste("As code is running in parallel, the progressBar will only show 100% upon completion of each component.",sep=""))

    } else {
        parallel = FALSE
        cl = NULL
    }
    #-- end checking --#
    #------------------#
    
   
    #---------------------------------------------------------------------------#
    #-- logration + multilevel approach ----------------------------------------#
    # we can do logratio and multilevel on the whole data as these transformation are done per sample
    X = logratio.transfo(X = X, logratio = logratio)
    
    if (!is.null(multilevel) & logratio == "none") # if no logratio, we can do multilevel on the whole data; otherwise it needs to be done after each logratio inside the CV
    {

        Xw = withinVariation(X, design = multilevel)
        X = Xw

        #-- Need to set Y variable for 1 or 2 factors
        Y = multilevel[, -1, drop=FALSE]
        if (ncol(Y) >= 1)
        Y = apply(Y, 1, paste, collapse = ".")  #  paste is to combine in the case we have 2 levels
        
        Y = as.factor(Y)
    }
    #-- multilevel approach ----------------------------------------------------#
    #---------------------------------------------------------------------------#


    #-- cross-validation approach  ---------------------------------------------#
    #---------------------------------------------------------------------------#
    
    test.keepX = sort(test.keepX) #sort test.keepX so as to be sure to chose the smallest in case of several minimum
    names(test.keepX) = test.keepX
    # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
    if ((!is.null(already.tested.X)))
    {
        comp.real = (length(already.tested.X) + 1):ncomp
        #check and match already.tested.X to X
        if(constraint == TRUE)
        {
            already.tested.X = get.keepA.and.keepA.constraint (X = list(X=X), keepX.constraint = list(X=already.tested.X), ncomp = length(already.tested.X))$keepA.constraint$X
            #transform already.tested.X to characters
            already.tested.X = relist(colnames(X), skeleton = already.tested.X)
        }

    } else {
        comp.real = 1:ncomp
    }
    
    
    choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist")
    dist = match.arg(dist, choices, several.ok = TRUE)
    
    mat.error = matrix(nrow = length(test.keepX), ncol = nrepeat,
    dimnames = list(test.keepX,c(paste('repeat', 1:nrepeat))))
    rownames(mat.error) = test.keepX
    
    mat.error.rate = list()
    error.per.class = list()

    mat.sd.error = matrix(0,nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real, sep=''))))
    mat.mean.error = matrix(nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real, sep=''))))

    error.per.class.mean = matrix(nrow = nlevels(Y), ncol = ncomp-length(already.tested.X),
        dimnames = list(c(levels(Y)), c(paste('comp', comp.real, sep=''))))
    error.per.class.sd = matrix(0,nrow = nlevels(Y), ncol = ncomp-length(already.tested.X),
        dimnames = list(c(levels(Y)), c(paste('comp', comp.real, sep=''))))
        
   
    # first: near zero var on the whole data set
    if(near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
            X = X[, -nzv$Position, drop=TRUE]
            
            if(ncol(X)==0)
            stop("No more predictors after Near Zero Var has been applied!")
            
        }
    }

    if (is.null(multilevel) | (!is.null(multilevel) && ncol(multilevel) == 2))
    {
        if(light.output == FALSE)
        prediction.all = class.all = list()
        
        if(auc)
        {
            auc.mean.sd=list()
            if(light.output == FALSE)
            auc.all=list()
        }
        
        error.per.class.keepX.opt = list()
        error.per.class.keepX.opt.mean = matrix(0, nrow = nlevels(Y), ncol = length(comp.real),
        dimnames = list(c(levels(Y)), c(paste('comp', comp.real, sep=''))))
        # successively tune the components until ncomp: comp1, then comp2, ...
        for(comp in 1:length(comp.real))
        {

            if (progressBar == TRUE)
            cat("\ncomp",comp.real[comp], "\n")
            
            result = MCVfold.splsda (X, Y, multilevel = multilevel, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = 1 + length(already.tested.X),
            choice.keepX = if(constraint){NULL}else{already.tested.X},
            choice.keepX.constraint = if(constraint){already.tested.X}else{NULL},
            test.keepX = test.keepX, measure = measure, dist = dist,
            near.zero.var = near.zero.var, progressBar = progressBar, class.object = "splsda", max.iter = max.iter, auc = auc, cl = cl)
            
            # in the following, there is [[1]] because 'tune' is working with only 1 distance and 'MCVfold.splsda' can work with multiple distances
            mat.error.rate[[comp]] = result[[measure]]$mat.error.rate[[1]]
            mat.mean.error[, comp]=result[[measure]]$error.rate.mean[[1]]
            if (!is.null(result[[measure]]$error.rate.sd[[1]]))
            mat.sd.error[, comp]=result[[measure]]$error.rate.sd[[1]]
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt[[comp]]=result[[measure]]$confusion[[1]]
            error.per.class.keepX.opt.mean[, comp]=apply(result[[measure]]$confusion[[1]], 1, mean)


            # best keepX
            if(!constraint)
            {
                already.tested.X = c(already.tested.X, result[[measure]]$keepX.opt[[1]])
            } else {
                fit = mixOmics::splsda(X, Y, ncomp = 1 + length(already.tested.X),
                keepX.constraint = already.tested.X, keepX = result[[measure]]$keepX.opt[[1]], near.zero.var = near.zero.var, mode = "regression")
                
                already.tested.X[[paste("comp",comp.real[comp],sep="")]] = selectVar(fit, comp = 1 + length(already.tested.X))$name
            }
            
            if(light.output == FALSE)
            {
                #prediction of each samples for each fold and each repeat, on each comp
                class.all[[comp]] = result$class.comp[[1]]
                prediction.all[[comp]] = result$prediction.comp
            }
            
            if(auc)
            {
                auc.mean.sd[[comp]] = result$auc
                if(light.output == FALSE)
                auc.all[[comp]] = result$auc.all
            }
            
        } # end comp
        if (parallel == TRUE)
        stopCluster(cl)
        
        names(mat.error.rate) = c(paste('comp', comp.real, sep=''))
        names(error.per.class.keepX.opt) = c(paste('comp', comp.real, sep=''))
        names(already.tested.X) = c(paste('comp', 1:ncomp, sep=''))
        
        if (progressBar == TRUE)
        cat('\n')
        
        #save(list=ls(),file="temp.Rdata")
        # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
        if(nrepeat > 2 & length(comp.real) >1)
        {
            keepX = if(constraint){lapply(already.tested.X, length)}else{already.tested.X}
            error.keepX = NULL
            for(comp in 1:length(comp.real))
            {
                ind.row = match(keepX[[comp.real[comp]]],test.keepX)
                error.keepX = cbind(error.keepX, mat.error.rate[[comp]][ind.row,])
            }
            colnames(error.keepX) = c(paste('comp', comp.real, sep=''))
            
            max = length(comp.real) #number max of components included
            pval = NULL
            opt = 1 #initialise the first optimal number of genes
            alpha = 0.01
            for(j in 2:max)
            {
                pval[j] = NA
                temp = try(t.test(error.keepX[,opt],error.keepX[,j],alternative="greater")$p.value, silent=T) #t.test of "is adding X comp improves the overall results"
                if(any(class(temp) == "try-error") || is.na(temp)) # temp can be NaN when error.keepX is constant
                {
                    ncomp_opt = NULL
                    break
                } else {
                    pval[j] = temp #had to be two steps (temp then pval =temp) otherwise the class is lost
                }
                if( (pval[j]< (alpha)))
                opt=j #if the p-value is lower than 0.05, the optimal number of comp is updated
            }
            if(all(class(temp) != "try-error") & !is.na(temp))
            ncomp_opt = comp.real[opt]
        } else {
            ncomp_opt = error.keepX = NULL
        }
        
        
        result = list(
        error.rate = mat.mean.error,
        error.rate.sd = mat.sd.error,
        error.rate.all = mat.error.rate,
        choice.keepX = if(constraint){sapply(already.tested.X, length)}else{already.tested.X},
        choice.keepX.constraint = if(constraint){already.tested.X}else{NULL},
        choice.ncomp = list(ncomp = ncomp_opt, values = error.keepX),
        error.rate.class = error.per.class.keepX.opt.mean,
        error.rate.class.all = error.per.class.keepX.opt)
        
        if(light.output == FALSE)
        {
            names(class.all) = names(prediction.all) = c(paste('comp', comp.real, sep=''))
            result$predict = prediction.all
            result$class = class.all
        }
        if(auc)
        {
            names(auc.mean.sd) = c(paste('comp', comp.real, sep=''))
            result$auc = auc.mean.sd
            if(light.output == FALSE)
            {
                names(auc.all) = c(paste('comp', comp.real, sep=''))
                result$auc.all =auc.all
            }
        }
        result$measure = measure
        result$call = match.call()

        class(result) = "tune.splsda"

        return(result)
    } else {
        # if multilevel with 2 factors, we can not do as before because withinvariation depends on the factors, we maximase a correlation
        message("For a two-factor analysis, the tuning criterion is based on the maximisation of the correlation between the components on the whole data set")

        cor.value = vector(length = length(test.keepX))
        names(cor.value) = test.keepX

        for (i in 1:length(test.keepX))
        {
            spls.train = mixOmics::splsda(X, Y, ncomp = ncomp, keepX = c(already.tested.X, test.keepX[i]), logratio = logratio, near.zero.var = FALSE, mode = "regression")
            
            # Note: this is performed on the full data set
            # (could be done with resampling (bootstrap) (option 1) and/or prediction (option 2))
            cor.value[i] = cor(spls.train$variates$X[, ncomp], spls.train$variates$Y[, ncomp])
            #
        }
        return(list(cor.value = cor.value))

    }
}
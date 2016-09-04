#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2015
# last modified: 24-08-2016
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
# tune.splsda: chose the optimal number of parameters per component on a splsda method
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 1.
# test.keepX: grid of keepX among which to chose the optimal one
# already.tested.X: a vector giving keepX on the components that were already tuned
# validation: Mfold or loo cross validation
# folds: if validation = Mfold, how many folds?
# dist: distance to classify samples. see predict
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# progressBar: show progress,
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# nrepeat: number of replication of the Mfold process
# logratio = c('none','CLR'). see splsda
# verbose: if TRUE, shows component and nrepeat being tested.


get.confusion_matrix = function(Y.learn,Y.test,pred)
{
    ClassifResult = array(0,c(nlevels(factor(Y.learn)),nlevels(factor(Y.learn))))
    rownames(ClassifResult) = levels(factor(Y.learn))
    colnames(ClassifResult) = paste("predicted.as.",levels(factor(Y.learn)),sep = "")
    #--------record of the classification accuracy for each level of Y
    for(i in 1:nlevels(factor(Y.learn)))
    {
        ind.i = which(Y.test == levels(factor(Y.learn))[i])
        for(ij in 1:nlevels(factor(Y.learn)))
        {
            ClassifResult[i,ij] = sum(pred[ind.i] == levels(Y.learn)[ij])
            
        }
    }
    ClassifResult
}

get.BER = function(X)
{
    if(!is.numeric(X)| !is.matrix(X) | length(dim(X)) != 2 | nrow(X)!=ncol(X))
    stop("'X' must be a square numeric matrix")
    
    nlev = nrow(X)
    #calculation of the BER
    ClassifResult.temp = X
    diag(ClassifResult.temp) = 0
    BER = sum(apply(ClassifResult.temp,1,sum,na.rm = TRUE)/apply(X,1,sum,na.rm = TRUE),na.rm = TRUE)/nlev
    return(BER)
}


stratified.subsampling = function(Y, folds = 10)
{
    stop = 0
    for(i in 1:nlevels(Y))
    {
        ai=sample(which(Y==levels(Y)[i]),replace=FALSE) # random sampling of the samples from level i
        aai=suppressWarnings(split(ai,factor(1:min(folds,length(ai)))))                       # split of the samples in k-folds
        if(length(ai)<folds)                                                # if one level doesn't have at least k samples, the list is completed with "integer(0)"
        {
            for(j in (length(ai)+1):folds)
            aai[[j]]=integer(0)
            stop = stop +1
        }
        assign(paste("aa",i,sep="_"),sample(aai,replace=FALSE))         # the `sample(aai)' is to avoid the first group to have a lot more data than the rest
    }
    
    # combination of the different split aa_i into SAMPLE
    SAMPLE=list()
    for(j in 1:folds)
    {
        SAMPLE[[j]]=integer(0)
        for(i in 1:nlevels(Y))
        {
            SAMPLE[[j]]=c(SAMPLE[[j]],get(paste("aa",i,sep="_"))[[j]])
        }
    }# SAMPLE is a list of k splits
    
    ind0 = sapply(SAMPLE, length)
    if(any(ind0 == 0))
    {
        SAMPLE = SAMPLE [-which(ind0 == 0)]
        message("Because of a too high number of 'folds' required, ",length(which(ind0 == 0))," folds were randomly assigned no data: the number of 'folds' is reduced to ", length(SAMPLE))
    }
    
    return(list(SAMPLE = SAMPLE, stop = stop))
}


MCVfold.splsda = function(
X,
Y,
multilevel = NULL, # repeated measurement only
validation,
folds,
nrepeat = 1,
ncomp,
choice.keepX = NULL, #either choice.keepX or choice.keepX.constraint, not both
choice.keepX.constraint = NULL,
test.keepX, # can be either a vector of names (keepX.constraint) or a value(keepX). In case of a value, there needs to be names(test.keepX)
measure = c("overall"), # one of c("overall","BER")
dist = "max.dist",
auc = FALSE,
max.iter = 100,
near.zero.var = FALSE,
progressBar = TRUE,
class.object = NULL
)
{    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    #-- set up a progress bar --#
    if (progressBar ==  TRUE)
    {
        pb = txtProgressBar(style = 3)
        nBar = 1
    } else {
        pb = FALSE
    }
    
    M = length(folds)
    features = NULL
    auc.all = prediction.comp = class.comp = list()
    for(ijk in dist)
    class.comp[[ijk]] = array(0, c(nrow(X), nrepeat, length(test.keepX)))# prediction of all samples for each test.keepX and  nrep at comp fixed
    folds.input = folds
    for(nrep in 1:nrepeat)
    {
        prediction.comp[[nrep]] = array(0, c(nrow(X), nlevels(Y), length(test.keepX)), dimnames = list(rownames(X), levels(Y), names(test.keepX)))
        rownames(prediction.comp[[nrep]]) = rownames(X)
        colnames(prediction.comp[[nrep]]) = levels(Y)
        
        if(nlevels(Y)>2)
        {
            auc.all[[nrep]] = array(0, c(nlevels(Y),2, length(test.keepX)), dimnames = list(paste(levels(Y), "vs Other(s)"), c("AUC","p-value"), names(test.keepX)))
        }else{
            auc.all[[nrep]] = array(0, c(1,2, length(test.keepX)), dimnames = list(paste(levels(Y)[1], levels(Y)[2], sep = " vs "), c("AUC","p-value"), names(test.keepX)))
        }
        
        n = nrow(X)
        repeated.measure = 1:n
        if (!is.null(multilevel))
        {
            repeated.measure = multilevel[,1]
            n = length(unique(repeated.measure)) # unique observation: we put every observation of the same "sample" in the either the training or test set
        }
       
        
        #-- define the folds --#
        if (validation ==  "Mfold")
        {
            
            if (nrep > 1) # reinitialise the folds
            folds = folds.input
            
            if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
            {
                stop("Invalid number of folds.")
            } else {
                M = round(folds)
                if (is.null(multilevel))
                {
                    temp = stratified.subsampling(Y, folds = M)
                    folds = temp$SAMPLE
                    if(temp$stop > 0 & nrep == 1) # to show only once
                    warning("At least one class is not represented in one fold, which may unbalance the error rate.\n  Consider a number of folds lower than the minimum in table(Y): ", min(table(Y)))
                } else {
                    folds = split(sample(1:n), rep(1:M, length = n)) # needs to have all repeated samples in the same fold
                }
            }
        } else if (validation ==  "loo") {
            folds = split(1:n, rep(1:n, length = n))
            M = n
        }
        
        M = length(folds)
        
        error.sw = matrix(0,nrow = M, ncol = length(test.keepX))
        rownames(error.sw) = paste0("fold",1:M)
        colnames(error.sw) = names(test.keepX)
        # for the last keepX (i) tested, prediction combined for all M folds so as to extract the error rate per class
        # prediction.all = vector(length = nrow(X))
        # in case the test set only includes one sample, it is better to advise the user to
        # perform loocv
        stop.user = FALSE
        for (j in 1:M)
        {
            if (progressBar ==  TRUE)
            setTxtProgressBar(pb, (M*(nrep-1)+j-1)/(M*nrepeat))
            
            #print(j)
            #set up leave out samples.
            omit = which(repeated.measure %in% folds[[j]] == TRUE)

            # get training and test set
            X.train = X[-omit, ]
            Y.train = Y[-omit]
            X.test = X[omit, , drop = FALSE]#matrix(X[omit, ], nrow = length(omit)) #removed to keep the colnames in X.test
            Y.test = Y[omit]
   
            #---------------------------------------#
            #-- near.zero.var ----------------------#
            
            # first remove variables with no variance
            var.train = apply(X.train, 2, var)
            ind.var = which(var.train == 0)
            if (length(ind.var) > 0)
            {
                X.train = X.train[, -c(ind.var),drop = FALSE]
                X.test = X.test[, -c(ind.var),drop = FALSE]
            }
            
            
            if(near.zero.var == TRUE)
            {
                remove.zero = nearZeroVar(X.train)$Position
                
                if (length(remove.zero) > 0)
                {
                    X.train = X.train[, -c(remove.zero),drop = FALSE]
                    X.test = X.test[, -c(remove.zero),drop = FALSE]
                }
            }
            
            #-- near.zero.var ----------------------#
            #---------------------------------------#
            for (i in 1:length(test.keepX))
            {
                if (progressBar ==  TRUE)
                setTxtProgressBar(pb, (M*(nrep-1)+j-1)/(M*nrepeat) + (i-1)/length(test.keepX)/(M*nrepeat))
                
                # depending on whether it is a constraint and whether it is from tune or perf, keepX and keepX.constraint differ:
                # if it's from perf, then it's only either keepX or keepX.constraint
                # if it's from tune, then it's either keepX, or a combination of keepX.constraint and keepX
                # we know if it's perf+constraint or tune+constraint depending on the test.keepX that is either a vector or a list
                object.res = splsda(X.train, Y.train, ncomp = ncomp,
                keepX = if(is.null(choice.keepX.constraint) & !is.list(test.keepX)){c(choice.keepX, test.keepX[i])}else if(!is.list(test.keepX)){test.keepX[i]} else {NULL} ,
                keepX.constraint = if(is.null(choice.keepX.constraint)& !is.list(test.keepX)){NULL}else if(!is.list(test.keepX)){choice.keepX.constraint} else {c(choice.keepX.constraint, test.keepX)},
                logratio = "none", near.zero.var = FALSE, mode = "regression", max.iter = max.iter)
                  
                # added: record selected features
                if (any(class.object %in% c("splsda")) & length(test.keepX) ==  1) # only done if splsda and if only one test.keepX as not used if more so far
                # note: if plsda, 'features' includes everything: to optimise computational time, we don't evaluate for plsda object
                features = c(features, selectVar(object.res, comp = ncomp)$name)
                
                test.predict.sw <- predict(object.res, newdata = X.test, method = dist)
                prediction.comp[[nrep]][omit, , i] =  test.predict.sw$predict[, , ncomp]
                
                for(ijk in dist)
                class.comp[[ijk]][omit,nrep,i] =  test.predict.sw$class[[ijk]][, ncomp] #levels(Y)[test.predict.sw$class[[ijk]][, ncomp]]
            } # end i
            
        } # end j 1:M (M folds)
        
        if (progressBar ==  TRUE)
        setTxtProgressBar(pb, (M*nrep)/(M*nrepeat))
        
        if(auc)
        {
            data=list()
            for (i in 1:length(test.keepX))
            {
                data$outcome = Y
                data$data = prediction.comp[[nrep]][, , i]
                auc.all[[nrep]][, , i] = as.matrix(statauc(data))
            }
        }

    } #end nrep 1:nrepeat
    names(prediction.comp) = names (auc.all) = paste0("nrep.", 1:nrepeat)
    # class.comp[[ijk]] is a matrix containing all prediction for test.keepX, all nrepeat and all distance, at comp fixed
    
    save(list=ls(),file="temp.Rdata")
    # average auc over the nrepeat, for each test.keepX
    if(auc)
    {
        
        if(nlevels(Y)>2)
        {
            auc.mean.sd =  array(0, c(nlevels(Y),2, length(test.keepX)), dimnames = list(rownames(auc.all[[1]]), c("AUC.mean","AUC.sd"), names(test.keepX)))
        }else{
            auc.mean.sd =  array(0, c(1,2, length(test.keepX)), dimnames = list(rownames(auc.all[[1]]), c("AUC.mean","AUC.sd"), names(test.keepX)))
        }
        
        for(i in 1:length(test.keepX))
        {
            temp = NULL
            for(nrep in 1:nrepeat)
            {
                temp = cbind(temp, auc.all[[nrep]][, 1, i])
            }
            auc.mean.sd[, 1, i] = apply(temp,1,mean)
            auc.mean.sd[, 2, i] = apply(temp,1,sd)
        }
    } else {
        auc.mean.sd = auc.all = NULL
        
    }
    
    result = list()
    error.mean = error.sd = error.per.class.keepX.opt.comp = keepX.opt = test.keepX.out = mat.error.final = choice.keepX.out = list()
    
    if (any(measure == "overall"))
    {
        for(ijk in dist)
        {
            rownames(class.comp[[ijk]]) = rownames(X)
            colnames(class.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            dimnames(class.comp[[ijk]])[[3]] = paste0("test.keepX.",names(test.keepX))
            
            #finding the best keepX depending on the error measure: overall or BER
            # classification error for each nrep and each test.keepX: summing over all samples
            error = apply(class.comp[[ijk]],c(3,2),function(x)
            {
                sum(as.character(Y) != x)
            })
            rownames(error) = names(test.keepX)
            colnames(error) = paste0("nrep.",1:nrepeat)
            
            # we want to average the error per keepX over nrepeat and choose the minimum error
            error.mean[[ijk]] = apply(error,1,mean)/length(Y)
            if (!nrepeat ==  1)
            error.sd[[ijk]] = apply(error,1,sd)/length(Y)
            
            mat.error.final[[ijk]] = error/length(Y)  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
            
            keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt.comp[[ijk]] = apply(class.comp[[ijk]][, , keepX.opt[[ijk]], drop = FALSE], 2, function(x)
            {
                conf = get.confusion_matrix(Y.learn = factor(Y), Y.test = factor(Y), pred = x)
                out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
            })
            
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            colnames(error.per.class.keepX.opt.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            
            
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            if(is.null(choice.keepX))
            {
                choice.keepX.out[[ijk]] = c(lapply(choice.keepX.constraint,length), test.keepX.out)
            }else{
                choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
            }
            result$"overall"$error.rate.mean = error.mean
            if (!nrepeat ==  1)
            result$"overall"$error.rate.sd = error.sd
            
            result$"overall"$confusion = error.per.class.keepX.opt.comp
            result$"overall"$mat.error.rate = mat.error.final
            result$"overall"$keepX.opt = test.keepX.out
        }
    }
    
    if (any(measure ==  "BER"))
    {
        for(ijk in dist)
        {
            rownames(class.comp[[ijk]]) = rownames(X)
            colnames(class.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            dimnames(class.comp[[ijk]])[[3]] = paste0("test.keepX.",names(test.keepX))
            
            error = apply(class.comp[[ijk]],c(3,2),function(x)
            {
                conf = get.confusion_matrix(Y.learn = factor(Y),Y.test = factor(Y),pred = x)
                get.BER(conf)
            })
            rownames(error) = names(test.keepX)
            colnames(error) = paste0("nrep.",1:nrepeat)
            
            # average BER over the nrepeat
            error.mean[[ijk]] = apply(error,1,mean)
            if (!nrepeat ==  1)
            error.sd[[ijk]] = apply(error,1,sd)
            
            mat.error.final[[ijk]] = error  # BER for each test.keepX (rows) and each nrepeat (columns)
            
            keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1]
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt.comp[[ijk]] = apply(class.comp[[ijk]][, , keepX.opt[[ijk]], drop = FALSE], 2, function(x)
            {
                conf = get.confusion_matrix(Y.learn = factor(Y), Y.test = factor(Y), pred = x)
                out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
            })
            
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            colnames(error.per.class.keepX.opt.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            if(is.null(choice.keepX))
            {
                choice.keepX.out[[ijk]] = c(lapply(choice.keepX.constraint,length), test.keepX.out)
            }else{
                choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
            }
            result$"BER"$error.rate.mean = error.mean
            if (!nrepeat ==  1)
            result$"BER"$error.rate.sd = error.sd
            
            result$"BER"$confusion = error.per.class.keepX.opt.comp
            result$"BER"$mat.error.rate = mat.error.final
            result$"BER"$keepX.opt = test.keepX.out
            
        }
        
        
    }
    
    
    result$prediction.comp = prediction.comp
    result$auc = auc.mean.sd
    result$auc.all = auc.all
    result$class.comp = class.comp
    result$features$stable = sort(table(as.factor(features))/M/nrepeat, decreasing = TRUE)
    return(result)
}




#############################################################################################################
# Authors:
#   Amrit Singh, University of British Columbia, Vancouver.
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 01-04-2015
# last modified: 27-05-2016
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


# ----------------------------------------------------------------------------------------------------------
# perf.sgccda - Function to evaluate the performance of the fitted PLS (cross-validation)
#   inputs: object - object obtain from running sgccda or splsda
#           dist - to evaluate the classification performance
#           validation - type of validation
#           folds - number of folds if validation = "Mfold"
# ----------------------------------------------------------------------------------------------------------

perf.sgccda = function (object,
dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
validation = c("Mfold", "loo"),
folds = 10,
...)
{
    
    ### Start: Initialization parameters
    X = object$X; level.Y = object$names$colnames$Y;
    J = length(X)
    Y = object$Y#ind.mat; Y = map(Y); Y = factor(Y, labels = level.Y)
    n = nrow(X[[1]]);
    indY = object$indY
    max.iter = object$max.iter

    if (any(dist == "all"))
    {
        dist.select = c("max.dist", "centroids.dist", "mahalanobis.dist")
    } else {
        dist.select = dist
    }
    ### End: Initialization parameters
    
    dist = match.arg(dist, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
    
    ### Start: Check parameter validation / set up sample
    
    
    if (length(validation) > 1 )
    validation = validation [1]
    
    if (!(validation %in% c("Mfold", "loo")))
    stop("Choose 'validation' among the two following possibilities: 'Mfold' or 'loo'")

    
    if (validation == "Mfold")
    {
        
        if (!(abs(folds - round(folds)) < .Machine$double.eps) || is.null(folds) || folds < 2 || folds > n)
        {
            stop(paste("Invalid number of folds.", "folds must be an integer contained between", 2, "and", n))
        } else {
            M = folds
            folds = split(sample(1:n), rep(1:M, length = n))
            folds = lapply(folds, sort)
        }
    } else if (validation == "loo") {
        M = n
        folds = split(1:n, rep(1:n, length = n))
    } else {
        stop("validation can be only 'Mfold' or 'loo'")
    }
    ### Start: Check parameter validation / set up sample
    
    ### Start: Training samples (X.training and Y.training) and Test samples (X.test / Y.test)
    X.training = lapply(folds, function(x){out=lapply(1:J, function(y) {X[[y]][-x, ]}); names(out) = names(X); out}) #need to name the block for prediction
    Y.training = lapply(folds, function(x) {Y[-x]});
    
    X.test = lapply(folds, function(x){out=lapply(1:J, function(y) {X[[y]][x, , drop = FALSE]}); names(out) = names(X); out})#need to name the block for prediction
    Y.test = lapply(folds, function(x) {Y[x]})
    ### End: Training samples (X.training and Y.training) and Test samples (X.test / Y.test)
    
    ### Estimation models
    model = lapply(1 : M, function(x) {block.splsda(X = X.training[[x]], Y = Y.training[[x]], ncomp = max(object$ncomp[-indY]), keepX = object$keepX,
            design = object$design, max.iter = object$max.iter, tol = object$tol, init = object$init, scheme = object$scheme,
            bias = object$bias, mode = object$mode)})
    
    ### Retrieve selected variables per component
    features = lapply(1 : J, function(x)
    {
        lapply(1 : object$ncomp[x], ### List level / ncomp level
        function(y)
        {
            unlist(lapply(1 : M,  ### Validation level
            function(z)
            {
                if (is.null(colnames(X[[x]])))
                {
                    paste0("X", which(model[[z]]$loadings[[x]][, y] != 0))
                } else {
                    if (!is.null(colnames(X[[x]])))
                    colnames(X[[x]])[model[[z]]$loadings[[x]][, y] != 0]
                }
            }
            ))
        })
    })
    
    ### Start: Analysis feature selection
    # Statistics: stability
    list.features = lapply(1 : J, function(x){lapply(features[[x]], function(y){(sort(table(factor(y))/M, decreasing = TRUE))})})
    
    # Statistics: original model (object)
    final.features = lapply(1 : J, function(x)
    {
        lapply(1 : object$ncomp[x],function(y)
        {
            temp = as.data.frame(object$loadings[[x]][which(object$loadings[[x]][, y] != 0), y, drop = FALSE])
            if (is.null(colnames(X[[x]])))
            {
                row.names(temp) = paste0("X", which(object$loadings[[x]][, y] != 0))
            } else {
                if (!is.null(colnames(X[[x]])))
                row.names(temp) = colnames(X[[x]])[which(object$loadings[[x]][, y] != 0)]
            }
            names(temp) = "value.var"
            return(temp[sort(abs(temp[, 1]), index.return = TRUE, decreasing = TRUE)$ix, 1, drop = FALSE])
        })
    })
    
    list.features = lapply(1 : J, function(x){ names(list.features[[x]]) = paste("comp", 1 : object$ncomp[x])
        return(list.features[[x]])})
    
    final.features = lapply(1 : J, function(x){ names(final.features[[x]]) = paste("comp", 1 : object$ncomp[x])
        return(final.features[[x]])})
    ### End: Analysis feature selection
    
    ### Warning: no near.zero.var applies with sgcca
    
    ### Start: Prediction (score / class) sample test
    # Prediction model on test dataset
    Y.all = lapply(1 : M, function(x) {predict(model[[x]], X.test[[x]], method = "all")})
    
    # Retrieve class prediction
    Y.predict = lapply(1 : M, function(x) {Y.all[[x]]$class})
    
    ## Start: retrieve score for each component
    # Keep score values
    Y.all = lapply(1 : M, function(x) {Y.all[[x]]$predict})
    
    # Reorganization list Y.all data / ncomp / folds
    Y.all = lapply(1 : J, function(x)
    {
        lapply(1 : object$ncomp[x], function(y)
        {
            lapply(1 : M, function(z)
            {
                Y.all[[z]][[x]][, , y]
            })
        })
    })
    # Merge score
    Y.all = lapply(1 : J, function(x)
    {
        lapply(1 : object$ncomp[x], function(y)
        {
            do.call(rbind, Y.all[[x]][[y]])
        })
    })
    
    # Define row.names
    Y.all = lapply(1 : J, function(x)
    {
        lapply(1 : object$ncomp[x], function(y)
        {
            row.names(Y.all[[x]][[y]]) = row.names(X[[x]])[unlist(folds)]
            return(Y.all[[x]][[y]])
        })
    })
    
    # Define colnames
    Y.all = lapply(1 : J, function(x)
    {
        lapply(1 : object$ncomp[x], function(y)
        {
            colnames(Y.all[[x]][[y]]) = levels(Y)
            return(Y.all[[x]][[y]])
        })
        names(Y.all[[x]])=paste("comp", 1:object$ncomp[x])
        return(Y.all[[x]])
    })
    ## End: retrieve score for each component
    
    ## Start: retrieve class for each component
    # Reorganization input data / folds / dist.select
    Y.predict = lapply(1 : J, function(x)
    {
        lapply(1 : M, function(y)
        {
            lapply(which(c("max.dist", "centroids.dist", "mahalanobis.dist") %in% dist.select), function(z)
            {
                Y.predict[[y]][[z]][[x]]
            })
        })
    })
    
    # Define row.names
    Y.predict = lapply(1 : J, function(x)
    {
        lapply(1 : M, function(y)
        {
            lapply(1 : length(dist.select), function(z)
            {
                if (!is.null(row.names(X[[x]])[folds[[y]]]))
                {
                    row.names(Y.predict[[x]][[y]][[z]]) = row.names(X[[x]])[folds[[y]]]
                } else {
                    row.names(Y.predict[[x]][[y]][[z]]) = paste("Ind", unlist(folds))
                }
                return(Y.predict[[x]][[y]][[z]])
            })
        })
    })
    
    # Reorgainzation list input / dist.select / folds
    Y.predict = lapply(1 : J, function(x)
    {
        lapply(1 : length(dist.select), function(y)
        {
            lapply(1 : M, function(z)
            {
                Y.predict[[x]][[z]][[y]]
            })
        })
    })
    # Merge Score
    Y.predict = lapply(1 : J, function(x)
    {
        lapply(1 : length(dist.select), function(y)
        {
            do.call(rbind, Y.predict[[x]][[y]])
        })
    })
    
    # Define names
    Y.predict = lapply(1 : J, function(x)
    {
        names(Y.predict[[x]]) = dist.select
        return(Y.predict[[x]])
    })
    ## End: retrieve class for each component
    ### End: Prediction (score / class) sample test
    
    ### Start: Estimation error rate
    ## Start: Estimation overall error rate
    #Statistics overall error rate
    error.mat = lapply(1 : J, function(x)
    {
        lapply(dist.select , function(y)
        {
            apply(Y.predict[[x]][[y]], 2, function(z)
            {
                1 - sum(diag(table(factor(z, levels = levels(Y)), Y[unlist(folds)])))/length(Y)
            })
        })
    })
    
    # Merge error rate according to dist
    error.mat = lapply(1 : J, function(x)
    {
        do.call(cbind, error.mat[[x]])
    })
    
    # Define name
    error.mat = lapply(1 : J, function(x)
    {
        colnames(error.mat[[x]]) = dist.select
        return(error.mat[[x]])
    })
    ## End: Estimation overall error rate
    
    ## Start: Estimation error rate per class
    # Statistics error rate per class
    error.mat.class = lapply(1 : J, function(x)
    {
        lapply(dist.select , function(y)
        {
            apply(Y.predict[[x]][[y]], 2, function(z)
            {
                temp = diag(table(factor(z, levels = levels(Y)), Y[unlist(folds)]))
                1 - c(temp/summary(Y))#, sum(temp)/length(Y))
            })
        })
    })
    
    if(FALSE)
    {
     # Transpose data
    error.mat.class = lapply(1 : J, function(x)
    {
        lapply(1 : length(dist.select), function(y)
        {
            t(error.mat.class[[x]][[y]])
        })
    })
    
    # Define colnames
    error.mat.class = lapply(1 : J, function(x)
    {
        lapply(1 : length(dist.select), function(y)
        {
            colnames(error.mat.class[[x]][[y]])[nlevels(Y) + 1] = "Overall"
            return(error.mat.class[[x]][[y]])
        })
    })
    }
    
    # Define names
    error.mat.class = lapply(1 : J, function(x)
    {
        names(error.mat.class[[x]]) = dist.select
        return(error.mat.class[[x]])
    })
    
    ## End: Estimation error rate per class
    ### End: Estimation error rate
    
    if (!is.null(names(X)))
    {
        names(error.mat) = names(X); names(error.mat.class) = names(X)
        names(list.features) = names(X); names(final.features) = names(X);
        names(Y.all) = names(X); names(Y.predict) = names(X)
    }
    
    ### Start: Supplementary analysis for sgcca
    if (length(X) > 1)
    {
        
        ### Start: Average prediction
        # if ncomp[[X]] < max(ncomp), copy the last prediction
        Y.mean = lapply(1 : max(object$ncomp[-indY]), function(x)
        {
            lapply(1 : J, function(y)
            {
                Y.all[[y]][[min(x, object$ncomp[y])]]
            })
        })
        
        # Sort matrix
        Y.mean = lapply(1 : max(object$ncomp[-indY]), function(x)
        {
            lapply(1 : J, function(y)
            {
                Y.mean[[x]][[y]][sort(unlist(folds), index.return = TRUE)$ix, , drop = FALSE]
            })
        })
        
        # Average score
        Y.mean = lapply(1 : max(object$ncomp[-indY]), function(x)
        {
            Reduce("+", Y.mean[[x]])/length(X)
        })
        
        # Determine index of the max
        Y.mean = sapply(1 : max(object$ncomp[-indY]), function(x)
        {
            apply(Y.mean[[x]], 1, function(y)
            {
                levels(Y)[which.max(y)]
            })
        })
        
        # Define colnames
        colnames(Y.mean) = paste("comp", 1 : max(object$ncomp[-(J + 1)]))
        
        # Estimation error.rate
        #Y.mean.res = sapply(1:max(object$ncomp[-(J + 1)]), function(x){temp = diag(table(factor(Y.mean[, x], levels = c(1:nlevels(Y))), Y))
        #                                                              c(temp/summary(Y), sum(temp)/length(Y))})
        Y.mean.res = sapply(1:max(object$ncomp[-indY]), function(x)
        {
            mat = table(factor(Y.mean[, x], levels = levels(Y)), Y)
            mat2 <- mat
            diag(mat2) <- 0
            err = c(c(colSums(mat2)/summary(Y), sum(mat2)/length(Y)), mean(colSums(mat2)/colSums(mat)))
        })
        
        #Y.mean.res = t(Y.mean.res)
        colnames(Y.mean.res) = paste("comp", 1:max(object$ncomp[-indY]))
        row.names(Y.mean.res) = c(levels(Y), "Overall.ER", "Overall.BER")
        ### End: Average prediction
        
        ### Start: Vote on the dataset
        # if ncomp[[X]] < max(ncomp), copy the last prediction
        Y.vote = lapply(1 : J, function(x)
        {
            lapply(dist.select, function(y)
            {
                if(ncol(Y.predict[[x]][[y]]) < max(object$ncomp[-indY]))
                {
                    Y.predict[[x]][[y]] = cbind(Y.predict[[x]][[y]], matrix(rep(Y.predict[[x]][[y]][, object$ncomp[x]], max(object$ncomp[-indY]) - ncol(Y.predict[[x]][[y]])), ncol = (max(object$ncomp[-indY]) - ncol(Y.predict[[x]][[y]])), nrow = nrow(X[[x]])))
                } else {
                    Y.predict[[x]][[y]]
                }
            })
        })
        
        # Sort matrix
        Y.vote = lapply(1 : J, function(x)
        {
            lapply(1 : length(dist.select), function(y)
            {
                Y.vote[[x]][[y]][sort(unlist(folds), index.return = TRUE)$ix, , drop = FALSE]
            })
        })
        
        # Reorganization dist.select / component / Data
        Y.vote = lapply(1 : length(dist.select), function(x)
        {
            lapply(1 : max(object$ncomp[-indY]), function(y){
                lapply(1 : J, function(z){
                    Y.vote[[z]][[x]][, y, drop = FALSE]
                })
            })
        })
        
        # Merge dist.select
        Y.vote = lapply(1 : length(dist.select), function(x)
        {
            lapply(1 : max(object$ncomp[-indY]), function(y){
                do.call(cbind, Y.vote[[x]][[y]])
            })
        })
        
        # Estimation Majority Vote
        Y.vote = lapply(1 : length(dist.select), function(x)
        {
            lapply(1 : max(object$ncomp[-indY]), function(y)
            {
                apply(Y.vote[[x]][[y]], 1, function(z)
                {
                    temp = table(z)
                    if (length(names(temp)[temp == max(temp)]) > 1)
                    {
                        NA
                    } else {
                        names(temp)[temp == max(temp)]
                    }
                })
            })
        })
        
        Y.vote = lapply(1 : length(dist.select), function(x)
        {
            do.call(cbind, Y.vote[[x]])
        })
        

        ## subjects with NA are considered false
        Y.vote.res = lapply(1 : length(dist.select), function(x)
        {
            apply(Y.vote[[x]], 2, function(y)
            {
                y[is.na(y)] <- nlevels(Y)+5   ## adding a new level for unsure subjects (replacing NA with this level)
                temp=table(factor(y, levels = c(levels(Y), nlevels(Y)+5)), Y)
                diag(temp) <- 0
                err = c(colSums(temp)/summary(Y), sum(temp)/length(Y), mean(colSums(temp)/summary(Y)))
                return(err=err)
            })
        })
        
        Y.vote = lapply(1 : length(dist.select), function(x)
        {
            colnames(Y.vote[[x]]) = paste("comp", 1:max(object$ncomp[-(J + 1)]))
            return(Y.vote[[x]])
        })
        
        Y.vote.res = lapply(1 : length(dist.select), function(x)
        {
            colnames(Y.vote.res[[x]]) = paste("comp", 1:max(object$ncomp[-(J + 1)]))
            row.names(Y.vote.res[[x]]) = c(levels(Y), "Overall.ER", "Overall.BER")
            return((Y.vote.res[[x]]))
        })
        names(Y.vote) = dist.select; names(Y.vote.res) = dist.select
        
        ### End: Vote on the dataset
    }
    ### End: Supplementary analysis for sgcca
    
    result = list()
    result$error.rate = error.mat
    result$error.rate.class = error.mat.class
    result$predict = Y.all
    result$class = Y.predict
    
    result$features$stable = list.features
    #result$features$final = final.features
    
    if (length(X) > 1)
    {
        result$AveragePredict.class = Y.mean
        result$AveragePredict.error.rate = Y.mean.res
        result$MajorityClass = Y.vote
        result$MajorityClass.error.rate = Y.vote.res
    }
    
    method = "sgccda.mthd"
    result$meth = "sgccda.mthd"
    class(result) = "perf.sgccda.mthd"
    result$call = match.call()

    return(invisible(result))
}
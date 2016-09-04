#############################################################################################################
# Authors:
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, University of Queensland Diamantina Institute, Brisbane, Australia
#
# created: 23-08-2016
# last modified:  23-08-2016
#
# Copyright (C) 2016
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

auroc = function(object, ...)
UseMethod("auroc")


# PLSDA object
# ----------------------
auroc.plsda = auroc.splsda = function(
object,
newdata = object$X,
outcome.test = as.factor(object$Y),
multilevel = NULL,
plot = TRUE,
roc.comp = 1,
...)
{
    if(dim(newdata)[[1]]!=length(outcome.test))
    stop("Factor outcome.test must be a factor with ",dim(newdata)[[1]], " elements.",call. = FALSE)
    
    data=list()
    statauc=list()
    data$outcome=factor(outcome.test)
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict  =  predict(object, newdata = newdata, dist = "max.dist", multilevel = multilevel)$predict
    
    for (i in 1:object$ncomp)
    {
        data$data=res.predict[,,i]
        title=paste("ROC Curve Comp",i)
        statauc[[paste("Comp", i, sep = "")]]=statauc(data, plot = ifelse(i%in%roc.comp,plot,FALSE), title = title)
    }
    return(statauc)
}


# MINT object
# ----------------------
auroc.mint.plsda = auroc.mint.splsda = function(
object,
newdata = object$X,
outcome.test = as.factor(object$Y),
study.test = object$study,
multilevel = NULL,
plot = TRUE,
roc.comp = 1,
...)
{
    
    if(dim(newdata)[[1]]!=length(outcome.test))
    stop("Factor outcome.test must be a factor with ",dim(newdata)[[1]], " elements.",call. = FALSE)
    
    if(dim(newdata)[[1]]!=length(study.test))
    stop("Factor study.test must be a factor with ",dim(newdata)[[1]], " elements.",call. = FALSE)
    study.test=factor(study.test)
    
    data=list()
    statauc=list()
    data$outcome=factor(outcome.test)
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict  =  predict(object, newdata = newdata, dist = "max.dist", multilevel = multilevel, study.test = study.test)$predict
    
    for (i in 1:object$ncomp)
    {
        data$data=res.predict[,,i]
        title=paste("ROC Curve Comp",i)
        statauc[[paste("Comp", i, sep = "")]]=statauc(data, plot = ifelse(i%in%roc.comp,plot,FALSE), title = title)
    }
    return(statauc)
}


# block.splsda object
# ----------------------
auroc.sgccda = function(
object,
newdata = object$X,
outcome.test = as.factor(object$Y),
multilevel = NULL,
plot = TRUE,
roc.block = 1,
roc.comp = 1,
...)
{
    
    data=list()
    auc.mean=list()
    data$outcome=factor(outcome.test)
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict  =  predict(object, newdata = newdata, dist = "max.dist", multilevel = multilevel)$predict
    block.all = names(res.predict)
    block.temp = names(res.predict[roc.block])
    
    for(j in 1:length(res.predict))
    {
        for (i in 1:object$ncomp[j])
        {
            data$data=res.predict[[j]][,,i]
            title=paste("ROC Curve\nBlock: ", names(res.predict)[j], ", comp: ",i, sep="")
            
            plot.temp = ifelse(i%in%roc.comp && names(res.predict)[j]%in%block.temp, plot, FALSE)
            auc.mean[[names(res.predict)[j]]][[paste("comp",i,sep = "")]] = statauc(data, plot = plot.temp, title = title)
            
        }
    }
    return(auc.mean)
}

# mint.block.splsda object
# ----------------------
auroc.mint.block.splsda=auroc.mint.block.plsda = function(
object,
newdata = object$X,

study.test = object$study,
outcome.test = as.factor(object$Y),
multilevel = NULL,
plot = TRUE,
roc.block = 1,
roc.comp = 1,
...)
{
    
    data=list()
    auc.mean=list()
    data$outcome=factor(outcome.test)
    study.test=factor(study.test)
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict  =  predict(object, newdata = newdata, study.test=study.test,dist = "max.dist", multilevel = multilevel)$predict
    block.all = names(res.predict)
    block.temp = names(res.predict[roc.block])
    
    for(j in 1:length(res.predict))
    {
        for (i in 1:object$ncomp[j])
        {
            data$data=res.predict[[j]][,,i]
            title=paste("ROC Curve\nBlock: ", names(res.predict)[j], ", comp: ",i, sep="")
            
            plot.temp = ifelse(i%in%roc.comp && names(res.predict)[j]%in%block.temp, plot, FALSE)
            auc.mean[[names(res.predict)[j]]][[paste("comp",i,sep = "")]] = statauc(data, plot = plot.temp, title = title)
            
        }
    }
    return(auc.mean)
}


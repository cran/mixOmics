#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2013
# last modified: 04-2015
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

## note: this should probably be set up as an S3 function.
tune.multilevel <- function (X,
Y,
multilevel,
ncomp = 1,
test.keepX = c(5, 10, 15),
test.keepY = NULL,
already.tested.X = NULL,
already.tested.Y = NULL,
method, #splsda or spls
mode = "regression",
validation = "Mfold",
folds = 10,
dist = "max.dist",
measure = "BER", # one of c("overall","BER")
auc = FALSE,
progressBar = TRUE,
near.zero.var = FALSE,
logratio = "none",
nrepeat=1,
light.output = TRUE
)
{
    
    X = as.matrix(X)
    if (length(dim(X)) != 2 || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
    
    if (missing(multilevel) | is.null(multilevel))
    stop("the 'multilevel' design matrix is missing.", call. = FALSE)
    
    multilevel = as.data.frame(multilevel)

    if ((nrow(multilevel) != nrow(X)))
    stop("unequal number of rows in 'X' and 'multilevel'.", call. = FALSE)
    
    # added condition for the spls case (no need to have the 2n and 3rd column in design)
    if ((ncol(multilevel) != 1) & (method == 'splsda'))
    stop("'multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.",
    call. = FALSE)
    
   
    if (is.factor(multilevel[, 1]))
    {
        design[, 1] = as.numeric(design[, 1])
        warning("the vector sample was converted into a numeric vector", call. = FALSE)
    }
    
    if (length(summary(as.factor(multilevel[, 1]))) == nrow(X))
    stop("Check that the vector sample reflects the repeated measurements")
    
    if (any(table(as.factor(multilevel[, 1])) == "1"))
    {
        cat("The vector sample includes the values: ", as.vector(summary(as.factor(multilevel[, 1]))), "\n")
        stop("sample vector", call. = FALSE)
    }
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    if (is.null(method))
    stop("Input method missing, should be set to 'splsda' or 'spls'", call. = FALSE)
    
    if (!method %in% c("splsda", "spls"))
    stop("'method' should be set to 'splsda' or 'spls'", call. = FALSE)


    if ((method == "splsda") && (!is.null(already.tested.Y)))
    warning("Unecessary already.tested.Y parameter")
    
    if (method == "splsda")
    {
        result = tune.splsda(X = X, Y = Y,
        multilevel = multilevel,
        ncomp = ncomp, test.keepX = test.keepX, dist = dist,
        already.tested.X = already.tested.X, validation = validation, folds = folds,
        measure = measure, auc = auc, progressBar = progressBar, near.zero.var = near.zero.var,
        logratio = logratio, nrepeat = nrepeat, light.output = light.output)
        
    } else {
        if(missing(test.keepY))
        test.keepY = rep(ncol(Y), ncomp)
        
        result = tune.splslevel(X = X, Y = Y,
        multilevel = multilevel,
        mode = mode,
        ncomp = ncomp, test.keepX = test.keepX, test.keepY = test.keepY,
        already.tested.X = already.tested.X, already.tested.Y = already.tested.Y)
    }
    result$call = match.call()
    class(result) = paste("tune",method,sep=".")

    return(result)
}

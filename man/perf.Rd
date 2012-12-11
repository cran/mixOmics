\name{perf}
\encoding{latin1}
\alias{perf}
\alias{perf.pls}
\alias{perf.spls}
\alias{perf.plsda}
\alias{perf.splsda}

\title{Compute evaluation criteria for PLS, sPLS, PLS-DA and sPLS-DA}

\description{
Function to evaluate the performance of the fitted PLS, sparse PLS,
PLS-DA and sparse PLS-DA models using various criteria.
}

\usage{
           
\method{perf}{pls}(object, 
            criterion = c("all", "MSEP", "R2", "Q2"),  
            validation = c("Mfold", "loo"),
           folds = 10, progressBar = TRUE, \ldots)	

\method{perf}{spls}(object, 
          criterion = c("all", "MSEP", "R2", "Q2"), 
          validation = c("Mfold", "loo"), 
          folds = 10, progressBar = TRUE, \ldots)

\method{perf}{plsda}(object, 
          method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
          validation = c("Mfold", "loo"), 
          folds = 10, progressBar = TRUE, \ldots)	

\method{perf}{splsda}(object, 
          method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
          validation = c("Mfold", "loo"),  
          folds = 10, progressBar = TRUE, \ldots)		  
}

\arguments{
  \item{object}{object of class inheriting from \code{"pls"}, \code{"plsda"}, 
    \code{"spls"} or \code{"splsda"}. The function will retrieve some key parameters stored in that object.}
    \item{criterion}{only applies to object inheriting from \code{"pls"}, and \code{"spls"}, the criteria measures to be calculated (see details). Can be set to either \code{"all"}, \code{"MSEP"}, \code{"R2"}, \code{"Q2"}. By default set to \code{"all"}}
  \item{method.predict}{only applies to an object inheriting from \code{"plsda"} or \code{"splsda"} to evaluate the classification performance of the model. Should be a subset of \code{"max.dist"}, \code{"centroids.dist"}, \code{"mahalanobis.dist"}. Default is \code{"all"}. See \code{\link{predict}}.}
  \item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or 
    \code{"loo"} (see below). Default is \code{"Mfold"}.}
  \item{folds}{the folds in the Mfold cross-validation. See Details.}
  \item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
  \item{...}{arguments to pass to \code{\link{nearZeroVar}}.}
}

\details{
For fitted PLS and sPLS regression models, \code{perf} estimates the 
mean squared error of prediction (MSEP), \eqn{R^2}, and \eqn{Q^2} to assess the predictive 
perfity of the model using M-fold or leave-one-out cross-validation. Note that only the \code{classic}, \code{regression} and  \code{invariant} modes can be applied.

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{folds}.
The folds also can be supplied as a list of vectors containing the indexes defining each 
fold as produced by \code{split}.
If \code{validation = "loo"}, leave-one-out cross-validation is performed.

For fitted PLS-DA and sPLS-DA models, \code{perf} estimates the classification error rate 
using cross-validation. How many folds to generate is selected such that there is at least 1 
sample for each class in the test set.

For the sparse approaches (sPLS and sLDA), note that the \code{perf} function will retrieve the \code{keepX} and \code{keepY} inputs from the previously run object. The sPLS or sPLS-DA functions will be run again on several and different subsets of data (the cross-folds) and certainly on different subset of selected features. For sPLS, the MSEP, \eqn{R^2}, and \eqn{Q^2} criteria are averaged across all folds. For sPLS-DA, the classification erro rate is averaged across all folds.
A feature stability measure is output for the user to assess how often the variables are selected across all folds.
}

\value{
For PLS and sPLS models, \code{perf} produces a list with the following components: 
  \item{MSEP}{Mean Square Error Prediction for each \eqn{Y} variable, only applies to object inherited from \code{"pls"}, and \code{"spls"}.}
  \item{R2}{a matrix of \eqn{R^2} values of the \eqn{Y}-variables for models 
    with \eqn{1, \ldots ,}\code{ncomp} components, only applies to object inherited from \code{"pls"}, and \code{"spls"}.}
  \item{Q2}{if \eqn{Y} containts one variable, a vector of \eqn{Q^2} values else a list with 
    a matrix of \eqn{Q^2} values for each \eqn{Y}-variable. Note that in the specific case of an sPLS model, it is better to have a look at the Q2.total criterion, only applies to object inherited from \code{"pls"}, and \code{"spls"}}	
    \item{Q2.total}{a vector of \eqn{Q^2}-total values for models with \eqn{1, \ldots ,}\code{ncomp} components, only applies to object inherited from \code{"pls"}, and \code{"spls"}}
    \item{features}{a list of features selected across the folds (\code{$stable.X} and \code{$stable.Y}) or on the whole data set (\code{$final}) for the \code{keepX} and \code{keepY} parameters from the input object.}
    \item{error.rate}{
For PLS-DA and sPLS-DA models, \code{perf} produces a matrix of classification error rate estimation. 
The dimensions correspond to the components in the model and to the prediction method used, respectively. Note that error rates reported in any component include the performance of the model in earlier components for the specified \code{keepX} parameters (e.g. error rate reported for component 3 for \code{keepX = 20} already includes the fitted model on components 1 and 2 for  \code{keepX = 20}). For more advanced usage of the \code{perf} function, see \url{http://perso.math.univ-toulouse.fr/mixomics/methods/spls-da/} and consider using the \code{predict} function.}
}

\references{
Tenenhaus, M. (1998). \emph{La r\'egression PLS: th\'eorie et pratique}. Paris: Editions Technic.  

L\^e Cao, K. A., Rossouw D., Robert-Grani\'e, C. and Besse, P. (2008). A sparse PLS for variable 
selection when integrating Omics data. \emph{Statistical Applications in Genetics and Molecular 
Biology} \bold{7}, article 35.

Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction (MSEP) Estimates for Principal Component 
Regression (PCR) and Partial Least Squares Regression (PLSR). \emph{Journal of Chemometrics} \bold{18}(9), 422-429.
}

\author{S\'ebastien D\'ejean, Ignacio Gonz`alez, Amrit Singh, Kim-Anh L\^e Cao.}

\seealso{\code{\link{predict}}, \code{\link{nipals}}, \code{\link{plot.perf}} and \url{http://perso.math.univ-toulouse.fr/mixomics/} for more details.}

\examples{
## validation for objects of class 'pls' (regression)
# ----------------------------------------
\dontrun{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic


# try tune the number of component to choose
# ---------------------
# first learn the full model
liver.pls <- pls(X, Y, ncomp = 10)

# with 5-fold cross validation: we use the same parameters as in model above
# but we perform cross validation to compute the MSEP, Q2 and R2 criteria
# ---------------------------
liver.val <- perf(liver.pls, validation = "Mfold", folds = 5)

# Q2 total should decrease until it reaches a threshold
liver.val$Q2.total

# ncomp = 3 is enough
plot(liver.val$Q2.total, type = 'l', col = 'red', ylim = c(-0.1, 0.5), 
	xlab = 'PLS components', ylab = 'Q2 total')
abline(h = 0.0975, col = 'darkgreen')
legend('topright', col = c('red', 'darkgreen'), legend = c('Q2 total', 'threshold 0.0975')
	, lty = 1)
title('Liver toxicity PLS 5-fold, Q2 values')

#have a look at the other criteria
# ----------------------
# R2
liver.val$R2
matplot(t(liver.val$R2), type = 'l', xlab = 'PLS components', ylab = 'R2 for each variable')
title('Liver toxicity PLS 5-fold, R2 values')

# MSEP
liver.val$MSEP
matplot(t(liver.val$MSEP), type = 'l', xlab = 'PLS components', ylab = 'MSEP for each variable')
title('Liver toxicity PLS 5-fold, MSEP values')


## validation for objects of class 'spls' (regression)
# ----------------------------------------
ncomp = 7
# first, learn the model on the whole data set
model.spls = spls(X, Y, ncomp = ncomp, mode = 'regression',
	 keepX = c(rep(5, ncomp)), keepY = c(rep(2,ncomp)))


# with leave-one-out cross validation
set.seed(45)
model.spls.loo.val <- perf(model.spls, validation = "loo")

#Q2 total
model.spls.loo.val$Q2.total

# R2:we can see how the performance degrades when ncomp increases
# results are similar to 5-fold
model.spls.loo.val$R2


## validation for objects of class 'splsda' (classification)
# ----------------------------------------
data(srbct)
X <- srbct$gene
Y <- srbct$class  

ncomp = 5

srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))  

# with Mfold
# ---------
set.seed(45)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8, 
               method.predict = "all")

plot(error, type = "l")
}

}

\keyword{regression}
\keyword{multivariate}
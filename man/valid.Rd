\name{valid}
\encoding{latin1}
\alias{valid}

\title{Compute validation criterion for PLS, sPLS, PLS-DA and sPLS-DA}

\description{
Function to estimate measures of the prediction error for fitted PLS, sparse PLS,
PLS-DA and sparse PLS-DA models. M-fold and leave-one-out cross-validation are implemented.
}

\usage{
valid(X, Y, ncomp = min(6, ncol(X)), 
      method = c("pls", "spls", "plsda", "splsda"),
      mode = c("regression", "invariant", "classic"),
      pred.method = c("all", "max.dist", 
                      "centroids.dist", "mahalanobis.dist"),
      criterion = c("all", "MSEP", "R2", "Q2"),
      keepX = NULL, keepY = NULL, 
      validation = c("Mfold", "loo"), M = 10,
      max.iter = 500, tol = 1e-06, ...)
}

\arguments{
  \item{X}{numeric matrix of predictors.}
  \item{Y}{numeric vector or matrix of responses (for multi-response models).}
  \item{ncomp}{the number of components to include in the model. Default is \code{1:min(6, ncol(X))}.}
  \item{method}{Choose between \code{pls}, \code{spls}, \code{plsda} or \code{splsda}.}
  \item{mode}{character string. What type of algorithm to use, matching one of \code{"classic"}, \code{"invariant"} 
    or \code{"regression"}.}
  \item{pred.method}{prediction method to be applied for \code{plsda} or \code{splsda}. Should be a subset 
    of \code{"max.dist"}, \code{"centroids.dist"}, \code{"mahalanobis.dist"}. 
	Default is \code{"all"}. See \code{\link{predict}}.}
  \item{criterion}{what type of validation criterion to be used for \code{pls} or \code{spls}. Should be a subset 
    of \code{"MSEP"}, \code{"R2"} or \code{"Q2"}. Default is \code{"all"}.}
  \item{keepX}{if \code{method="spls"} numeric vector of length \code{ncomp}, the number of variables
    weights to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
  \item{keepY}{if \code{method="spls"} numeric vector of length \code{ncomp}, the number of variables
    weights to keep in \eqn{Y}-loadings. By default all variables are kept in the model.}
  \item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or 
    \code{"loo"} (see below). Default is \code{"Mfold"}.}
  \item{M}{the number of folds in the Mfold cross-validation.}
  \item{max.iter}{integer, the maximum number of iterations.}
  \item{tol}{a not negative real, the tolerance used in the iterative algorithm.}
  \item{...}{arguments to pass to \code{\link{nearZeroVar}}.}
}

\details{
For fitted PLS and sPLS regression models, \code{valid} estimates the 
mean squared error of prediction (MSEP), \eqn{R^2}, and \eqn{Q^2} to assess the predictive 
validity of the model using M-fold or 
leave-one-out cross-validation. Note that only the \code{classic}, \code{regression} and 
\code{invariant} modes can be applied.

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{M}.
If \code{validation = "loo"}, leave-one-out cross-validation is performed.

For fitted PLS-DA and sPLS-DA models, \code{valid} estimates the classification error rate 
using cross-validation. How many folds to generate is selected such that there is at least 1 
sample for each class in the test set.
}

\value{
For PLS and sPLS models, \code{valid} produces a list with the following components: 
  \item{MSEP}{Mean Square Error Prediction for each \eqn{Y} variable.}
  \item{R2}{a matrix of \eqn{R^2} values of the \eqn{Y}-variables for models 
    with \eqn{1, \ldots ,}\code{ncomp} components.}
  \item{Q2}{if \eqn{Y} containts one variable, a vector of \eqn{Q^2} values else a list with 
    a matrix of \eqn{Q^2} values for each \eqn{Y}-variable and a vector of \eqn{Q^2}-total values for models 
    with \eqn{1, \ldots ,}\code{ncomp} components.}	
	
For PLS-DA and sPLS-DA models, \code{valid} produces a matrix of classification error rate estimation. 
The dimensions correspond to the components in the model and to the prediction method used, respectively.
}

\references{
Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.  

Lê Cao, K. A., Rossouw D., Robert-Granié, C. and Besse, P. (2008). A sparse PLS for variable 
selection when integrating Omics data. \emph{Statistical Applications in Genetics and Molecular 
Biology} \bold{7}, article 35.

Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction (MSEP) Estimates for Principal Component 
Regression (PCR) and Partial Least Squares Regression (PLSR). \emph{Journal of Chemometrics} \bold{18}(9), 422-429.
}

\author{Sébastien Déjean, Ignacio González and Kim-Anh Lê Cao.}

\seealso{\code{\link{predict}}, \code{\link{nipals}}, \code{\link{plot.valid}} and http://www.math.univ-toulouse.fr/~biostat/mixOmics/ for more details.}

\examples{
## validation for objects of class 'pls' or 'spls'
\dontrun{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

liver.val <- valid(X, Y, ncomp = 5, mode = "regression", 
                   method = "pls", validation = "Mfold")
				   
plot(liver.val, criterion = "R2", type = "l", layout = c(2, 2))
}

## validation for objects of class 'plsda' or 'splsda'
data(srbct)
X <- srbct$gene
Y <- srbct$class  
  
error <- valid(X, Y, ncomp = 10, method = "splsda", M = 3, 
               keepX = rep(10, 10), pred.method = "all")

plot(error, type = "l")
}

\keyword{regression}
\keyword{multivariate}

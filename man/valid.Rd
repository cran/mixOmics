\name{valid}
\encoding{latin1}
\alias{valid}

\title{Compute validation criterion for PLS and sparse PLS}

\description{
Function to estimate the mean squared error of prediction (MSEP), 
root mean squared error of prediction (RMSEP) and \eqn{R^2} for fitted PLS and 
sPLS regression models. M-fold and leave-one-out cross-validation are implemented.
}

\usage{
valid(X, Y, ncomp = min(6, ncol(X)), 
      mode = c("regression", "invariant", "classic"),
      method = c("pls", "spls"),
      keepX = if(method == "pls") NULL else c(rep(ncol(X), ncomp)),
      keepY = if(method == "pls") NULL else c(rep(ncol(Y), ncomp)),
      validation = c("loo", "Mfold"),
      M = if(validation == "Mfold") 10 else nrow(X),
      max.iter = 500, 
      tol = 1e-06,
      na.action = c("omit", "predict"),
      predict.par = NULL)
}

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{numeric vector or matrix of responses (for multi-response models). \code{NA}s are allowed.}
  \item{ncomp}{the number of components to include in the model. Default is from one to \code{min(6, ncol(X)}.}
  \item{mode}{character string. What type of algorithm to use, matching one of \code{"classic"}, \code{"invariant"} 
    or \code{"regression"}.}
  \item{method}{Choose between \code{pls} or \code{spls}.}
  \item{keepX}{if \code{method="spls"} numeric vector of length \code{ncomp}, the number of variables
    weights to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
  \item{keepY}{if \code{method="spls"} numeric vector of length \code{ncomp}, the number of variables
    weights to keep in \eqn{Y}-loadings. By default all variables are kept in the model.}
  \item{validation}{character.  What kind of (internal) validation to use.  See below.}
  \item{M}{the number of folds in the Mfold cross-validation.}
  \item{max.iter}{integer, the maximum number of iterations.}
  \item{tol}{a not negative real, the tolerance used in the iterative algorithm.}
  \item{na.action}{action determining what should be done with missing values in \code{X}. 
    One of \code{"predict"} or \code{"omit"} (see Details).}
  \item{predict.par}{further arguments sent to \code{\link{nipals}} function.}
}

\details{
If \code{na.action = "predict"} the estimation of the missing values is performed 
by the reconstitution of the data matrix using the \code{nipals} function. Otherwise, missing 
values are handled by deletion of incomplete cases. 

The validation criterion \code{"MSEP"}, \code{"RMSEP"} or \code{"R2"} allows one to assess the 
predictive validity of the model using M-fold or leave-one-out cross-validation.
Note that only the \code{classic}, \code{regression} and \code{invariant} modes can be applied.

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{M}.
If \code{validation = "loo"}, leave-one-out cross-validation is performed.
}

\value{
\code{valid} produces a list with the following components: 
  \item{msep}{Mean Square Error Prediction for each Y variable.}
  \item{rmsep}{Root Mean Square Error Prediction for each Y variable.}
  \item{r2}{a matrix of \eqn{R^2} values of the \eqn{Y}-variables for models 
    with \eqn{1, ... ,\code{ncomp}} components.}
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

\seealso{\code{\link{predict}}, \code{\link{nipals}}, code{\link{plot.valid}}.}

\examples{
\dontrun{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

liver.val <- valid(X, Y, ncomp = 5, mode = "regression", 
                   method = "pls", validation = "loo")
}
}
\keyword{regression}
\keyword{multivariate}

\name{valid}
\encoding{latin1}
\alias{valid}

\title{Compute validation criterion for PLS and Sparse PLS}

\description{
Function to estimate the root mean squared error of prediction (RMSEP) and the Q2 criterion for PLS (classic, regression and invariant modes) and sPLS (regression).  Cross-validation or leave-one-out cross-validation are implemented.
}

\usage{
valid(X, Y, ncomp = 3,
      mode = c("regression", "invariant", "classic"),
      max.iter = 500, tol = 1e-06, criterion = c("rmsep", "q2"),
      method = c("pls", "spls"),
      keepX = if(method == "pls") NULL else c(rep(ncol(X), ncomp)),
      keepY = if(method == "pls") NULL else c(rep(ncol(Y), ncomp)),
      scaleY = TRUE,
      validation = c("loo", "Mfold"),
      M = if(validation == 'Mfold') 10 else nrow(X))
}

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{numeric vector or matrix of responses (for multi-response models). \code{NA}s are allowed.}
  \item{ncomp}{the number of components to include in the model.  Default is from one to the rank of \code{X}.}
  \item{mode}{character string. What type of algorithm to use, matching one of \code{"regression"}, \code{"invariant"} or \code{"classic"}.}
  \item{max.iter}{integer, the maximum number of iterations.}
  \item{tol}{a not negative real, the tolerance used in the iterative algorithm.}
  \item{criterion}{character string. What type of validation criterion to use, see details.}
  \item{method}{ character. \code{pls} or \code{spls} methods.}
  \item{keepX}{if \code{method="spls"} numeric vector of length \code{ncomp}, the number of variables
    weights to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
  \item{keepY}{if \code{method="spls"} numeric vector of length \code{ncomp}, the number of variables
    weights to keep in \eqn{Y}-loadings. By default all variables are kept in the model.}
  \item{scaleY}{should the Y data be scaled ? In the case of a 'discriminant' version of the (s)PLS
where the Y data are of discrete type, this should be set to \code{FALSE}.}
  \item{validation}{character.  What kind of (internal) validation to use.  See below.}
  \item{M}{the number of folds in the Mfold cross-validation.}


}

\details{

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{M}.
If \code{validation = "loo"}, leave-one-out cross-validation is performed.

The validation criterion \code{"rmsep"} allows one to assess the predictive validity of the model (using loo or cross-validation). It produces the estimated error obtained by evaluating the PLS or the sPLS models. \code{"q2"} helps choosing the number of (s)PLS dimensions. \code{rmsep}. Note that only the classic, regression and invariant modes can be applied.

What follows is the definition of these criteria: 

Let \eqn{n} the number of individuals (experimetals units).
The fraction of the variation of a variable \eqn{y_{k}} that can be predicted 
by a component, as estimated by cross-validation, is computed as: 

\deqn{Q_{kh}^2 = 1-\frac{PRESS_{kh}}{RSS_{k(h-1)}}} 

where 

\deqn{PRESS_{kh} = \sum_{i=1}^{n}(y_{ik} - \hat{y}_{(-i)k}^h)^2}

is the PRediction Error Sum of Squares and

\deqn{RSS_{kh} = \sum_{i=1}^{n}(y_{ik} - \hat{y}_{ik}^h)^2}

is the Residual Sum of Squares for the variable \eqn{k}, (\eqn{k=1, \ldots ,q}) 
and the PLS variate \eqn{h}, (\eqn{h=1, \ldots ,H}).
For \eqn{h=0}, \eqn{RSS_{kh} = n-1}. 

The fraction of the total variation of \eqn{Y} that can be predicted by a component, 
as estimated by cross-validation, is computed as: 

\deqn{Q_h^2 = 1-\frac{\sum_{k=1}^{q}PRESS_{kh}}{\sum_{k=1}^{q}RSS_{k(h-1)}}} 

The cumulative \eqn{(Q_{cum}^2)_{kh}} of a variable is computed as:

\deqn{(Q_{cum}^2)_{kh} = 1-\prod_{j=1}^h\frac{PRESS_{kj}}{RSS_{k(j-1)}}}

and the cumulative \eqn{(Q_{cum}^2)_h} for the extracted components is computed as:

\deqn{(Q_{cum}^2)_h = 1-\prod_{j=1}^h\frac{\sum_{k=1}^{q}PRESS_{kj}}{\sum_{k=1}^{q}RSS_{k(j-1)}}} 
}


\value{
\code{valid} produces a list with the following components: 
  \item{Y.hat}{the predicted values using cross-validation}
  \item{fold}{indicates which folds the samples belong to wen using k-fold cross-validation}

  \item{rmsep}{if \code{validation="rmsep"} Root Mean Square Error Prediction for each Y variable}
  \item{RSS}{if \code{validation="q2"} a matrix of RSS values of the \eqn{Y}-variables for models 
    with \eqn{1, \ldots ,\code{ncomp}} components.}
  \item{PRESS}{if \code{validation="q2"} prediction error sum of squares of the \eqn{Y}-variables.
    A matrix of PRESS values for models with \eqn{1, \ldots ,\code{ncomp}} components.}
  \item{q2}{if \code{validation="q2"} vector of \eqn{Q^2} values for the extracted components.}


}

\references{
Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.  

Lê Cao, K. A., Rossouw D., Robert-Granié, C. and Besse, P. (2008). A sparse PLS for variable 
selection when integrating Omics data. \emph{Statistical Applications in Genetics and Molecular 
Biology} \bold{7}, article 35.

}

\author{Sébastien Déjean, Ignacio González and Kim-Anh Lê Cao.}

\seealso{\code{\link{predict}}.}

\examples{
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological

## computing the RMSEP with 10-fold CV with pls
error <- valid(X, Y, mode = "regression", ncomp = 3, method = "pls", 
               validation = "Mfold", criterion = "rmsep")
error$rmsep 
}

\keyword{regression}
\keyword{multivariate}

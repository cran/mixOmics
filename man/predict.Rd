\name{predict}
\encoding{latin1}
\alias{predict.splsda}
\alias{predict.plsda}
\alias{predict.pls}
\alias{predict.spls}

\title{Predict Method for PLS, sparse PLS, PLSDA Regression or Sparse PLSDA}

\description{
Predicted values based on PLS, sparse PLS, PLSDA Regression or Sparse PLSDA models. New responses and 
variates are predicted using a fitted model and a new matrix of observations.
}

\usage{
\method{predict}{pls}(object, newdata, ...)
\method{predict}{spls}(object, newdata, ...)
\method{predict}{splsda}(object, newdata, method = c("class.dist", "centroids.dist", "Sr.dist", "max.dist"), ...)
\method{predict}{plsda}(object, newdata, method = c("class.dist", "centroids.dist", "Sr.dist", "max.dist"), ...)
}

\arguments{
  \item{object}{object of class inheriting from \code{"splsda"} or \code{"plsda"}.}
  \item{newdata}{data matrix in which to look for for explanatory variables to be used for prediction.}
  \item{method}{method to be applied for \code{"splsda"} or \code{"plsda"} to predict the class of new data. It can be either "max.dist", "centroids.dist", "class.dist", or "Sr.dist", see details.}
  \item{...}{other arguments to be passed \code{predict.default}.}
}

\details{
\code{predict} produces predicted values, obtained by evaluating the PLS/sPLS/PLSDA or sparse PLSDA  
model returned by \code{pls}, \code{spls}, \code{plsda} or \code{splsda} in the frame \code{newdata}. 
Variates for \code{newdata} are also returned. 
Different prediction methods are proposed:
-max.dist is the naive method to predict the class. It is based on the predicted Y matrix (object$predict) which can be seen as a probability matrix to assign each  test data to a class.
}

\value{
\code{predict} produces a list with the following components: 
  \item{predict}{a three dimensional array of predicted response values. The dimensions 
	correspond to the observations, the response variables and the model dimension, respectively.}
  \item{variates}{matrix of predicted variates.}
  \item{B.hat}{matrix of regression coefficients (without the intercept).}
  \item{class}{predicted class of the newdata observation for each PLS dimension.}
  \item{centroid}{}
}

\references{
Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.  
}

\author{Sébastien Déjean, Ignacio González, Kim-Anh Lê Cao and Pierre Monget}

\seealso{\code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}}, \code{\link{splsda}}.}

\examples{
## First example with plsda
data(liver.toxicity)
X = as.matrix(liver.toxicity$gene)
Y = as.factor(liver.toxicity$treatment[,4] )

# if training is perfomed on 4/5th of the original data
samp = sample(1:5, nrow(X), replace = TRUE)  
test = which(samp == 1)   # testing on the first fold
train = setdiff(1:nrow(X), test)

plsda.train = plsda(X[train,], Y[train], ncomp = 1,  mode = 'regression')
test.predict = predict(plsda.train, X[test,], method = "class.dist")
test.predict$cl

## Second example with splsda
data(liver.toxicity)
X = as.matrix(liver.toxicity$gene)
Y = as.factor(liver.toxicity$treatment[,4])  # time points

# if training is perfomed on 4/5th of the original data
samp = sample(1:5, nrow(X), replace = TRUE)  
test = which(samp == 1)   # testing on the first fold
train = setdiff(1:nrow(X), test)


splsda.train = splsda(X[train,], Y[train], ncomp = 1, keepX = 20, mode = 'regression')
test.predict = predict(splsda.train, X[test,], method = "class.dist")
test.predict$cl
}

\keyword{regression}
\keyword{multivariate}

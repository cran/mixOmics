\name{predict}
\encoding{latin1}
\alias{predict.splsda}
\alias{predict.plsda}
\alias{predict.pls}
\alias{predict.spls}

\title{Predict Method for PLS, sPLS, PLS-DA or sPLS-DA}

\description{
Predicted values based on PLS, sparse PLS, PLS-DA or sparse PLS-DA models. New responses and 
variates are predicted using a fitted model and a new matrix of observations.
}

\usage{
\method{predict}{pls}(object, newdata, ...)

\method{predict}{spls}(object, newdata, ...)

\method{predict}{plsda}(object, newdata, method = c("all", "max.dist", 
        "centroids.dist", "mahalanobis.dist"), ...)

\method{predict}{splsda}(object, newdata, method = c("all", "max.dist", 
        "centroids.dist", "mahalanobis.dist"), ...)
}

\arguments{
  \item{object}{object of class inheriting from \code{"pls"}, \code{"spls"}, \code{"plsda"} or \code{"splsda"}.}
  \item{newdata}{data matrix in which to look for for explanatory variables to be used for prediction.}
  \item{method}{method to be applied for \code{plsda} or \code{splsda} to predict the class of new data, 
    should be a subset of \code{"centroids.dist"}, \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details). 
	Defaults to \code{"all"}.}
  \item{...}{not used currently.}
}

\details{
\code{predict} produces predicted values, obtained by evaluating the PLS, sparse PLS, PLSDA or sparse PLSDA  
model returned by \code{pls}, \code{spls}, \code{plsda} or \code{splsda} in the frame \code{newdata}. 
Variates for \code{newdata} are also returned. 

Different class prediction methods are proposed for \code{plsda} or \code{splsda}: \code{"max.dist"} 
is the naive method to predict the class. It is based on the predicted matrix (\code{object$predict}) 
which can be seen as a probability matrix to assign each test data to a class. The class with the largest 
class value is the predicted class. \code{"centroids.dist"} allocates the individual \eqn{x} to the class of 
\eqn{Y} minimizing \eqn{dist(\code{x-variate}, G_l)}, where \eqn{G_l}, \eqn{l = 1,...,L} are the centroids of 
the classes calculated on the \eqn{X}-variates of the model. \code{"mahalanobis.dist"} allocates the individual 
\eqn{x} to the class of \eqn{Y} as in \code{"centroids.dist"} but by using the Mahalanobis metric 
in the calculation of the distance.
}

\value{
\code{predict} produces a list with the following components: 
  \item{predict}{a three dimensional array of predicted response values. The dimensions 
	correspond to the observations, the response variables and the model dimension, respectively.}
  \item{variates}{matrix of predicted variates.}
  \item{B.hat}{matrix of regression coefficients (without the intercept).}
  \item{class}{vector or matrix of predicted class by using \eqn{1,...,}\code{ncomp} 
    (sparse)PLS-DA components.}
  \item{centroid}{matrix of coordinates for centroids.}
}

\references{
Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.  
}

\author{Sébastien Déjean, Ignacio González, Kim-Anh Lê Cao and Pierre Monget}

\seealso{\code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}}, \code{\link{splsda}} and http://www.math.univ-toulouse.fr/~biostat/mixOmics/ for more details.}

\examples{
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, ncomp = 2, mode = "classic")

indiv1 <- c(200, 40, 60)
indiv2 <- c(190, 45, 45)
newdata <- rbind(indiv1, indiv2)
colnames(newdata) <- colnames(X)
newdata

pred <- predict(linn.pls, newdata)

plotIndiv(linn.pls, comp = 1:2, rep.space = "X-variate")
points(pred$variates[, 1], pred$variates[, 2], pch = 19, cex = 1.2)
text(pred$variates[, 1], pred$variates[, 2], 
     c("new ind.1", "new ind.2"), pos = 3)
	 
## First example with plsda
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- as.factor(liver.toxicity$treatment[, 4])

## if training is perfomed on 4/5th of the original data
samp <- sample(1:5, nrow(X), replace = TRUE)  
test <- which(samp == 1)   # testing on the first fold
train <- setdiff(1:nrow(X), test)

plsda.train <- plsda(X[train, ], Y[train], ncomp = 2)
test.predict <- predict(plsda.train, X[test, ], method = "max.dist")
Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
cbind(Y = as.character(Y[test]), Prediction)

\dontrun{
## Second example with splsda
splsda.train <- splsda(X[train, ], Y[train], ncomp = 2, keepX = c(30, 30))
test.predict <- predict(splsda.train, X[test, ], method = "max.dist")
Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
cbind(Y = as.character(Y[test]), Prediction)
}
}

\keyword{regression}
\keyword{multivariate}

\name{predict}
\encoding{latin1}
\alias{predict.pls}
\alias{predict.spls}
\alias{predict.splsda}
\alias{predict.plsda}

\title{Predict Method for PLS Regression and Sparse PLS}

\description{
Predicted values based on PLS regression or sparse PLS models. New responses and 
variates are predicted using a fitted model and a new matrix of observations.
}

\usage{
\method{predict}{pls}(object, newdata, \ldots)

\method{predict}{spls}(object, newdata, \ldots)
}

\arguments{
  \item{object}{object of class inheriting from \code{"pls"} or \code{"spls"}.}
  \item{newdata}{data matrix in which to look for for explanatory variables 
    to be used for prediction.}
  \item{\ldots}{not used currently.}
}

\details{
\code{predict} produces predicted values, obtained by evaluating the PLS  
model returned by \code{pls} or \code{spls} in the frame \code{newdata}. 
Variates for \code{newdata} are also returned. 
}

\value{
\code{predict} produces a list with the following components: 
  \item{predict}{a three dimensional array of predicted response values. The dimensions 
	correspond to the observations, the response variables and the model dimension, respectively.}
  \item{variates}{matrix of predicted variates.}
  \item{B.hat}{matrix of regression coefficients (without the intercept).}
}

\references{
Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.  
}

\author{Sébastien Déjean, Ignacio González and Kim-Anh Lê Cao}

\seealso{\code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}}, \code{\link{splsda}}.}

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

plotIndiv(linn.pls, comp= 1:2, rep.space = "X-variate")
points(pred$variates[, 1], pred$variates[, 2], pch = 19, cex = 1.2)
text(pred$variates[, 1], pred$variates[, 2], 
     c("new ind.1", "new ind.2"), pos = 3)
}

\keyword{regression}
\keyword{multivariate}

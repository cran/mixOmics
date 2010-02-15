\name{pls}
\encoding{latin1}
\alias{pls}

\title{Partial Least Squares (PLS) Regression}

\description{Functions to perform Partial Least Squares (PLS) regression.

}

\usage{
pls(X, Y, ncomp = 3, 
    mode = c("regression", "canonical", "invariant", "classic"), 
    max.iter = 500, tol = 1e-06, scaleY=TRUE)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{numeric vector or matrix of responses (for multi-response models). 
    \code{NA}s are allowed.}
  \item{ncomp}{the number of components to include in the model.  
    Default is from one to the rank of \code{X}.}
  \item{mode}{character string. What type of algorithm to use, (partially) matching 
    one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}. 
	See Details.}
  \item{max.iter}{integer, the maximum number of iterations.}
  \item{tol}{a not negative real, the tolerance used in the iterative algorithm.}
  \item{scaleY}{should the Y data be scaled ? In the case of a 'discriminant' version of the PLS
where the Y data are of discrete type, this should be set to \code{FALSE}.}
  \item{...}{not used currently.}

}

\details{
\code{pls} function fit PLS models with \eqn{1, \ldots ,}\code{ncomp} components. 
Multi-response models are fully supported. The \code{X} and \code{Y} datasets
can contain missing values.

The type of algorithm to use is specified with the \code{mode} argument. Four PLS 
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis 
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS 
algorithm \code{("classic")} (see References).

The number of components to fit is specified with the argument \code{ncomp}. 
It this is not supplied, the rank of \code{X} is used. The rank is compute by
using the \code{\link{mat.rank}} function. 
}

\value{
\code{pls} returns an object of class \code{"pls"}, a list 
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized original response vector or matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{mode}{the algoritthm used to fit the model.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the \code{X} and \code{Y} variates.}
  \item{loadings}{list containing the estimated loadings for the variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
}

\references{
Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.

Wold H. (1966). Estimation of principal components and related models by iterative least squares. 
In: Krishnaiah, P. R. (editors), \emph{Multivariate Analysis}. Academic Press, N.Y., 391-420.
}

\author{Sébastien Déjean and Ignacio González and Kim-Anh Lê Cao.}

\seealso{\code{\link{spls}}, \code{\link{summary}}, \code{\link{mat.rank}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}.}

\examples{
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, mode = "classic")

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.pls <- pls(X, Y, ncomp = 3)
}

\keyword{regression}
\keyword{multivariate}

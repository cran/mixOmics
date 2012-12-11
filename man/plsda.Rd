\name{plsda}
\encoding{latin1}
\alias{plsda}

\title{Partial Least Squares Discriminate Analysis (PLS-DA).}

\description{Function to perform standard Partial Least Squares regression
to classify samples.
}

\usage{
plsda(X, Y, ncomp = 2, max.iter = 500, tol = 1e-06, ...)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{a factor or a class vector for the discrete outcome.}
  \item{ncomp}{the number of components to include in the model. Default to 2.}
  \item{max.iter}{integer, the maximum number of iterations.}  
  \item{tol}{a not negative real, the tolerance used in the iterative algorithm.}
  \item{...}{arguments to pass to \code{\link{nearZeroVar}}.}
}

\details{
\code{plsda} function fit PLS models with \eqn{1,...,}\code{ncomp} components
to the factor or class vector \code{Y}. The appropriate indicator 
matrix is created.
}

\value{
\code{plsda} returns an object of class \code{"plsda"}, a list 
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized indicator response vector or matrix.}
  \item{ind.mat}{the indicator matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the \code{X} and \code{Y} variates.}
  \item{loadings}{list containing the estimated loadings for the variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}
}

\references{
Pérez-Enciso, M. and Tenenhaus, M. (2003). Prediction of clinical outcome with microarray data: 
a partial least squares discriminant analysis (PLS-DA) approach. \emph{Human Genetics} 
\bold{112}, 581-592.

Nguyen, D. V. and Rocke, D. M. (2002). Tumor classification by partial
least squares using microarray gene expression data. \emph{Bioinformatics}
\bold{18}, 39-50.

Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.
}

\author{Ignacio González, Kim-Anh Lê Cao and Pierre Monget.}

\seealso{\code{\link{splsda}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}}, \code{\link{predict}}, \code{\link{valid}} and http://www.math.univ-toulouse.fr/~biostat/mixOmics/ for more details.}

\examples{
## First example
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)
palette(c("red", "blue"))
col.breast <- as.numeric(as.factor(Y))
plotIndiv(plsda.breast, ind.names = TRUE, col = col.breast)
legend('bottomleft', c("After", "Before"), pch = c(16, 16), 
       col = unique(col.breast), cex = 1, pt.cex = c(1.2, 1.2), 
       title = "Treatment")
palette("default")

## Second example
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$treatment[, 4]

plsda.liver <- plsda(X, Y, ncomp = 2)
col.rat <- as.numeric(as.factor(Y))
plotIndiv(plsda.liver, col = col.rat, ind.names = Y)
}

\keyword{regression}
\keyword{multivariate}

\name{splsda}
\encoding{latin1}
\alias{splsda}

\title{Sparse Partial Least Squares Discriminate Analysis (sPLS-DA)}

\description{Function to perform sparse Partial Least Squares to classify samples. The sPLS-DA approach
embeds variable selection for this purpose.
}

\usage{
splsda(X, Y, ncomp = 2, keepX = rep(ncol(X), ncomp),
       max.iter = 500, tol = 1e-06, ...)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{a factor or a class vector for the discrete outcome.}
  \item{ncomp}{the number of components to include in the model (see Details).}
  \item{keepX}{numeric vector of length \code{ncomp}, the number of variables
    to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
  \item{max.iter}{integer, the maximum number of iterations.}   
  \item{tol}{a positive real, the tolerance used in the iterative algorithm.} 
  \item{...}{arguments to pass to \code{\link{nearZeroVar}}.}
}

\details{
\code{splsda} function fit sPLS models with \eqn{1, \ldots ,}\code{ncomp} components
to the factor or class vector \code{Y}. The appropriate indicator 
matrix is created.
}

\value{
\code{splsda} returns an object of class \code{"splsda"}, a list 
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized indicator response vector or matrix.}
  \item{ind.mat}{the indicator matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{keepX}{number of \eqn{X} variables kept in the model on each component.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the variates.}
  \item{loadings}{list containing the estimated loadings for the \code{X} and 
	\code{Y} variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}
}

\references{
On sPLS-DA:
Lê Cao, K.-A., Boitard, S. and Besse, P. (2011). Sparse PLS Discriminant Analysis: biologically relevant feature selection and graphical displays for multiclass problems. \emph{BMC Bioinformatics} \bold{12}:253.
}

\author{Ignacio González, Kim-Anh Lê Cao and Pierre Monget.}

\seealso{\code{\link{spls}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}},
\code{\link{cim}}, \code{\link{network}}, \code{\link{predict}}, \code{\link{valid}} and http://www.math.univ-toulouse.fr/~biostat/mixOmics/ for more details.}

\examples{
## First example
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))
palette(c("red", "blue"))
col.breast <- as.numeric(as.factor(Y))
plotIndiv(res, ind.names = TRUE, col = col.breast)
legend('bottomleft', c("After", "Before"), pch = c(16, 16), 
       col = unique(col.breast), cex = 1, pt.cex = c(1.2, 1.2), 
       title = "Treatment")
palette("default")

## Second example
\dontrun{
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- liver.toxicity$treatment[, 4]

splsda.liver = splsda(X, Y, ncomp = 2, keepX = c(20, 20))
col.rat <- as.numeric(as.factor(Y))
plotIndiv(splsda.liver, col = col.rat, ind.names = Y)
}
}

\keyword{regression}
\keyword{multivariate}

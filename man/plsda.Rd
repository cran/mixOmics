\name{plsda}
\encoding{latin1}
\alias{plsda}

\title{Partial Least Squares Discriminate Analysis (PLS-DA).}

\description{Function to perform standard Partial Least Squares regression
to classify samples.
}

\usage{
plsda(X, Y, ncomp = 2, max.iter = 500, tol = 1e-06)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{a factor or indicator matrix for the discrete outcome. If a matrix, the 
	entries must be either 0 (not in class) or 1 (in class) for each sample (rows must sum to one).}
  \item{ncomp}{the number of components to include in the model. Default to 2.}
  \item{max.iter}{integer, the maximum number of iterations.}
  \item{tol}{a not negative real, the tolerance used in the iterative algorithm.}
}

\details{
\code{plsda} function fit PLS models with \eqn{1, \ldots ,}\code{ncomp} components
to the indicator matrix \code{Y}. If a factor is supplied, the appropriate indicator 
matrix is created.
}

\value{
\code{plsda} returns an object of class \code{"pls"}, a list 
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the indicator matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the \code{X} and \code{Y} variates.}
  \item{loadings}{list containing the estimated loadings for the variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
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

\author{Ignacio González and Kim-Anh Lê Cao.}

\seealso{\code{\link{splsda}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}}.}

\examples{
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- as.factor(breast.tumors$sample$treatment)

res <- plsda(X, Y, ncomp = 2)
palette(c("red", "blue"))
plotIndiv(res, ind.names = TRUE, col = as.numeric(Y))
legend(18, -9, c("After", "Before"), pch = c(16, 16), 
       col = c("red", "blue"), cex = 1, pt.cex = c(1.2, 1.2), 
       title = "Treatment")
palette("default")
}

\keyword{regression}
\keyword{multivariate}

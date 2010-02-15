\name{splsda}
\encoding{latin1}
\alias{splsda}

\title{Sparse Partial Least Squares Discriminate Analysis (sPLS-DA)}

\description{Function to perform sparse Partial Least Squares to classify samples. The sPLS-DA approach
embeds variable selection for this purpose.
}

\usage{
splsda(X, Y, ncomp = 2, keepX = c(rep(ncol(X), ncomp)), max.iter = 500, 
       mode = 'regression', tol = 1e-06, ...)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{a factor or a class vector for the discrete outcome.}
  \item{ncomp}{the number of components to include in the model (see Details).}
  \item{keepX}{numeric vector of length \code{ncomp}, the number of variables
    to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
  \item{max.iter}{integer, the maximum number of iterations.}   
  \item{mode}{indicates the mode of analysis.
  It can be either 'regression', 'classic' or 'invariant'}
  \item{tol}{a positive real, the tolerance used in the iterative algorithm.} 
  \item{...}{arguments to pass to \code{\link{nearZeroVar}}.}
}

\details{
\code{splsda} function fit sPLS models with \eqn{1,...,}\code{ncomp} components
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
  \item{keepX}{number of X variables kept in the model on each component.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the variates.}
  \item{loadings}{list containing the estimated loadings for the \code{X} and 
	\code{Y} variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}
}

\references{
Lê Cao, K.-A., Martin, P.G.P., Robert-Granié, C. and Besse, P. (2009). Sparse canonical 
methods for biological data integration: application to a cross-platform study. 
\emph{BMC Bioinformatics} \bold{10}:34.

Lê Cao, K.-A., Rossouw, D., Robert-Granié, C. and Besse, P. (2008). A sparse PLS for variable 
selection when integrating Omics data. \emph{Statistical Applications in Genetics and Molecular 
Biology} \bold{7}, article 35.

Shen, H. and Huang, J. Z. (2008). Sparse principal component analysis via regularized 
low rank matrix approximation. \emph{Journal of Multivariate Analysis} \bold{99}, 1015-1034.    

Pérez-Enciso, M. and Tenenhaus, M. (2003). Prediction of clinical outcome with microarray data: 
a partial least squares discriminant analysis (PLS-DA) approach. \emph{Human Genetics} 
\bold{112}, 581-592.

Nguyen, D. V. and Rocke, D. M. (2002). Tumor classification by partial
least squares using microarray gene expression data. \emph{Bioinformatics}
\bold{18}, 39-50.

Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.

Wold H. (1966). Estimation of principal components and related models by iterative least squares. 
In: Krishnaiah, P. R. (editors), \emph{Multivariate Analysis}. Academic Press, N.Y., 391-420.
}

\author{Ignacio González, Kim-Anh Lê Cao and Pierre Monget.}

\seealso{\code{\link{spls}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}},
\code{\link{cim}}, \code{\link{network}}.}

\examples{
# First example
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))
palette(c("red", "blue"))
plotIndiv(res, ind.names = TRUE, col = as.numeric(Y))
legend(-0.35, -0.19, c("After", "Before"), pch = c(16, 16), 
       col = c("red", "blue"), cex = 1, pt.cex = c(1.2, 1.2), 
       title = "Treatment")
palette("default")


# Second example
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- liver.toxicity$treatment[, 4]

splsda.liver = splsda(X, Y, ncomp = 2, keepX = c(20, 20))
col = as.numeric(Y)
plotIndiv(splsda.liver, col = col)
}

\keyword{regression}
\keyword{multivariate}

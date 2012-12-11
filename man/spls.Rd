\name{spls}
\encoding{latin1}
\alias{spls}

\title{Sparse Partial Least Squares (sPLS)}

\description{Function to perform sparse Partial Least Squares (sPLS). The sPLS approach
combines both integration and variable selection simultaneously on two data sets
in a one-step strategy.
}

\usage{
spls(X, Y, ncomp = 2, mode = c("regression", "canonical"),
     max.iter = 500, tol = 1e-06, keepX = rep(ncol(X), ncomp), 
     keepY = rep(ncol(Y), ncomp), near.zero.var = TRUE,...)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{numeric vector or matrix of responses (for multi-response models). 
    \code{NA}s are allowed.}
  \item{ncomp}{the number of components to include in the model (see Details). 
    Default is set to from one to the rank of \code{X}.}
  \item{mode}{character string. What type of algorithm to use, (partially) matching 
    one of \code{"regression"} or \code{"canonical"}. See Details.}
  \item{max.iter}{integer, the maximum number of iterations.}
  \item{tol}{a positive real, the tolerance used in the iterative algorithm.}
  \item{keepX}{numeric vector of length \code{ncomp}, the number of variables
    to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
  \item{keepY}{numeric vector of length \code{ncomp}, the number of variables
    to keep in \eqn{Y}-loadings. By default all variables are kept in the model.}
      \item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations.}

  \item{...}{arguments to pass to \code{\link{nearZeroVar}}.}
}

\details{
\code{spls} function fit sPLS models with \eqn{1, \ldots ,}\code{ncomp} components. 
Multi-response models are fully supported. The \code{X} and \code{Y} datasets
can contain missing values.

The type of algorithm to use is specified with the \code{mode} argument. Two sPLS 
algorithms are available: sPLS regression \code{("regression")} and sPLS canonical analysis 
\code{("canonical")} (see References). 

The estimation of the missing values can be performed 
by the reconstitution of the data matrix using the \code{nipals} function. Otherwise, missing 
values are handled by casewise deletion in the \code{spls} function without having to 
delete the rows with missing data.
}

\value{
\code{spls} returns an object of class \code{"spls"}, a list 
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized original response vector or matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{mode}{the algorithm used to fit the model.}
  \item{keepX}{number of \eqn{X} variables kept in the model on each component.}
  \item{keepY}{number of \eqn{Y} variables kept in the model on each component.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the variates.}
  \item{loadings}{list containing the estimated loadings for the \eqn{X} and 
	\eqn{Y} variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}
  \item{tol}{the tolerance used in the iterative algorithm, used for subsequent S3 methods}
  \item{max.iter}{the maximum number of iterations, used for subsequent S3 methods}

}

\references{

Lê Cao, K.-A., Martin, P.G.P., Robert-Granié, C. and Besse, P. (2009). Sparse canonical methods for biological data integration: application to a cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.

Lê Cao, K.-A., Rossouw, D., Robert-Granié, C. and Besse, P. (2008). A sparse PLS for variable 
selection when integrating Omics data. \emph{Statistical Applications in Genetics and Molecular 
Biology} \bold{7}, article 35.

Shen, H. and Huang, J. Z. (2008). Sparse principal component analysis via regularized 
low rank matrix approximation. \emph{Journal of Multivariate Analysis} \bold{99}, 1015-1034.    

Tenenhaus, M. (1998). \emph{La régression PLS: théorie et pratique}. Paris: Editions Technic.

Wold H. (1966). Estimation of principal components and related models by iterative least squares. 
In: Krishnaiah, P. R. (editors), \emph{Multivariate Analysis}. Academic Press, N.Y., 391-420.
}

\author{Sébastien Déjean, Ignacio González and Kim-Anh Lê Cao.}

\seealso{\code{\link{pls}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}},
\code{\link{cim}}, \code{\link{network}}, \code{\link{predict}}, \code{\link{perf}} and http://www.math.univ-toulouse.fr/~biostat/mixOmics/ for more details.}

\examples{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.spls <- spls(X, Y, ncomp = 2, keepX = c(50, 50), 
                      keepY = c(10, 10))
}

\keyword{regression}
\keyword{multivariate}

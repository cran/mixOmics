\name{rcc}
\encoding{latin1}
\alias{rcc}	
\alias{rcc.default}

\title{Regularized Canonical Correlation Analysis}

\description{
The function performs the regularized extension of the Canonical Correlation Analysis
to seek correlations between two data matrices.
}

\usage{
rcc(X, Y, lambda1 = 0, lambda2 = 0, ...)
}

\arguments{
  \item{X}{numeric matrix or data frame \eqn{(n \times p)}, the observations on the \eqn{X} variables.
    \code{NA}s are allowed.}
  \item{Y}{numeric matrix or data frame \eqn{(n \times q)}, the observations on the \eqn{Y} variables.
    \code{NA}s are allowed.}
  \item{lambda1, lambda2}{a not negative real. The regularization parameter for the \emph{X} and
    \emph{Y} data. Defaults to \code{lambda1=lambda2=0}.}
  \item{...}{not used currently.}	
}

\details{
The main purpose of Canonical Correlations Analysis (CCA) is the exploration of sample 
correlations between two sets of variables \eqn{X} and \eqn{Y}   
observed on the same individuals (experimental units) 
whose roles in the analysis are strictly symmetric.
 
The \code{cancor} function performs the core of computations 
but additional tools are required to deal with data sets highly 
correlated (nearly collinear), data sets with more variables 
than units by example.

The \code{rcc} function, the regularized version of CCA, 
is one way to deal with this problem by
including a regularization step in the computations of CCA.
Such a regularization in this context  
was first proposed by Vinod (1976), then developped by Leurgans \emph{et al.} (1993). 
It consists in the regularization of the empirical covariances matrices of  
\eqn{X} and \eqn{Y} by adding a multiple of the matrix identity, that is,
Cov\eqn{(X)+ \lambda_1 I} and Cov\eqn{(Y)+ \lambda_2 I}.

When \code{lambda1=0} and \code{lambda2=0}, \code{rcc} perform classic
CCA, if posible.

\code{rcc} was implemented to allows
the use of data missing.
}

\value{
\code{rcc} returns a object of class \code{"rcc"}, a list that 
contains the following components:
  \item{X}{the original \eqn{X} data.}
  \item{Y}{the original \eqn{Y} data.}
  \item{lambda}{a vector containing the regularization parameters.}
  \item{cor}{a vector containing the canonical correlations.}
  \item{loadings}{list containing the estimated loadings for the \eqn{X} and 
	\eqn{Y} canonical variates.}
  \item{variates}{list containing the canonical variates.}	
  \item{names}{list containing the names to be used for individuals and variables.} 
}

\references{
Leurgans, S. E., Moyeed, R. A. and Silverman, B. W. (1993). 
Canonical correlation analysis when the data are curves. 
\emph{Journal of the Royal Statistical Society. Series B} \bold{55}, 725-740.

Vinod, H. D. (1976). Canonical ridge and econometrics of joint production. 
\emph{Journal of Econometrics} 6, 129-137.
}

\author{Sébastien Déjean and Ignacio González.}

\seealso{\code{\link{summary}}, \code{\link{estim.regul}}, 
\code{\link{plot.rcc}}, \code{\link{plotIndiv}}, 
\code{\link{plotVar}}, \code{\link{cim}}}

\examples{
## Classic CCA
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.res <- rcc(X, Y)

## Regularized CCA
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, lambda1 = 0.064, lambda2 = 0.008)
}

\keyword{multivariate}

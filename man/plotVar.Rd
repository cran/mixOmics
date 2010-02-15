\name{plotVar}
\encoding{latin1}
\alias{plotVar}
\alias{plotVar.rcc}
\alias{plotVar.pls}
\alias{plotVar.spls}

\title{Plot of Variables}

\description{
This function provides variables representation for (regularized) CCA and (sparse) PLS regression.
}

\usage{
\method{plotVar}{rcc}(object, comp = 1:2, rad.in = 0.5, cutoff = NULL,
        X.label = FALSE, Y.label = FALSE,
        pch = NULL, cex = NULL, col = NULL, font = NULL, \ldots)

\method{plotVar}{pls}(object, comp = 1:2, rad.in = 0.5, 
        X.label = FALSE, Y.label = FALSE, keep.var = FALSE, 
        pch = NULL, cex = NULL, col = NULL, font = NULL, \ldots)	

\method{plotVar}{spls}(object, comp = 1:2, rad.in = 0.5, 
        X.label = FALSE, Y.label = FALSE, keep.var = FALSE, 
        pch = NULL, cex = NULL, col = NULL, font = NULL, \ldots)	 
}

\arguments{
  \item{object}{object of class inheriting from \code{"rcc"}, \code{"pls"} or \code{"spls"}.}
  \item{comp}{integer vector of length two. The components that will be used on the horizontal 
	and the vertical axis respectively to project the variables.}
  \item{rad.in}{numeric between 0 and 1, the radius of the inner circle. Defaults to \code{0.5}.}
  \item{cutoff}{numeric between 0 and 1. Variables with correlations below this cutoff 
    in absolute value are not plotted (see Details).}
  \item{X.label, Y.label}{either a character vector of names for the \eqn{X}- and  
    \eqn{Y}-variables or \code{FALSE} for no names. If \code{TRUE}, the columns names 
	of the matrice are used as labels.}	
  \item{col}{character or integer vector of colors for plotted character and symbols. 
	See Details.}
  \item{pch}{plot character. A vector of single characters 
	or integers. See \code{\link{points}} for all alternatives.}
  \item{cex}{numeric vector of character expansion sizes for the plotted character and symbols.}
  \item{font}{numeric vector of font to be used. See \code{\link{par}} for details.}
  \item{keep.var}{boolean. If \code{TRUE} only the variables with loadings not zero are plotted
    (as selected by \code{spls}). Defaults to \code{FALSE}.}
  \item{\ldots}{not used currently.}
}

\details{
\code{plotVar} produce a "correlation circle", i.e. the correlations 
between each variable and the selected components are plotted as scatter plot,
with concentric circles of radius one et radius given by \code{rad.in}. Each point 
corresponds to a variable. For (regularized) CCA the 
components correspond to the equiangular vector between \eqn{X}- and \eqn{Y}-variates. 
For (sparse) PLS regression mode the components correspond to the \eqn{X}-variates. If mode is
canonical, the components for \eqn{X} and \eqn{Y} variables correspond to
the \eqn{X}- and \eqn{Y}-variates respectively. 

The arguments \code{col}, \code{pch}, \code{cex} and \code{font} can be either vectors of
length two or a list with two vector components of length \eqn{p} and \eqn{q} respectively, 
where \eqn{p} is the number of \eqn{X}-variables and \eqn{q}
is the number of \eqn{Y}-variables. In the first case, the first and second component of the
vector determine the graphics attributes for the \eqn{X}- and \eqn{Y}-variables respectively. 
Otherwise, multiple arguments values can be specified so that each point (variable)
can be given its own graphic attributes. In this case, the first component of the list 
correspond to the \eqn{X} attributs and the second component correspond to 
the \eqn{Y} attributs. Default values exist for this arguments.  
}

\author{Sébastien Déjean and Ignacio González.}

\seealso{\code{\link{plot3dVar}}, \code{\link{cim}}, \code{\link{network}}, \code{\link{par}}.}

\examples{
## variable representation for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

plotVar(nutri.res) #(default)

plotVar(nutri.res, comp = 1:2, cutoff = 0.5, 
        X.label = TRUE, Y.label = TRUE)

## variable representation for objects of class 'pls' or 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
	
plotVar(toxicity.spls, keep.var = TRUE, Y.label = TRUE)	
}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}

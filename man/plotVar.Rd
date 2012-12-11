\name{plotVar}
\encoding{latin1}
\alias{plotVar}
\alias{plotVar.rcc}
\alias{plotVar.pls}
\alias{plotVar.spls}
\alias{plotVar.plsda}
\alias{plotVar.splsda}
\alias{plotVar.pca}
\alias{plotVar.spca}

\title{Plot of Variables}

\description{
This function provides variables representation for (regularized) CCA, (sparse) PLS regression and PCA.
}

\usage{
\method{plotVar}{rcc}(object, comp = 1:2, rad.in = 0.5, cutoff = NULL,
        X.label = FALSE, Y.label = FALSE,
        pch = NULL, cex = NULL, col = NULL, font = NULL, \ldots)

\method{plotVar}{pls}(object, comp = 1:2, rad.in = 0.5, 
        X.label = FALSE, Y.label = FALSE, pch = NULL, cex = NULL, 
        col = NULL, font = NULL, \ldots)

\method{plotVar}{plsda}(object, comp = 1:2, rad.in = 0.5, 
        var.label = FALSE, pch = NULL, cex = NULL, col = NULL, 
        font = NULL, \ldots)		

\method{plotVar}{spls}(object, comp = 1:2, rad.in = 0.5, 
        X.label = FALSE, Y.label = FALSE, pch = NULL, cex = NULL, 
        col = NULL, font = NULL, \ldots)	 

\method{plotVar}{splsda}(object, comp = 1:2, rad.in = 0.5, 
        var.label = FALSE, pch = NULL, cex = NULL, col = NULL, 
        font = NULL, \ldots)

\method{plotVar}{pca}(object, comp = 1:2, rad.in = 0.5, 
        var.label = FALSE, \ldots)

\method{plotVar}{spca}(object, comp = 1:2, rad.in = 0.5, 
        var.label = FALSE, pch = NULL, cex = NULL, col = NULL, 
        font = NULL, \ldots)	
}

\arguments{
  \item{object}{object of class inheriting from \code{"rcc"}, \code{"pls"}, \code{"plsda"}, 
    \code{"spls"}, \code{"splsda"}, \code{"pca"} or \code{"spca"}.}
  \item{comp}{integer vector of length two. The components that will be used on the horizontal 
	and the vertical axis respectively to project the variables.}
  \item{rad.in}{numeric between 0 and 1, the radius of the inner circle. Defaults to \code{0.5}.}
  \item{cutoff}{numeric between 0 and 1. Variables with correlations below this cutoff 
    in absolute value are not plotted (see Details).}
  \item{X.label, Y.label, var.label}{either a character vector of names for the variables 
    or \code{FALSE} for no names. If \code{TRUE}, the columns names 
	of the matrice are used as labels.}	
  \item{col}{character or integer vector of colors for plotted character and symbols. 
	See Details.}
  \item{pch}{plot character. A vector of single characters 
	or integers. See \code{\link{points}} for all alternatives.}
  \item{cex}{numeric vector of character expansion sizes for the plotted character and symbols.}
  \item{font}{numeric vector of font to be used. See \code{\link{par}} for details.}
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

For \code{plsda} and \code{splsda} objects, only the \eqn{X} variables are represented.

For \code{spls} and \code{splsda} objects, only the \eqn{X} and \eqn{Y} variables selected on 
dimensions \code{comp} are represented.

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

\value{
A list containing the following components:
  \item{coord.X}{matrix of \eqn{X}-variables coordinates.}
  \item{coord.Y}{matrix of \eqn{Y}-variables coordinates.}
}

\author{Sébastien Déjean and Ignacio González.}

\seealso{\code{\link{plot3dVar}}, \code{\link{cim}}, \code{\link{network}}, \code{\link{par}} and http://www.math.univ-toulouse.fr/~biostat/mixOmics/ for more details.}

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
	
plotVar(toxicity.spls, keep.var = TRUE, Y.label = TRUE, cex = c(1,0.8))	

## variable representation for objects of class 'splsda'
\dontrun{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- as.factor(liver.toxicity$treatment[, 4])

ncomp <- 2
keepX <- rep(20, ncomp)

splsda.liver <- splsda(X, Y, ncomp = ncomp, keepX = keepX)
plotVar(splsda.liver, var.label = FALSE)
}
}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}

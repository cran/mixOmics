\name{plot3dVar}
\encoding{latin1}
\alias{plot3dVar}
\alias{plot3dVar.rcc}
\alias{plot3dVar.pls}
\alias{plot3dVar.spls}

\title{Plot of Variables in three dimensions}

\description{
This function provides variables representation in three dimensions for 
(regularized) CCA and (sparse) PLS regression.
}

\usage{
\method{plot3dVar}{rcc}(object, comp = 1:3, rad.in = 0.5, 
          cutoff = NULL, X.label = FALSE, Y.label = FALSE, 
          pch = NULL, cex = NULL, col = NULL, font = NULL,
          axes.box = "all", label.axes.box = "both",		
          xlab = NULL, ylab = NULL, zlab = NULL, \ldots)

\method{plot3dVar}{pls}(object, comp = 1:3, rad.in = 0.5, 
          keep.var = FALSE, X.label = FALSE, Y.label = FALSE,
          pch = NULL, cex = NULL, col = NULL, font = NULL, 
          axes.box = "all", label.axes.box = "both",		
          xlab = NULL, ylab = NULL, zlab = NULL, \ldots)

\method{plot3dVar}{spls}(object, comp = 1:3, rad.in = 0.5, 
          keep.var = FALSE, X.label = FALSE, Y.label = FALSE,
          pch = NULL, cex = NULL, col = NULL, font = NULL, 
          axes.box = "all", label.axes.box = "both",		
          xlab = NULL, ylab = NULL, zlab = NULL, \ldots)	 
}

\arguments{
  \item{object}{object of class inheriting from \code{"rcc"}, \code{"pls"} or \code{"spls"}.}
  \item{comp}{an integer vector of length 3, the components that will be used on the \eqn{x}-axis,  
	the \eqn{y}-axis and \eqn{z}-axis respectively to project the variables.}
  \item{rad.in}{numeric between 0 and 1, the radius of the small sphere. Defaults to \code{0.5}.}
  \item{cutoff}{numeric between 0 and 1. Variables with correlations below this cutoff 
    in absolute value are not plotted (see Details).}
  \item{X.label, Y.label}{either a character vector of names for the \eqn{X}- and  
    \eqn{Y}-variables or \code{FALSE} for no names. If \code{TRUE}, the columns names 
	of the matrice are used as labels.}	
  \item{col}{character or integer vector of colors for plotted character and items. 
	See Details.}
  \item{pch}{character indicating the type of item to plot, possibly vector (see Details).}
  \item{cex}{numeric character (or item) expansion (see Details).}
  \item{font}{numeric font to be used. See \code{\link[rgl]{text3d}} for details.}
  \item{keep.var}{boolean. If \code{TRUE} only the variables with loadings not zero are plotted
    (as selected by \code{spls}). Defaults to \code{FALSE}.}
  \item{axes.box}{character string or vector specifying the axes box type. Should be a subset 
    of \code{c("axes", "box", "bbox", "all")}.}
  \item{label.axes.box}{character string indicating whether to label the axes and/or boxes.
    Should be a of \code{c("axes", "box", "both")}}	
  \item{xlab, ylab, zlab}{\eqn{x}-, \eqn{y}- and \eqn{z}-axis titles.}	
  \item{\ldots}{further title parameters are passed to \code{\link[rgl]{title3d}}.}
}

\details{
\code{plot3dVar} gives an extension of the "correlation circle" to three dimensions, i.e. the correlations 
between each variable and the selected components are plotted as scatter plot in three dimensions,
with a sphere of radius given by \code{rad.in} centred in the origin. If \code{cutoff} is not \code{NULL} 
the \code{rad.in} is ignored and the radius of the sphere is equal to \code{cutoff}. For (sparse) PLS 
regression the components correspond to
the \eqn{X}-variates components and for (regularized) CCA the components correspond 
to the bisector vector components between \eqn{X}- and \eqn{Y}-variates. 

The arguments \code{cex}, \code{pch} and \code{font} are vectors of
length two. The firts and second component value determines the graphical attribute for 
\eqn{X}- and \eqn{Y}-variables respectively. Default values exist for this arguments.
Supported types for \code{pch} are: \code{"s"} for spheres, \code{"t"} for tetrahedrons, \code{"c"} for cubes, 
\code{"o"} for octahedrons, \code{"i"} for icosahedrons and \code{"d"} dodecahedrons. Partial matching is used. 

The argument \code{col} can be either a vector of
length two or a list with two vector components of length \eqn{p} and \eqn{q} respectively, 
where \eqn{p} is the number of \eqn{X}-variables and \eqn{q}
is the number of \eqn{Y}-variables. In the first case, the first and second component of the
vector determine the color for the \eqn{X}- and \eqn{Y}-variables respectively. 
Otherwise, multiple colors can be specified so that each item (variable)
can be given its own color. In this case, the first component of the list 
correspond to the \eqn{X}-color attribut and the second component correspond to 
the \eqn{Y}-color attribut.  
}

\author{Ignacio González.}

\seealso{\code{\link{plotVar}}, \code{\link{plot3dIndiv}}, \code{\link[rgl]{text3d}}, 
\code{\link[rgl]{title3d}} and \code{\link[rgl]{rgl.postscript}} to save the screen shot to a file in PostScript 
or other vector graphics format.}

\examples{
require(rgl)

\dontrun{
## 3D variable representation for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, 0.064, 0.008)

# default
plot3dVar(nutri.res)

# cutoff active, labeling the variables
plot3dVar(nutri.res, cutoff = 0.7, X.label = TRUE, cex = c(0.8, 0.8))

## 3D variable representation for objects of class 'pls' or 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls.1 <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                        keepY = c(10, 10, 10))
        
plot3dVar(toxicity.spls.1, rad.in = 0.5, keep.var = TRUE, 
    cex = c(1, 0.8), main = "Variables 3D representation") 

toxicity.spls.2 <- spls(X, Y, ncomp = 3, keepX = c(10, 10, 10), 
                        keepY = c(10, 10, 10))
					  
plot3dVar(toxicity.spls.2, rad.in = 0.5, keep.var = TRUE, 
    Y.label = TRUE, main = "Variables 3D representation", 
    label.axes.box = "axes")
}
}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}

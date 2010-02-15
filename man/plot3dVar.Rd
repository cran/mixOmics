\name{plot3dVar}
\encoding{latin1}
\alias{plot3dVar}
\alias{plot3dVar.rcc}
\alias{plot3dVar.pls}
\alias{plot3dVar.spls}
\alias{plot3dVar.plsda}
\alias{plot3dVar.splsda}
\alias{plot3dVar.pca}

\title{Plot of Variables in three dimensions}

\description{
This function provides variables representation in three dimensions for 
(regularized) CCA, (sparse) PLS regression and PCA.
}

\usage{
\method{plot3dVar}{rcc}(object, comp = 1:3, rad.in = 1, 
          cutoff = NULL, X.label = FALSE, Y.label = FALSE, 
          pch = NULL, cex = NULL, col = NULL, font = NULL,
          axes.box = "all", label.axes.box = "both",		
          xlab = NULL, ylab = NULL, zlab = NULL, \ldots)

\method{plot3dVar}{pls}(object, comp = 1:3, rad.in = 1, 
          X.label = FALSE, Y.label = FALSE,
          pch = NULL, cex = NULL, col = NULL, font = NULL, 
          axes.box = "all", label.axes.box = "both",		
          xlab = NULL, ylab = NULL, zlab = NULL, \ldots)

\method{plot3dVar}{spls}(object, comp = 1:3, rad.in = 1, 
          X.label = FALSE, Y.label = FALSE,
          pch = NULL, cex = NULL, col = NULL, font = NULL, 
          axes.box = "all", label.axes.box = "both",		
          xlab = NULL, ylab = NULL, zlab = NULL, \ldots)

\method{plot3dVar}{plsda}(object, comp = 1:3, rad.in = 1, 
          X.label = FALSE, 
          pch = NULL, cex = NULL, col = NULL, font = NULL, 
          axes.box = "all", label.axes.box = "both",		
          xlab = NULL, ylab = NULL, zlab = NULL, \ldots)

\method{plot3dVar}{splsda}(object, comp = 1:3, rad.in = 1, 
          X.label = FALSE, 
          pch = NULL, cex = NULL, col = NULL, font = NULL, 
          axes.box = "all", label.axes.box = "both",		
          xlab = NULL, ylab = NULL, zlab = NULL, \ldots)	 
		  
\method{plot3dVar}{pca}(object, comp = 1:3, var.label = FALSE,
          pch = NULL, cex = NULL, col = NULL, font = NULL, 
          axes.box = "all", label.axes.box = "both",		
          xlab = NULL, ylab = NULL, zlab = NULL, \ldots)
}

\arguments{
  \item{object}{object of class inheriting from \code{"rcc"}, \code{"pls"}, \code{"spls"} or \code{"pca"}.}
  \item{comp}{an integer vector of length 3, the components that will be used on the \eqn{x}-axis,  
	the \eqn{y}-axis and \eqn{z}-axis respectively to project the variables.}
  \item{rad.in}{numeric between 0 and 1, the radius of the small sphere. Defaults to \code{1}.}
  \item{cutoff}{numeric between 0 and 1. Variables with correlations below this cutoff 
    in absolute value are not plotted (see Details).}
  \item{X.label, Y.label, var.label}{either a character vector of names for the variables 
    or \code{FALSE} for no names. If \code{TRUE}, the columns names 
	of the matrice are used as labels.}	
  \item{col}{character or integer vector of colors for plotted character and items. 
	See Details.}
  \item{pch}{character indicating the type of item to plot, possibly vector (see Details).}
  \item{cex}{numeric character (or item) expansion (see Details).}
  \item{font}{numeric font to be used. See \code{\link[rgl]{text3d}} for details.}
  \item{axes.box}{character string or vector specifying the axes box type. Should be a subset 
    of \code{c("axes", "box", "bbox", "all")}.}
  \item{label.axes.box}{character string indicating whether to label the axes and/or boxes.
    Should be a of \code{c("axes", "box", "both")}}	
  \item{xlab, ylab, zlab}{\eqn{x}-, \eqn{y}- and \eqn{z}-axis titles.}	
  \item{\ldots}{further title parameters are passed to \code{\link[rgl]{title3d}}.}
}

\details{
\code{plot3dVar} gives an extension of the "correlation circle" to three dimensions, i.e. the correlations 
between each variable and the selected components are plotted as scatter plot in three dimensions.

For the objects of class \code{rcc}, \code{pls} and \code{spls} a sphere of radius given by \code{rad.in} 
centred in the origin is displayed. If \code{cutoff} is not \code{NULL} 
the \code{rad.in} is ignored and the radius of the sphere is equal to \code{cutoff}.


For \code{plsda} and \code{splsda} objects, only the X variables are represented.

For \code{spls} and \code{splsda} objects, only the X and Y variables selected on dimensions \code{comp} are represented.


The arguments \code{cex}, \code{pch} and \code{font} are vectors of
length two for the objects of class \code{rcc}, \code{pls} and \code{spls}. 
The firts and second component value determines the graphical attribute for 
\eqn{X}- and \eqn{Y}-variables respectively. Default values exist for this arguments.
Supported types for \code{pch} are: \code{"s"} for spheres, \code{"t"} for tetrahedrons, \code{"c"} for cubes, 
\code{"o"} for octahedrons, \code{"i"} for icosahedrons and \code{"d"} dodecahedrons. 

The argument \code{col} can be either a vector of
length two or a list with two vector components of length \eqn{p} and \eqn{q} respectively, 
where \eqn{p} is the number of \eqn{X}-variables and \eqn{q}
is the number of \eqn{Y}-variables. In the first case, the first and second component of the
vector determine the color for the \eqn{X}- and \eqn{Y}-variables respectively. 
Otherwise, multiple colors can be specified so that each item (variable)
can be given its own color. In this case, the first component of the list 
correspond to the \eqn{X}-color attribut and the second component correspond to 
the \eqn{Y}-color attribute.  

The pointing device of your graphics user-interface can also be used to set the viewpoint interactively. With the pointing device the buttons are by default set as follows:
- left adjust viewpoint position
- middle adjust field of view angle
- right or wheel adjust zoom factor 
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
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

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
	
## 3D variable representation for objects of class 'pca'	   
data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)

plot3dVar(pca.res)

## variable representation for objects of class 'splsda'
data(liver.toxicity)
X = as.matrix(liver.toxicity$gene)
Y = as.factor(liver.toxicity$treatment[,4])

keepX = 20
ncomp = 3

splsda.liver = splsda(X, Y, ncomp = ncomp, keepX = c(rep(keepX, ncomp)), mode = 'regression')
plot3dVar(splsda.liver, X.label = FALSE, Y.label = TRUE, keep.var = TRUE)

}
}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}

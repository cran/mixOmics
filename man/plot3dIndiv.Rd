\name{plot3dIndiv}
\encoding{latin1}
\alias{plot3dIndiv}
\alias{plot3dIndiv.rcc}
\alias{plot3dIndiv.pls}
\alias{plot3dIndiv.spls}
\alias{plot3dIndiv.plsda}
\alias{plot3dIndiv.splsda}
\alias{plot3dIndiv.pca}
\alias{plot3dIndiv.spca}

\title{Plot of Individuals (Experimental Units) in three dimensions}

\description{
This function provides scatter plots of individuals (experimental units)
in three dimensions for (r)CCA, (s)PLS regression and PCA.
}

\usage{
\method{plot3dIndiv}{rcc}(object, comp = 1:3, ind.names = FALSE,
            rep.space = "XY-variate",
            xlab = NULL, ylab = NULL, zlab = NULL,
            col = "blue", cex = 1, pch = "s", font = 1,
            axes.box = "box", \ldots)
	  
\method{plot3dIndiv}{pls}(object, comp = 1:3, ind.names = FALSE,
            rep.space = "X-variate",
            xlab = NULL, ylab = NULL, zlab = NULL,
            col = "blue", cex = 1, pch = "s", font = 1,
            axes.box = "box", \ldots)

\method{plot3dIndiv}{splsda}(object, comp = 1:3, ind.names = FALSE,
            rep.space = "X-variate",
            xlab = NULL, ylab = NULL, zlab = NULL,
            col = "blue", cex = 1, pch = "s", font = 1,
            axes.box = "box", \ldots)

\method{plot3dIndiv}{plsda}(object, comp = 1:3, ind.names = FALSE,
            rep.space = "X-variate",
            xlab = NULL, ylab = NULL, zlab = NULL,
            col = "blue", cex = 1, pch = "s", font = 1,
            axes.box = "box", \ldots)

\method{plot3dIndiv}{spls}(object, comp = 1:3, ind.names = FALSE,
            rep.space = "X-variate",
            xlab = NULL, ylab = NULL, zlab = NULL,
            col = "blue", cex = 1, pch = "s", font = 1,
            axes.box = "box", \ldots)

\method{plot3dIndiv}{pca}(object, comp = 1:3, ind.names = FALSE,
            xlab = NULL, ylab = NULL, zlab = NULL,
            col = "blue", cex = 1, pch = "s", font = 1,
            axes.box = "box", \ldots)

\method{plot3dIndiv}{spca}(object, comp = 1:3, ind.names = FALSE,
            xlab = NULL, ylab = NULL, zlab = NULL,
            col = "blue", cex = 1, pch = "s", font = 1,
            axes.box = "box", \ldots)			
}

\arguments{
  \item{object}{object of class inheriting from \code{"rcc"}, \code{"pls"}, \code{"plsda"}, 
    \code{"spls"}, \code{"plsda"} or \code{"pca"}.}
  \item{comp}{an integer vector of length 3, the components that will be used on the \eqn{x}-axis,  
	the \eqn{y}-axis and \eqn{z}-axis respectively to project the individuals.}
  \item{ind.names}{either a character vector of names for the individuals to be plotted, 
    or \code{FALSE} for no names. If \code{TRUE}, the row names of the first (or second) 
	data matrix is used as names (see Details).}	
  \item{rep.space}{character string, (partially) matching one of \code{"X-variate"}, 
    \code{"Y-variate"} or \code{"XY-variate"}, determining the subspace to project the 
	  individuals. Defaults to \code{"XY-variate"} and \code{"X-variate"} for (r)CCA and (s)PLS 
      respectively.}
  \item{xlab, ylab, zlab}{\eqn{x}-, \eqn{y}- and \eqn{z}-axis titles.}
  \item{col}{character (or item) color to be used, possibly vector (see Details).}
  \item{cex}{numeric character (or item) expansion, possibly vector (see Details).}
  \item{pch}{character indicating the type of item to plot, possibly vector (see Details).}
  \item{font}{character font to be used, possibly vector (see Details).}
  \item{axes.box}{character string specifying the axes box type. Should be a 
    of \code{c("box",} \code{"bbox",} \code{"both")}.}
  \item{\dots}{further title parameters are passed to \code{\link[rgl]{title3d}}.}
}

\details{
\code{plot3dIndiv} method makes scatter plot for individuals representation in three dimensions. 
Each item corresponds to an individual.

If \code{ind.names=TRUE} and row names is \code{NULL}, then \code{ind.names=1:n}, where
\code{n} is the number of individuals. 

The arguments \code{col}, \code{cex}, \code{pch} and \code{font} can be atomic vectors or vectors of
length \code{n}. If atomic, this argument value determines the graphical attribute for all the individuals. 
In the last case, multiple arguments values can be specified so that each item (individual)
can be given its own graphic attributes. Default values exist for this arguments. 
Supported types for \code{pch} are: \code{"s"} for spheres, \code{"t"} for tetrahedrons, \code{"c"} for cubes, 
\code{"o"} for octahedrons, \code{"i"} for icosahedrons and \code{"d"} dodecahedrons.    

The pointing device of your graphics user-interface can also be used to set the viewpoint interactively. With the pointing device the buttons are by default set as follows:
- left adjust viewpoint position
- middle adjust field of view angle
- right or wheel adjust zoom factor 
}

\author{Ignacio González}

\seealso{\code{\link{plotIndiv}}, \code{\link{plot3dVar}}, \code{\link[rgl]{text3d}}, 
\code{\link[rgl]{title3d}} and \code{\link[rgl]{rgl.postscript}} to save the screen shot 
to a file in PostScript or other vector graphics format.}

\examples{
require(rgl)

\dontrun{
## plot3d of individuals for objects of class 'rcc' 
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

col = nutrimouse$diet
font = c(rep(1, 20), rep(3, 20))
plot3dIndiv(nutri.res, ind.names = nutrimouse$diet, 
                axes.box = "box", font = font, col = col)
				
pch = c(rep("s", 20), rep("t", 20))
plot3dIndiv(nutri.res, ind.names = FALSE, axes.box = "both", 
                col = col, cex = 1.5, pch = pch)

## plot3d of individuals for objects of class 'pls' or 'spls'	   
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
		  	
Time.Group = liver.toxicity$treatment[, "Time.Group"]				  
col <- rep(c("blue", "red", "darkgreen", "darkviolet"), rep(16, 4))
plot3dIndiv(toxicity.spls, ind.names = Time.Group, 
                 col = col, cex = 0.8)		  
		  				  
col <- rainbow(48)[Time.Group]
plot3dIndiv(toxicity.spls, ind.names = FALSE, 
                 col = col, cex = 0.3, axes.box = "both")	

## plot3d of individuals for objects of class 'pca'	   
data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)

palette(rainbow(9))
col = as.numeric(as.factor(multidrug$cell.line$Class))
plot3dIndiv(pca.res, cex = 0.25, col = col)
palette("default")				 
}
}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}

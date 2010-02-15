\name{plotIndiv}
\encoding{latin1}
\alias{plotIndiv}
\alias{plotIndiv.rcc}
\alias{plotIndiv.pls}
\alias{plotIndiv.spls}
\alias{plotIndiv.plsda}
\alias{plotIndiv.splsda}

\title{Plot of Individuals (Experimental Units)}

\description{
This function provides scatter plots for individuals (experimental units)
representation in (r)CCA and (s)PLS regression.
}

\usage{
\method{plotIndiv}{rcc}(object, comp = 1:2, ind.names = TRUE,
          rep.space = "XY-variate",
          x.label = NULL, y.label = NULL,
          col = "black", cex = 1, pch = 1, \ldots)
	  
\method{plotIndiv}{pls}(object, comp = 1:2, ind.names = TRUE,
          rep.space = "X-variate",
          x.label = NULL, y.label = NULL,  
          col = "black", cex = 1, pch = 1, \ldots)	

\method{plotIndiv}{spls}(object, comp = 1:2, ind.names = TRUE,
          rep.space = "X-variate",
          x.label = NULL, y.label = NULL,  
          col = "black", cex = 1, pch = 1, \ldots)		  
}

\arguments{
  \item{object}{object of class inheriting from \code{"rcc"}, \code{"pls"} or \code{"spls"}.}
  \item{comp}{integer vector of length two. The components that will be used on the horizontal 
	and the vertical axis respectively to project the individuals.}
  \item{ind.names}{either a character vector of names for the individuals to be plotted, 
    or \code{FALSE} for no names. If \code{TRUE}, the row names of the first (or second) 
	data matrix is used as names (see Details).}
  \item{rep.space}{character string, (partially) matching one of \code{"X-variate"}, 
    \code{"Y-variate"} or \code{"XY-variate"}, determining the subspace to project the 
	  individuals. Defaults to \code{"XY-variate"} and \code{"X-variate"} for (r)CCA and (s)PLS 
      respectivvely.}
  \item{x.label, y.label}{\eqn{x}- and \eqn{y}-axis titles.}
  \item{col}{character (or symbol) color to be used, possibly vector.}
  \item{cex}{numeric character (or symbol) expansion, possibly vector.}
  \item{pch}{plot character. A character string or a vector of single characters 
	or integers. See \code{\link{points}} for all alternatives.}
  \item{\dots}{further graphical parameters are passed 
	to \code{\link{text}}.}
}

\details{
\code{plotIndiv} method makes scatter plot for individuals representation  
depending on the subspace of projection. Each point corresponds to an individual.

If \code{ind.names=TRUE} and row names is \code{NULL}, then \code{ind.names=1:n}, where
\code{n} is the number of individuals. 

The arguments \code{col}, \code{cex} and \code{pch} can be atomic vectors or vectors of
length \code{n}. If atomic, this argument value determines the graphical attribute for all the individuals. 
In the last case, multiple arguments values can be specified so that each point (individual)
can be given its own graphic attributes (see \code{\link{par}}). Default values exist for this arguments. 
}

\author{Sébastien Déjean and Ignacio González.}

\seealso{\code{\link{plot3dIndiv}}, \code{\link{text}}, \code{\link{points}}.}


\examples{
## plot of individuals for objects of class 'rcc' 
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

plotIndiv(nutri.res) #(default) 

col <- rep(c("blue", "red"), c(20, 20))
plotIndiv(nutri.res, ind.names = nutrimouse$diet, col = col)
legend(-2.2, -1.1, c("WT", "PPAR"), pch = c(16, 16), 
       col = c("blue", "red"), text.col = c("blue", "red"),
       cex = 1, pt.cex = c(1.2, 1.2))

## plot of individuals for objects of class 'pls' or 'spls'	   
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
		  
col <- rep(c("blue", "red", "darkgreen", "darkviolet"), c(16, 16, 16, 16))
cex <- rep(c(1, 1.2, 1, 1.4), c(16, 16, 16, 16))
pch <- rep(c(15, 16, 17, 18), c(16, 16, 16, 16))
plotIndiv(toxicity.spls, comp = 1:2, ind.names = FALSE,
          rep.space = "X-variate", col = col, cex = cex, pch = pch)
legend(0.22, -0.24, c("50 mg/kg", "150 mg/kg", "1500 mg/kg", "2000 mg/kg"), 
       col = c("blue", "red", "darkgreen", "darkviolet"), 
       pch = c(15, 16, 17, 18), pt.cex = c(1, 1.2, 1, 1.4), 
       title = "Treatment")
}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}

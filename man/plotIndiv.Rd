\name{plotIndiv}
\encoding{latin1}
\alias{plotIndiv}
\alias{plotIndiv.rcc}
\alias{plotIndiv.pls}
\alias{plotIndiv.spls}
\alias{plotIndiv.plsda}
\alias{plotIndiv.splsda}
\alias{plotIndiv.pca}
\alias{plotIndiv.spca}
\alias{plotIndiv.rgcca}
\alias{plotIndiv.sgcca}

\title{Plot of Individuals (Experimental Units)}

\description{
This function provides scatter plots for individuals (experimental units)
representation in (regularized) CCA, (sparse) PLS regression, (sparse)RGCCA and PCA.
}

\usage{
\method{plotIndiv}{rcc}(object, comp = 1:2, ind.names = TRUE,
          rep.space = "XY-variate",
          X.label = NULL, Y.label = NULL,
          col = "black", cex = 1, pch = 1, 
          abline.line = FALSE, \ldots)
	  
\method{plotIndiv}{pls}(object, comp = 1:2, 
			ind.names = TRUE,
          rep.space = "X-variate",
          X.label = NULL, Y.label = NULL,  
          col = "black", cex = 1, pch = 1, 
          abline.line = FALSE,\ldots)	

\method{plotIndiv}{spls}(object, comp = 1:2, ind.names = TRUE,
          rep.space = "X-variate",
          X.label = NULL, Y.label = NULL,  
          col = "black", cex = 1, pch = 1, 
          abline.line = FALSE,\ldots)

\method{plotIndiv}{pca}(object, comp = 1:2, ind.names = TRUE,
          X.label = NULL, Y.label = NULL, 
          abline.line = FALSE,\ldots)

\method{plotIndiv}{spca}(object, comp = 1:2, ind.names = TRUE,
          X.label = NULL, Y.label = NULL, 
          abline.line = FALSE,\ldots)

\method{plotIndiv}{rgcca}(object, comp = 1:2, ind.names = TRUE, rep.space = 1,
    X.label = NULL, Y.label = NULL,  col = "black", cex = 1, pch = 1, 
    abline.line = FALSE,\ldots)

\method{plotIndiv}{sgcca}(object, comp = 1:2, ind.names = TRUE, rep.space = 1,
    X.label = NULL, Y.label = NULL,  col = "black", cex = 1, pch = 1, 
    abline.line = FALSE, \ldots)

		  
}

\arguments{
  \item{object}{object of class inheriting from \code{"rcc"}, \code{"pls"}, \code{"plsda"}, 
    \code{"spls"}, \code{"plsda"}, \code{"pca"} or \code{"spca"}.}
  \item{comp}{integer vector of length two. The components that will be used on the horizontal 
	and the vertical axis respectively to project the individuals.}
  \item{ind.names}{either a character vector of names for the individuals to be plotted, 
    or \code{FALSE} for no names. If \code{TRUE}, the row names of the first (or second) 
	data matrix is used as names (see Details).}
  \item{rep.space}{For objects of class \code{"rcc"}, \code{"pls"}, \code{"spls"}, character string, (partially) matching one of \code{"X-variate"}, 
    \code{"Y-variate"} or \code{"XY-variate"}, determining the subspace to project the 
	  individuals. Defaults to \code{"XY-variate"} \code{"rcc"} object for  and \code{"X-variate"} \code{"pls"} and \code{"spls"} objects.
For objects of class \code{"rgcca"} and \code{"sgcca"}, numerical value indicating the block data set form which to represent the individuals.}
  \item{X.label, Y.label}{\eqn{x}- and \eqn{y}-axis titles.}
  \item{col}{character (or symbol) color to be used, possibly vector.}
  \item{cex}{numeric character (or symbol) expansion, possibly vector.}
  \item{pch}{plot character. A character string or a vector of single characters 
	or integers. See \code{\link{points}} for all alternatives.}
\item{abline.line}{should the vertical and horizontal line through the center be plotted?default set to \code{FALSE}}
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

\author{S?bastien D?jean and Ignacio Gonz?lez.}

\seealso{\code{\link{plot3dIndiv}}, \code{\link{text}}, \code{\link{points}} and http://www.math.univ-toulouse.fr/~biostat/mixOmics/ for more details.}


\examples{
## plot of individuals for objects of class 'rcc' 
# ----------------------------------------------------
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
# ----------------------------------------------------   
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
		  
col <- rep(c("blue", "red", "darkgreen", "darkviolet"), rep(16, 4))
cex <- rep(c(1, 1.2, 1, 1.4), c(16, 16, 16, 16))
pch <- rep(c(15, 16, 17, 18), c(16, 16, 16, 16))
plotIndiv(toxicity.spls, comp = 1:2, ind.names = FALSE,
          rep.space = "X-variate", col = col, cex = cex, pch = pch)
legend("topright", c("50 mg/kg", "150 mg/kg", "1500 mg/kg", "2000 mg/kg"), 
       col = c("blue", "red", "darkgreen", "darkviolet"), 
       pch = c(15, 16, 17, 18), pt.cex = c(1, 1.2, 1, 1.4), 
       title = "Treatment")
       
## variable representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------
data(nutrimouse)

# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

#note: the penalty parameters will need to be tuned
wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)


# on the first data set
plotIndiv(wrap.result.sgcca, rep.space = 1, ind.names = TRUE, 
          col = as.numeric(nutrimouse$diet), cex = .6)
# on the second data set
plotIndiv(wrap.result.sgcca, rep.space = 2, ind.names = TRUE, 
          col = as.numeric(nutrimouse$diet), cex = .6)


}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}

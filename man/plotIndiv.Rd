\name{plotIndiv}
\encoding{latin1}
\alias{plotIndiv}

\title{Plot of Individuals (Experimental Units)}

\description{
This function provides scatter plots for individuals (experimental units)
representation in (sparse)(I)PCA, (regularized)CCA, (sparse)PLS(DA) and (sparse)(R)GCCA(DA).
}

\usage{
plotIndiv(object,
           comp = c(1, 2),
           ind.names = TRUE,
           rep.space = "X-variate",
           blocks = NULL, 
           X.label = NULL,
           Y.label = NULL,
           abline.line = FALSE,
           col.per.group,
           cex,
           pch,
           plot.ellipse = FALSE, 
           ellipse.level = 0.95,
           group = NULL, 
           main="plotIndiv",
           add.legend=FALSE,
           style="ggplot2", 
           \ldots)
}

\arguments{
  \item{object}{object of class inheriting from any \pkg{mixOmics}:  \code{PLS, sPLS, PLS-DA, SPLS-DA, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA, sGCCDA}}, 
  \item{comp}{integer vector of length two. The components that will be used on the horizontal 
  and the vertical axis respectively to project the individuals.}
    \item{ind.names}{either a character vector of names for the individuals to be plotted, 
    or \code{FALSE} for no names. If \code{TRUE}, the row names of the first (or second) 
  data matrix is used as names (see Details).}
  \item{rep.space}{For objects of class \code{"rcc"}, \code{"pls"}, \code{"spls"}, character string, (partially) matching one of \code{"X-variate"}, 
    \code{"Y-variate"} or \code{"XY-variate"}, determining the subspace to project the individuals. Defaults to \code{"XY-variate"} \code{"rcc"} object for  and \code{"X-variate"} \code{"pls"} and \code{"spls"} objects. For objects of class \code{"rgcca"} and \code{"sgcca"}, numerical value indicating the block data set form which to represent the individuals.}

  \item{blocks}{integer value of name of a block to be plotted using the GCCA module. See examples.}  
  \item{X.label}{x axis titles.}
  \item{Y.label}{y axis titles}
  \item{abline.line}{should the vertical and horizontal line through the center be plotted? Default set to \code{FALSE}}
  \item{col.per.group}{character (or symbol) color to be used. If a supervised method is used, a color vector indicating the color per outcome group. }
  \item{cex}{numeric character (or symbol) expansion, possibly vector.}
  \item{pch}{plot character. A character string or a vector of single characters 
  or integers. See \code{\link{points}} for all alternatives.}
  \item{plot.ellipse}{boolean indicating if ellipse plots should be plotted. In the \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA, sGCCA} non supervised objects, by default will include all data points unless the argument \code{group} is set up. In the \code{PLS-DA, SPLS-DA,sGCCDA} supervised object, by default the ellipse will be plotted accoding to the outcome \code{Y}.} 
    \item{ellipse.level}{Numerical value indicating the confidence level of ellipse being plotted when \code{plot.ellipse =TRUE} (i.e. the size of the ellipse). The default is set to 0.95, for a 95\% region.}  
  \item{group}{factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the supervised methods \code{PLS-DA, SPLS-DA,sGCCDA}, but needs to be input for the unsupervised methods \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA, sGCCA}. By default it will assign the \code{color.mixo}.}    
  \item{main}{character indicating the title plot.}
  \item{add.legend}{boolean. Whether the legend should be added. Default is TRUE.}
  \item{style}{argument to be set to either \code{'graphics'}, \code{'lattice'} or \code{'ggplot2'} for a style of plotting.} 
  \item{\dots}{further graphical parameters are passed 
	to \code{\link{text}}.}
}

\details{
\code{plotIndiv} method makes scatter plot for individuals representation  
depending on the subspace of projection. Each point corresponds to an individual.

If \code{ind.names=TRUE} and row names is \code{NULL}, then \code{ind.names=1:n}, where
\code{n} is the number of individuals. 

The arguments \code{col.per.group}, \code{cex} and \code{pch} can be atomic vectors or vectors of
length \code{n}. If atomic, this argument value determines the graphical attribute for all the individuals. 
In the last case, multiple arguments values can be specified so that each point (individual)
can be given its own graphic attributes (see \code{\link{par}}). Default values exist for this arguments. 

}

\author{Ignacio Gonzalez, Benoit Gautier, Francois Bartolo, Florian Rohart}

\seealso{\code{\link{plot3dIndiv}}, \code{\link{text}}, \code{\link{points}} and http://mixOmics.org/graphics for more details.}


\examples{
## plot of individuals for objects of class 'rcc' 
# ----------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# default, only in the X space
plotIndiv(nutri.res) 

# ellipse with respect to genotype in the XY space, 
# names also indicate genotype
plotIndiv(nutri.res, rep.space= 'XY-variate', 
  plot.ellipse = TRUE, ellipse.level = 0.9, 
  group = nutrimouse$genotype, ind.names = nutrimouse$genotype)

# ellipse with respect to genotype in the XY space, with legend
plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype, 
  add.legend = TRUE)


# lattice style
plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype, 
  add.legend = TRUE, style = 'lattice')

# classic style, in the Y space
plotIndiv(nutri.res, rep.space= 'Y-variate', group = nutrimouse$genotype, 
  add.legend = TRUE, style = 'graphics')


## plot of individuals for objects of class 'pls' or 'spls'  
# ----------------------------------------------------   
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

#default
plotIndiv(toxicity.spls)

# in the Y space, colors indicate time of necropsy, text is the dose
plotIndiv(toxicity.spls, rep.space= 'Y-variate', 
  group = liver.toxicity$treatment[, 'Time.Group'], 
  ind.names = liver.toxicity$treatment[, 'Dose.Group'], 
  add.legend = TRUE)


## plot of individuals for objects of class 'plsda' or 'splsda'  
# ----------------------------------------------------   
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

splsda.breast <- splsda(X, Y,keepX=c(10,10),ncomp=2)

# default option: note the outcome color is included by default!
plotIndiv(splsda.breast)

# default option with no ind name: pch and color are set automatically
plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2))

# default option with no ind name: pch and color are set automatically, with legend
plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2), add.legend = TRUE)

# trying the different styles
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), plot.indiv = FALSE, 
  plot.ellipse = TRUE, style = "ggplot2", cex = c(1, 1))
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), plot.indiv = FALSE, 
  plot.ellipse = TRUE, style = "lattice", cex = c(1, 1))


## variable representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------
data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
nutrimouse.sgcca <- wrapper.sgcca(blocks = data,
                                  design = design1,
                                  penalty = c(0.3, 0.5, 1),
                                  ncomp = c(2, 2, 3),
                                  scheme = "centroid",
                                  verbose = FALSE, 
                                  bias = FALSE)

# default style: one panel for each block
plotIndiv(nutrimouse.sgcca)

# for the block 'lipid' with ellipse plots and legend, different styles
plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE, 
  plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", main = 'my plot')
plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet, 
  add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", 
  main = 'my plot')
plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, 
  add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", 
  main = 'my plot')


## variable representation for objects of class 'sgccda' 
# ----------------------------------------------------
\dontrun{
# Note: the code differs from above as we use a 'supervised' GCCA analysis
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)

nutrimouse.sgccda1 <- wrapper.sgccda(blocks = data,
                                     Y = Y,
                                     design = design1,
                                     ncomp = c(2, 2),
                                     keep = list(c(10,10), c(15,15)),
                                     scheme = "centroid",
                                     verbose = FALSE,
                                     bias = FALSE)


# plotIndiv
# ----------
# displaying all blocks. bu default colors correspond to outcome Y
plotIndiv(nutrimouse.sgccda1)

# displaying only 2 blocks
plotIndiv(nutrimouse.sgccda1, blocks = c(1,2), group = nutrimouse$diet)

# with some ellipse, legend and title
plotIndiv(nutrimouse.sgccda1, blocks = c(1,2), group = nutrimouse$diet, 
  plot.ellipse = TRUE, add.legend = TRUE, main = 'my sample plot')

}
}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}

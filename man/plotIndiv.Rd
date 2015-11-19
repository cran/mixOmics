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
           comp = NULL,
           rep.space = NULL,
           blocks = NULL, # to choose which block data to plot with GCCA module
           ind.names = TRUE,
           group,  
           col.per.group,
           style="ggplot2", # choose between graphics,3d, lattice or ggplot2
           plot.ellipse = FALSE,
           ellipse.level = 0.95,
           plot.centroid=FALSE,
           plot.star=FALSE,
           main=NULL,
           add.legend=FALSE,
           X.label = NULL,
           Y.label = NULL,
           Z.label = NULL,
           abline.line = FALSE,
           xlim = NULL,
           ylim = NULL,
           col,
           cex,
           pch,
           alpha=0.2,
           axes.box = "box",
           ...
)
}

\arguments{
\item{object}{object of class inheriting from any \pkg{mixOmics}:  \code{PLS, sPLS, PLS-DA, SPLS-DA, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA, sGCCDA}},
\item{comp}{integer vector of length two (or three to 3d). The components that will be used on the horizontal
and the vertical axis respectively to project the individuals.}
\item{rep.space}{For objects of class \code{"rcc"}, \code{"pls"}, \code{"spls"}, character string, (partially) matching one of \code{"X-variate"}, \code{"Y-variate"} ,or \code{"XY-variate"}, determining the subspace to project the individuals. Defaults to \code{"X-variate"} \code{"pca"} object and for \code{"plsda"} objects. For objects of class \code{"pls"} and \code{"rcc"},defaults, the tree subspaces  represent the individuals. For objects of class \code{"rgcca"} and \code{"sgcca"}, numerical value indicating the block data set form which to represent the individuals.}

\item{blocks}{integer value of name of a block to be plotted using the GCCA module. See examples.}
\item{ind.names}{either a character vector of names for the individuals to be plotted,
or \code{FALSE} for no names. If \code{TRUE}, the row names of the first (or second)
data matrix is used as names (see Details).}
\item{group}{factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the supervised methods \code{PLS-DA, SPLS-DA,sGCCDA}, but needs to be input for the unsupervised methods \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA, sGCCA}}
\item{col.per.group}{character (or symbol) color to be used when 'group' is defined. Vector of the same length than the number of groups.}
\item{style}{argument to be set to either \code{'graphics'}, \code{'lattice'}, \code{'ggplot2'} or \code{'3d'} for a style of plotting. Default set to 'ggplot2'. See details.}
\item{plot.ellipse}{boolean indicating if ellipse plots should be plotted. In the non supervised objects \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA, sGCCA} ellipse plot is only be plotted if the argument \code{group} is provided. In the \code{PLS-DA, SPLS-DA,sGCCDA} supervised object, by default the ellipse will be plotted accoding to the outcome \code{Y}.}
\item{ellipse.level}{Numerical value indicating the confidence level of ellipse being plotted when \code{plot.ellipse =TRUE} (i.e. the size of the ellipse). The default is set to 0.95, for a 95\% region.}
\item{plot.centroid}{boolean indicating whether centroid points should be plotted. In the non supervised objects \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA, sGCCA} the centroid will  only be plotted if the argument \code{group} is provided. The centroid will be calculated based on the group categories. In the supervised objects \code{PLS-DA, SPLS-DA,sGCCDA} the centroid will be calculated according to the outcome \code{Y}.}
\item{plot.star}{boolean indicating whether a star plot should be plotted, with arrows starting from the centroid (see argument 
\code{plot.centroid}, and ending for each sample belonging to each group or outcome. In the non supervised objects \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA, sGCCA} star plot is only be plotted if the argument \code{group} is provided. In the supervised objects \code{PLS-DA, SPLS-DA,sGCCDA} the star plot is plotted according to the outcome \code{Y}.}
\item{main}{set of characters indicating the title plot.}
\item{add.legend}{boolean. Whether the legend should be added. Default is FALSE.}
\item{X.label}{x axis titles.}
\item{Y.label}{y axis titles.}
\item{Z.label}{z axis titles (when style = '3d').}
\item{abline.line}{should the vertical and horizontal line through the center be plotted? Default set to \code{FALSE}}
\item{xlim,ylim}{numeric list of vectors of length 2 and length =length(blocks), giving the x and y coordinates ranges.}
\item{col}{character (or symbol) color to be used, possibly vector.}
\item{cex}{numeric character (or symbol) expansion, possibly vector.}
\item{pch}{plot character. A character string or a vector of single characters
or integers. See \code{\link{points}} for all alternatives.}
\item{alpha}{Semi-transparent colors (0 < \code{'alpha'} < 1) }
\item{axes.box}{for style '3d', argument to be set to either \code{'axes'}, \code{'box'}, \code{'bbox'} or \code{'all'}, defining the shape of the box.}
\item{...}{external arguments or type par can be added with \code{style = 'graphics'}}
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

For customized plots (i.e. adding points, text), use the style = 'graphics' (default is ggplot2).

Note: the ellipse options were borrowed from the \pkg{ellipse}.
}

\author{Ignacio Gonzalez, Benoit Gautier, Francois Bartolo, Florian Rohart}

\seealso{\code{\link{text}}, \code{\link{points}} and http://mixOmics.org/graphics for more details.}


\examples{
## plot of individuals for objects of class 'rcc'
# ----------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# default, only in the X space
plotIndiv(nutri.res)

\dontrun{
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
}

## plot of individuals for objects of class 'pls' or 'spls'
# ----------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
keepY = c(10, 10, 10))

#default
plotIndiv(toxicity.spls)


# indicating the centroid
plotIndiv(toxicity.spls, rep.space= 'X-variate', ind.names = FALSE,
          group = liver.toxicity$treatment[, 'Time.Group'], plot.centroid = TRUE)

# indicating the star and centroid
plotIndiv(toxicity.spls, rep.space= 'X-variate', ind.names = FALSE,
          group = liver.toxicity$treatment[, 'Time.Group'], plot.centroid = TRUE, plot.star = TRUE)


# indicating the star and ellipse
plotIndiv(toxicity.spls, rep.space= 'X-variate', ind.names = FALSE,
          group = liver.toxicity$treatment[, 'Time.Group'], plot.centroid = TRUE,
            plot.star = TRUE, plot.ellipse = TRUE)



\dontrun{
    # in the Y space, colors indicate time of necropsy, text is the dose
    plotIndiv(toxicity.spls, rep.space= 'Y-variate',
    group = liver.toxicity$treatment[, 'Time.Group'],
    ind.names = liver.toxicity$treatment[, 'Dose.Group'],
    add.legend = TRUE)
}

## plot of individuals for objects of class 'plsda' or 'splsda'
# ----------------------------------------------------
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

splsda.breast <- splsda(X, Y,keepX=c(10,10),ncomp=2)

# default option: note the outcome color is included by default!
plotIndiv(splsda.breast)

\dontrun{
    # default option with no ind name: pch and color are set automatically
    plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2))

    # default option with no ind name: pch and color are set automatically, with legend
    plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2), add.legend = TRUE)

    # trying the different styles
    plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2),
    plot.ellipse = TRUE, style = "ggplot2", cex = c(1, 1))
    plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2),
    plot.ellipse = TRUE, style = "lattice", cex = c(1, 1))
}

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

\dontrun{
    # for the block 'lipid' with ellipse plots and legend, different styles
    plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE,
    plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", main = 'my plot')
    plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet,
    add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid",
    main = 'my plot')
    plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet,
    add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid",
    main = 'my plot')
}

## variable representation for objects of class 'sgccda'
# ----------------------------------------------------
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

\dontrun{
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
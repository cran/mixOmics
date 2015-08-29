\name{plotVar}
\encoding{latin1}
\alias{plotVar}

\title{Plot of Variables}

\description{
This function provides variable representation through Correlation Circle Plots for sparse)(I)PCA, (regularized)CCA, (sparse)PLS(DA) and (sparse)(R)GCCA(DA).
}

\usage{
plotVar(object,
           comp = c(1, 2),
           comp.select = NULL,
           var.names = NULL,
           blocks = NULL, 
           X.label = NULL,
           Y.label = NULL,
           abline.line = TRUE,
           col,
           cex,
           pch,
           font,
           cutoff = 0,
           rad.in = 0.5,
           main="Correlation Circle Plots",
           style="ggplot2", 
           overlap = TRUE,
           ...)
         
}

\arguments{
  \item{object}{object of class inheriting from \code{"rcc"}, \code{"pls"}, \code{"plsda"}, 
    \code{"spls"}, \code{"splsda"}, \code{"pca"} or \code{"spca"}.}
    
   \item{comp}{integer vector of length two. The components that will be used on the horizontal 
  and the vertical axis respectively to project the variables.}
   
     \item{comp.select}{for the sparse versions, an input vector indicating the components on which the variables were selected. Only those selected variables are displayed. See examples with \code{sgcca}}

  \item{var.names}{a list of character vectors indicating alternative names to be displayed for all variables in all data sets.  Each vector should be the length of the the total number of variables in each data set. By default set to \code{NULL} to display the original variable names from the input data (i.e. columns names 
  of the matrices). If you do not wish to show the variable name, then specify the argument \code{pch}.}	

  \item{blocks}{for an object of class \code{"rgcca"} or \code{"sgcca"}, a numerical vector indicating the blocks of variables to display.}
  
  \item{X.label}{x axis titles, by default set to 'Component'}
  \item{Y.label}{y axis titles, by default set to 'Component'}

  \item{abline.line}{should the vertical and horizontal line through the center be plotted? Default set to \code{TRUE}}

  \item{col}{character or integer vector of colors for plotted character and symbols, can be of length the number of data sets that are integrated, or a list of vectors associated to each data set, where each vector has for length the number of variables in that data set.}

  \item{cex}{numeric vector of character expansion sizes for the plotted character and symbols, can be of length the number of data sets that are integrated, or a list of vectors associated to each data set, where each vector has for length the number of variables in that data set. See \code{\link{par}} for details.}

  \item{pch}{numeric vector for symbols, can be of length the number of data sets that are integrated, or a list of vectors associated to each data set, where each vector has for length the number of variables in that data set. See \code{\link{points}} for all alternatives.}

  \item{font}{numeric vector for font, can be of length the number of data sets that are integrated, or a list of vectors associated to each data set, where each vector has for length the number of variables in that data set. }

  \item{cutoff}{numeric between 0 and 1. Variables with correlations below this cutoff in absolute value are not plotted on any of the components specified in \code{comp}.}

  \item{rad.in}{numeric between 0 and 1, the radius of the inner circle. Defaults to \code{0.5}.}
  
    \item{main}{character indicating the title plot.}
    
      \item{style}{argument to be set to either \code{'graphics'}, \code{'lattice'} or \code{'ggplot2'} for a style of plotting.} 
      
      
  \item{overlap}{boolean indicating whether the correlation circle plots should be idsplayed in an overlap fashion from each data set. By default set to \code{TRUE}.}	

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
You can store the function in an object which will return a matrix containing the \eqn{x}- and \eqn{y}- coordinates of the plotted variables, along with the various grpahical parameters used to generate the plot.
}

\author{Benoit Gautier, Ignacio Gonzalez, Francois Bartolo, Florian Rohart, Kim-Anh Le Cao.}

\references{Gonzalez I., Le Cao K-A., Davis, M.J. and Dejean, S. (2012). Visualising associations between paired 'omics data sets. J. Data Mining 5:19. \url{http://www.biodatamining.org/content/5/1/19/abstract}}


\seealso{\code{\link{plot3dVar}}, \code{\link{cim}}, \code{\link{network}}, \code{\link{par}} and http://www.mixOmics.org for more details.}

\examples{
## variable representation for objects of class 'rcc'
# ----------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

plotVar(nutri.res) #(default)

# playing with the style
plotVar(nutri.res, style = 'lattice') #(default)

# changing x and y labels
plotVar(nutri.res, comp = c(1,3), cutoff = 0.5, 
        X.label = 'PC1', Y.label = 'PC3')

# one correlation circle plot per data set
plotVar(nutri.res, comp = c(1,2), cutoff = 0.5, 
        overlap = FALSE)


# with pch symbols
plotVar(nutri.res, comp = c(1,2), pch = c(16,2))


## variable representation for objects of class 'pls' or 'spls'
# ----------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

# default shows only the variables selected on the plotted components
plotVar(toxicity.spls)

# shows only the variables selected on the plotted components
plotVar(toxicity.spls, comp = c(1,3))

# shows only the variables selected on the selected components
plotVar(toxicity.spls, comp.select = c(1:3))


# change variable names
new.names = list(paste('gene', 1:ncol(X)), paste('clinic', 1:ncol(Y)))
plotVar(toxicity.spls, overlap = FALSE, var.names = new.names)

# prefilter even further and use of pch
plotVar(toxicity.spls, comp.select = c(1:3), cutoff = 0.8, pch = c(15,16))

# change colors
plotVar(toxicity.spls, col = color.mixo(3:4))

my.col = list(c(rep(1, ncol(X))), c(rep(3,ncol(Y))))
plotVar(toxicity.spls, col = my.col)


## variable representation for objects of class 'splsda'
# ----------------------------------------------------
\dontrun{
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- as.factor(liver.toxicity$treatment[, 4])
  
  ncomp <- 2
  keepX <- rep(20, ncomp)
  
  splsda.liver <- splsda(X, Y, ncomp = ncomp, keepX = keepX)
  # use of pch symbols
  plotVar(splsda.liver, pch = 16, col = 3)

}

## variable representation for objects of class 'sgcca' 
# ----------------------------------------------------

## see example in ??wrapper.sgcca
data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
# set up the data as list
data = list(nutrimouse$gene, nutrimouse$lipid,Y)

# set up the design matrix:
# with this design, gene expression and lipids are connected to the diet factor
# design = matrix(c(0,0,1,
#                   0,0,1,
#                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


#note: the penalty parameters will need to be tuned
wrap.result.sgcca = wrapper.sgcca(blocks = data, design = design, penalty = c(.3,.3, 1),
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
#wrap.result.sgcca


# showing 2 blocks, with variables selected on comp 1 for block 1 and comp 1 for block 2
plotVar(wrap.result.sgcca, comp = c(1,2), 
        blocks = c(1,2), comp.select = c(1,1), 
        overlap = FALSE,
        main = 'Variables selected on component 1 only')


# displaying variables selected on comp 2 for block 1 and comp 2 for block 2
plotVar(wrap.result.sgcca, comp = c(1,2), blocks = c(1,2), comp.select = c(2,2), 
        main = 'Variables selected on component 2 only')


## variable representation for objects of class 'rgcca'
# ----------------------------------------------------
\dontrun{
data(nutrimouse)
# need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
Y = unmap(nutrimouse$diet)

data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, 
                byrow = TRUE, dimnames = list(names(data), names(data)))

nutrimouse.rgcca <- wrapper.rgcca(blocks = data,
                                  design = design,
                                  tau = "optimal",
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid",
                                  verbose = FALSE)

# changing cex
plotVar(nutrimouse.rgcca, comp = c(1,2), blocks = c(1,2), cex = c(1.5, 1.5))
# changing font
plotVar(nutrimouse.rgcca, comp = c(1,2), blocks = c(1,2), font = c(1,3))


# set up the data as list
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
# with this design, gene expression and lipids are connected to the diet factor
# design = matrix(c(0,0,1,
#                   0,0,1,
#                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#note: the tau parameter is the regularization parameter
wrap.result.rgcca = wrapper.rgcca(blocks = data, design = design, tau = c(1, 1, 0),
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
#wrap.result.rgcca
plotVar(wrap.result.rgcca, comp = c(1,2), blocks = c(1,2))
}
}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}

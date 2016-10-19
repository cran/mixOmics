\name{circosPlot}
\encoding{latin1}
\alias{circosPlot}

\title{circosPlot for DIABLO}

\description{Displays variable correlation among different blocks
}

\usage{
circosPlot(object,
comp = 1 : min(object$ncomp),
cutoff,
var.names = NULL,
showIntraLinks = FALSE,
line = TRUE,
size.legend = 0.8,
ncol.legend = 1,
size.variables = 0.25,
size.labels = 1)
}

\arguments{
\item{object}{An object of class inheriting from \code{"block.splsda"}.}
\item{comp}{Numeric vector indicating which component to plot. Default to all}
\item{cutoff}{Only shows links with a correlation higher than \code{cutoff}}
\item{var.names}{Optional parameter. A list of length the number of blocks in \code{object$X}, containing the names of the variables of each block.
If \code{NULL}, the colnames of the data matrix are used.}
\item{showIntraLinks}{if TRUE, shows the correlation higher than the threshold inside each block.}
\item{line}{if TRUE, shows the overall expression of the selected variables. see examples.}
\item{size.legend}{size of the legend}
\item{ncol.legend}{number of columns for the legend}
\item{size.variables}{size of the variable labels}
\item{size.labels}{size of the block labels}
}


\details{
\code{circosPlot} function depicts correlations of variables selected with \code{block.splsda} among different blocks, using a generalisation of the method presented in Gonzalez et al 2012.
If \code{ncomp} is specified, then only the variables selected on that component are displayed.
}




\references{
Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Le Cao K.A. (2016).
DIABLO - multi omics integration for biomarker discovery.  BioRxiv \url{http://dx.doi.org/10.1101/067611}

Gonzalez I., Le Cao K.A., Davis M.J., Dejean S. (2012). Visualising associations between
paired 'omics' data sets. \emph{BioData Mining}; \bold{5}(1).

}

\author{Michael Vacher, Amrit Singh, Florian Rohart, Kim-Anh Le Cao}

\seealso{\code{\link{block.splsda}} and http://www.mixOmics.org/mixDIABLO for more details.}

\examples{
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- wrapper.sgccda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,
scheme = "horst",
verbose = FALSE,
bias = FALSE)

circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2, size.legend = 1.1)

par(mfrow=c(2,2))
circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2,
    size.legend = 1.1)
circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2,
    size.legend = 1.1, showIntraLinks = TRUE)
circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 1,
    size.legend = 1.1, showIntraLinks = TRUE)
circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2,
    size.legend = 1.1, showIntraLinks = TRUE, line = FALSE, size.variables = 0.5)

}

\keyword{regression}
\keyword{multivariate}

\name{circosPlot}
\encoding{latin1}
\alias{circosPlot}

\title{circosPlot for DIABLO}

\description{Displays variable correlation among different blocks
}

\usage{
circosPlot(object,
corThreshold,
showIntraLinks = FALSE,
line = TRUE,
size.legend = 0.8,
ncol.legend = 1,
size.labels = 1)
}

\arguments{
\item{object}{An object of class inheriting from \code{"block.splsda"}.}
\item{corThreshold}{Only shows links with a correlation higher than \code{corThreshold}}
\item{showIntraLinks}{if TRUE, shows the correlation higher than the threshold inside each block.}
\item{line}{if TRUE, shows the overall expression of the selected variables. see examples.}
\item{size.legend}{size of the legend}
\item{ncol.legend}{number of columns for the legend}
\item{size.labels}{size of the block labels}
}


\details{
\code{circosPlot} function shows the correlation of variables selected with \code{block.splsda} among different blocks.}




\references{
Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Le Cao K.A. (2016).
DIABLO - multi omics integration for biomarker discovery.
}

\author{Michale Vacher, Amrit Singh, Florian Rohart, Kim-Anh Le Cao}

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
scheme = "centroid",
verbose = FALSE,
bias = FALSE)

circosPlot(nutrimouse.sgccda, corThreshold = 0.7, ncol.legend = 2, size.legend = 1.1)

par(mfrow=c(2,2))
circosPlot(nutrimouse.sgccda, corThreshold = 0.7, ncol.legend = 2,
    size.legend = 1.1)
circosPlot(nutrimouse.sgccda, corThreshold = 0.7, ncol.legend = 2,
    size.legend = 1.1, showIntraLinks = TRUE)
circosPlot(nutrimouse.sgccda, corThreshold = 0.7, ncol.legend = 1,
    size.legend = 1.1, showIntraLinks = TRUE)
circosPlot(nutrimouse.sgccda, corThreshold = 0.7, ncol.legend = 2,
    size.legend = 1.1, showIntraLinks = TRUE, line = FALSE)

}

\keyword{regression}
\keyword{multivariate}

\name{s.match}
\encoding{latin1}
\alias{s.match}
\title{Plot of Paired Coordinates}
\description{
performs the scatter diagram for a paired coordinates. 
}

\usage{
s.match(df1xy, df2xy, xax = 1, yax = 2, pch = 20, cpoint = 1, 
        label = row.names(df1xy), clabel = 1, edge = TRUE, 
        xlim = NULL, ylim = NULL, grid = FALSE, addaxes = TRUE, 
        cgrid = 1, include.origin = TRUE, origin = c(0,0), 
        sub = "", csub = 1.25, possub = "bottomleft", 
        pixmap = NULL, contour = NULL, area = NULL, 
        add.plot = FALSE, col, lty = 1)
}

\arguments{
  \item{df1xy}{a data frame containing two columns from the first system.}
  \item{df2xy}{a data frame containing two columns from the second system.}
  \item{xax}{the column number for the \eqn{x}-axis of both the two systems.}
  \item{yax}{the column number for the \eqn{y}-axis of both the two systems.}
  \item{pch}{if \code{cpoint} \eqn{> 0}, an integer specifying the symbol or the 
	single character to be used in plotting points.}
  \item{cpoint}{a character size for plotting the points, used with \code{par("cex")*cpoint}. 
	If zero, no points are drawn.}
  \item{label}{a vector of strings of characters for the couple labels.}
  \item{clabel}{if not \code{NULL}, a character size for the labels, used with \code{par("cex")*clabel}.}
  \item{edge}{if \code{TRUE} the arrows are plotted, otherwise only the segments are drawn.}
  \item{xlim}{the ranges to be encompassed by the \eqn{x} axis, if \code{NULL} they are computed.}
  \item{ylim}{the ranges to be encompassed by the \eqn{y} axis, if \code{NULL} they are computed.}
  \item{grid}{a logical value indicating whether a grid in the background of the plot should be drawn.}
  \item{addaxes}{a logical value indicating whether the axes should be plotted.}
  \item{cgrid}{a character size, parameter used with \code{par("cex")*cgrid} to indicate the mesh of the grid.}
  \item{include.origin}{a logical value indicating whether the point "origin" should be belonged to the graph space.}
  \item{origin}{the fixed point in the graph space, for example \code{c(0,0)} the origin axes.}
  \item{sub}{a string of characters to be inserted as legend.}
  \item{csub}{a character size for the legend, used with \code{par("cex")*csub}.}
  \item{possub}{character string indicating the sub-title position 
	\code{("topleft", "topright", "bottomleft", "bottomright")}.}
  \item{pixmap}{an object \code{pixmap} displayed in the map background.}
  \item{contour}{a data frame with 4 columns to plot the contour of the map : each row gives a segment 
	\code{(x1,y1,x2,y2)}.}
  \item{area}{a data frame of class 'area' to plot a set of surface units in contour.}
  \item{add.plot}{if \code{TRUE} uses the current graphics window.}
  \item{col}{vector of size \code{n} (the number of samples) to indicate the biological conditions of each sample.}
  \item{lty}{set by default to 1.}
}

\details{
Graphical of the samples (individuals) is displayed in a superimposed manner where each sample will be 
indicated using an arrow. The start of the arrow indicates the location of the sample in $X$ in one plot, 
and the tip the location of the sample in $Y$ in the other plot. Short arrows indicate if both data sets 
strongly agree and long arrows a disagreement between the two data sets.
}

\value{
The matched call.
}

\author{Daniel Chessel from \pkg{ade4 R} package, modified by Kim-Anh Lê Cao.}

\references{
Lê Cao, K.-A., Martin, P.G.P., Robert-Granié, C. and Besse, P. (2009). Sparse canonical methods 
for biological data integration: application to a cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.
}

\examples{

## relevant only for canonical mode
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, mode = "canonical", ncomp = 3, 
                      keepX = c(50, 50, 50), keepY = c(10, 10, 10))

color.toxicity <- as.numeric(liver.toxicity$treatment[, 2])
label.toxicity <- liver.toxicity$treatment[, 1]
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)], 
        clabel = 0.5, label = label.toxicity, col = color.toxicity)

}
\keyword{multivariate}
\keyword{hplot}

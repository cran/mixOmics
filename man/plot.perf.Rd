\encoding{latin1}
\name{plot.perf}
\alias{plot.perf}

\title{Plot for model performance}

\description{
Function to plot performance criteria, such as MSEP, RMSEP, \eqn{R^2} or \eqn{Q^2}, as a
function of the number of components.
}

\usage{
\method{plot}{perf}(x, criterion = c("MSEP", "RMSEP", "R2", "Q2"),
     pred.method = "all",
     xlab = "number of components", ylab = NULL,
     LimQ2 = 0.0975, LimQ2.col = "darkgrey",
     cTicks = NULL, layout = NULL, \dots)
}

\arguments{
  \item{x}{an \code{perf} object.}
  \item{criterion}{character string. What type of validation criterion to plot for 
    \code{pls} or \code{spls}. One of \code{"MSEP"}, \code{"RMSEP"}, \code{"R2"} or \code{"Q2"}. 
	See \code{\link{perf}}.}
  \item{pred.method}{prediction method applied in \code{perf} for \code{plsda} or \code{splsda}. 
	See \code{\link{perf}}.}
  \item{xlab, ylab}{titles for \eqn{x} and \eqn{y} axes.  Typically character strings, 
    but can be expressions (e.g., \code{expression(R^2)}).}
  \item{LimQ2}{numeric value. Signification limit for the components in the model. 
    Default is \code{LimQ2 = 0.0975}.}
  \item{LimQ2.col}{character string specifying the color for the \code{LimQ2} line to be plotted.
  If \code{"none"} the line will not be plotted.
  }	
  \item{cTicks}{integer vector. Axis tickmark locations for the used number of components.
    Default is \code{1:ncomp} (see \code{\link{perf}}).}
  \item{layout}{numeric vector of length two giving the number of rows and columns in a 
    multi panel display. If not specified, \code{plot.perf} tries to be intelligent.}
  \item{\dots}{Further arguments sent to \code{\link{xyplot}} function.}
}

\details{
\code{plot.perf} creates one plot for each response variable in the
model, laid out in a multi panel display.  It uses \code{\link{xyplot}} for
performing the actual plotting.
}

\author{Ignacio Gonzalez, Florian Rohart.}

\seealso{\code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}}, \code{\link{splsda}}, 
\code{\link{perf}}.}

\examples{
require(lattice)

## validation for objects of class 'pls' or 'spls'
\dontrun{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

liver.pls <- pls(X, Y, ncomp = 3)
liver.perf <- perf(liver.pls, validation = "Mfold")
				   
plot(liver.perf, criterion = "R2", type = "l", layout = c(2, 2))
}
}

\keyword{regression}
\keyword{multivariate}
\keyword{hplot}

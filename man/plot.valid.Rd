\encoding{latin1}
\name{plot.valid}
\alias{plot.valid}

\title{Validation Plot}

\description{
Function to plot validation statistics, such as MSEP, RMSEP or \eqn{R^2}, as a
function of the number of components.
}

\usage{
\method{plot}{valid}(x, criterion = c("MSEP", "RMSEP", "R2"),
     xlab = "number of components", ylab = NULL,
     cTicks = 1:ncol(x[[1]]), layout = NULL, \dots)
}

\arguments{
  \item{x}{an \code{valid} object.}
  \item{criterion}{character string. What type of validation criterion to plot. One of 
    \code{"MSEP"}, \code{"RMSEP"} or \code{"R2"}.}
  \item{xlab, ylab}{titles for \eqn{x} and \eqn{y} axes.  Typically character strings, 
    but can be expressions (e.g., \code{expression(R^2)}.}  
  \item{cTicks}{integer vector. Axis tickmark locations for the used number of components.
    Default is from one to \code{ncomp} (see \code{\link{valid}}).}
  \item{layout}{numeric vector of length two giving the number of rows and columns in a 
    multi panel display. If not specified, \code{plot.valid} tries to be intelligent.}
  \item{\dots}{Further arguments sent to \code{\link{xyplot}} function.}
}

\details{
\code{plot.valid} creates one plot for each response variable in the
model, laid out in a multi panel display.  It uses \code{\link{xyplot}} for
performing the actual plotting.
}

\author{Ignacio González.}

\seealso{\code{\link{pls}}, \code{\link{spls}}, \code{\link{predict}}.}

\examples{
require(lattice)

\dontrun{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

liver.val <- valid(X, Y, ncomp = 5, mode = "regression", 
                   method = "pls", validation = "loo")
				   
plot.valid(liver.val, criterion = "RMSEP", type = "l", layout = c(2, 2))
}
}

\keyword{regression}
\keyword{multivariate}
\keyword{hplot}

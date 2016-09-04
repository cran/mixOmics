\encoding{latin1}
\name{plot.tune.splsda}
\alias{plot.tune.splsda}
\alias{plot.tune}

\title{Plot for model performance}

\description{
Function to plot performance criteria, such as classification error rate or balanced error rate on a tune.splsda result.
}

\usage{
\method{plot}{tune.splsda}(x, optimal = TRUE, sd = TRUE, \dots)
}

\arguments{
  \item{x}{an \code{tune.splsda} object.}
\item{optimal}{If TRUE, highlights the optimal keepX per component}
\item{sd}{If 'nrepeat' was used in the call to 'tune.splsda', error bar shows the standard deviation if sd=TRUE}

  \item{\dots}{Further arguments sent to \code{\link{xyplot}} function.}
}

\details{
\code{plot.tune.splsda} plots the classification error rate or the balanced error rate from x$error.rate, for each component of the model.
A circle highlights the optimal number of variables oneach component.

}

\author{Kim-Anh Le Cao, Florian Rohart, Francois Bartolo.}

\seealso{ \code{\link{tune.mint.splsda}}, \code{\link{tune.splsda}} and http://www.mixOmics.org for more details.}

\examples{

## validation for objects of class 'splsda'
\dontrun{
data(breast.tumors)
X = breast.tumors$gene.exp
Y = as.factor(breast.tumors$sample$treatment)
out = tune.splsda(X, Y, ncomp = 3, nrepeat = 2, logratio = "none",
test.keepX = c(5, 10, 15), folds = 10, dist = "max.dist",
progressBar = TRUE)


plot(out)

}

## validation for objects of class 'mint.splsda'
\dontrun{
data(stemcells)
data = stemcells$gene
type.id = stemcells$celltype
exp = stemcells$study

out = tune(method="mint.splsda", X=data,Y=type.id, ncomp=2, study=exp, test.keepX=seq(1,10,1))
out$choice.keepX.constraint

plot(out)

}
}

\keyword{regression}
\keyword{multivariate}
\keyword{hplot}

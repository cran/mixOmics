\name{mint.pls}
\encoding{latin1}
\alias{mint.pls}

\title{Vertical integration}


\description{Function to integrate and combine multiple independent studies using variants of multi-group PLS.
}

\usage{
mint.pls(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
scale = TRUE,
tol = 1e-06,
max.iter = 500,
near.zero.var = FALSE)

}

\arguments{
\item{X}{numeric matrix of predictors combining multiple independent studies on the same set of predictors. \code{NA}s are allowed.}
\item{Y}{Matrix or vector response for a multivariate regression framework. Data should be continuous variables (see \code{mint.plsda} for supervised classification and factor reponse)}
\item{ncomp}{Number of components to include in the model (see Details). Default to 2}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"} or \code{"canonical"}. See Details.}
\item{study}{factor indicating the membership of each sample to each of the studies being combined}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances. Default = \code{TRUE}.}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default = \code{FALSE}.}

}


\details{
\code{mint.pls} fits a vertical PLS-DA models with \code{ncomp} components in which several independent studies measured on the same variables are integrated. The aim is to explain the continuous outcome \code{Y}. The \code{study} factor indicates the membership of each sample in each study. We advise to only combine studies with more than 3 samples as the function performs internal scaling per study.

Multi (continuous)response are supported. \code{X} and \code{Y} can contain missing values. Missing values are handled by being disregarded during the cross product computations in the algorithm \code{block.pls} without having to delete rows with missing data. Alternatively, missing data can be imputed prior using the  \code{nipals} function.

The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References).

Useful graphical outputs are available, e.g. \code{\link{plotIndiv}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}.
}

\value{
\code{mint.pls} returns an object of class \code{"mint.pls", "pls"}, a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{Y}{the centered and standardized original response vector or matrix.}
\item{ncomp}{the number of components included in the model.}
\item{study}{The study grouping factor}
\item{mode}{the algorithm used to fit the model.}
\item{variates}{list containing the variates of X - global variates.}
\item{loadings}{list containing the estimated loadings for the variates - global loadings.}
\item{variates.partial}{list containing the variates of X relative to each study - partial variates.}
\item{loadings.partial}{list containing the estimated loadings for the partial variates - partial loadings.}
\item{names}{list containing the names to be used for individuals and variables.}
\item{nzv}{list containing the zero- or near-zero predictors information.}
\item{iter}{Number of iterations of the algorthm for each component}
\item{explained_variance}{Percentage of explained variance for each component and each study}

}

\references{
F Rohart, A Eslami, S Bougeard and K-A. Lê Cao. MINT: A multivariate integrative approach to identify a reproducible biomarker signature across multiple experiments and platforms. Submitted.

Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2014). Algorithms for multi-group PLS. J. Chemometrics, 28(3), 192-201.
}

\author{Florian Rohart, Kim-Anh Le Cao}

\seealso{\code{\link{spls}}, \code{\link{summary}},
\code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}}, \code{\link{mint.spls}}, \code{\link{mint.plsda}}, \code{\link{mint.splsda}} and http://www.mixOmics.org/mixMINT for more details.}

\examples{
# we will soon provide more examples on our website (data too large to be included in the package)
}

\keyword{regression}
\keyword{multivariate}

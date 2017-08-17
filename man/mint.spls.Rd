\name{mint.spls}
\encoding{UTF-8}
\alias{mint.spls}

\title{P-integration with variable selection}


\description{Function to integrate and combine multiple independent studies measured on the same variables or predictors (P-integration) using variants of multi-group sparse PLS for variable selection (unsupervised analysis).
}


\usage{
mint.spls(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
keepX.constraint = NULL,
keepY.constraint = NULL,
keepX = rep(ncol(X), ncomp),
keepY = rep(ncol(Y), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE)

}	

\arguments{
\item{X}{numeric matrix of predictors combining multiple independent studies on the same set of predictors. \code{NA}s are allowed.}
\item{Y}{Matrix or vector response for a multivariate regression framework. Data should be continuous variables (see \code{mint.splsda} for supervised classification and factor reponse)}
\item{ncomp}{Number of components to include in the model. Default to 2}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details.}
\item{study}{grouping factor indicating which samples are from the same study}
\item{keepX.constraint}{A list of vectors. Each vector entry indicates specific variables to be selected in \code{X} for each previous PLS component. See details. }
\item{keepY.constraint}{A list of vectors. Each vector entry indicates specific variables to be selected  in \code{Y} for each previous PLS component. See details. }
\item{keepX}{numeric vector indicating the number of variables to select in \code{X} on each component, or remaining components should \code{keepX.constraint} be provided. By default all variables are kept in the model.}
\item{keepY}{numeric vector indicating the number of variables to select in \code{Y} on each component, or remaining components should \code{keepY.constraint} be provided. By default all variables are kept in the model.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances. Default = \code{TRUE}.}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default = \code{FALSE}.}

}


\details{
\code{mint.spls} fits a vertical sparse PLS-DA models with \code{ncomp} components in which several independent studies measured on the same variables are integrated. The aim is to explain the continuous outcome \code{Y} and selecting correlated features between both data sets \code{X} and \code{Y}. The \code{study} factor indicates the membership of each sample in each study. We advise to only combine studies with more than 3 samples as the function performs internal scaling per study.

Multi (continuous)response are supported. \code{X} and \code{Y} can contain missing values. Missing values are handled by being disregarded during the cross product computations in the algorithm \code{mint.spls} without having to delete rows with missing data. Alternatively, missing data can be imputed prior using the  \code{nipals} function.

The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References and more details in \code{?pls}).

Constraint model (\code{constraint=TRUE)}. The model is run on the optimal list of selected features \code{keepX.constraint} from the previous components, as opposed to considering only a specified number of features \code{keepX}.


Variable selection is performed on each component for each block of \code{X}, and for \code{Y} if specified, via input parameter \code{keepX} and \code{keepY}, and \code{keepX.constraint} and \code{keepY.constraint} for a constraint model. \code{keepX.constraint} and \code{keepY.constraint} specify variable names to be selected in the previous components of the model.

Note that both \code{keepX.constraint} and \code{keepY.constraint} must have same length.

Useful graphical outputs are available, e.g. \code{\link{plotIndiv}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}.
}


\value{
\code{mint.spls} returns an object of class \code{"mint.spls","spls"}, a list
that contains the following components:

\item{X}{numeric matrix of predictors combining multiple independent studies on the same set of predictors. \code{NA}s are allowed.}
\item{Y}{the centered and standardized original response vector or matrix.}
\item{ncomp}{the number of components included in the model.}
\item{study}{The study grouping factor}
\item{mode}{the algorithm used to fit the model.}
\item{keepX}{Number of variables used to build each component of X}
\item{keepY}{Number of variables used to build each component of Y}
\item{keepX.constraint}{list indicating which variables where constrained to be used to build the first component of X}
\item{keepY.constraint}{list indicating which variables where constrained to be used to build the first component of Y}
\item{variates}{list containing the variates of X - global variates.}
\item{loadings}{list containing the estimated loadings for the variates - global loadings.}
\item{variates.partial}{list containing the variates of X relative to each study - partial variates.}
\item{loadings.partial}{list containing the estimated loadings for the partial variates - partial loadings.}
\item{names}{list containing the names to be used for individuals and variables.}
\item{nzv}{list containing the zero- or near-zero predictors information.}
\item{iter}{Number of iterations of the algorthm for each component}
\item{explained_variance}{Percentage of explained variance for each component and each study (note that contrary to PCA, this amount may not decrease as the aim of the method is not to maximise the variance, but the covariance between data sets).}

}


\references{
Rohart F, Eslami A, Matigian, N, Bougeard S, Le Cao K-A (2017). MINT: A multivariate integrative approach to identify a reproducible biomarker signature across multiple experiments and platforms. BMC Bioinformatics 18:128.

Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2014). Algorithms for multi-group PLS. J. Chemometrics, 28(3), 192-201.
}

\author{Florian Rohart, Kim-Anh Le Cao}

\seealso{\code{\link{spls}}, \code{\link{summary}},
\code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}}, \code{\link{mint.pls}}, \code{\link{mint.plsda}}, \code{\link{mint.splsda}} and http://www.mixOmics.org/mixMINT for more details.}



\examples{
# we will soon provide more examples on our website (data too large to be included in the package)
}

\keyword{regression}
\keyword{multivariate}

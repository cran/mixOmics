\name{mint.splsda}
\encoding{latin1}
\alias{mint.splsda}

\title{Vertical Discriminant Analysis integration with variable selection}


\description{Function to combine multiple independent studies (vertical integration) using variants of multi-group sparse PLS-DA for supervised classification with variable selection.
}


\usage{
mint.splsda(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
keepX.constraint = NULL,
keepX = rep(ncol(X), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE)

}

\arguments{
\item{X}{numeric matrix of predictors combining multiple independent studies on the same set of predictors. \code{NA}s are allowed.}
\item{Y}{A factor or a class vector indicating the discrete outcome of each sample.}
\item{ncomp}{Number of components to include in the model (see Details). Default to 2}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"} or \code{"canonical"}. See Details.}
\item{study}{factor indicating the membership of each sample to each of the studies being combined}
\item{keepX.constraint}{A list of vectors. Each vector entry indicates specific variables to be selected in \code{X} for each previous PLS component. }
\item{keepX}{numeric vector indicating the number of variables to select in \code{X} on each component, or remaining components should \code{keepX.constraint} be provided. By default all variables are kept in the model.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances. Default = \code{TRUE}.}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default = \code{FALSE}.}

}



\details{
\code{mint.splsda} function fits a vertical sparse PLS-DA models with \code{ncomp} components in which several independent studies measured on the same variables are integrated. The aim is to classify the discrete outcome \code{Y} and select variables that explain the outcome. The \code{study} factor indicates the membership of each sample in each study. We advise to only combine studies with more than 3 samples as the function performs internal scaling per study, and where all outcome categories are represented.

\code{X} can contain missing values. Missing values are handled by being disregarded during the cross product computations in the algorithm \code{mint.splsda} without having to delete rows with missing data. Alternatively, missing data can be imputed prior using the  \code{nipals} function.

The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References).


Variable selection is performed on each component for \code{X} via input parameter \code{keepX} and \code{keepX.constraint} afor a constraint model. \code{keepX.constraint} specifies variable names to be selected in the previous components of the model, see example.

Useful graphical outputs are available, e.g. \code{\link{plotIndiv}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}.
}


\value{
\code{mint.splsda} returns an object of class \code{"mint.splsda", "splsda"}, a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{Y}{the centered and standardized original response vector or matrix.}
\item{ind.mat}{the centered and standardized original response vector or matrix.}
\item{ncomp}{the number of components included in the model.}
\item{study}{The study grouping factor}
\item{mode}{the algorithm used to fit the model.}
\item{keepX}{Number of variables used to build each component of X}
\item{keepX.constraint}{list indicating which variables where constrained to be used to build the first component of X}
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
\code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}}, \code{\link{mint.pls}}, \code{\link{mint.plsda}}, \code{\link{mint.plsda}} and http://www.mixOmics.org/mixMINT for more details.}



\examples{
data(stemcells)

# -- feature selection with keepX
res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3, keepX = c(10, 5, 15),
study = stemcells$study)

plotIndiv(res)
#plot study-specific outputs for all studies
plotIndiv(res, study = "all.partial")

#plot study-specific outputs for study "2"
plotIndiv(res, study = "2")

# -- feature selection with keepX.constraint and keepX
# we force mint.spls to build the PLS-component 1 and 2 with only specific genes
keepX.constraint = list(comp1 = c("ENSG00000138756", "ENSG00000101470"),
comp2 = c("ENSG00000120049", "ENSG00000138685", "ENSG00000204248", "ENSG00000103260"))

# we select 10 and 15 variables on component 3 and 4
res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 4, study = stemcells$study,
keepX.constraint = keepX.constraint, keepX = c(10, 15))

#plot study-specific outputs for study "2", "3" and "4"
plotIndiv(res, study = c(2, 3, 4))
# poor results as the genes constraint on comp 1 and 2 were not relevant

}

\keyword{regression}
\keyword{multivariate}

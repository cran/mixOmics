\name{block.pls}
\encoding{latin1}
\alias{block.pls}

\title{Horizontal Partial Least Squares (PLS) integration}

\description{Integration of multiple data sets measured on the same samples, ie. Horizontal integration.}

\usage{
block.pls(X,
Y,
indY,
ncomp = 2,
design,
scheme,
mode,
scale = TRUE,
bias,
init ,
tol = 1e-06,
verbose,
max.iter = 500,
near.zero.var = FALSE)
}

\arguments{
\item{X}{A list of data sets (called 'blocks') measured on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets.}
\item{Y}{Matrix or vector response for a multivariate regression framework. Data should be continuous variables (see \code{block.splsda} for supervised classification and factor reponse)}
\item{indY}{To be supplied if Y is missing, indicates the position of the matrix / vector response in the list \code{X}}
\item{ncomp}{the number of components to include in the model. Default to 2. Applies to all blocks.}
\item{design}{numeric matrix of size (number of blocks in X) x (number of blocks in X) with 0 or 1 values. A value of 1 (0) indicates a relationship (no relationship) between the blocks to be modelled. If \code{Y} is provided instead of \code{indY}, the \code{design} matrix is changed to include relationships to \code{Y}. }
\item{scheme}{Either "horst", "factorial" or "centroid". Default = \code{centroid}, see reference.}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details. Default = \code{regression}.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances. Default = \code{TRUE}.}
\item{bias}{boleean. A logical value for biaised or unbiaised estimator of the var/cov. Default = \code{FALSE}.}
\item{init}{Mode of initialization use in the algorithm, either by Singular Value Decompostion of the product of each block of X with Y ("svd") or each block independently ("svd.single"). Default = \code{svd}.}
\item{tol}{Convergence stopping value.}
\item{verbose}{if set to \code{TRUE}, reports progress on computing.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default = \code{FALSE}.}
}


\details{
\code{block.spls} function fits a horizontal integration PLS model with a specified number of components per block).
An outcome needs to be provided, either by \code{Y} or by its position \code{indY} in the list of blocks \code{X}.
Multi (continuous)response are supported. \code{X} and \code{Y} can contain missing values. Missing values are handled by being disregarded during the cross product computations in the algorithm \code{block.pls} without having to delete rows with missing data. Alternatively, missing data can be imputed prior using the  \code{nipals} function.


The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References).
}

\value{
\code{block.pls} returns an object of class \code{"block.pls"}, a list
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{indY}{the position of the outcome Y in the output list X.}
  \item{ncomp}{the number of components included in the model for each block.}
  \item{mode}{the algorithm used to fit the model.}
  \item{variates}{list containing the variates of each block of X.}
  \item{loadings}{list containing the estimated loadings for the variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}
    \item{iter}{Number of iterations of the algorthm for each component}
\item{explained_variance}{Percentage of explained variance for each component and each block}

}

\references{
Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris: Editions Technic.

Wold H. (1966). Estimation of principal components and related models by iterative least squares.
In: Krishnaiah, P. R. (editors), \emph{Multivariate Analysis}. Academic Press, N.Y., 391-420.

Tenenhaus A. and Tenenhaus M., (2011), Regularized
Generalized Canonical Correlation Analysis,
Psychometrika, Vol. 76, Nr 2, pp 257-284.
}

\author{Florian Rohart, Benoit Gautier, Kim-Anh Le Cao}

\seealso{\code{\link{plotIndiv}}, \code{\link{plotArrow}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}}, \code{\link{selectVar}}, \code{\link{block.spls}}, \code{\link{block.plsda}} and http://www.mixOmics.org for more details.}

\examples{

# Example with TCGA multi omics study
# -----------------------------------
data("breast.TCGA")
# this is the X data as a list of mRNA and miRNA; the Y data set is a single data set of proteins
data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna)
# set up a full design where every block is connected
design = matrix(1, ncol = length(data), nrow = length(data),
dimnames = list(names(data), names(data)))
diag(design) =  0
design
# set number of component per data set
ncomp = c(2)

TCGA.block.pls = block.pls(X = data, Y = breast.TCGA$data.train$protein, ncomp = ncomp,
design = design)
TCGA.block.pls
# in plotindiv we color the samples per breast subtype group but the method is unsupervised!
# here Y is the protein data set
plotIndiv(TCGA.block.pls, group =  breast.TCGA$data.train$subtype, ind.names = FALSE)


}

\keyword{regression}
\keyword{multivariate}

\name{mixOmics}
\encoding{latin1}
\alias{mixOmics}

\title{PLS-derived methods: one function to rule them all}

\description{This function performs one of the PLS derived methods included in the mixOmics package that is the most appropriate for your input data, one of (mint).(block).(s)pls(da) depending on your input data (single data, list of data, discrete outcome, \dots)
}

\usage{
mixOmics(X,
Y,
indY,
study,
ncomp,
keepX.constraint,
keepY.constraint,
keepX,
keepY,
design,
tau = NULL,# rgcca, number between 0,1 or "optimal"
scheme, 
mode,
scale,
bias,
init,
tol =  1e-06,
verbose = FALSE,
max.iter = 500,
near.zero.var = FALSE)

}

\arguments{
\item{X}{Input data. Either a matrix or a list of data sets (called 'blocks') matching on the same samples. Data should be arranged in samples x variables, with samples order matching in all data sets.}
\item{Y}{Outcome. Either a numeric matrix of responses or a factor or a class vector for the discrete outcome.}
\item{indY}{To supply if Y is missing, indicates the position of the outcome in the list X}
\item{study}{grouping factor indicating which samples are from the same study}
\item{ncomp}{If \code{X} is a data matrix, \code{ncomp} is a single value. If \code{X} is a list of data sets, \code{ncomp} is a numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block).}
\item{keepX.constraint}{A list containing which variables of X are to be kept on each of the first PLS-components}
\item{keepY.constraint}{A list containing which variables of Y are to be kept on each of the first PLS-components}
\item{keepX}{Number of variables to keep in the \eqn{X}-loadings}
\item{keepY}{Number of variables to keep in the \eqn{Y}-loadings}
\item{design}{numeric matrix of size (number of blocks) x (number of blocks) with only 0 or 1 values. A value of 1 (0) indicates a relationship (no relationship) between the blocks to be modelled. If \code{Y} is provided instead of \code{indY}, the \code{design} matrix is changed to include relationships to \code{Y}. }
\item{tau}{numeric vector of length the number of blocks in \code{X}. Each regularization parameter will be applied on each block and takes the value between 0 (no regularisation) and 1. If tau = "optimal" the shrinkage paramaters are estimated for each block and
each dimension using the Schafer and Strimmer (2005)
analytical formula.}
\item{scheme}{Either "horst", "factorial" or "centroid" (Default: "centroid"), see reference paper.}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
\item{bias}{boleean. A logical value for biaised or unbiaised estimator of the var/cov (defaults to FALSE).}
\item{init}{Mode of initialization use in the algorithm, either by Singular Value Decompostion of the product of each block of X with Y ("svd") or each block independently ("svd.single") . Default to "svd".}
\item{tol}{Convergence stopping value.}
\item{verbose}{if set to \code{TRUE}, reports progress on computing.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE}
}

\details{
This function performs one of the PLS derived methods included in the mixOmics package that is the most appropriate for your input data, one of (mint).(block).(s)pls(da).

If your input data \code{X} is a matrix, then the algorithm is directed towards one of (mint).(s)pls(da) depending on your input data \code{Y} (factor for the discrete outcome directs the algorithm to DA analysis) and whether you input a \code{study} parameter (MINT analysis) or a \code{keepX}/\code{keepX.constraint} parameter (sparse analysis).

If your input data \code{X} is a list of matrices, then the algorithm is directed towards one of (mint).block.(s)pls(da) depending on your input data \code{Y} (factor for the discrete outcome directs the algorithm to DA analysis) and whether you input a \code{study} parameter (MINT analysis) or a \code{keepX}/\code{keepX.constraint} parameter (sparse analysis).


}

\author{Florian Rohart}

\seealso{\code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}}, \code{\link{splsda}},
\code{\link{mint.pls}}, \code{\link{mint.spls}}, \code{\link{mint.plsda}}, \code{\link{mint.splsda}},
\code{\link{block.pls}}, \code{\link{block.spls}}, \code{\link{block.plsda}}, \code{\link{block.splsda}},
\code{\link{mint.block.pls}}, \code{\link{mint.block.spls}}, \code{\link{mint.block.plsda}}, \code{\link{mint.block.splsda}}}

\references{
Rohart F. et al (2016, submitted). MINT: A multivariate integrative approach to identify a reproducible biomarker signature across multiple experiments and platforms.

Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2013). Multi-group PLS Regression: Application to Epidemiology.
In New Perspectives in Partial Least Squares and Related Methods, pages 243-255. Springer.

Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Le Cao K.A. (2016).
DIABLO - multi omics integration for biomarker discovery.

Le Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009). Sparse canonical methods
for biological data integration: application to a cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.

Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris: Editions Technic.

Wold H. (1966). Estimation of principal components and related models by iterative least squares.
In: Krishnaiah, P. R. (editors), \emph{Multivariate Analysis}. Academic Press, N.Y., 391-420.

}

\examples{


## -- directed towards PLS framework because X is a matrix and the study argument is missing
# ----------------------------------------------------
data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic
Y.factor = as.factor(liver.toxicity$treatment[, 4])

# directed towards PLS
out = mixOmics(X, Y, ncomp = 2)

# directed towards sPLS because of keepX and/or keepY
out = mixOmics(X, Y, ncomp = 2, keepX = c(50, 50), keepY = c(10, 10))

# directed towards PLS-DA because Y is a factor
out = mixOmics(X, Y.factor, ncomp = 2)

# directed towards sPLS-DA because Y is a factor and there is a keepX
out = mixOmics(X, Y.factor, ncomp = 2, keepX = c(20, 20))


## -- directed towards block.pls framework because X is a list
# ----------------------------------------------------
data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)

# directed towards block PLS
out = mixOmics(X = data, Y = Y,ncomp = 3)

# directed towards block sPLS because of keepX and/or keepY
out = mixOmics(X = data, Y = Y,ncomp = 3,
keepX = list(gene = c(10,10), lipid = c(15,15)))

# directed towards block PLS-DA because Y is a factor
out = mixOmics(X = data, Y = nutrimouse$diet, ncomp = 3)

# directed towards block sPLS-DA because Y is a factor and there is a keepX
out = mixOmics(X = data, Y = nutrimouse$diet, ncomp = 3,
keepX = list(gene = c(10,10), lipid = c(15,15)))


## -- directed towards mint.pls framework because of the study factor
# ----------------------------------------------------
data(stemcells)
# directed towards PLS
out = mixOmics(X = stemcells$gene, Y = unmap(stemcells$celltype), ncomp = 2)

# directed towards mint.PLS
out = mixOmics(X = stemcells$gene, Y = unmap(stemcells$celltype),
ncomp = 2, study = stemcells$study)

# directed towards mint.sPLS because of keepX and/or keepY
out = mixOmics(X = stemcells$gene, Y = unmap(stemcells$celltype),
ncomp = 2, study = stemcells$study, keepX = c(10, 5, 15))

# directed towards mint.PLS-DA because Y is a factor
out = mixOmics(X = stemcells$gene, Y = stemcells$celltype, ncomp = 2,
study = stemcells$study)

# directed towards mint.sPLS-DA because Y is a factor and there is a keepX
out = mixOmics(X = stemcells$gene, Y = stemcells$celltype, ncomp = 2,
study = stemcells$study, keepX = c(10, 5, 15))
}

\keyword{multivariate}
\keyword{hplot}
\keyword{dplot}
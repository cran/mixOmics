\name{wrapper.rgcca}
\encoding{latin1}
\alias{wrapper.rgcca}

\title{mixOmics wrapper for Regularised Generalised Canonical Correlation Analysis (rgcca)}

\description{Wrapper function to perform Regularized Generalised Canonical Correlation Analysis (rGCCA), a generalised approach for the integration of multiple datasets. For more details, see the \code{help(rgcca)} from the \pkg{RGCCA} package.
}

\usage{	
wrapper.rgcca(X,
    design = 1 - diag(length(X)),
    tau = rep(1, length(X)),
    ncomp = 1,
    keepX.constraint,
    keepX,
    scheme = "horst",
    scale = TRUE,
    bias = FALSE,
    init = "svd.single",
    tol = .Machine$double.eps,
    verbose = FALSE,
    max.iter=1000,
    near.zero.var = FALSE)
}

\arguments{
  \item{X}{a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables. \code{NA}s are not allowed.}
  \item{design}{numeric matrix of size (number of blocks) x (number of blocks) with only 0 or 1 values. A value of 1 (0) indicates a relationship (no relationship) between the blocks to be modelled using sGCCA.}
  \item{tau}{numeric vector of length the number of blocks in \code{X}. Each regularization parameter will be applied on each block and takes the value between 0 (no regularisation) and 1. If tau = "optimal" the shrinkage paramaters are estimated for each block and
  each dimension using the Schafer and Strimmer (2005)
  analytical formula.}
\item{ncomp}{the number of components to include in the model. Default to 1.}
\item{keepX.constraint}{A list of same length as X. Each entry keepX.constraint[[i]] is a list containing which variables of X[[i]] are to be kept on each of the first PLS-components}
\item{keepX}{ A vector of same length as X.  Each entry keepX[i] is the number of X[[i]]-variables kept in the model on the last components (once all keepX.constraint[[i]] are used).}
  \item{scheme}{Either "horst", "factorial" or "centroid"
  (Default: "horst").}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
\item{bias}{boleean. A logical value for biaised or unbiaised estimator of the var/cov (defaults to FALSE).}
\item{init}{Mode of initialization use in the algorithm, either by Singular Value Decompostion of the product of each block of X with Y ("svd") or each block independently ("svd.single") . Default to "svd.single".}
\item{tol}{Convergence stopping value.}
\item{verbose}{if set to \code{TRUE}, reports progress on computing.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE}
}


\details{
	This wrapper function performs rGCCA (see \pkg{RGCCA}) with \eqn{1, \ldots ,}\code{ncomp} components on each block data set.
	A supervised or unsupervised model can be run. For a supervised model, the \code{\link{unmap}} function should be used as an input data set.
	More details can be found on the package \pkg{RGCCA}.
}

\value{
\code{wrapper.rgcca} returns an object of class \code{"rgcca"}, a list 
that contains the following components:

  \item{data}{the input data set (as a list).}
  \item{design}{the input design.}
  \item{variates}{the sgcca components.}
  \item{loadings}{the loadings for each block data set (outer wieght vector).}
  \item{loadings.star}{the laodings, standardised.}
  \item{tau}{the input tau parameter.}
  \item{scheme}{the input schme.}
\item{ncomp}{the number of components included in the model for each block.}
  \item{crit}{the convergence criterion.}
  \item{AVE}{Indicators of model quality based on the
  Average Variance Explained (AVE): AVE(for one block),
  AVE(outer model), AVE(inner model)..}
  \item{names}{list containing the names to be used for individuals and variables.}
  More details can be found in the references.
}

\references{

  Tenenhaus A. and Tenenhaus M., (2011), Regularized
  Generalized Canonical Correlation Analysis,
  Psychometrika, Vol. 76, Nr 2, pp 257-284.
  
  Schafer J. and Strimmer K., (2005), A shrinkage approach
  to large-scale covariance matrix estimation and
  implications for functional genomics. Statist. Appl.
  Genet. Mol. Biol. 4:32.
}

\author{Arthur Tenenhaus, Vincent Guillemot and Kim-Anh Le Cao.}

\seealso{\code{\link{wrapper.rgcca}}, \code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{wrapper.sgcca}} and \url{http://www.mixOmics.org} for more details.}

\examples{
data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, gene expression and lipids are connected to the diet factor
# design = matrix(c(0,0,1,
#                   0,0,1,
#                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#note: the tau parameter is the regularization parameter
wrap.result.rgcca = wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0),
                                  ncomp = 2,
                                  scheme = "centroid", verbose = FALSE)
#wrap.result.rgcca
}


\keyword{multivariate}


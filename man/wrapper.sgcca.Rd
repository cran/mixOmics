\name{wrapper.sgcca}
\encoding{latin1}
\alias{wrapper.sgcca}

\title{mixOmics wrapper for Sparse Generalised Canonical Correlation Analysis (sgcca)}

\description{Wrapper function to perform Sparse Generalised Canonical Correlation Analysis (sGCCA), a generalised approach for the integration of multiple datasets. For more details, see the \code{help(sgcca)} from the \pkg{RGCCA} package.
}

\usage{	
wrapper.sgcca(X,
design = 1 - diag(length(X)),
penalty = NULL,
ncomp = 1,
keepX.constraint,
keepX,
scheme = "centroid",
mode="canonical",
scale = TRUE,
bias = TRUE,
init = "svd",
tol = .Machine$double.eps,
verbose = FALSE,
max.iter=1000,
near.zero.var = FALSE)
}

\arguments{
  \item{X}{a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables. \code{NA}s are not allowed.}
  \item{design}{numeric matrix of size (number of blocks) x (number of blocks) with only 0 or 1 values. A value of 1 (0) indicates a relationship (no relationship) between the blocks to be modelled using sGCCA.}
  \item{penalty}{numeric vector of length the number of blocks in \code{X}. Each penalty parameter will be applied on each block and takes the value between 0 (no variable selected) and 1 (all variables included).}
\item{ncomp}{the number of components to include in the model. Default to 1.}
\item{keepX.constraint}{A list of same length as X. Each entry keepX.constraint[[i]] is a list containing which variables of X[[i]] are to be kept on each of the first PLS-components}
\item{keepX}{ A vector of same length as X.  Each entry keepX[i] is the number of X[[i]]-variables kept in the model on the last components (once all keepX.constraint[[i]] are used).}
\item{scheme}{Either "horst", "factorial" or "centroid"
  (Default: "centroid").}
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
	This wrapper function performs sGCCA (see \pkg{RGCCA}) with \eqn{1, \ldots ,}\code{ncomp} components on each block data set.
	A supervised or unsupervised model can be run. For a supervised model, the \code{\link{unmap}} function should be used as an input data set.
	More details can be found on the package \pkg{RGCCA}.

Note that this function is the same as \code{\link{block.spls}} with different default arguments.
}

\value{
\code{wrapper.sgcca} returns an object of class \code{"sgcca"}, a list 
that contains the following components:

  \item{data}{the input data set (as a list).}
  \item{design}{the input design.}
  \item{variates}{the sgcca components.}
  \item{loadings}{the loadings for each block data set (outer wieght vector).}
  \item{loadings.star}{the laodings, standardised.}
  \item{penalty}{the input penalty parameter.}
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
  
  Tenenhaus A., Phillipe C., Guillemot, V., Le Cao K-A., Grill J., Frouin, V. Variable Selection For Generalized Canonical Correlation Analysis. 2013. (in revision)
}

\author{Arthur Tenenhaus, Vincent Guillemot and Kim-Anh Le Cao.}

\seealso{\code{\link{wrapper.sgcca}}, \code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{wrapper.rgcca}} and \url{http://www.mixOmics.org} for more details.}

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

#note: the penalty parameters will need to be tuned
wrap.result.sgcca = wrapper.sgcca(X = data, design = design, penalty = c(.3,.5, 1),
                                  ncomp = 2,
                                  scheme = "centroid", verbose = FALSE)
wrap.result.sgcca
#did the algo converge?
wrap.result.sgcca$crit  # yes
}


\keyword{multivariate}

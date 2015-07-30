\name{tau.estimate}
\alias{tau.estimate}
\title{Optimal shrinkage intensity parameters.}
\usage{
  tau.estimate(x)
}
\arguments{
  \item{x}{A list of matched data matrices on which the covariance matrix will be estimated.}
}
\value{
  \item{tau}{Optimal shrinkage intensity parameter}
}

\description{
  Estimation of the optimal shrinkage parameters as
  described in [1,2] and implemented in
  \pkg{RGCCA} [2].
  
  The shrinkage estimates can be used in the \code{\link{rcc}} function instead of the \code{\link{tune.rcc}} which can be long and costly to compute with very large data sets. Note that both functions \code{\link{tune.rcc}}  (which uses cross-validation) and \code{\link{tau.estimate}} (which uses the formula from [1]) may output different results. 
  
  Note: we advise to center and scale the input matrices in \code{\link{rcc}} since the intput data for \code{tau.estimate} are centered and scaled to estimate the shrinkage parameters (as in [1]).
  
  The shrinkage parameters are also estimated using that function when running \code{\link{wrapper.rgcca}} when \code{tau = 'optimal'}. This esxplains the name of this function (tau.estimate), which may not fit with our current nomenclature with \code{\link{rcc}}.
  
}

\references{
  [1] Schaefer J. and Strimmer K., 2005. A shrinkage
  approach to large-scale covariance matrix estimation and
  implications for functional genomics. Statist. Appl.
  Genet. Mol. Biol. 4:32.

  [2] Jelizarow M., Guillemot V., Tenenhaus A., Strimmer
  K., Boulesteix A.-L., 2010. Over-optimism in
  bioinformatics: an illustration. Bioinformatics
  26:1990-1998.
  
  [3] Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76(2), pp 257-284.
}

\seealso{\code{\link{tune.rcc}}, \code{\link{rcc}}, \code{\link{wrapper.rgcca}}}

\examples{
## tuning the regularization parameters for 'rcc'
# ----------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

# create a list as an input
A = list(X, Y)

param = sapply(A, tau.estimate)
# regularisation parameters for X (first) and then Y
param

# note: we advise scaling the data as the rcc input in that case
rcc.result = rcc(scale(X), scale(Y), ncomp = 2, lambda1 = param[1], lambda2 = param[2])

# we observe quite a good deparation of the diets
plotIndiv(rcc.result, ind.names = nutrimouse$diet, group = nutrimouse$genotype)
}

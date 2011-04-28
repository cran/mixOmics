\name{spca}
\encoding{latin1}
\alias{spca}

\title{Sparse Principal Components Analysis}

\description{ 
Performs a sparse principal components analysis to perform variable selection by using singular value decomposition.}

\usage{
spca(X, ncomp = 3, center = TRUE, scale. = TRUE, 
     keepX = rep(ncol(X),ncomp), iter.max = 500,
     tol = 1e-06)
}

\arguments{
  \item{X}{a numeric matrix (or data frame) which provides the 
    data for the sparse principal components analysis.}
  \item{ncomp}{integer, the number of components to keep.} 	
  \item{center}{a logical value indicating whether the variables should be shifted to be zero centered. 
    Alternatively, a vector of length equal the number of columns of \code{X} can be supplied. 
    The value is passed to \code{\link{scale}}.}  
  \item{scale.}{a logical value indicating whether the variables should be scaled to have 
    unit variance before the analysis takes place. The default is \code{TRUE}. See details.}   
  \item{iter.max}{integer, the maximum number of iterations to check convergence in each component.}
  \item{tol}{a positive real, the tolerance used in the iterative algorithm.} 
  \item{keepX}{numeric vector of length ncomp, the number of variables to keep in loading vectors. By default all 
   variables are kept in the model. See details.}
}

\details{
The calculation employs singular value decomposition of the (centered and scaled) data matrix and LASSO to generate
sparsity on the loading vectors.

\code{scale.= TRUE} is highly recommended as it will help obtaining orthogonal sparse loading vectors. 

\code{keepX} is the number of variables to keep in loading vectors. The difference between number of columns 
of \code{X} and \code{keepX} is the degree of sparsity, which refers to the number of zeros in each loading vector.

Note that \code{spca} does not apply to the data matrix with missing values. The biplot function for \code{spca} is not available.
}    

\value{
\code{spca} returns a list with class \code{"spca"} containing the following components: 
  \item{ncomp}{the number of components to keep in the calculation.} 
  \item{varX}{the adjusted cumulative percentage of variances explained.}
  \item{keepX}{the number of variables kept in each loading vector.}
  \item{iter}{the number of iterations needed to reach convergence for each component.}
  \item{rotation}{the matrix containing the sparse loading vectors.}
  \item{x}{the matrix containing the principal components.}  
}

\references{
Shen, H. and Huang, J. Z. (2008). Sparse principal component analysis via regularized 
low rank matrix approximation. \emph{Journal of Multivariate Analysis} \bold{99}, 1015-1034. 
}

\author{Kim-Anh Lê Cao, Fangzhou Yao, Leigh Coonan}

\seealso{\code{\link{pca}}} 

\examples{
data(liver.toxicity)
spca.rat <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))
spca.rat

## variable representation
plotVar(spca.rat, X.label = TRUE, cex = 0.5)
plot3dVar(spca.rat)

## samples representation
plotIndiv(spca.rat, ind.names = liver.toxicity$treatment[, 3], cex = 0.5, 
          col = as.numeric(liver.toxicity$treatment[, 3]))
plot3dIndiv(spca.rat, cex = 0.01, 
            col = as.numeric(liver.toxicity$treatment[, 3]))
}

\keyword{algebra}

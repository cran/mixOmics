\name{pca}
\encoding{latin1}
\alias{pca}

\title{Principal Components Analysis}

\description{ 
Performs a principal components analysis on the given data matrix that can contain missing values
by using the nipals algorithm. 
}

\usage{
pca(x, ncomp = NULL, retx = TRUE, 
    center = TRUE, scale. = FALSE, 
    comp.tol = NULL, max.iter = 500, 
    tol = 1e-09)
}

\arguments{
  \item{x}{a numeric matrix (or data frame) which provides the 
    data for the principal components analysis. It can contain missing values.}
  \item{ncomp}{integer, the number of components to keep to perform the 
	reconstitution of the data using the nipals algorithm. If \code{NULL}, 
    a criterion of stop is given.}	
  \item{retx}{a logical value indicating whether the rotated variables should be returned.} 
  \item{center}{a logical value indicating whether the variables should be shifted to be zero centered. 
    Alternately, a vector of length equal the number of columns of \code{x} can be supplied. 
    The value is passed to \code{\link{scale}}.} 
  \item{scale.}{a logical value indicating whether the variables should be scaled to have 
    unit variance before the analysis takes place. The default is FALSE for consistency with \code{prcomp} 
    function, but in general scaling is advisable. Alternatively, a vector of length equal the number of 
    columns of \code{x} can be supplied. The value is passed to \code{\link{scale}}.} 
  \item{comp.tol}{a value indicating the magnitude below which components should be omitted.} 
  \item{max.iter}{integer, the maximum number of iterations in the nipals algorithm.}
  \item{tol}{a positive real, the tolerance used in the nipals algorithm.}
}

\details{
The calculation is done by a singular value decomposition of the (centered and possibly scaled) data matrix, 
by using the nipals algorithm. Unlike \code{\link{princomp}}, the print method for these objects prints the 
results in a nice format and the \code{plot} method produces a scree plot.

Note that \code{scale. = TRUE} cannot be used if there are zero or constant (for \code{center = TRUE}) variables. 
 
Components are omitted if their standard deviations are less than or equal to \code{comp.tol} times 
the standard deviation of the first component. With the default null setting, no components are omitted. 
Other settings for \code{comp.tol} could be \code{comp.tol = 0} or \code{comp.tol = sqrt(.Machine$double.eps)}, 
which would omit essentially constant components.
}    

\value{
\code{pca} returns a list with class "pca" and "prcomp" containing the following components: 
  \item{sdev}{the standard deviations of the principal components (i.e., the square roots of the 
    eigenvalues of the covariance/correlation matrix, though the calculation is actually 
    done with the singular values of the data matrix).} 
  \item{rotation}{the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors).}
  \item{x}{if \code{retx} is true the value of the rotated data (the centred (and scaled if requested) data 
    multiplied by the rotation matrix) is returned.}  
  \item{center, scale}{the centering and scaling used, or \code{FALSE}.} 
}

\author{Ignacio González.}

\seealso{\code{\link{nipals}}, \code{\link{prcomp}}, \code{\link{biplot}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}},
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}}.}

\examples{
data(multidrug)

## this data set contains missing values, therefore 
## the 'prcomp' function cannot be applied
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
plot(pca.res)
summary(pca.res)
biplot(pca.res, xlabs = multidrug$cell.line$Class, cex = 0.7)
}

\keyword{algebra}
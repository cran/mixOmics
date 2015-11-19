\name{pca}
\encoding{latin1}
\alias{pca}

\title{Principal Components Analysis}

\description{ 
Performs a principal components analysis on the given data matrix that can 
contain missing values. If data are complete 'pca' uses Singular Value 
Decomposition, if there are some missing values, it uses the NIPALS algorithm. 
}

\usage{
pca(X,
ncomp = 2,
center = TRUE,
scale = FALSE,
max.iter = 500,
tol = 1e-09)
}

\arguments{
  \item{X}{a numeric matrix (or data frame) which provides the
    data for the principal components analysis. It can contain missing values.}
  \item{ncomp}{integer, if data is complete \code{ncomp} decides the number of components and associated 
    eigenvalues to display from the \code{pcasvd} algorithm and if the data has missing values, 
	\code{ncomp} gives the number of components to keep to perform the reconstitution of the data using 
	the NIPALS algorithm. If \code{NULL}, function sets \code{ncomp = min(nrow(X), ncol(X))}}	
  \item{center}{a logical value indicating whether the variables should be shifted to be zero centered. 
    Alternately, a vector of length equal the number of columns of \code{X} can be supplied. 
    The value is passed to \code{\link{scale}}.} 
  \item{scale}{a logical value indicating whether the variables should be scaled to have 
    unit variance before the analysis takes place. The default is \code{FALSE} for consistency with \code{prcomp} 
    function, but in general scaling is advisable. Alternatively, a vector of length equal the number of 
    columns of \code{X} can be supplied. The value is passed to \code{\link{scale}}.} 
  \item{max.iter}{integer, the maximum number of iterations in the NIPALS algorithm.}
  \item{tol}{a positive real, the tolerance used in the NIPALS algorithm.}
}

\details{
The calculation is done either by a singular value decomposition of the (possibly centered and scaled) 
data matrix, if the data is complete or by using the NIPALS algorithm if there is data missing. Unlike 
\code{\link{princomp}}, the print method for these objects prints the results in a nice format and the 
\code{plot} method produces a bar plot of the percentage of variance explaned by the principal 
components (PCs).

Note that \code{scale= TRUE} cannot be used if there are zero or constant (for \code{center = TRUE}) variables. 
 
Components are omitted if their standard deviations are less than or equal to \code{comp.tol} times 
the standard deviation of the first component. With the default null setting, no components are omitted. 
Other settings for \code{comp.tol} could be \code{comp.tol = sqrt(.Machine$double.eps)}, 
which would omit essentially constant components, or \code{comp.tol = 0}.
}    

\value{
\code{pca} returns a list with class \code{"pca"} and \code{"prcomp"} containing the following components: 
  \item{ncomp}{the number of principal components used.}
  \item{sdev}{the eigenvalues of the covariance/correlation matrix, though the calculation is actually 
    done with the singular values of the data matrix or by using NIPALS.} 
  \item{rotation}{the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors).}
  \item{loadings}{same as 'rotation' to keep the mixOmics spirit}
  \item{x}{the value of the rotated data (the centred (and scaled if requested) data 
    multiplied by the rotation/loadings matrix), also called the principal components.}  
      \item{variates}{same as 'x' to keep the mixOmics spirit}
  \item{center, scale}{the centering and scaling used, or \code{FALSE}.} 
}



\author{Ignacio Gonzalez.}

\seealso{\code{\link{nipals}}, \code{\link{prcomp}}, \code{\link{biplot}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}
and http://www.mixOmics.org for more details.}

\examples{
data(multidrug)

## this data set contains missing values, therefore 
## the 'prcomp' function cannot be applied
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
plot(pca.res)
print(pca.res)
biplot(pca.res, xlabs = multidrug$cell.line$Class, cex = 0.7)

# samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, 
          group = as.numeric(as.factor(multidrug$cell.line$Class)))
\dontrun{
plotIndiv(pca.res, cex = 0.2,
            col = as.numeric(as.factor(multidrug$cell.line$Class)),style="3d")
}
# variables representation
plotVar(pca.res)
\dontrun{
plotVar(pca.res, rad.in = 0.5, cex = 0.5,style="3d")
}
}

\keyword{algebra}

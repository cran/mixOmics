\name{pcatune}
\encoding{latin1}
\alias{pcatune}

\title{Tune the number of principal components in PCA}

\description{ 
\code{pcatune} can be used to quickly visualise the proportion of explained variance 
 for a large number of principal components in PCA. 
}

\usage{
pcatune(X, ncomp = NULL, center = TRUE, scale. = FALSE,
        max.iter = 500, tol = 1e-09)
}

\arguments{
  \item{X}{a numeric matrix (or data frame) which provides the 
    data for the principal components analysis. It can contain missing values.}
  \item{ncomp}{integer, the number of components to initially analyse in \code{pcatune} to choose a final 
    \code{ncomp} for \code{pca}. If \code{NULL}, 
     function sets \code{ncomp = min(nrow(X), ncol(X))}}	
  \item{center}{a logical value indicating whether the variables should be shifted to be zero centered. 
    Alternately, a vector of length equal the number of columns of \code{X} can be supplied. 
    The value is passed to \code{\link{scale}}.} 
  \item{scale.}{a logical value indicating whether the variables should be scaled to have 
    unit variance before the analysis takes place. The default is \code{FALSE} for consistency with \code{prcomp} 
    function, but in general scaling is advisable. Alternatively, a vector of length equal the number of 
    columns of \code{X} can be supplied. The value is passed to \code{\link{scale}}.} 
  \item{max.iter}{integer, the maximum number of iterations for the NIPALS algorithm.}
  \item{tol}{a positive real, the tolerance used for the NIPALS algorithm.}
}

\details{
The calculation is done either by a singular value decomposition of the (possibly centered and scaled) 
data matrix, if the data is complete or by using the NIPALS algorithm if there is data missing. Unlike 
\code{\link{princomp}}, the print method for these objects prints the results in a nice format and the 
\code{plot} method produces a bar plot of the percentage of variance explained by the principal 
components (PCs).

When using NIPALS (missing values), we make the assumption that the first (\code{min(ncol(X),} \code{nrow(X)}) 
principal components will account for 100 \% of the explained variance.

Note that \code{scale.= TRUE} cannot be used if there are zero or constant (for \code{center = TRUE}) variables.
}    

\value{
\code{pcatune} returns a list with class \code{"pcatune"} containing the following components: 
  \item{var}{the eigenvalues of the covariance/correlation matrix, though the calculation is 
actually  done with the singular values of the data matrix).} 
  \item{prop.var}{the proportion of explained variance accounted for by each 
	principal component is calculated using the eigenvalues}
  \item{cum.var}{the cumulative proportion of explained variance 
	accounted for by the sequential accumulation of principal components is calculated using the sum of the
	proportion of explained variance}
}

\author{Ignacio González and Leigh Coonan}

\seealso{\code{\link{nipals}}, \code{\link{biplot}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}},
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}}.}

\examples{
data(liver.toxicity)
tune <- pcatune(liver.toxicity$gene, center = TRUE, scale. = TRUE)
}

\keyword{algebra}

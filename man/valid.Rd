\name{valid}
\encoding{latin1}
\alias{valid}

\title{Compute validation criterion for PLS, sPLS, PLS-DA and sPLS-DA}

\description{
The valid function has been supersided by the perf function to evaluate the performance of the PLS, sPLS, PLS-DA and sPLS-DA methods. 
     See \code{help(perf)} for more details}

\usage{
valid(object, \ldots)	
}

\arguments{
  \item{object}{object of class inheriting from \code{"pls"}, \code{"plsda"}, 
    \code{"spls"} or \code{"splsda"}.}
  \item{...}{arguments to pass to \code{\link{nearZeroVar}}.}
}


\references{
Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris: Editions Technic.  

Le Cao, K. A., Rossouw D., Robert-Granie, C. and Besse, P. (2008). A sparse PLS for variable 
selection when integrating Omics data. \emph{Statistical Applications in Genetics and Molecular 
Biology} \bold{7}, article 35.

Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction (MSEP) Estimates for Principal Component 
Regression (PCR) and Partial Least Squares Regression (PLSR). \emph{Journal of Chemometrics} \bold{18}(9), 422-429.
}

\author{Sebastien Dejean, Ignacio Gonzalez and Kim-Anh Le Cao.}

\seealso{\code{\link{perf}}, \code{\link{predict}} and http://www.mixOmics.org for more details.}

\keyword{regression}
\keyword{multivariate}

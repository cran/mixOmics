\name{map}
\encoding{latin1}
\alias{map}
\title{
  Converts an indicator matrix to class vector
}
\description{
  Converts a matrix in which each row sums to \emph{1} 
  into the nearest matrix of \emph{(0,1)} indicator variables. Borrowed from the 
\pkg{mclust} package, used as internal function.
}
\usage{
map(Y, \dots)
}
\arguments{
  \item{Y}{
    An indicator (dummy) matrix.
  }  
  \item{\dots }{
    Provided to allow lists with elements other than the arguments can
    be passed in indirect or list calls with \code{do.call}.
  }
}
\value{
  A integer vector with one entry for each row of Y,
  in which the \emph{i}-th value is the column index at which the
  \emph{i}-th row of \code{Y} attains a maximum.
  
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{unmap}}
}
\keyword{cluster}
% docclass is function

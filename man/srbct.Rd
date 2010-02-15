\name{srbct}
\encoding{latin1}
\docType{data}
\alias{srbct}


\title{Small version of the small round blue cell tumors of childhood data}

\description{This data set from Khan et al. 2001 gives the expression measure of 2308 genes measured on 63 samples.} 


\usage{data(srbct)}


\format{
A list containing the following components:
  \describe{
    \item{\code{gene}}{data frame with 63 rows and 2308 columns. 
		The expression measure of 2308 genes for the 63 subjects.}
    \item{\code{class}}{A class vector containing the class tumour of each case (4 classes in total).}
    \item{\code{gene.name}}{data frame with 2308 rows and 2 columns 
		containing further information on the genes.}
}
}

\source{\url{http://research.nhgri.nih.gov/microarray/Supplement}}
\references{Khan et al. (2001) Classification and diagnostic prediction of cancers using gene expression profiling and artificial neural networks. Nature Medicine, 7, Number 6, June.}

\keyword{datasets} % at least one, from doc/KEYWORDS
%\keyword{microarray} a mettre mais j'ai pas la liste des kwords



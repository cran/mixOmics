\name{prostate}
\encoding{UTF-8}
\alias{prostate}

\docType{data}

\title{Human Prostate Tumors Data}

\description{
This data set contains the expression of 6,033 genes derived from 52 prostate tumors and from 50 non tumor prostate samples (referred to as normal) using oligonucleotide microarrays.} 

\usage{data(prostate)}

\format{
A list containing the following components:
  \describe{
    \item{\code{data}}{data matrix with 102 rows and 6033 columns. Each row represents 
      an experimental sample, and each column a single gene.}
	\item{\code{type}}{a vector containing the clinical status (normal or tumor).}
  }
}

\details{
This study investigated whether gene expression differences could distinguished
common clinical and pathological features of prostate cancer. Expression profiles were derived from 52
prostate tumors and from 50 non tumour prostate samples (referred to as normal) using oligonucleotide
microarrays containing probes for approximately 12,600 genes and ESTs. After preprocessing remains the
expression of 6033 genes.}


\references{
Singh D, Febbo P, Ross K, Jackson D, Manola J, Ladd C, Tamayo P, Renshaw A, D'Amico A, Richie J, et al. (2002). Gene expression correlates of clinical prostate cancer behavior. \emph{Cancer cell}, \bold{1}(2):203–209.
}

\keyword{datasets}

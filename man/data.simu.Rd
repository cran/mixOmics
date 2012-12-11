\name{data.simu}
\encoding{latin1}
\alias{data.simu}

\docType{data}

\title{Simulation study for multilevel analysis}

\description{
Simulation study to illustrate the use of the multilevel analysis for one and two-factor analysis with sPLS-DA.
This data set contains the expression measure of 1000 genes.
}

\usage{data(data.simu)}

\format{
A list containing the following components:
\describe{
    \item{\code{X}}{data frame with 48 samples and 1000 genes.}
    \item{\code{stimu}}{a factor indicating the conditions on each sample (stimulations)}
    \item{\code{sample}}{a vector indicating the repeated measurements on each unique subject. See details.}
}
}

\details{
In this cross-over design, repeated measurements are performed 12 experiments units (or unique subjects) for each of the 4 stimulations.

The simulation study was based on a mixed effects model (see reference for details). Ten clsuters of 100 genes were generated. Amongt those, 4 clusters of genes discriminate the 4 stimulations (denoted LIPO5, GAG+, GAG- and NS) as follows: \\
-2 gene clusters discriminate (LIPO5, GAG+) versus (GAG-, NS) \\
-2 gene clusters discriminate LIPO5 versus GAG+, while GAG+ and  NS have the same effect \\
- gene clusters discriminate GAG- versus NS, while LIPO5 and GAG+ have the same effect \\
-the 4 remaining clusters represent noisy signal (no stimulation effect)
}


\references{
Liquet, B., Lê Cao, K.-A., Hocini, H. and Thiebaut, R. A novel approach for biomarker selection and the
integration of repeated measures experiments from two platforms. \emph{Submitted}.
}

\keyword{datasets}

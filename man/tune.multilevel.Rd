\name{tune.multilevel}
\encoding{latin1}
\alias{tune.multilevel}
\alias{tune.splsdalevel1}
\alias{tune.splsdalevel2}
\alias{tune.splslevel}

\title{Tuning functions for multilevel analyses}

\description{
These functions were implemented to help tuning the variable selection parameters in the multilevel analyses.
}

\usage{
tune.multilevel(X, Y = NULL, 
                    cond = NULL,
                    sample, 
                    ncomp=1,
                    test.keepX=c(5, 10, 15), 
                    test.keepY=NULL, 
                    already.tested.X = NULL, 
                    already.tested.Y = NULL, 
                    method = NULL, 
                    dist)  
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{\code{if(method = 'spls')} numeric vector or matrix of continuous responses (for multi-response models) \code{NA}s are allowed.}
  \item{cond}{a factor or a class vector for one-factor discrete outcome, a numeric matrix of 2 columns for two-factor discrete outcome. See Details}
  \item{sample}{a vector indicating the repeated measured on each individual.}
  \item{ncomp}{the number of components to include in the model.}
  \item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
  \item{test.keepY}{If \code{method = 'spls'}, numeric vector for the different number of variables to test from the \eqn{Y} data set}
  \item{already.tested.X}{if \code{ncomp > 1} numeric vector indicating the number of variables to select rom the \eqn{X} data set on the previous ncomp-1 components}
   \item{already.tested.Y}{if \code{method = 'spls'} and \code{if(ncomp > 1)} numeric vector indicating the number of variables to select rom the \eqn{Y} data set on the previous ncomp-1 components}
  \item{method}{character string. Which multivariate method and type of analysis to choose, matching one of 'splsda' (Discriminant Analysis) or 'spls' (unsupervised integrative analysis). See Details.}
  \item{dist}{distance metric to use for \code{splsda} to estimate the classification error rate, 
    should be a subset of \code{"centroids.dist"}, \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details).}
}

\details{
This tuning function should be used to tune the parameters in the \code{multilevel} function.

If \code{method = 'splsda'}, a distance metric must be used, see \code{help(predict.splsda)} for details about the distances.

For a sPLS-DA one-factor analysis, leave-one-out cross-validation is performed, internally the training data is decomposed into within-subject variation.

For a sPLS-DA two-factor analysis, the correlation between components from the within-subject variation of X and the \code{cond} matrix is computed on the whole data set.


For a sPLS-DA two-factor analysis, the correlation between components from the within-subject variation of X and Y is computed on the whole data set.

}

\value{
Depending on the type of analysis performed, a list that contains:
  \item{error}{leave-one-out cross-validation is performed for one-factor sPLS-DA analysis.}
  \item{cor.value}{compute the correlation between latent variables for two-factor sPLS-DA analysis or sPLS.}

}

\references{
On multilevel analysis:
Liquet, B., Lê Cao, K.-A., Hocini, H. and Thiebaut, R. A novel approach for biomarker selection and the
integration of repeated measures experiments from two platforms. \emph{Submitted}.

Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K. (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA. \emph{Metabolomics}, \bold{6}(1), 119-128.

}

\author{Benoit Liquet, Kim-Anh Lê Cao.}

\seealso{\code{\link{multilevel}} and http://www.math.univ-toulouse.fr/~biostat/mixOmics/ for more details..}

\examples{
## First example: one-factor analysis with sPLS-DA
\dontrun{
data(data.simu) # simulated data
result.ex1 = tune.multilevel(data.simu$X,
                               cond = data.simu$stimu,
                               sample = data.simu$sample, 
                               ncomp=2,
                               test.keepX=c(5, 10, 15), 
                               already.tested.X = c(50),
                               method = 'splsda',
                               dist = 'mahalanobis.dist') 

result.ex1
}

## Second example: two-factor analysis with sPLS-DA
\dontrun{
data(liver.toxicity)
dose = liver.toxicity$treatment$Dose.Group
time = liver.toxicity$treatment$Time.Group
dose.time = cbind(dose, time)
repeat.indiv = c(1,2, 1,  2,  1,  2,  1,  2,  3,  3,  4,  3,  4,  3,  4,  4,  5,  6,  5,  5,  6,  5,  6,  7,  7,  8,  6,  7,  8,  7,  8,  8,  9, 10,  9, 10, 11, 9,  9, 10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14, 13, 14, 15, 16, 15, 16, 15, 16, 15, 16)


result.ex2 = tune.multilevel (liver.toxicity$gene,
                               cond = dose.time,
                               sample = repeat.indiv, 
                               ncomp=2,
                               test.keepX=c(5, 10, 15), 
                               already.tested.X = c(50),
                               method = 'splsda',
			       dist = 'mahalanobis.dist') 
result.ex2
}

## Third example: one-factor integrative analysis with sPLS
\dontrun{
result.ex3 = tune.multilevel (liver.toxicity$gene, liver.toxicity$clinic, 
                               cond = dose,
                               sample = repeat.indiv, 
                               ncomp=2,
                               test.keepX=c(5, 10, 15), 
                               test.keepY=c(2,3), 
                               already.tested.X = c(50), already.tested.Y = c(5),
                               method = 'spls') 

result.ex3
}
}

\keyword{regression}
\keyword{multivariate}

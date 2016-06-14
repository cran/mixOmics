\name{tune.multilevel}
\encoding{latin1}
\alias{tune.multilevel}
\alias{tune.splslevel}

\title{Tuning functions for multilevel analyses}

\description{
These functions were implemented to help tuning the variable selection parameters in the multilevel analyses.
}

\usage{
tune.multilevel(X, 
Y,
multilevel , 
ncomp = 1, 
test.keepX = c(5, 10, 15), 
test.keepY = NULL, 
already.tested.X = NULL, 
already.tested.Y = NULL, 
method,
mode = "regression",
validation = "Mfold",
folds = 10,
dist = "max.dist",
measure = "BER",
progressBar = TRUE,
near.zero.var = FALSE,
logratio = "none",
nrepeat=1,
light.output = TRUE)

}

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{\code{if(method = 'spls')} numeric vector or matrix of continuous responses (for multi-response models) \code{NA}s are allowed.}
\item{multilevel}{a numeric matrix that indicates the repeated measures on each individual, i.e. the individuals ID.}
\item{ncomp}{the number of components to include in the model.}
\item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
\item{test.keepY}{If \code{method = 'spls'}, numeric vector for the different number of variables to test from the \eqn{Y} data set}
\item{already.tested.X}{if \code{ncomp > 1} numeric vector indicating the number of variables to select rom the \eqn{X} data set on the pfirst components}
\item{already.tested.Y}{if \code{method = 'spls'} and \code{if(ncomp > 1)} numeric vector indicating the number of variables to select from the \eqn{Y} data set on the first components}
\item{method}{character string. Which multivariate method and type of analysis to choose, matching one of 'splsda' (Discriminant Analysis) or 'spls' (unsupervised integrative analysis). See Details.}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details.}

\item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or \code{"loo"} (see below). Default is \code{"Mfold"}.}
\item{folds}{the folds in the Mfold cross-validation. See Details.}
\item{dist}{distance metric to use for \code{splsda} to estimate the classification error rate,
should be a subset of \code{"centroids.dist"}, \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details).}

\item{measure}{Two misclassification measure are available: overall misclassification error \code{overall} or the Balanced Error Rate \code{BER}}
\item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}
\item{logratio}{one of ('none','CLR'). Default to 'none'}
\item{nrepeat}{Number of times the Cross-Validation process is repeated.}
\item{light.output}{if set to FALSE, the prediction/classification of each sample for each of \code{test.keepX} and each comp is returned.}
}

\details{
This tuning function should be used to tune the parameters in the \code{multilevel} function.

If \code{method = 'splsda'}, a distance metric must be used, see \code{help(predict.splsda)} for details about the distances.

For a sPLS-DA one-factor analysis, M-fold cross-validation is performed, internally the training data is decomposed into within-subject variation.

For a sPLS-DA two-factor analysis, the correlation between components from the within-subject variation of X and the \code{cond} matrix is computed on the whole data set. The reason why we cannot obtain a corss-validation error rate as for the spls-DA one-factor analysis is because of the dififculty to decompose and predict the within matrices within each fold.

For a sPLS two-factor analysis a sPLS canonical mode is run, and the correlation between components from the within-subject variation of X and Y is computed on the whole data set.

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{folds}.
The folds also can be supplied as a list of vectors containing the indexes defining each 
fold as produced by \code{split}.

If \code{validation = "loo"}, leave-one-out cross-validation is performed. By default \code{folds} is set to the number of unique individuals.


}

\value{
Depending on the type of analysis performed, a list that contains:
\item{mat.mean.error,}{returns the prediction error for each \code{test.keepX} on each component}
\item{choice.keepX.constraint}{returns the variables selected by the optimal keepX on each component.}
\item{error.per.class}{returns the error rate for each level of \code{Y} and for each component computed with the optimal keepX}

\item{predict}{Prediction values for each sample, each \code{test.keepX}, each comp and each repeat. Only if light.output=FALSE}
\item{class}{Predicted class for each sample, each \code{test.keepX}, each comp and each repeat. Only if light.output=FALSE}

\item{cor.value}{compute the correlation between latent variables for two-factor sPLS-DA analysis.}

}

\references{
On multilevel analysis:
Liquet, B., Le Cao, K.-A., Hocini, H. and Thiebaut, R. (2012) A novel approach for
biomarker selection and the integration of repeated measures experiments from two
platforms. \emph{BMC Bioinformatics} \bold{13}:325.

Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K. (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA. \emph{Metabolomics}, \bold{6}(1), 119-128.

}

\author{Benoit Liquet, Kim-Anh Le Cao, Benoit Gautier.}

\seealso{\code{\link{multilevel}} and http://www.mixOmics.org for more details.}

\examples{
## First example: one-factor analysis with sPLS-DA
\dontrun{
  data(vac18.simulated) # simulated data
  design <- data.frame(sample = vac18.simulated$sample,
                       stimu = vac18.simulated$stimulation)
  
    result.ex1 = tune.multilevel(vac18.simulated$genes,
                               multilevel = design,
                               ncomp=2,
                               test.keepX=c(5, 10, 15), 
                               already.tested.X = c(50),
                               method = 'splsda',
                               dist = 'mahalanobis.dist',
                               validation = 'loo') 
  
  # error rate for the tested parameters est.keepX=c(5, 10, 15)
  result.ex1$error
  # prediction for ncomp = 2 and keepX = c(50, 15) (15 is the last tested parameter)
  result.ex1$prediction.all
  table(vac18.simulated$stimulation, result.ex1$prediction.all)
}



## Second example: two-factor analysis with sPLS-DA
\dontrun{
  data(liver.toxicity)
  dose <- as.factor(liver.toxicity$treatment$Dose.Group)
  time <- as.factor(liver.toxicity$treatment$Time.Group)
  # note: we made up those data, pretending they are repeated measurements
  repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                    6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                    10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                    13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
  summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
  
  design <- data.frame(sample = repeat.indiv,
                       dose = dose,
                       time = time)
  
  result.ex2 = tune.multilevel(liver.toxicity$gene,
                                multilevel = design,
                                ncomp=2,
                                test.keepX=c(5, 10, 15), 
                                already.tested.X = c(50),
                                method = 'splsda',
                                dist = 'mahalanobis.dist') 
  result.ex2
}

## Third example: one-factor integrative analysis with sPLS
\dontrun{
  data(liver.toxicity)
  # note: we made up those data, pretending they are repeated measurements
  repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                    6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                    10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                    13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
  summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
  
  # here we are only interested in a one level variation split since spls is an unsupervised method
  design <- data.frame(sample = repeat.indiv)
  
  result.ex3 = tune.multilevel(X = liver.toxicity$gene, Y = liver.toxicity$clinic, 
                                multilevel = design,
                                mode = 'canonical',
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

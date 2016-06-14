\name{tune.splsda}
\encoding{latin1}
\alias{tune.splsda}


\title{Tuning functions for sPLS-DA method}

\description{
Computes M-fold or Leave-One-Out Cross-Validation scores on a user-input
grid to determine optimal values for the sparsity parameters in \code{splsda}.
}

\usage{
tune.splsda(X, Y, ncomp = 1,
test.keepX = c(5, 10, 15), already.tested.X = NULL, validation = "Mfold", folds = 10,
dist = "max.dist", measure = "BER", progressBar = TRUE, near.zero.var = FALSE,
nrepeat = 1, logratio = c('none','CLR'), multilevel = NULL, light.output = TRUE)

}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{\code{if(method = 'spls')} numeric vector or matrix of continuous responses (for multi-response models) \code{NA}s are allowed.}
  \item{ncomp}{the number of components to include in the model.}
  \item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
  \item{already.tested.X}{if \code{ncomp > 1} numeric vector indicating the number of variables to select from the \eqn{X} data set on the firsts components}
\item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or
\code{"loo"} (see below). Default is \code{"Mfold"}.}
\item{folds}{the folds in the Mfold cross-validation. See Details.}
\item{dist}{distance metric to use for \code{splsda} to estimate the classification error rate,
should be a subset of \code{"centroids.dist"}, \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details).}
\item{measure}{Two misclassification measure are available: overall misclassification error \code{overall} or the Balanced Error Rate \code{BER}}
\item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}
\item{nrepeat}{Number of times the Cross-Validation process is repeated.}
\item{logratio}{one of ('none','CLR'). Default to 'none'}
\item{multilevel}{Design matrix for multilevel anaylis (for repeated measurements) that indicates the repeated measures on each individual, i.e. the individuals ID. See Details.}

\item{light.output}{if set to FALSE, the prediction/classification of each sample for each of \code{test.keepX} and each comp is returned.}

}

\details{

mention stratified subsamplings

This tuning function should be used to tune the parameters in the \code{splsda} function.

For a sPLS-DA, M-fold or LOO cross-validation is performed with stratified subsampling where all classes are represented in each fold.

For sPLS-DA multilevel one-factor analysis, M-fold or LOO cross-validation is performed where all repeated measurements of one sample are in the same fold. Note that logratio transform and the multilevel analysis are performed internally and independently on the training and test set.

For a sPLS-DA multilevel two-factor analysis, the correlation between components from the within-subject variation of X and the \code{cond} matrix is computed on the whole data set. The reason why we cannot obtain a cross-validation error rate as for the spls-DA one-factor analysis is because of the dififculty to decompose and predict the within matrices within each fold.

For a sPLS two-factor analysis a sPLS canonical mode is run, and the correlation between components from the within-subject variation of X and Y is computed on the whole data set.

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{folds}.

If \code{validation = "loo"}, leave-one-out cross-validation is performed. By default \code{folds} is set to the number of unique individuals.


}

\value{
Depending on the type of analysis performed, a list that contains:
\item{mat.mean.error,}{returns the prediction error for each \code{test.keepX} on each component}
\item{choice.keepX}{returns the variables selected by the optimal keepX on each component.}
\item{error.per.class}{returns the error rate for each level of \code{Y} and for each component computed with the optimal keepX}

\item{predict}{Prediction values for each sample, each \code{test.keepX}, each comp and each repeat. Only if light.output=FALSE}
\item{class}{Predicted class for each sample, each \code{test.keepX}, each comp and each repeat. Only if light.output=FALSE}

\item{cor.value}{compute the correlation between latent variables for two-factor sPLS-DA analysis.}

}


\author{Kim-Anh Le Cao, Benoit Gautier, Francois Bartolo, Florian Rohart.}

\seealso{\code{\link{splsda}} and http://www.mixOmics.org for more details.}

\examples{
## First example: analysis with sPLS-DA
\dontrun{
data(breast.tumors)
X = breast.tumors$gene.exp
Y = as.factor(breast.tumors$sample$treatment)
tune = tune.splsda(X, Y, ncomp = 1, nrepeat = 10, logratio = "none",
    test.keepX = c(5, 10, 15), folds = 10, dist = "max.dist",
    already.tested.X = NULL, progressBar = TRUE)
}

## only tune component 3 and 4
# keeping 5 and 10 variables on the first two components respectively
\dontrun{
tune = tune.splsda(X = X,Y = Y, ncomp = 4,
already.tested.X = c(5,10),
test.keepX = seq(1,10,2), progressBar = TRUE)

}



## Second example: multilevel one-factor analysis with sPLS-DA
\dontrun{
data(vac18)
X = vac18$genes
Y = vac18$stimulation
# sample indicates the repeated measurements
design = data.frame(sample = vac18$sample)

tune = tune.splsda(X, Y = Y, ncomp = 3, nrepeat = 10, logratio = "none",
    test.keepX = c(5,50,100),folds = 10, dist = "max.dist", multilevel = design)

}





}

\keyword{regression}
\keyword{multivariate}

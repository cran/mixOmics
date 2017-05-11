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
test.keepX = c(5, 10, 15), already.tested.X, constraint = FALSE, validation = "Mfold",
folds = 10, dist = "max.dist", measure = "BER", auc = FALSE,
progressBar = TRUE, max.iter = 100, near.zero.var = FALSE, nrepeat = 1,
logratio = c('none','CLR'), multilevel = NULL, light.output = TRUE, cpus)

}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{\code{if(method = 'spls')} numeric vector or matrix of continuous responses (for multi-response models) \code{NA}s are allowed.}
  \item{ncomp}{the number of components to include in the model.}
  \item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
\item{already.tested.X}{Optional, if \code{ncomp > 1} A list containing which variables of X are to be kept on each of the first PLS-components (if \code{constraint=TRUE}) or a numeric vector indicating the number of variables to select from the \eqn{X} data set on the firsts components (if \code{constraint = FALSE})}
\item{constraint}{Indicate whether the performance of the model should be evaluated w.r.t the number of variables selected (keepX, default constraint = FALSE) or w.r.t the specific variables selected by the model (keepX.constraint, constraint = TRUE), see details and examples. }
\item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or
\code{"loo"} (see below). Default is \code{"Mfold"}.}
\item{folds}{the folds in the Mfold cross-validation. See Details.}
\item{dist}{distance metric to use for \code{splsda} to estimate the classification error rate,
should be a subset of \code{"centroids.dist"}, \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details).}
\item{measure}{Two misclassification measure are available: overall misclassification error \code{overall} or the Balanced Error Rate \code{BER}}
\item{auc}{if \code{TRUE} calculate the Area Under the Curve (AUC) performance of the model.}
\item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}
\item{nrepeat}{Number of times the Cross-Validation process is repeated.}
\item{logratio}{one of ('none','CLR'). Default to 'none'}
\item{multilevel}{Design matrix for multilevel analysis (for repeated measurements) that indicates the repeated measures on each individual, i.e. the individuals ID. See Details.}

\item{light.output}{if set to FALSE, the prediction/classification of each sample for each of \code{test.keepX} and each comp is returned.}
\item{cpus}{Number of cpus to use when running the code in parallel.}

}

\details{

This tuning function should be used to tune the parameters in the \code{splsda} function (number of components and number of variables in \code{keepX} to select).

For a sPLS-DA, M-fold or LOO cross-validation is performed with stratified subsampling where all classes are represented in each fold.

If \code{validation = "loo"}, leave-one-out cross-validation is performed. By default \code{folds} is set to the number of unique individuals.

The function outputs the optimal number of components that achieve the best performance based on the overall error rate or BER. The assessment is data-driven and similar to the process detailed in (Rohart et al., 2016), where one-sided t-tests assess whether there is a gain in performance when adding a component to the model. Our experience has shown that in most case, the optimal number of components is the number of categories in \code{Y} - 1.


For sPLS-DA multilevel one-factor analysis, M-fold or LOO cross-validation is performed where all repeated measurements of one sample are in the same fold. Note that logratio transform and the multilevel analysis are performed internally and independently on the training and test set.

For a sPLS-DA multilevel two-factor analysis, the correlation between components from the within-subject variation of X and the \code{cond} matrix is computed on the whole data set. The reason why we cannot obtain a cross-validation error rate as for the spls-DA one-factor analysis is because of the dififculty to decompose and predict the within matrices within each fold.

For a sPLS two-factor analysis a sPLS canonical mode is run, and the correlation between components from the within-subject variation of X and Y is computed on the whole data set.

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{folds}.


If \code{auc = TRUE} and there are more than 2 categories in \code{Y}, the Area Under the Curve is averaged using one-vs-all comparison. Note however that the AUC criteria may not be particularly insightful as the prediction threshold we use in sPLS-DA differs from an AUC threshold (sPLS-DA relies on prediction distances for predictions, see \code{?predict.splsda}).

Constraint model. To improve performance, a constraint model is now implemented (\code{constraint=TRUE)}. The performance is assessed on same list of selected features from the previous components in the different CV folds, as opposed to considering only the number of features. Our experience has shown that the constraint model can improve the performance of the methods.

}

\value{
Depending on the type of analysis performed, a list that contains:
\item{error.rate}{returns the prediction error for each \code{test.keepX} on each component, averaged across all repeats and subsampling folds. Standard deviation is also output. All error rates are also available as a list.}
\item{choice.keepX}{returns the number of variables selected (optimal keepX) on each component.}
\item{choice.keepX.constraint}{returns the variables selected by the optimal keepX on each component. only if constraint=TRUE.}
\item{choice.nconp}{returns the optimal number of components in the sPLS-DA model}
\item{error.rate.class}{returns the error rate for each level of \code{Y} and for each component computed with the optimal keepX}

\item{predict}{Prediction values for each sample, each \code{test.keepX}, each comp and each repeat. Only if light.output=FALSE}
\item{class}{Predicted class for each sample, each \code{test.keepX}, each comp and each repeat. Only if light.output=FALSE}

\item{auc}{AUC mean and standard deviation if the number of categories in \code{Y} is greater than 2, see details above. Only if auc = TRUE}

\item{cor.value}{only if multilevel analysis with 2 factors: correlation between latent variables.}

}

\references{
mixOmics manuscript: Rohart F, Gautier B, Singh A, Le Cao K-A. mixOmics: an R package for 'omics feature selection and multiple data integration. BioRxiv available here: \url{http://biorxiv.org/content/early/2017/02/14/108597}

}



\author{Kim-Anh Le Cao, Benoit Gautier, Francois Bartolo, Florian Rohart.}

\seealso{\code{\link{splsda}}, \code{\link{predict.splsda}} and http://www.mixOmics.org for more details.}

\examples{
## First example: analysis with sPLS-DA
\dontrun{
data(breast.tumors)
X = breast.tumors$gene.exp
Y = as.factor(breast.tumors$sample$treatment)
tune = tune.splsda(X, Y, ncomp = 1, nrepeat = 10, logratio = "none",
    test.keepX = c(5, 10, 15), folds = 10, dist = "max.dist",
    progressBar = TRUE)

# 2 components, optimising 'keepX' with 'constraint = FALSE'
tune = tune.splsda(X, Y, ncomp = 2, test.keepX = c(5, 10, 15),
folds = 10, dist = "max.dist", constraint = FALSE,
progressBar = TRUE)
tune$choice.keepX
tune$choice.keepX.constraint

# 2 components, optimising 'keepX.constraint' with 'constraint = TRUE'
tune = tune.splsda(X, Y, ncomp = 2, test.keepX = c(5, 10, 15),
folds = 10, dist = "max.dist", constraint = TRUE,
progressBar = TRUE)
tune$choice.keepX
tune$choice.keepX.constraint

}

## only tune component 3 and 4
# keeping 5 and 10 variables on the first two components respectively
\dontrun{
tune = tune.splsda(X = X,Y = Y, ncomp = 4,
already.tested.X = c(5,10),
test.keepX = seq(1,10,2), progressBar = TRUE)

}

## only tune component 3 and 4
# keeping specific variables on the first two components respectively
\dontrun{
keepX.constraint = list(comp1 = c("818", "791", "911", "592"),
comp2 = c("593", "743", "309")
)
# constraint is set to TRUE because we use a keepX.constraint on the first 2 PLS-components
tune = tune.splsda(X = X,Y = Y, ncomp = 4,
already.tested.X = keepX.constraint,
test.keepX = seq(1,10,2), constraint = TRUE, progressBar = TRUE)

tune$choice.keepX.constraint

plot(tune)

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

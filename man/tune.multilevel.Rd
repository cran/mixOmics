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
constraint = FALSE,
method,
mode = "regression",
validation = "Mfold",
folds = 10,
dist = "max.dist",
measure = "BER",
auc = FALSE,
progressBar = TRUE,
near.zero.var = FALSE,
logratio = "none",
nrepeat=1,
light.output = TRUE)

}

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{\code{if(method = 'spls')} numeric vector or matrix of continuous responses (for multi-response models) \code{NA}s are allowed.}
\item{multilevel}{Design matrix for multilevel analysis (for repeated measurements). A numeric matrix or data frame. For a one level factor decomposition, the input is a vector indicating the repeated measures on each individual, i.e. the individuals ID. For a two level decomposition with splsda models, the two factors are included in Y. Finally for a two level decomposition with spls models, 2nd AND 3rd columns in design indicate those factors (see example in \code{?splsda} and \code{?spls}).}
\item{ncomp}{the number of components to include in the model.}
\item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
\item{test.keepY}{If \code{method = 'spls'}, numeric vector for the different number of variables to test from the \eqn{Y} data set}
\item{already.tested.X}{if \code{ncomp > 1} numeric vector indicating the number of variables to select rom the \eqn{X} data set on the pfirst components}
\item{already.tested.Y}{if \code{method = 'spls'} and \code{if(ncomp > 1)} numeric vector indicating the number of variables to select from the \eqn{Y} data set on the first components}
\item{constraint}{Indicate whether the performance of the model should be evaluated w.r.t the number of variables selected (keepX, default constraint = FALSE) or w.r.t the specific variables selected by the model (keepX.constraint, constraint = TRUE), see details and examples. }
\item{method}{character string. Which multivariate method and type of analysis to choose, matching one of 'splsda' (Discriminant Analysis) or 'spls' (unsupervised integrative analysis). See Details.}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details.}

\item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or \code{"loo"} (see below). Default is \code{"Mfold"}.}
\item{folds}{the folds in the Mfold cross-validation. See Details.}
\item{dist}{distance metric to use for \code{splsda} to estimate the classification error rate,
should be a subset of \code{"centroids.dist"}, \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details).}

\item{measure}{Two misclassification measure are available: overall misclassification error \code{overall} or the Balanced Error Rate \code{BER}}
\item{auc}{if \code{TRUE} calculate the Area Under the Curve (AUC) performance of the model. Only used for \code{method='splsda'}}
\item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}
\item{logratio}{one of ('none','CLR'). Default to 'none'}
\item{nrepeat}{Number of times the Cross-Validation process is repeated.}
\item{light.output}{if set to FALSE, the prediction/classification of each sample for each of \code{test.keepX} and each comp is returned.}
}

\details{
This tuning function should be used to tune the parameters when using a variance decomposition ('multilevel') with a repeated measurement design, see also details in \code{?predict.splsda}.

If \code{method = 'splsda'}, a distance metric must be used, see \code{help(predict.splsda)} for details about the distances.

The function outputs the optimal number of components that achieve the best performance based on the overall error rate or BER. The assessment is data-driven and similar to the process detailed in (Rohart et al., 2016), where one-sided t-tests assess whether there is a gain in performance when adding a component to the model. Our experience has shown that in most case, the optimal number of components is the number of categories in \code{Y} - 1.

For a sPLS-DA one-factor analysis, M-fold cross-validation is performed, internally the training data is decomposed into within-subject variation.

For a sPLS-DA two-factor analysis, the correlation between components from the within-subject variation of X and a matrix including the two factors (design[,-1]) is computed on the whole data set. We cannot obtain a cross-validation error rate as for the spls-DA one-factor analysis because of the difficulty to decompose and predict the within matrices within each fold.

For a sPLS two-factor analysis a sPLS canonical mode is run, and the correlation between components from the within-subject variation of X and Y is computed on the whole data set.

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{folds}.
The folds also can be supplied as a list of vectors containing the indexes defining each 
fold as produced by \code{split}.

If \code{validation = "loo"}, leave-one-out cross-validation is performed. By default \code{folds} is set to the number of unique individuals.

If \code{auc = TRUE} and there are more than 2 categories in \code{Y}, the Area Under the Curve is averaged using one-vs-all comparison. Note however that the AUC criteria may not be particularly insightful as the prediction threshold we use in sPLS-DA differs from an AUC threshold (sPLS-DA relies on prediction distances for predictions, see \code{?predict.splsda}).

Constraint model. To improve performance, a new variant in all tuning functions was recently added to fit a constraint model (\code{constraint=TRUE)}. The tuning is performed on the optimal list of selected features \code{keepX.constraint} from the previous components, as opposed to considering only the number of features \code{keepX}. Such strategy was implemented in the sister package bootPLS and successfully applied in our recent integrative study (Rohart et al., 2016). Our experience has shown that the constraint tuning and models can substantially improve the performance of the methods. 

}

\value{
Depending on the type of analysis performed, a list that may contain:
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
On multilevel analysis:
Liquet, B., Le Cao, K.-A., Hocini, H. and Thiebaut, R. (2012) A novel approach for
biomarker selection and the integration of repeated measures experiments from two
platforms. \emph{BMC Bioinformatics} \bold{13}:325.

Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K. (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA. \emph{Metabolomics}, \bold{6}(1), 119-128.

}

\author{Benoit Liquet, Kim-Anh Le Cao, Florian Rohart.}

\seealso{\code{\link{tune.splsda}} and http://www.mixOmics.org for more details.}

\examples{
## First example: one-factor analysis with sPLS-DA
\dontrun{
  data(vac18.simulated) # simulated data
  design <- data.frame(sample = vac18.simulated$sample)
  
    result.ex1 = tune.multilevel(X = vac18.simulated$genes,
                                Y = vac18.simulated$stimulation,
                               multilevel = design,
                               ncomp=2,
                               test.keepX=c(5, 10, 15), 
                               already.tested.X = c(50),
                               method = 'splsda',
                               dist = 'mahalanobis.dist',
                               validation = 'loo') 
  
  # overall error rate
  result.ex1$error.rate 
  # classification error rate per class after 2 components
  result.ex1$error.rate.class
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
  
  design <- data.frame(sample = repeat.indiv)

  result.ex2 = tune.multilevel(liver.toxicity$gene,
                                Y = data.frame(dose, time),
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

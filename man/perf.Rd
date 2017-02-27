\name{perf}
\encoding{latin1}
\alias{perf}
\alias{perf.pls}
\alias{perf.spls}
\alias{perf.plsda}
\alias{perf.splsda}
\alias{perf.mint.splsda}
\alias{perf.sgccda}

\title{Compute evaluation criteria for PLS, sPLS, PLS-DA, sPLS-DA, MINT and DIABLO}

\description{
Function to evaluate the performance of the fitted PLS, sparse PLS,
PLS-DA, sparse PLS-DA, MINT (mint.splsda) and DIABLO (block.splsda) models using various criteria.
}

\usage{
           
\method{perf}{pls}(object,validation = c("Mfold", "loo"),
           folds = 10, progressBar = TRUE, \ldots)	

\method{perf}{spls}(object,validation = c("Mfold", "loo"),
          folds = 10, progressBar = TRUE, \ldots)

\method{perf}{plsda}(object,
          dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
          constraint = FALSE, validation = c("Mfold", "loo"), 
          folds = 10, nrepeat =1, auc = FALSE, progressBar = TRUE, cpus, \ldots)

\method{perf}{splsda}(object, 
dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), constraint = FALSE,
validation = c("Mfold", "loo"),
folds = 10, nrepeat =1, auc = FALSE, progressBar = TRUE, cpus, \ldots)

\method{perf}{mint.splsda}(object,
            dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
            constraint = TRUE, auc = FALSE, progressBar = TRUE, \ldots)

\method{perf}{sgccda}(object,
            dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
            validation = c("Mfold", "loo"),
            folds = 10, nrepeat =1, cpus, \ldots)

}

\arguments{
  \item{object}{object of class inheriting from \code{"pls"}, \code{"plsda"},
    \code{"spls"}, \code{"splsda"} or  \code{"mint.splsda"}. The function will retrieve some key parameters stored in that object.}
\item{dist}{only applies to an object inheriting from \code{"plsda"}, \code{"splsda"} or \code{"mint.splsda"} to evaluate the classification performance of the model. Should be a subset of \code{"max.dist"}, \code{"centroids.dist"}, \code{"mahalanobis.dist"}. Default is \code{"all"}. See \code{\link{predict}}.}
\item{constraint}{Indicate whether the performance of the model should be evaluated w.r.t the number of variables selected (keepX, default constraint = FALSE) or w.r.t the specific variables selected by the model (keepX.constraint, constraint = TRUE). }
\item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or
    \code{"loo"} (see below). Default is \code{"Mfold"}.}
  \item{folds}{the folds in the Mfold cross-validation. See Details.}
\item{nrepeat}{Number of times the Cross-Validation process is repeated.}
\item{auc}{if \code{TRUE} calculate the Area Under the Curve (AUC) performance of the model.}
  \item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
\item{cpus}{Number of cpus to use when running the code in parallel.}
\item{...}{not used}
}

\details{
The process of evaluating the performance of a fitted model \code{object} is similar for all PLS-derived methods; a cross-validation approach is used to fit the method of \code{object} on \code{folds-1} subsets of the data and then to predict on the subset left out. Different measures of performance are available depending on the model. Parameters such as \code{logratio}, \code{multilevel}, \code{keepX} or \code{keepY} are retrieved from \code{object}.

If \code{validation = "Mfold"}, M-fold cross-validation is performed.
\code{folds} specifies the number of folds to generate.
The folds also can be supplied as a list of vectors containing the indexes defining each
fold as produced by \code{split}. When using \code{validation = "Mfold"}, make sure that you repeat the process several times (as the results will be highly dependent on the random splits and the sample size).

If \code{validation = "loo"}, leave-one-out cross-validation is performed (in that case, there is no need to repeat the process).


For fitted PLS and sPLS regression models, \code{perf} estimates the 
mean squared error of prediction (MSEP), \eqn{R^2}, and \eqn{Q^2} to assess the predictive 
perfity of the model using M-fold or leave-one-out cross-validation. Note that only the \code{classic}, \code{regression} and  \code{invariant} modes can be applied.
For sPLS, the MSEP, \eqn{R^2}, and \eqn{Q^2} criteria are averaged across all folds.
Note that for PLS and sPLS objects, perf is performed on the pre-processed data after log ratio transform and multilevel analysis, if any.


The sPLS, sPLS-DA and sgccda functions are run on several and different subsets of data (the cross-folds) and will certainly lead to different subset of selected features. Those are summarised in the output \code{features$stable} (see output Value below) to assess how often the variables are selected across all folds.
Note that for PLS-DA and sPLS-DA objects, perf is performed on the original data, i.e. before the pre-processing step of the log ratio transform and multilevel analysis, if any. In addition for these methods, the classification error rate is averaged across all folds.

The mint.sPLS-DA function estimates errors based on Leave-one-group-out cross validation (where each levels of object$study is left out (and predicted) once) and provides study-specific outputs (\code{study.specific.error}) as well as global outputs (\code{global.error}).

For PLS-DA, sPLS-DA, mint.PLS-DA and mint.sPLS-DA methods: if \code{auc=TRUE}, Area Under the Curve (AUC) values are calculated from the predicted values obtained from the cross-validation process, either for all samples or for study-specific samples (for mint models). See \code{\link{auroc}} for more details.

}

\value{
For PLS and sPLS models, \code{perf} produces a list with the following components: 
  \item{MSEP}{Mean Square Error Prediction for each \eqn{Y} variable, only applies to object inherited from \code{"pls"}, and \code{"spls"}.}
  \item{R2}{a matrix of \eqn{R^2} values of the \eqn{Y}-variables for models 
    with \eqn{1, \ldots ,}\code{ncomp} components, only applies to object inherited from \code{"pls"}, and \code{"spls"}.}
  \item{Q2}{if \eqn{Y} containts one variable, a vector of \eqn{Q^2} values else a list with 
    a matrix of \eqn{Q^2} values for each \eqn{Y}-variable. Note that in the specific case of an sPLS model, it is better to have a look at the Q2.total criterion, only applies to object inherited from \code{"pls"}, and \code{"spls"}}	
    \item{Q2.total}{a vector of \eqn{Q^2}-total values for models with \eqn{1, \ldots ,}\code{ncomp} components, only applies to object inherited from \code{"pls"}, and \code{"spls"}}
    \item{features}{a list of features selected across the folds (\code{$stable.X} and \code{$stable.Y}) for the \code{keepX} and \code{keepY} parameters from the input object.}
    \item{error.rate}{
For PLS-DA and sPLS-DA models, \code{perf} produces a matrix of classification error rate estimation. 
The dimensions correspond to the components in the model and to the prediction method used, respectively. Note that error rates reported in any component include the performance of the model in earlier components for the specified \code{keepX} parameters (e.g. error rate reported for component 3 for \code{keepX = 20} already includes the fitted model on components 1 and 2 for  \code{keepX = 20}). For more advanced usage of the \code{perf} function, see \url{www.mixomics.org/methods/spls-da/} and consider using the \code{predict} function.}
\item{auc}{Averaged AUC values over the \code{nrepeat}}

For mint.splsda models, \code{perf} produces the following outputs:
\item{study.specific.error}{A list that gives BER, overall error rate and error rate per class, for each study}
\item{global.error}{A list that gives BER, overall error rate and error rate per class for all samples}
\item{predict}{A list of length \code{ncomp} that produces the predicted values of each sample for each class}
\item{class}{A list which gives the predicted class of each sample for each \code{dist} and each of the \code{ncomp} components. Directly obtained from the \code{predict} output.}
\item{auc}{AUC values}
\item{auc.study}{AUC values for each study}



For sgccda models, \code{perf} produces the following outputs:
\item{error.rate}{Prediction error rate for each block of \code{object$X} and each \code{dist}}
\item{error.rate.per.class}{Prediction error rate for each block of \code{object$X}, each \code{dist} and each class}
\item{predict}{Predicted values of each sample for each class, each block and each component}
\item{class}{Predicted class of each sample for each block, each \code{dist}, each component and each nrepeat}
\item{features}{a list of features selected across the folds (\code{$stable.X} and \code{$stable.Y}) for the \code{keepX} and \code{keepY} parameters from the input object.}
\item{AveragedPredict.class}{if more than one block, returns the average predicted class over the blocks (averaged of the \code{Predict} output and prediction using the \code{max.dist} distance)}
\item{AveragedPredict.error.rate}{if more than one block, returns the average predicted error rate over the blocks (using the \code{AveragedPredict.class} output)}
\item{WeightedPredict.class}{if more than one block, returns the weighted predicted class over the blocks (weighted average of the \code{Predict} output and prediction using the \code{max.dist} distance)}
\item{WeightedPredict.error.rate}{if more than one block, returns the weighted average predicted error rate over the blocks (using the \code{WeightedPredict.class} output)}
\item{MajorityVote}{if more than one block, returns the majority class over the blocks. NA for a sample means that there is no consensus on the predicted class for this particular sample over the blocks.}
\item{MajorityVote.error.rate}{if more than one block, returns the error rate of the \code{MajorityVote} output}
\item{WeightedVote}{if more than one block, returns the weighted majority class over the blocks. NA for a sample means that there is no consensus on the predicted class for this particular sample over the blocks.}
\item{WeightedVote.error.rate}{if more than one block, returns the error rate of the \code{WeightedVote} output}
\item{weights}{Returns the weights of each block used for the weighted predictions, for each nrepeat and each fold}
}

\references{
Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Le Cao K.A. (2016).
DIABLO - multi omics integration for biomarker discovery.

Rohart F. et al (2016, submitted). MINT: A multivariate integrative approach to identify a reproducible biomarker signature across multiple experiments and platforms.

Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris: Editions Technic.

Chavent, Marie and Patouille, Brigitte (2003). Calcul des coefficients de r{e}gression et du PRESS en r{e}gression PLS1. \emph{Modulad n}, \bold{30} 1-11. (this is the formula we use to calculate the Q2 in perf.pls and perf.spls)

Le Cao, K. A., Rossouw D., Robert-Granie, C. and Besse, P. (2008). A sparse PLS for variable 
selection when integrating Omics data. \emph{Statistical Applications in Genetics and Molecular 
Biology} \bold{7}, article 35.

Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction (MSEP) Estimates for Principal Component 
Regression (PCR) and Partial Least Squares Regression (PLSR). \emph{Journal of Chemometrics} \bold{18}(9), 422-429.

}

\author{Ignacio Gonzalez, Amrit Singh, Kim-Anh Le Cao, Benoit Gautier, Florian Rohart.}

\seealso{\code{\link{predict}}, \code{\link{nipals}}, \code{\link{plot.perf}}, \code{\link{auroc}} and \url{www.mixOmics.org} for more details.}

\examples{
\dontrun{
## validation for objects of class 'pls' (regression)
# ----------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic


# try tune the number of component to choose
# ---------------------
# first learn the full model
liver.pls <- pls(X, Y, ncomp = 10)

# with 5-fold cross validation: we use the same parameters as in model above
# but we perform cross validation to compute the MSEP, Q2 and R2 criteria
# ---------------------------
liver.val <- perf(liver.pls, validation = "Mfold", folds = 5)

# Q2 total should decrease until it reaches a threshold
liver.val$Q2.total

# ncomp = 2 is enough
plot(liver.val$Q2.total, type = 'l', col = 'red', ylim = c(-0.5, 0.5),
xlab = 'PLS components', ylab = 'Q2 total')
abline(h = 0.0975, col = 'darkgreen')
legend('topright', col = c('red', 'darkgreen'),
legend = c('Q2 total', 'threshold 0.0975'), lty = 1)
title('Liver toxicity PLS 5-fold, Q2 total values')

#have a look at the other criteria
# ----------------------
# R2
liver.val$R2
matplot(t(liver.val$R2), type = 'l', xlab = 'PLS components', ylab = 'R2 for each variable')
title('Liver toxicity PLS 5-fold, R2 values')

# MSEP
liver.val$MSEP
matplot(t(liver.val$MSEP), type = 'l', xlab = 'PLS components', ylab = 'MSEP for each variable')
title('Liver toxicity PLS 5-fold, MSEP values')


## validation for objects of class 'spls' (regression)
# ----------------------------------------
ncomp = 7
# first, learn the model on the whole data set
model.spls = spls(X, Y, ncomp = ncomp, mode = 'regression',
keepX = c(rep(10, ncomp)), keepY = c(rep(4,ncomp)))


# with leave-one-out cross validation
##set.seed(45)
model.spls.val <- perf(model.spls, validation = "Mfold", folds = 5 )#validation = "loo")

#Q2 total
model.spls.val$Q2.total

# R2:we can see how the performance degrades when ncomp increases
model.spls.val$R2
plot(model.spls.val, criterion="R2", type = 'l')
plot(model.spls.val, criterion="Q2", type = 'l')


## validation for objects of class 'splsda' (classification)
# ----------------------------------------
data(srbct)
X <- srbct$gene
Y <- srbct$class

ncomp = 5

srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))

# with Mfold
# ---------
set.seed(45)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
dist = "all", auc = TRUE)
error
error$auc

plot(error, type = "l")

# parallel code
set.seed(45)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
dist = "all", auc = TRUE, cpus =2)


## validation for objects of class 'mint.splsda' (classification)
# ----------------------------------------

data(stemcells)
res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3, keepX = c(10, 5, 15),
study = stemcells$study, auc = TRUE)

out = perf(res)
out

out$auc
out$auc.study

## validation for objects of class 'sgccda' (classification)
# ----------------------------------------

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- block.splsda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,
scheme = "horst",
verbose = FALSE,
bias = FALSE)

perf = perf(nutrimouse.sgccda)
perf
perf$WeightedVote.error.rate
}

}

\keyword{regression}
\keyword{multivariate}

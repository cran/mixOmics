\name{tune.mint.splsda}
\encoding{latin1}
\alias{tune.mint.splsda}

\title{Estimate the parameters of mint.splsda method}

\description{
Computes Leave-One-Group-Out-Cross-Validation (LOGOCV) scores on a user-input
grid to determine optimal values for the sparsity parameters in \code{mint.splsda}.
}

\usage{
tune.mint.splsda(X, Y, ncomp = 1, study, test.keepX = c(5, 10, 15), already.tested.X,
constraint = TRUE, dist = "max.dist", measure = "BER", auc = FALSE, progressBar = TRUE,
scale = TRUE, tol = 1e-06, max.iter = 100, near.zero.var = FALSE, light.output = TRUE )
}

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{Outcome. Numeric vector or matrix of responses (for multi-response models)}
\item{ncomp}{Number of components to include in the model (see Details). Default to 1}
\item{study}{grouping factor indicating which samples are from the same study}
\item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
\item{already.tested.X}{if \code{ncomp > 1} A list containing which variables of X are to be kept on each of the first PLS-components (if \code{constraint=TRUE}) or a numeric vector indicating the number of variables to select from the \eqn{X} data set on the firsts components (if \code{constraint = FALSE})}
\item{constraint}{Indicate whether the performance of the model should be evaluated w.r.t the number of variables selected (keepX, default constraint = FALSE) or w.r.t the specific variables selected by the model (keepX.constraint, constraint = TRUE), see examples. }
\item{dist}{only applies to an object inheriting from \code{"plsda"} or \code{"splsda"} to evaluate the classification performance of the model. Should be a subset of \code{"max.dist"}, \code{"centroids.dist"}, \code{"mahalanobis.dist"}. Default is \code{"all"}. See \code{\link{predict}}.}
\item{measure}{Two misclassification measure are available: overall misclassification error \code{overall} or the Balanced Error Rate \code{BER}}
\item{auc}{if \code{TRUE} calculate the Area Under the Curve (AUC) performance of the model.}
\item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}
\item{light.output}{if set to FALSE, the prediction/classification of each sample for each of \code{test.keepX} and each comp is returned.}
}


\details{
This function performs a Leave-One-Group-Out-Cross-Validation (LOGOCV), where each of \code{study} is left out once. It returns a list of variables of \code{X} that were selected on each of the \code{ncomp} components.
Then, a \code{\link{mint.splsda}} can be performed with \code{keepX.constraint} set as the output \code{choice.keepX.constraint}.

All component \eqn{1:\code{ncomp}} are tuned, except the first ones for which a \code{already.tested.X} is provided. See examples below.


}

\value{
The returned value is a list with components: 
\item{error.rate}{returns the prediction error for each \code{test.keepX} on each component}
\item{choice.keepX.constraint}{returns the variables selected by the optimal keepX on each component.}
\item{error.rate.class}{returns the error rate for each level of \code{Y} and for each component computed with the optimal keepX.constraint}

  \item{predict}{Prediction values for each sample, each \code{test.keepX} and each comp.}
  \item{class}{Predicted class for each sample, each \code{test.keepX} and each comp.}
}

\references{
Rohart F. et al (2016, submitted). MINT: A multivariate integrative approach to identify a reproducible biomarker signature across multiple experiments and platforms.
}


\author{Florian Rohart}

\seealso{\code{\link{mint.splsda}} and http://www.mixOmics.org for more details.}

\examples{
data(stemcells)
data = stemcells$gene
type.id = stemcells$celltype
exp = stemcells$study

res = mint.splsda(X=data,Y=type.id,ncomp=3,keepX=c(10,5,15),study=exp)

\dontrun{
out = tune.mint.splsda(X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,
    study=exp,test.keepX=seq(1,10,1), constraint = TRUE)
out$choice.keepX.constraint

## only tune component 2 and keeping
# "ENSG00000159199" "ENSG00000106012" "ENSG00000129317" on the first component
out = tune.mint.splsda(X=data,Y=type.id,ncomp=2, study=exp,
already.tested.X = list(comp1 = c("ENSG00000159199", "ENSG00000106012", "ENSG00000129317")),
test.keepX=seq(1,10,1), constraint = TRUE)
out$choice.keepX.constraint
out$choice.keepX
# gives the number of variables selected (which comes from out$choice.keepX.constraint)

}

}

\keyword{multivariate}
\keyword{dplot}

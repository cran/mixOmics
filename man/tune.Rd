\name{tune}
\encoding{latin1}
\alias{tune}


\title{Overall tuning function that can be used to tune several methods}

\description{
Wrapper of all the tuning functions.
}

\usage{
tune(method,
X,
Y,
multilevel,
ncomp,
study, # mint.splsda
test.keepX = c(5, 10, 15), # all but pca, rcc
test.keepY = NULL, # rcc, multilevel
already.tested.X, # all but pca, rcc
already.tested.Y, #multilevel
constraint,
mode = "regression", # multilevel
nrepeat = 1, #multilevel, splsda
grid1 = seq(0.001, 1, length = 5), # rcc
grid2 = seq(0.001, 1, length = 5), # rcc
validation = "Mfold", # all but pca
folds = 10, # all but pca
dist = "max.dist", # all but pca, rcc
measure = c("BER"), # all but pca, rcc
auc = FALSE,
progressBar = TRUE, # all but pca, rcc
near.zero.var = FALSE, # all but pca, rcc
logratio = "none", # all but pca, rcc
center = TRUE, # pca
scale = TRUE, # mint, splsda
max.iter = 100, #pca
tol = 1e-09,
light.output = TRUE # mint, splsda
)

}	

\arguments{
\item{method}{This parameter is used to pass all other argument to the suitable function. \code{method} has to be one of the following:
"spls", "splsda", "mint.splsda", "rcc", "pca".
}
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{Either a factor or a class vector for the discrete outcome, or a numeric vector or matrix of continuous responses (for multi-response models).}
\item{multilevel}{Design matrix for multilevel anaylis (for repeated measurements) that indicates the repeated measures on each individual, i.e. the individuals ID. See Details.}
\item{ncomp}{the number of components to include in the model.}
\item{study}{grouping factor indicating which samples are from the same study}
\item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
\item{test.keepY}{If \code{method = 'spls'}, numeric vector for the different number of variables to test from the \eqn{Y} data set}
\item{already.tested.X}{if \code{ncomp > 1} A list containing which variables of X are to be kept on each of the first PLS-components (if \code{constraint=TRUE}) or a numeric vector indicating the number of variables to select from the \eqn{X} data set on the firsts components (if \code{constraint = FALSE})}
\item{already.tested.Y}{if \code{method = 'spls'} and \code{if(ncomp > 1)} numeric vector indicating the number of variables to select from the \eqn{Y} data set on the first components}
\item{constraint}{Indicate whether the performance of the model is evaluated in regards to the number of variables selected (keepX, constraint = FALSE) or which specific variables are selected (keepX.constraint, constraint = TRUE)}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details.}
\item{nrepeat}{Number of times the Cross-Validation process is repeated.}
\item{grid1, grid2}{vector numeric defining the values of \code{lambda1} and \code{lambda2}
at which cross-validation score should be computed. Defaults to
\code{grid1=grid2=seq(0.001, 1, length=5)}.}
\item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or
\code{"loo"} (see below). Default is \code{"Mfold"}.}
\item{folds}{the folds in the Mfold cross-validation. See Details.}
\item{dist}{distance metric to use for \code{splsda} to estimate the classification error rate,
should be a subset of \code{"centroids.dist"}, \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details).}
\item{measure}{Two misclassification measure are available: overall misclassification error \code{overall} or the Balanced Error Rate \code{BER}}
\item{auc}{if \code{TRUE} calculate the Area Under the Curve (AUC) performance of the model.}
\item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}
\item{logratio}{one of ('none','CLR'). Default to 'none'}

\item{center}{a logical value indicating whether the variables should be shifted to be zero centered.
Alternately, a vector of length equal the number of columns of \code{X} can be supplied.
The value is passed to \code{\link{scale}}.}
\item{scale}{a logical value indicating whether the variables should be scaled to have
unit variance before the analysis takes place. The default is \code{FALSE} for consistency with \code{prcomp}
function, but in general scaling is advisable. Alternatively, a vector of length equal the number of
columns of \code{X} can be supplied. The value is passed to \code{\link{scale}}.}
\item{max.iter}{integer, the maximum number of iterations for the NIPALS algorithm.}
\item{tol}{a positive real, the tolerance used for the NIPALS algorithm.}
\item{light.output}{if set to FALSE, the prediction/classification of each sample for each of \code{test.keepX} and each comp is returned.}

}

\details{
For more details, please see the help file corresponding to your \code{method}.
Note that only the arguments used in the tune function corresponding to \code{method} are passed on.
}

\value{
Depending on the type of analysis performed, a list that contains:
  \item{error}{cross-validation overall error rate when one-factor sPLS-DA analysis is performed.}
  \item{predict}{cross-validation prediction for all samples and for the LAST keepX tested parameter when one-factor sPLS-DA analysis is performed.}
  \item{cor.value}{compute the correlation between latent variables for two-factor sPLS-DA analysis or sPLS.}

}


\author{Florian Rohart}

\seealso{\code{\link{tune.rcc}}, \code{\link{tune.mint.splsda}},
\code{\link{tune.pca}}, \code{\link{tune.multilevel}},
\code{\link{tune.splsda}}, \code{\link{tune.splslevel}} and http://www.mixOmics.org for more details.}

\examples{
## sPLS-DA
\dontrun{
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- as.factor(breast.tumors$sample$treatment)
tune= tune(method = "splsda", X, Y, ncomp=1, nrepeat=10, logratio="none",
    test.keepX = c(5, 10, 15), folds=10, dist="max.dist", progressBar = TRUE)

plot(tune)
}


## mint.splsda
\dontrun{
data(stemcells)
data = stemcells$gene
type.id = stemcells$celltype
exp = stemcells$study

out = tune(method="mint.splsda", X=data,Y=type.id, ncomp=2, study=exp, test.keepX=seq(1,10,1))
out$choice.keepX.constraint

plot(out)

}



}

\keyword{regression}
\keyword{multivariate}

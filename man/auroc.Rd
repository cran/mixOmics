\name{auroc}
\encoding{latin1}
\alias{auroc}
\alias{auroc.plsda}
\alias{auroc.splsda}
\alias{auroc.mint.plsda}
\alias{auroc.mint.splsda}
\alias{auroc.sgccda}


\title{Area Under the Curve (AUC) and Receiver Operating Characteristic (ROC) curves for supervised classification}

\description{
Calculates the AUC and plots ROC for supervised objects from s/plsda, mint.s/plsda and block.plsda, block.splsda or wrapper.sgccda.
}

\usage{
\method{auroc}{plsda}(object, newdata = object$X, outcome.test = as.factor(object$Y),
multilevel = NULL, plot = TRUE, roc.comp = 1, ...)

\method{auroc}{splsda}(object, newdata = object$X, outcome.test = as.factor(object$Y),
multilevel = NULL, plot = TRUE, roc.comp = 1, ...)

\method{auroc}{mint.plsda}(object, newdata = object$X, outcome.test = as.factor(object$Y),
study.test = object$study, multilevel = NULL,plot = TRUE, roc.comp = 1, ...)

\method{auroc}{mint.splsda}(object, newdata = object$X, outcome.test = as.factor(object$Y),
study.test = object$study, multilevel = NULL,plot = TRUE, roc.comp = 1, ...)

\method{auroc}{sgccda}(object, newdata = object$X, outcome.test = as.factor(object$Y),
multilevel = NULL,plot = TRUE, roc.block = 1, roc.comp = 1, ...)

}	

\arguments{
\item{object}{Object from one of the following supervised analysis class: "plsda", "splsda", "mint.plsda", "mint.splsda", "block.splsda" or "wrapper.sgccda"}
\item{newdata}{numeric matrix of predictors, by default set to the training data set (see details).}
\item{outcome.test}{Either a factor or a class vector for the discrete outcome, by default set to the outcome vector from the training set (see details). }
\item{study.test}{For MINT objects, grouping factor indicating which samples of \code{newdata} are from the same study. Overlap with \code{object$study} are allowed.}
\item{multilevel}{Sample information when a newdata matrix is input and when multilevel decomposition for repeated measurements is required. A numeric matrix or data frame indicating the repeated measures on each individual, i.e. the individuals ID. See examples in \code{splsda}.}
\item{plot}{Whether the ROC curves should be plotted, by default set to TRUE (see details).}
\item{roc.comp}{Specify the component (integer) for which the ROC will be plotted from the multivariate model, default to 1.}
\item{roc.block}{Specify the block number (integer) or the name of the block (set of characters) for which the ROC will be plotted for a block.plsda or block.splsda object, default to 1.}
\item{...}{external optional arguments for plotting}
}

\details{
For more than two classes in the categorical outcome Y, the AUC is calculated as one class vs. the other and the ROC curves one class vs. the others are output.

The ROC and AUC are calculated based on the predicted scores obtained from the multivariate methods. Oour multivariate supervised methods already use a prediction threshold based on distances (see \code{predict}) that optimally determine class membership of the samples tested. As such AUC and ROC are not needed to estimate the performance of the model (see \code{perf}, \code{tune} that report classification error rates). We provide those outputs as complementary performance measures. 

The pvalue is from a Wilcoxon test between the predicted scores between one class vs the others. 

External independent data set could be input with the same number of predictors than in the training data set. In that case the multivariate model predicts the predicted scores on the test data (based on the trained model stored in object) and compares the prediction to the available and known outcome vector from the external independent data set
}

\value{
Depending on the type of object used, a list that contains:
The AUC and Wilcoxon test pvalue for each 'one vs other' classes comparison performed, either per component (splsda, plsda, mint.plsda, mint.splsda), or per block and per component (wrapper.sgccda, block.plsda, blocksplsda).
}


\author{Benoit Gautier, Francois Bartolo, Florian Rohart}

\seealso{\code{\link{tune}}, \code{\link{perf}}, and http://www.mixOmics.org for more details.}

\examples{
## example with PLSDA, 2 classes
# ----------------
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)
auc.plsda.breast = auroc(plsda.breast, ncomp = 1)

## example with sPLSDA
# -----------------
splsda.breast <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))
auroc(plsda.breast, plot = FALSE)

## example with sPLSDA with 4 classes
# -----------------
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))
auc.splsda.liver = auroc(splsda.liver, ncomp = 1)


## example with mint.plsda
# -----------------
data(stemcells)

res = mint.plsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3,
study = stemcells$study)
auc.mint.pslda = auroc(res, plot = FALSE)

## example with mint.splsda
# -----------------
res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3, keepX = c(10, 5, 15),
study = stemcells$study)
auc.mint.spslda = auroc(res, plot = TRUE, roc.comp = 3)


## example with block.plsda
# ------------------
data(nutrimouse)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
# with this design, all blocks are connected
design = matrix(c(0,1,1,0), ncol = 2, nrow = 2,
byrow = TRUE, dimnames = list(names(data), names(data)))

block.plsda.nutri = block.plsda(X = data, Y = nutrimouse$diet)
auc.block.plsda.nutri = auroc(block.plsda.nutri, block = 'lipid')

## example with block.splsda
# ---------------
list.keepX = list(gene = rep(10, 2), lipid = rep(5,2))
block.splsda.nutri = block.splsda(X = data, Y = nutrimouse$diet, keepX = list.keepX)
auc.block.splsda.nutri = auroc(block.splsda.nutri, block = 1)
}

\keyword{regression}
\keyword{multivariate}

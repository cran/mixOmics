\name{multilevel}
\encoding{latin1}
\alias{multilevel}
\alias{multilevel.splsda}
\alias{multilevel.spls}

\title{Multilevel analysis for repeated measurements 9cross-over design), either with a discriminant or an integrative analysis}

\description{
The analysis of repeated measurements is performed by combining a multilevel approach with multivariate approaches: sPLS-DA (Discriminant Analysis) or sPLS (Integrative analysis). Both approaches embbed variable selection.
}

\usage{
multilevel(X, Y = NULL, 
           cond = NULL,                  
           sample = NULL,
           ncomp = 1,
           keepX = rep(ncol(X), ncomp),
           keepY = NULL,  # for spls
           method = NULL,
           tab.prob.gene=NULL, 
           max.iter = 500, 
           tol = 1e-06,...)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{\code{if(method = 'spls')} numeric vector or matrix of continuous responses (for multi-response models) \code{NA}s are allowed.}
  \item{cond}{a factor or a class vector for one-factor discrete outcome, a numeric matrix of 2 columns for two-factor discrete outcome. See Details}
  \item{sample}{a vector indicating the repeated measured on each individual.}
  \item{ncomp}{the number of components to include in the model (see Details).}
  \item{keepX}{numeric vector of length \code{ncomp}, the number of variables
    to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
  \item{keepY}{\code{if(method = 'spls')}, numeric vector of length \code{ncomp}, the number of variables
    to keep in \eqn{Y}-loadings. By default all variables are kept in the model.}
  \item{method}{character string. Which multivariate method and type of analysis to choose, matching one of 'splsda' (Discriminant Analysis) or 'spls' (unsupervised integrative analysis). See Details.}
    \item{tab.prob.gene}{matrix of dimension \code{ncol(X)}x2 indicating the probes names (column 1) and corresponding gene names (column 2).}
  \item{max.iter}{integer, the maximum number of iterations.}   
  \item{tol}{a positive real, the tolerance used in the iterative algorithm.} 
  \item{...}{other internal arguments.}
}

\details{
\code{multilevel} function first decomposes the variance in the data sets X (and Y) and applies either sPLS-DA (\code{method = 'splsda'}) or sPLS ((\code{method = 'splsda'})) on the within-subject deviation. 

One or two-factor analyses are available for (\code{method = 'splsda'}).

A sPLS-DA or sPLS model is performed with \eqn{1, \ldots ,}\code{ncomp} components
to the factor or class vector \code{cond}. The appropriate indicator matrix is created for sPLS-DA.

Multilevel sPLS-DA enables the selection of discriminant variables between the conditions in \code{cond}.

Multilevel sPLS enables the integration of data measured on two different platforms on the same individuals. This approach differs from multilevel sPLS-DA as the aim is to select subsets variables from both data sets that are highly correlated (positively or negatively) across the samples. The approach is unsupervised: no prior knowledge about the samples groups is included.
}

\value{
For sPLS-DA analysis, \code{multilevel} returns either an object of class \code{"splsda1fact"} for one-factor analysis and \code{"splsda2fact"}, a list that contains the following sPLS-DA outputs:
  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized indicator response vector or matrix.}
  \item{ind.mat}{the indicator matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{keepX}{number of \eqn{X} variables kept in the model on each component.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the variates.}
  \item{loadings}{list containing the estimated loadings for the \code{X} and 
	\code{Y} variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}
  
  
  For sPLS analysis, \code{multilevel} returns either an object of class \code{"splslevel"} for one-factor analysis, a list that contains the following sPLS-DA components:
    \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized original response vector or matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{mode}{the algorithm used to fit the model.}
  \item{keepX}{number of \eqn{X} variables kept in the model on each component.}
  \item{keepY}{number of \eqn{Y} variables kept in the model on each component.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the variates.}
  \item{loadings}{list containing the estimated loadings for the \eqn{X} and 
  \eqn{Y} variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}

}

\references{
On multilevel analysis:
Liquet, B., Lê Cao, K.-A., Hocini, H. and Thiebaut, R. A novel approach for biomarker selection and the
integration of repeated measures experiments from two platforms. \emph{Submitted}.

Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K. (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA. \emph{Metabolomics}, \bold{6}(1), 119-128.

On sPLS-DA:
Lê Cao, K.-A., Boitard, S. and Besse, P. (2011). Sparse PLS Discriminant Analysis: biologically relevant feature selection and graphical displays for multiclass problems. \emph{BMC Bioinformatics} \bold{12}:253.

On sPLS:
Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009). Sparse canonical 
methods for biological data integration: application to a cross-platform study. 
\emph{BMC Bioinformatics} \bold{10}:34.

Lê Cao, K.-A., Rossouw, D., Robert-Granie, C. and Besse, P. (2008). A sparse PLS for variable 
selection when integrating Omics data. \emph{Statistical Applications in Genetics and Molecular 
Biology} \bold{7}, article 35.

}

\author{Benoit Liquet, Kim-Anh Lê Cao.}

\seealso{\code{\link{spls}}, \code{\link{splsda}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}},
\code{\link{cim}}, \code{\link{network}}.}

\examples{
## First example: one-factor analysis with sPLS-DA
data(data.simu) # simulated data
res.1level <- multilevel(data.simu$X, cond = data.simu$stimu, sample=data.simu$sample, ncomp=3,
                         keepX=c(200,200,200),tab.prob.gene=NULL, method = 'splsda')


plot3dIndiv(res.1level, col = as.numeric(data.simu$stimu), cex = 0.6)


## Second example: two-factor analysis with sPLS-DA
data(data.simu) # simulated data

time = factor(rep(c(rep('t1', 6), rep('t2', 6)), 4))
stimu.time = data.frame(cbind(as.character(data.simu$stimu), as.character(time)))
repeat.simu2 = rep(c(1:6), 8)

res.2level <- multilevel(data.simu$X, cond = stimu.time,sample=repeat.simu2, 
                        ncomp=2,keepX=c(200, 200),tab.prob.gene=NULL, method = 'splsda')

# color for plotIndiv
col.stimu = as.numeric(data.simu$stimu)
# pch for plots
pch.time = rep(20, 48)
pch.time[time == 't2'] = 4


plotIndiv(res.2level, col = col.stimu, pch = pch.time, ind.names = FALSE)
legend('bottomright', col = unique(col.stimu), legend = levels(data.simu$stimu),  pch = 20, cex = 0.8)
legend('topright', col = 'black', legend = levels(time),  pch = unique(pch.time), cex = 0.8)


## Third example: one-factor integrative analysis with sPLS
data(liver.toxicity)
repeat.indiv = c(1,2, 1,  2,  1,  2,  1,  2,  3,  3,  4,  3,  4,  3,  4,  4,  5,  6,  5,  5,  6,  5,  6,  7,  7,  8,  6,  7,  8,  7,  8,  8,  9, 10,  9, 10, 11, 9,  9, 10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14, 13, 14, 15, 16, 15, 16, 15, 16, 15, 16)


result.rat <- multilevel(X = liver.toxicity$gene, Y=liver.toxicity$clinic, cond = liver.toxicity$treatment$Dose.Group, sample=repeat.indiv, ncomp = 2, keepX= c(50,50), keepY= c(5, 5), method = 'spls') 

# variable plots
plotVar(result.rat,comp=1:2,X.label=TRUE,Y.label=TRUE,cex=c(0.5,0.9)) 

CIM <- cim(result.rat , comp = 1:2, xlab = "genes", ylab = "clinic var", margins = c(5,6), zoom = FALSE)

\dontrun{
network(result.rat, comp = 1:2, threshold = 0.8, 
        Y.names = NULL, keep.var = TRUE,
        color.node = c( "lightcyan","mistyrose"),
        shape.node = c("circle", "rectangle"),
        color.edge = c("red", "green"),
        lty.edge = c("solid", "solid"), lwd.edge = c(1, 1), 
        show.edge.labels = FALSE, interactive = FALSE)
}
}

\keyword{regression}
\keyword{multivariate}

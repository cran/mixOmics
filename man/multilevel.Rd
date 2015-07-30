\name{multilevel}
\encoding{latin1}
\alias{multilevel}
\alias{multilevel.spls}
\alias{multilevel.splsda}

\title{Multilevel analysis for repeated measurements (cross-over design)}

\description{
The analysis of repeated measurements is performed by combining a multilevel 
approach with multivariate methods: sPLS-DA (Discriminant Analysis) or sPLS 
(Integrative analysis). Both approaches enable variable selection.
}

\usage{
multilevel(X, 
           Y = NULL,
           design,  
           ncomp = 2,
           keepX = NULL,
           keepY = NULL,  
           method = c("spls", "splsda"),  
           mode = c("regression", "canonical"), 
           max.iter = 500, 
           tol = 1e-06,
           near.zero.var = TRUE)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  
  \item{Y}{if \code{method = "spls"}, numeric vector or matrix of continuous responses
    (for multi-response models). \code{NA}s are allowed.}
    
  \item{design}{a numeric matrix or data frame. The first column 
    indicates the repeated measures on each individual, i.e. the individuals ID. If \code{  method = 'splsda'}, the 2nd and 3rd columns are factors. If \code{method = 'spls'} then you can choose to only input the repeated measures (column 1) or the 2nd AND 3rd columns to split the variation for a 2 level factor. See Details.}
    
  \item{ncomp}{the number of components to include in the model (see Details).}
  
  \item{keepX}{numeric vector of length \code{ncomp}, the number of variables
    to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
    
  \item{keepY}{if \code{method = "spls"}, numeric vector of length \code{ncomp}, 
    the number of variables to keep in \eqn{Y}-loadings. By default all variables 
    are kept in the model.}
    
  \item{method}{character string. Which multivariate method and type of analysis 
    to choose, matching \code{"spls"}
    (unsupervised integrative analysis) or \code{"splsda"} (discriminant analysis). See Details.}
    
  \item{mode}{character string. What type of algorithm to use, matching 
    \code{"regression"} or \code{"canonical"}. See detals in \code{?pls}.}
    
  \item{max.iter}{integer, the maximum number of iterations.}
  
  \item{tol}{a not negative real, the tolerance used in the iterative algorithm.}

    
    
    \item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations.}
}

\details{
\code{multilevel} function first decomposes the variance in the data sets \eqn{X} 
(and \eqn{Y}) and applies either sPLS-DA (\code{method = "splsda"}) or sPLS 
(\code{method = "spls"}) on the within-subject deviation. 

One or two-factor analyses are available for \code{method = "splsda"}.

A sPLS or sPLS-DA model is performed with 1,\ldots,\code{ncomp} components
to the factor in \code{design[, 2]} (or \code{design[, 2:3]} for two-factor in sPLS-DA).

Multilevel sPLS-DA enables the selection of discriminant variables between the factors 
in \code{design}.

Multilevel sPLS enables the integration of data measured on two different data sets 
on the same individuals. This approach differs from multilevel sPLS-DA as the aim is 
to select subsets of variables from both data sets that are highly positively or negatively correlated across samples. The approach is unsupervised, i.e. no prior knowledge about the sample groups is included.
}

\value{
\code{multilevel} returns either an object of class \code{"mlspls"} for sPLS analysis 
or an object of class \code{"mlsplsda"} for sPLS-DA analysis, a list that contains the
following components:
  \item{X}{the centered and standardized original predictor matrix.}
  
  \item{Y}{the centered and standardized original (or indicator) response vector 
    or matrix.}
    
  \item{Xw}{the within-subject \eqn{X}-deviation matrix.}
  
  \item{Yw}{the within-subject \eqn{Y}-deviation matrix if \code{method = "spls"}.}
  
  \item{design}{the design matrix.}
  
  \item{ind.mat}{the indicator matrix associated to \eqn{Y} if \code{method = "splsda"}.}
  
  \item{ncomp}{the number of components included in the model.}
  
  \item{keepX}{number of \eqn{X} variables kept in the model on each component.}
  
  \item{keepY}{number of \eqn{Y} variables kept in the model on each component \eqn{Y} 
    if \code{method = "spls"}.}
    
  \item{variates}{list containing the \eqn{X}- and \eqn{Y}-variates.}
  
  \item{loadings}{list containing the estimated loadings for the \code{X} and 
    \code{Y} variates.}
    
  \item{names}{list containing the names to be used for individuals and variables.}
  
  \item{mode}{the algorithm used to fit the model if \code{method = "spls"}.}
  
  \item{nzv}{list containing the zero- or near-zero predictors information.}
}

\references{
On multilevel analysis:

Liquet, B., Le Cao, K.-A., Hocini, H. and Thiebaut, R. (2012) A novel approach for
biomarker selection and the integration of repeated measures experiments from two
platforms. \emph{BMC Bioinformatics} \bold{13}:325.

Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K. (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA. \emph{Metabolomics},
\bold{6}(1), 119-128.

On sPLS-DA:

Le Cao, K.-A., Boitard, S. and Besse, P. (2011). Sparse PLS Discriminant Analysis:
biologically relevant feature selection and graphical displays for multiclass problems.
\emph{BMC Bioinformatics} \bold{12}:253.

On sPLS:

Le Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009). Sparse canonical 
methods for biological data integration: application to a cross-platform study. 
\emph{BMC Bioinformatics} \bold{10}:34.

Le Cao, K.-A., Rossouw, D., Robert-Granie, C. and Besse, P. (2008). A sparse PLS for
variable selection when integrating Omics data. \emph{Statistical Applications in
Genetics and Molecular Biology} \bold{7}, article 35.
}

\author{Benoit Liquet, Kim-Anh Le Cao, Benoit Gautier, Ignacio Gonzalez.}

\seealso{\code{\link{spls}}, \code{\link{splsda}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}},
\code{\link{cim}}, \code{\link{network}}.}

\examples{
## First example: one-factor analysis with sPLS-DA, selecting a subset of variables
# as in the paper Liquet et al.
#--------------------------------------------------------------
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)

# multilevel sPLS-DA model
res.1level <- multilevel(X, ncomp = 3, design = design,
                         method = "splsda", keepX = c(30, 137, 123))

# set up colors for plotIndiv
col.stim <- c("darkblue", "purple", "green4","red3")
plotIndiv(res.1level, ind.names = Y, col.per.group = col.stim)


## Second example: two-factor analysis with sPLS-DA, selecting a subset of variables
# as in the paper Liquet et al.
#--------------------------------------------------------------
\dontrun{
data(vac18.simulated) # simulated data

X <- vac18.simulated$genes
design <- data.frame(sample = vac18.simulated$sample,
                     stimu = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

res.2level <- multilevel(X, ncomp = 2, design = design,
                         keepX = c(200, 200), method = 'splsda')

plotIndiv(res.2level, group = design$stimu, ind.names = vac18.simulated$time,
    add.legend = TRUE, style = 'lattice')
}       

## Third example: one-factor analysis with sPLS, selecting a subset of variables
#--------------------------------------------------------------
\dontrun{
data(liver.toxicity)
# note: we made up those data, pretending they are repeated measurements
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                 6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                 10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                 13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each

# this is a spls (unsupervised analysis) so no need to mention any factor in design
# we only perform a one level variation split
design <- data.frame(sample = repeat.indiv) 
res.spls.1level <- multilevel(X = liver.toxicity$gene,
                                       Y=liver.toxicity$clinic,
                                       design = design,
                                       ncomp = 3,
                                       keepX = c(50, 50, 50), keepY = c(5, 5, 5),
                                       method = 'spls', mode = 'canonical')

# set up colors and pch for plotIndiv
col.stimu <- 1:nlevels(design$stimu)

plotIndiv(res.spls.1level, rep.space = 'X-variate', ind.names = FALSE, 
group = liver.toxicity$treatment$Dose.Group,
          pch = 20, main = 'Gene expression subspace',
          add.legend = TRUE)


plotIndiv(res.spls.1level, rep.space = 'Y-variate', ind.names = FALSE,
group = liver.toxicity$treatment$Dose.Group,
          pch = 20, main = 'Clinical measurements ssubpace',
          add.legend = TRUE)
          
plotIndiv(res.spls.1level, rep.space = 'XY-variate', ind.names = FALSE,
group = liver.toxicity$treatment$Dose.Group,
          pch = 20, main = 'Both Gene expression and Clinical subspaces',
          add.legend = TRUE)

}
}

\keyword{regression}
\keyword{multivariate}

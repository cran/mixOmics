\name{wrappers}
\encoding{latin1}
\alias{wrapper.rgcca}
\alias{wrapper.sgcca}
\alias{wrapper.sgccda}

\title{(Generalised Canonical Correlation Analysis }

\description{Wrapper function to perform (sparse) (regularized) Generalised Canonical Correlation Analysis (discriminant analysis), a generalised approach for the integration of multiple datasets.
}

\usage{	
wrapper.rgcca(  
                blocks, 
                design = NULL, 
                ncomp = rep(2, length(blocks)), 
                tau = "optimal",
                scheme = "centroid", 
                scale = TRUE, 
                bias = FALSE, 
                max.iter = 1000,
                tol = .Machine$double.eps, 
                verbose = FALSE, 
                 near.zero.var = FALSE
                 )
                 
wrapper.sgcca(
            blocks, 
            design = NULL, 
            penalty = NULL, 
            ncomp = rep(2, length(blocks)),
            keep = NULL, 
            scheme = "centroid",
            scale = TRUE, 
            bias = FALSE, 
            max.iter = 1000,
            tol = .Machine$double.eps, 
            verbose = FALSE, 
            near.zero.var = FALSE
           )
           
wrapper.sgccda(
               blocks, 
               Y, 
               design = NULL, 
               ncomp = rep(2, length(blocks)),                            
               keep = NULL, 
               scheme = "centroid", 
               scale = TRUE, 
               bias = FALSE, 
               max.iter = 1000,
               tol = .Machine$double.eps, 
               verbose = FALSE, 
               near.zero.var = FALSE
               )

}	

\arguments{
  \item{blocks}{a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets. \code{NA}s are not allowed.}
  \item{Y}{for \code{wrapper.sgccda} only: a factor or a class vector for the discrete outcome}
  \item{design}{numeric matrix of size (number of blocks) x (number of blocks) with only 0 or 1 values. A value of 1 (0) indicates a relationship (no relationship) between the blocks to be modelled. For \code{wrapper.sgccda} the Y outcome should not be added in the design matrix as it will be done internally.}
    \item{ncomp}{numeric vector of length the number of blocks in \code{blocks}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.}
  \item{tau}{for \code{wrapper.rgcca} only: numeric vector of length the number of blocks in \code{data}. Each regularization parameter will be applied on each block and takes the value between 0 (no regularisation) and 1. If tau = "optimal" the shrinkage paramaters are estimated for each block and each dimension using the Schafer and Strimmer (2005) analytical formula.}
    \item{penalty}{for \code{wrapper.sgcca} and \code{wrapper.sgccda} only: numeric vector of length the number of blocks in \code{data}. Each penalty parameter will be applied on each block and takes the value between 0 (no variable selected) and 1 (all variables included).}
    \item{keep}{for \code{wrapper.sgcca} \code{wrapper.sgccda} only: a list of integer values for each block specifying the number of variables to select on each specified component. Each element of the list corresponds to a block and is of length the number of components chosen for that particular block. This argument is an alternative to the argument \code{penalty} (similar to our keepX argument in the sPLS(DA) methods). See example.}
  \item{scheme}{Either "horst", "factorial" or "centroid" (Default: "centroid"), see reference paper.}
  \item{scale}{boleean. If scale = TRUE, each block is standardized
  to zero means and unit variances (default: TRUE)}
  \item{bias}{boleean. A logical value for biaised or unbiaised estimator of the var/cov (defaults to FALSE).}
  \item{max.iter}{integer, the maximum number of iterations.}
   \item{tol}{Convergence stopping value.}
  \item{verbose}{if set to \code{TRUE}, reports progress on computing.}
  \item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE}
}

\details{
	These wrapper functions are improved versions from the functions of the package \pkg{RGCCA}.  
  rGCCA is an unsupervised model is run, sGCCA is a sparse model and sGCC-DA is a supervised sparse model.
  In sGCC-DA the arguments \code{design}, \code{penalty}, \code{keep}, are specified for the blocks data in the input data \code{blocks} only.
}

\value{
\code{wrapper.rgcca}, \code{wrapper.sgcca} and \code{wrapper.sgccada} return an object of class \code{"rgcca"}, \code{"sgcca"} or \code{"sgccada"}, a list 
that contains the following components:

  \item{blocks}{the input data set (as a list).}
  \item{design}{the input design.}
  \item{variates}{the sgcca components.}
  \item{loadings}{the loadings for each block data set (outer weight vector).}
  \item{loadings.star}{the standardised loading vectors.}
  \item{tau}{the input tau parameter.}
  \item{scheme}{the input scheme.}
  \item{ncomp}{the number of components on each block.}
  \item{crit}{the convergence criterion.}
  \item{AVE}{Indicators of model quality based on the
  Average Variance Explained (AVE): AVE(for one block),
  AVE(outer model), AVE(inner model).}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{defl.matrix}{The deflated matrices at the end of the algorithm}
  More details can be found in the references.
              
}

\references{
  Tenenhaus A. and Tenenhaus M., (2011), Regularized
  Generalized Canonical Correlation Analysis,
  Psychometrika, Vol. 76, Nr 2, pp 257-284.
  
  Schafer J. and Strimmer K., (2005), A shrinkage approach
  to large-scale covariance matrix estimation and
  implications for functional genomics. Statist. Appl.
  Genet. Mol. Biol. 4:32.
  
  Tenenhaus A., Phillipe C., Guillemot, V., Le Cao K-A., Grill J., Frouin, V. (2014) Variable Selection For Generalized Canonical Correlation Analysis. Biostatistics, 15(3): 569-83.
  
  O. P. Gunther, H. Shin, R. T. Ng, W. R. McMaster, B. M. McManus, P. A. Keown, S.J. Tebbutt, K-A. Le Cao,  (2014) Novel multivariate methods for integration of genomics and proteomics data: Applications in a kidney transplant rejection study, OMICS: A journal of integrative biology, 18(11), 682-95.
    
    
}

\author{Benoit Gautier, Florian Rohart, Kim-Anh Le Cao and former RGCCA developers Arthur Tenenhaus, Vincent Guillemot.}

\seealso{\code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{wrapper.sgcca}}. \code{\link{wrapper.sgccda}}.}

\examples{
## RGCCA 
# --------------
data(nutrimouse)
# need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, 
                byrow = TRUE, dimnames = list(names(data), names(data)))

nutrimouse.rgcca <- wrapper.rgcca(blocks = data,
                                         design = design,
                                         tau = "optimal",
                                         ncomp = c(2, 2, 1),
                                         scheme = "centroid",
                                         verbose = FALSE)

# blocks should specify the block data set where the sample plot can be performed 
# (ideally when there are >= 2 components!)
# we indicate the diet variable colors.
plotIndiv(nutrimouse.rgcca, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)

# have a look at the looadings
head(nutrimouse.rgcca$loadings[[1]])
head(nutrimouse.rgcca$loadings[[2]])
head(nutrimouse.rgcca$loadings[[3]])


## sGCCA
# -------------
# same data as above but sparse approach

# version 1 using the penalisation penalty criterion
# ---
\dontrun{
nutrimouse.sgcca <- wrapper.sgcca(blocks = data,
                                   design = design,
                                   penalty = c(0.3, 0.5, 1),
                                   ncomp = c(2, 2, 1),
                                   scheme = "centroid",
                                   verbose = FALSE, 
                                   bias = FALSE)

# In plotIndiv we indicate the diet variable colors and the blocks to be plotted 
# (only blocks with comp  >=2!)
plotIndiv(nutrimouse.sgcca, blocks = c(1,2), group = nutrimouse$diet, 
  plot.ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgcca, comp = 1, block = 1)
selectVar(nutrimouse.sgcca, comp = 1, block = 2)

# variable plot on the selected variables
plotVar(nutrimouse.sgcca, col = color.mixo(1:2), cex = c(2,2))
}

# version 2 using the keep penalty criterion (number of variables to keep)
# it is a list per block and per component, need to specify all variables for the 
# Y 'outcome' here 
# (see below for sgccda code, which is more appropriate)
# ----
nutrimouse.sgcca <- wrapper.sgcca(blocks = data,
                                  design = design,
                                  ncomp = c(2, 2, 1),
                                  # for keep: each element of the list corresponds to a block 
                                  # and is of length the # comp per block
                                  keep = list(c(10,10), c(15,15), c(ncol(Y))),
                                  scheme = "centroid",
                                  verbose = FALSE, 
                                  bias = FALSE)


# In plotIndiv we indicate the diet variable colors and the blocks to be plotted 
# (only blocks with comp  >=2!)
plotIndiv(nutrimouse.sgcca, blocks = c(1,2), group = nutrimouse$diet, 
  plot.ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgcca, comp = 1, block = 1)
selectVar(nutrimouse.sgcca, comp = 1, block = 2)


## sGCC-DA
# -------------
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutrimouse.sgccda <- wrapper.sgccda(blocks = data,
                                    Y = Y,
                                    design = design,
                                    keep = list(c(10,10), c(15,15)),
                                    ncomp = c(2, 2, 1),
                                    scheme = "centroid",
                                    verbose = FALSE,
                                    bias = FALSE)

plotIndiv(nutrimouse.sgccda, blocks = c(1,2), group = nutrimouse$diet, 
  plot.ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgccda, comp = 1, block = 1)
selectVar(nutrimouse.sgccda, comp = 1, block = 2)

# variable plot on the selected variables
plotVar(nutrimouse.sgccda, col = color.mixo(1:2), cex = c(2,2))


}


\keyword{multivariate}


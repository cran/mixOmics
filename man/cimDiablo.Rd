\name{cimDiablo}
\encoding{latin1}
\alias{cimDiablo}

\title{Clustered Image Maps (CIMs) ("heat maps") for DIABLO}

\description{
This function generates color-coded Clustered Image Maps (CIMs) ("heat maps") to 
represent "high-dimensional" data sets analysed with DIABLO.
}

\usage{
cimDiablo(object,
ncomp=1,
margins = c(2, 15),
legend.position = "topright",
transpose = FALSE,
row.names = TRUE,
col.names = TRUE,
size.legend = 1.5)
}

\arguments{
  \item{object}{An object of class inheriting from \code{"block.splsda"}.}
\item{ncomp}{positive integer. The similarity matrix is computed based on the variables selected on those specified components. See example. Defaults to \code{comp = 1}.}
\item{margins}{numeric vector of length two containing the margins (see \code{\link{par}(mar)})
	for column and row names respectively.}
\item{legend.position}{position of the legend, one of "bottomright", "bottom", "bottomleft",
"left", "topleft", "top", "topright", "right" and "center". }

\item{transpose}{logical indicating if the matrix should be transposed for plotting.
Defaults to \code{FALSE}.}
\item{row.names, col.names}{logical, should the name of rows and/or columns of \code{mat} be shown? If \code{TRUE}
(defaults) \code{rownames(mat)} and/or \code{colnames(mat)} are used. Possible character vectors with
row and/or column labels can be used.}

\item{size.legend}{size of the legend}
}

\details{ 
This function is a small wrapper of \code{link{cim}} specific to the DIABLO framework.
}

\references{
Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Le Cao K.A. (2016).
DIABLO - multi omics integration for biomarker discovery. BioRxiv \url{http://dx.doi.org/10.1101/067611}

Eisen, M. B., Spellman, P. T., Brown, P. O. and Botstein, D. (1998). 
Cluster analysis and display of genome-wide expression patterns. 
\emph{Proceeding of the National Academy of Sciences of the USA} \bold{95}, 14863-14868.

Weinstein, J. N., Myers, T. G., O'Connor, P. M., Friend, S. H., Fornace Jr., A. J., 
Kohn, K. W., Fojo, T., Bates, S. E., Rubinstein, L. V., Anderson, N. L., 
Buolamwini, J. K., van Osdol, W. W., Monks, A. P., Scudiero, D. A., Sausville, E. A., 
Zaharevitz, D. W., Bunow, B., Viswanadhan, V. N., Johnson, G. S., Wittes, R. E. 
and Paull, K. D. (1997). An information-intensive approach to the molecular 
pharmacology of cancer. \emph{Science} \bold{275}, 343-349.

Gonzalez I., Le Cao K.A., Davis M.J., Dejean S. (2012). Visualising associations between
paired 'omics' data sets. \emph{BioData Mining}; \bold{5}(1).
}

\author{Amrit Singh, Florian Rohart}

\seealso{\code{\link{cim}}, \code{\link{heatmap}},
\code{\link{hclust}}, \code{\link{plotVar}}, 
\code{\link{network}} and 

\url{http://mixomics.org/mixDIABLO/} for more details on all options available.}

\examples{
## default method: shows cross correlation between 2 data sets
#------------------------------------------------------------------
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- block.splsda(X = data,
Y = Y,
design = design,
keepX = list(gene = c(10,10), lipid = c(15,15)),
ncomp = 2,
scheme = "centroid",
verbose = FALSE,
bias = FALSE)

cimDiablo(nutrimouse.sgccda)

}

\keyword{multivariate}
\keyword{iplot}
\keyword{hplot}
\keyword{graphs}
\keyword{cluster}


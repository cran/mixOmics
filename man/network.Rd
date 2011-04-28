\name{network}
\encoding{latin1}
\alias{network.default}
\alias{network.rcc}
\alias{network.pls}
\alias{network.spls}
\alias{network}

\title{Relevance Network for (r)CCA and (s)PLS regression}

\description{
Display relevance associations network for (regularized) canonical correlation analysis 
and (sparse) PLS regression. 
}

\usage{
\method{network}{default}(mat, threshold = 0.5, X.names = NULL, Y.names = NULL,
        color.node = c("white", "white"), 
        shape.node = c("circle", "rectangle"), 
        color.edge = c("blue", "red"), 
        lty.edge = c("solid", "solid"), 
        lwd.edge = c(1, 1), show.edge.labels = FALSE, 
        show.color.key = TRUE, symkey = TRUE, keysize = 1, 
        breaks, interactive = FALSE, alpha = 1, \ldots)

\method{network}{rcc}(object, comp = 1, X.names = NULL, Y.names = NULL, \ldots)

\method{network}{pls}(object, comp = 1, X.names = NULL, Y.names = NULL, 
        keep.var = TRUE, \ldots)

\method{network}{spls}(object, comp = 1, X.names = NULL, Y.names = NULL, 
        keep.var = TRUE, \ldots)		
}

\arguments{
  \item{mat}{numeric matrix of values to be represented.}
  \item{object}{object of class inheriting from \code{"rcc"}, \code{"pls"} or \code{"spls"}.}
  \item{comp}{atomic or vector of positive integers. The components to adequately account 
    for the data association. Defaults to \code{comp = 1}.}
  \item{threshold}{numeric value between \code{0} and \code{1}. The tuning 
    threshold for the relevant associations network (see Details).}
  \item{X.names, Y.names}{character vector containing the names of \eqn{X}- and  
    \eqn{Y}-variables.}
  \item{keep.var}{boolean. If \code{TRUE} only the variables with loadings not zero are plotted
    (as selected by \code{spls}). Defaults to \code{TRUE}.}
  \item{color.node}{vector of length two, the colors of the \eqn{X} and  
    \eqn{Y}	nodes (see Details).}
  \item{shape.node}{character vector of length two, the shape of the \eqn{X} and  
    \eqn{Y} nodes (see Details).}
  \item{color.edge}{vector of colors or character string specifying the colors function to using
     to color the edges (see Details).}
  \item{lty.edge}{character vector of length two, the line type for the edges (see Details).}
  \item{lwd.edge}{vector of length two, the line width of the edges (see Details).}
  \item{show.edge.labels}{logical. If \code{TRUE}, plot association values as edge 
	labels (defaults to \code{FALSE}).}
  \item{show.color.key}{boolean. If \code{TRUE} a color key should be plotted.}
  \item{symkey}{boolean indicating whether the color key should be made symmetric about 0. 
	Defaults to \code{TRUE}.}
  \item{keysize}{numeric value indicating the size of the color key.}
  \item{breaks}{(optional) either a numeric vector indicating the splitting points for 
    binning \code{mat} into colors, or a integer number of break points to be used, in 
	which case the break points will be spaced equally between \code{min(mat)} 
	and \code{max(mat)}.}
  \item{interactive}{logical. If \code{TRUE}, a scrollbar is created to change the
    threshold value interactively (defaults to \code{FALSE}). See Details.}
  \item{alpha}{numeric value. Tuning parameter to enhance the differences between edges corresponding to
    low and high correlations. Defaults to \code{alpha = 1}.}	
  \item{\ldots}{arguments passed to \code{network.default}.}
}

\details{
\code{network} allows to infer large-scale association networks
between the \eqn{X} and \eqn{Y} datasets in \code{rcc} or \code{spls}. The output is a graph where 
each \eqn{X}- and \eqn{Y}-variable corresponds to a node and the edges 
included in the graph portray associations between them.

In \code{rcc}, to identify \eqn{X}-\eqn{Y} pairs showing relevant associations, \code{network} calculate a 
similarity measure between \eqn{X} and \eqn{Y} variables in a pair-wise manner: the 
scalar product value between every pairs of vectors in dimension \code{length(comp)} representing 
the variables \eqn{X} and \eqn{Y} on the axis defined by \eqn{Z_i} with \eqn{i} in \code{comp}, 
where \eqn{Z_i} is the equiangular vector between the \eqn{i}-th \eqn{X} and \eqn{Y} canonical variate. 

In \code{spls}, if \code{object$mode} is \code{regression}, the similarity measure between \eqn{X} and 
\eqn{Y} variables is given by the scalar product value between every pairs of vectors in dimension 
\code{length(comp)} representing the variables \eqn{X} and \eqn{Y} on the axis defined by \eqn{U_i} with 
\eqn{i} in \code{comp}, where \eqn{U_i} is the \eqn{i}-th \eqn{X} variate. If \code{object$mode} is 
\code{canonical} then \eqn{X} and \eqn{Y} are represented on the axis defined by \eqn{U_i} and \eqn{V_i} 
respectively.

Variable pairs with a high similarity measure (in absolute value) are considered as relevant. 
By changing the threshold, one can tune the relevance of the associations to include or 
exclude relationships in the network.
   
\code{interactive=TRUE} open two device, one for association network, one for scrollbar,
and define an interactive process: by clicking either at each end 
(`\eqn{-}' or `\eqn{+}') of the scrollbar or 
at middle portion of this. The position of the slider indicate which is the `threshold' value
associated to the display network. 

The interactive process is terminated by clicking the second button and selecting `Stop' 
from the menu, or from the `Stop' menu on the graphics window. 

The \code{color.node} is a vector of length two, 
of any of the three kind of \code{R} colors, i.e., either a color name 
(an element of \code{colors()}), a hexadecimal string of the form \code{"#rrggbb"}, 
or an integer \code{i} meaning \code{palette()[i]}. \code{color.node[1]} and 
\code{color.node[2]} give the color for filled nodes of the \eqn{X}- and \eqn{Y}-variables
respectively. Defaults to \code{c("white", "white")}.
 
\code{color.edge} give the color to edges with colors corresponding to the values in 
\code{mat}. Defaults to \code{c("blue", "red")} for low and high weight respectively.

\code{shape.node[1]} and \code{shape.node[2]} provide the shape of the nodes associate
to \eqn{X}- and \eqn{Y}-variables respectively. Current acceptable values are 
\code{"circle"} and \code{"rectangle"}. Defaults to \code{c("circle", "rectangle")}.

\code{lty.edge[1]} and \code{lty.egde[2]} give the line type to edges with positive 
and negative weight respectively. Can be one of \code{"solid"}, \code{"dashed"}, 
\code{"dotted"}, \code{"dotdash"}, \code{"longdash"} and \code{"twodash"}. Defaults 
to \code{c("solid", "solid")}.

\code{lwd.edge[1]} and \code{lwd.edge[2]} provide the line width to edges with positive 
and negative weight respectively. This attribute is of type double with 
a default of \code{c(1, 1)}. 
}

\section{Warning}{
If the number of variables is high, the generation of the 
network can take some seconds.
}
 
\value{
\code{network} return a list containing the following components:
  \item{simMat}{the similarity (adjacent) matrix used by \code{network}.}
  \item{gR}{a \code{graph} object (see the \pkg{igraph R} package).}
}

\references{
Butte, A. J., Tamayo, P., Slonim, D., Golub, T. R. and Kohane, I. S. (2000). 
Discovering functional relationships between RNA expression and chemotherapeutic 
susceptibility using relevance networks. \emph{Proceedings of the National Academy 
of Sciences of the USA} \bold{97}, 12182-12186.

Moriyama, M., Hoshida, Y., Otsuka, M., Nishimura, S., Kato, N., Goto, T., Taniguchi, H., 
Shiratori, Y., Seki, N. and Omata, M. (2003). Relevance Network between Chemosensitivity 
and Transcriptome in Human Hepatoma Cells. \emph{Molecular Cancer Therapeutics} \bold{2}, 
199-205.
}

\author{Sébastien Déjean, Ignacio González and Kim-Anh Lê Cao.}

\seealso{\code{\link{plotVar}}, \code{\link{plot3dVar}}, \code{\link{cim}}.}

\examples{
require(igraph)

## network representation for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

\dontrun{   # may not work on the Linux version, use Windows instead
network(nutri.res, comp = 1:3, threshold = 0.6)
}

## Changing the attributes of the network
\dontrun{
network(nutri.res, comp = 1:3, threshold = 0.45,
        color.node = c("mistyrose", "lightcyan"),,
        shape.node = c("circle", "rectangle"), 
        color.edge = jet.colors(8),
        lty.edge = c("solid", "solid"), lwd.edge = c(2, 2), 
        show.edge.labels = FALSE)
}

## interactive 'threshold' 
\dontrun{
network(nutri.res, comp = 1:3, threshold = 0.55, interactive = TRUE)
## select the 'threshold' and "see" the new network
}

## network representation for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
\dontrun{
network(toxicity.spls, comp = 1:3, threshold = 0.8, 
        X.names = NULL, Y.names = NULL, keep.var = TRUE,
        color.node = c("mistyrose", "lightcyan"),
        shape.node = c("rectangle", "circle"),
        color.edge = c("red", "blue"),
        lty.edge = c("solid", "solid"), lwd.edge = c(1, 1), 
        show.edge.labels = FALSE, interactive = FALSE)
}
}

\keyword{multivariate}
\keyword{graphs}
\keyword{dplot}
\keyword{hplot}
\keyword{iplot}

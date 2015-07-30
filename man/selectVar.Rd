\name{selectVar}
\encoding{latin1}
\alias{selectVar}
\alias{selectVar.pls}
\alias{selectVar.plsda}
\alias{selectVar.pca}
\alias{selectVar.spls}
\alias{selectVar.splsda}
\alias{selectVar.spca}
\alias{selectVar.sipca}
\alias{selectVar.sgcca}
\alias{select.var}


\title{Output of selected variables}

\description{
This function outputs the selected variables on each component for the sparse versions of the approaches (can also be generalised to the non sparse versions). 
Note the change of name in version >= 5.0-4 to selectVar.
}

\usage{
\method{selectVar}{pls}(object, comp = 1, \ldots)
\method{selectVar}{plsda}(object, comp = 1, \ldots)
\method{selectVar}{pca}(object, comp = 1, \ldots)
\method{selectVar}{spls}(object, comp = 1, \ldots)
\method{selectVar}{splsda}(object, comp = 1, \ldots)
\method{selectVar}{spca}(object, comp = 1, \ldots)
\method{selectVar}{sipca}(object, comp = 1, \ldots)
\method{selectVar}{sgcca}(object, block = NULL, comp = 1,  \ldots)
}

\arguments{
  \item{object}{object of class inheriting from \code{"pls"}, \code{"spls"}, \code{"plsda"},\code{"splsda"}, \code{"pca"}, \code{"spca"}, \code{"sipca"}.}
  \item{comp}{integer value indicating the component of interest.}
  \item{block}{for an object of class \code{"sgcca"}, the block data sets can be specified as an input vector, for example \code{c(1,2)} for the first two blocks. Default to NULL (all block data sets)}
  \item{\ldots}{other arguments.}
}

\details{
\code{selectVar} provides the variables selected on a given component. \\
  \describe{
    \item{\code{name}}{outputs the name of the selected variables (provided that the input data have colnames) ranked in decreasing order of importance.}
    \item{\code{value}}{outputs the loading value for each selected variable, the loadings are ranked according to their absolute value.}
}
These functions are only implemented for the sparse versions.
}



\author{Kim-Anh Le Cao, Florian Rohart.}

\examples{
data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

# example with sPCA
# ------------------
liver.spca <- spca(X, ncomp = 1, keepX = 10)
selectVar(liver.spca, comp = 1)$name
selectVar(liver.spca, comp = 1)$value

#example with sIPCA
# -----------------
\dontrun{
liver.sipca <- sipca(X, ncomp = 3, keepX = rep(10, 3))
selectVar(liver.sipca, comp = 1)
}

# example with sPLS
# -----------------
\dontrun{
liver.spls = spls(X, Y, ncomp = 2, keepX = c(20, 40),keepY = c(5, 5))
selectVar(liver.spls, comp = 2)

# example with sPLS-DA
data(srbct)   # an example with no gene name in the data
X = srbct$gene
Y = srbct$class

srbct.splsda = splsda(X, Y, ncomp = 2, keepX = c(5, 10))
select = selectVar(srbct.splsda, comp = 2)
select
# this is a very specific case where a data set has no rownames. 
srbct$gene.name[substr(select$select, 2,5),]  
}

# example with sGCCA
# -----------------
\dontrun{
data(nutrimouse)

# ! need to unmap the Y factor
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
# in this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = T)
#note: the penalty parameters need to be tuned
wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.3, 1), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)

#variables selected on component 1 for the two blocs
selectVar(wrap.result.sgcca, comp = 1, block = c(1,2))$name.var

#variables selected on component 2 for each block 
selectVar(wrap.result.sgcca, comp = 2, block = c(1,2))$name.var

# loading value of the variables selected on the first block
selectVar(wrap.result.sgcca, comp = 1, block = 1)$value.var
}

}


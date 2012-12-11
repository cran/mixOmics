\name{select.var}
\encoding{latin1}
\alias{select.var}
\alias{select.var.spls}
\alias{select.var.splsda}
\alias{select.var.spca}
\alias{select.var.sipca}
\alias{select.var.sgcca}



\title{Output of selected variables}

\description{
This function outputs the selected variables on each component for the sparse versions of the approaches. 
}

\usage{
\method{select.var}{spls}(object, comp = 1, \ldots)
\method{select.var}{splsda}(object, comp = 1, \ldots)
\method{select.var}{spca}(object, comp = 1, \ldots)
\method{select.var}{sipca}(object, comp = 1, \ldots)
\method{select.var}{sgcca}(object, block = NULL, comp = 1,  \ldots)

	  
 
 }

\arguments{
  \item{object}{object of class inheriting from \code{"spls"}, \code{"splsda"}, \code{"spca"}, 
    \code{"sipca"}.}
  \item{comp}{integer value indicating the component of interest.}
  \item{block}{for an object of class \code{"sgcca"}, the block data sets can be specified as an input vector, for example \code{c(1,2)} for the first two blocks. Default to NULL (all block data sets)}
  \item{\ldots}{other arguments.}
}

\details{
\code{select.var} provides the variables selected on a given component. \\
  \describe{
    \item{\code{name}}{outputs the name of the selected variables (provided that the input data have colnames) ranked in decreasing order of importance.}
    \item{\code{value}}{outputs the loading value for each selected variable, the loadings are ranked according to their absolute value.}
}
These functions are only implemented for the sparse versions.
}



\author{Kim-Anh L\^e Cao.}

\examples{
data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

# example with sPCA
# ------------------
liver.spca <- spca(X, ncomp = 3, keepX = rep(10, 3))
select.var(liver.spca, comp = 1)$name
select.var(liver.spca, comp = 1)$value

#example with sIPCA
# -----------------
\dontrun{
liver.sipca <- sipca(X, ncomp = 3, keepX = rep(10, 3))
select.var(liver.sipca, comp = 1)
}

# example with sPLS
# -----------------
\dontrun{
liver.spls = spls(X, Y, ncomp = 2, keepX = c(20, 40),keepY = c(5, 5))
select.var(liver.spls, comp = 2)

# example with sPLS-DA
data(srbct)   # an example with no gene name in the data
X = srbct$gene
Y = srbct$class

srbct.splsda = splsda(X, Y, ncomp = 2, keepX = c(5, 10))
select = select.var(srbct.splsda, comp = 2)
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
wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)

select = select.var(wrap.result.sgcca, comp = 1)
# loading value of the variables selected on the first block
select$value.var[[1]]
}

}


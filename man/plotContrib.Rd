\name{plotContrib}
\encoding{latin1}
\alias{plotContrib}

\title{Contribution plot of variables}

\description{
This function provides a horizontal bar plot to visualise the highest or lowest mean/median value of the variables with color code corresponding to the outcome of interest. Only applies to (sparse)PLSDA (including multilevel).
}

\usage{
plotContrib(
                       object,
                       contrib=c("max"),  # choose between 'max" or "min"
                       method = c("mean", "median"), 
                       comp = 1,
                       ties = TRUE, 
                       ndisplay = NULL, 
                       cex.name = 0.7, 
                       cex.legend = 0.8,
                       name.var = NULL,
                       legend = TRUE,
                       legend.color = NULL, 
                       main = NULL,
                       legend.title = 'Outcome',
                       plot = TRUE 
                       
                       )
	  
}

\arguments{
  \item{object}{object of class inheriting from a \code{plsda} or \code{splsda} class}, 
  \item{contrib}{a character set to 'max' or 'min' indicating if the color of the bar should correspond to the group with the maximal or minimam expression levels / abundance. By default set to 'max'.}
  \item{method}{a character set to 'mean' or 'median' indicating the criterion to assess the contribution. We recommend using median in the case of count or skewed data.}
  \item{comp}{integer value indicating the component of interest from the object.}
  \item{ties}{boolean. If TRUE then tie groups appear in white color which will appear in the legend. Ties can happen when dealing with count data type. By default set to TRUE.}
  \item{ndisplay}{integer. }
  \item{cex.name}{A numerical value giving the amount by which plotting the variable name text should be magnified or reduced relative to the default.}
  \item{cex.legend}{A numerical value giving the amount by which plotting the legend text should be magnified or reduced relative to the default.}
  \item{name.var}{A character vector indicating the names of the variables. The names of the vector should match the names of the input data, see example.}
  \item{legend}{boolean indicating if the legend indicating the group outcomes should be added to the plot. Default value is TRUE.}
  \item{legend.color}{A color vector of length the number of group outcomes. See examples.}
  \item{main}{A set of characters to indicate the title of the plot. Default value is NULL.}
  \item{legend.title}{A set of characters to indicate the title of the legend. Default value is NULL.}
  \item{plot}{boolean indicating of the plot should be output. If set to FALSE the user can extract the contribution matrix, see example. Default value is TRUE.}
}

\details{

The contribution of each variable is represented in a barplot where eachbar length corresponds to the loading weight (importance) of the feature in (sparse)PLSDA for each component. The loading weight can be  positive or negative. The color corresponds to the group in which the feature is most 'abundant'.  Note that this type of graphical output is particularly insightful for count microbial data - in that latter case using the \code{method = 'median'} is advised.
}

\author{Kim-Anh Le Cao, Benoit Gautier, Florian Rohart}

\seealso{\code{\link{plotVar}}, \code{\link{plot3dVar}}, \code{\link{cim}}, for other variable plots and http://mixOmics.org/graphics for more details.}


\examples{
## object of class 'splsda' 
# --------------------------
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))

# contribution on comp 1, based on the median. 
# Colors indicate the group in which the median expression is maximal
plotContrib(splsda.liver, comp = 1, method = 'median')

# contribution on comp 2, based on median. 
#Colors indicate the group in which the median expression is maximal
plotContrib(splsda.liver, comp = 2, method = 'median')

# contribution on comp 2, based on median. 
# Colors indicate the group in which the median expression is minimal
plotContrib(splsda.liver, comp = 2, method = 'median', contrib = 'min')

# changing the name to gene names
# if the user input a name.var but names(name.var) is NULL, 
# then a warning will be output and assign names of name.var to colnames(X)
# this is to make sure we can match the name of the selected variables to the contribution plot.
name.var = liver.toxicity$gene.ID[, 'geneBank']
length(name.var)
plotContrib(splsda.liver, comp = 2, method = 'median', name.var = name.var)

# if names are provided: ok, even when NAs
name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)
plotContrib(splsda.liver, comp = 2, method = 'median', 
name.var = name.var, cex.name = 0.5)


# look at the contribution (median) for each variable
plot.contrib = plotContrib(splsda.liver, comp = 2, method = 'median', plot = FALSE)
head(plot.contrib$contrib)
# change the title of the legend and title name
plotContrib(splsda.liver, comp = 2, method = 'median', legend.title = 'Time', 
main = 'Contribution plot')

# no legend
plotContrib(splsda.liver, comp = 2, method = 'median', legend = FALSE)

# change the color of the legend
plotContrib(splsda.liver, comp = 2, method = 'median', legend.color = c(1:4))



# object 'splsda multilevel'
# -----------------
\dontrun{
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)

# multilevel sPLS-DA model
res.1level <- multilevel(X, ncomp = 3, design = design,
                         method = "splsda", keepX = c(30, 137, 123))


name.var = vac18$tab.prob.gene[, 'Gene']
names(name.var) = colnames(X)

plotContrib(res.1level, comp = 2, method = 'median', legend.title = 'Stimu', 
name.var = name.var, cex.name = 0.2)

# too many transcripts? only output the top ones
plotContrib(res.1level, comp = 2, method = 'median', legend.title = 'Stimu', 
name.var = name.var, cex.name = 0.5, ndisplay = 60)
}

# object 'plsda'
# ----------------
\dontrun{
# breast tumors
# ---
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)

name.var = as.character(breast.tumors$genes$name)
names(name.var) = colnames(X)

# with gene IDs, showing the top 60
plotContrib(plsda.breast, contrib = 'max', comp = 1, method = 'median', 
            ndisplay = 60, 
            name.var = name.var, cex.name = 0.6,
            legend.color = color.mixo(1:2))
}

# liver toxicity
# ---
\dontrun{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$treatment[, 4]

plsda.liver <- plsda(X, Y, ncomp = 2)
plotIndiv(plsda.liver, ind.names = Y, plot.ellipse = TRUE)


name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)

plotContrib(plsda.liver, contrib = 'max', comp = 1, method = 'median', ndisplay = 100, 
            name.var = name.var, cex.name = 0.4,
            legend.color = color.mixo(1:4))
}

}

\keyword{multivariate}


\name{pheatmap.multilevel}
\encoding{latin1}
\alias{pheatmap.multilevel}
\alias{pheatmap.multilevel.splsda1fact}
\alias{pheatmap.multilevel.splsda2fact}

\title{Clustered heatmap}

\description{
A function to draw clustered heatmaps from the pheatmap package for the sPLS-DA multilevel functions for one and two factors only.
}

\usage{
\method{pheatmap.multilevel}{splsda1fact}(
result, 
cluster = NULL, 
color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", 
"#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100), 
col_sample = NULL, col_stimulation = NULL, 
label_annotation = NULL, 
breaks = NA, border_color = "grey60", 
cellwidth = NA, 
cellheight = NA, 
scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, 
clustering_distance_rows = "euclidean", 
clustering_distance_cols = "euclidean", 
clustering_method = "complete", treeheight_row = ifelse(cluster_rows, 50, 0), 
treeheight_col = ifelse(cluster_cols, 50, 0), 
legend = TRUE, annotation = NA, annotation_colors = NA, annotation_legend = TRUE, 
show_rownames = TRUE, show_colnames = TRUE, fontsize = 10, 
fontsize_row = fontsize, fontsize_col = fontsize, filename = NA, 
width = NA, height = NA, order_sample = NULL, tab.prob.gene = NULL, ...)

\method{pheatmap.multilevel}{splsda2fact}(
result, cluster = NULL, 
color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", 
"#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100), 
col_sample = NULL, col_stimulation = NULL, col_time = NULL, 
label_color_stimulation = NULL, label_color_time = NULL, 
label_annotation = NULL, breaks = NA, border_color = "grey60", 
cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE, 
cluster_cols = TRUE, clustering_distance_rows = "euclidean", 
clustering_distance_cols = "euclidean", clustering_method = "complete", 
treeheight_row = ifelse(cluster_rows, 50, 0), 
treeheight_col = ifelse(cluster_cols, 50, 0),
legend = TRUE, annotation = NA, annotation_colors = NA, 
annotation_legend = TRUE, show_rownames = TRUE, show_colnames = TRUE, 
fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, 
filename = NA, width = NA, height = NA, order_sample = NULL, 
tab.prob.gene = NULL, ...)
}	

\arguments{
  \item{result}{a result from function \code{multilevel} with \code{method = 'splsda'} (can only apply to a splsda analysis type).}
  \item{cluster}{a vector that would indicate a specific ordering (clustering) on the variables. By default should be set to \code{NULL}}
  \item{col_sample}{vector of colors indicating the color of each individual (the length of that vector should be equal to the number of unique individuals)}
  \item{col_stimulation}{vector of colors indicating the color of each condition (the length of that vector should be equal to the number of conditions), i.e. the condition indicated in the 2nd colum of the \code{design} matrix in the multilevel function.}
  \item{col_time}{for a two-factor analysis, vector of colors indicating the color of the factors of the second level (the length of that vector should be equal to the number of conditions), i.e. the condition indicated in the 3rd colum of the \code{design} matrix in the multilevel function.}
  \item{label_color_stimulation}{character vector indicating the label of the first factor, see details}
  \item{label_color_time}{character vector indicating the label of the second factor, see details}
  \item{label_annotation}{set to \code{NULL} by default.}
  \item{order_sample}{vector indicatin the reordering of the samples, set to \code{NULL} by default. }
  \item{color}{vector of colors used in heatmap.}
  \item{breaks}{a sequence of numbers that covers the range of values in mat and is one element longer than color vector. Used for mapping values to colors. Useful, if needed to map certain values to certain colors, to certain values. If value is NA then the breaks are calculated automatically.}
  \item{border_color}{color of cell borders on heatmap, use NA if no border should be 
drawn.}
  \item{cellwidth}{individual cell width in points. If left as NA, then the values depend on the size of plotting window.}
  \item{cellheight}{individual cell height in points. If left as NA, then the values depend on the size of plotting window.}
\item{scale}{character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are \code{"row"}, \code{"column"} and \code{"none"}}
\item{cluster_rows}{boolean values determining if rows should be clustered,}
\item{cluster_cols}{boolean values determining if columns should be clustered.}
\item{clustering_distance_rows}{distance measure used in clustering rows. Possible values are \code{"correlation"} for Pearson correlation and all the distances supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none of the above it is assumed that a distance matrix is provided.}
\item{clustering_distance_cols}{distance measure used in clustering columns. Possible values the same as for clustering_distance_rows.}
\item{clustering_method}{clustering method used. Accepts the same values as \code{\link{hclust}}.}
\item{treeheight_row}{the height of a tree for rows, if these are clustered. 
Default value 50 points.}
\item{treeheight_col}{the height of a tree for columns, if these are clustered. 
Default value 50 points.}
\item{legend}{boolean value that determines if legend should be drawn or not.}
\item{annotation}{data frame that specifies the annotations shown on top of the columns. Each row defines the features for a specific column. The columns in the data and rows in the annotation are matched using corresponding row and column names. Note that color schemes takes into account if variable is continuous or discrete.}
\item{annotation_colors}{list for specifying annotation track colors manually. It is possible to define the colors for only some of the features. Check examples for details.}
\item{annotation_legend}{boolean value showing if the legend for annotation tracks 
should be drawn.}
\item{show_rownames}{boolean specifying if column names are be shown.}
\item{show_colnames}{boolean specifying if column names are be shown.}
\item{fontsize}{base fontsize for the plot}
\item{fontsize_row}{fontsize for rownames (Default: fontsize)}
\item{fontsize_col}{fontsize for colnames (Default: fontsize)}
\item{filename}{file path where to save the picture. Filetype is decided by 
the extension in the path. Currently following formats are supported: png, pdf, tiff,
bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is 
calculated so that the plot would fit there, unless specified otherwise.}
\item{width}{manual option for determining the output file width in}
\item{height}{manual option for determining the output file height in inches.}
\item{tab.prob.gene}{A character vector indicating the name of all genes.}
\item{\dots}{graphical parameters for the text used in plot. Parameters passed to \code{\link{grid.text}}, see \code{\link{gpar}}.}
}

\details{
This function has been borrowed from the \code{pheatmap} function of the pheatmap package. See \code{help(pheatmap)} for more details about the arguments of the function.

In the \code{multilevel} function, the factors indicated in the vector or the data frame \code{cond} must match the arguments \code{label_color_stimulation} and \code{label_color_time}, see example below.

}


\references{
On multilevel analysis:
Liquet, B., Le Cao, K.-A., Hocini, H. and Thiebaut, R. (2012) A novel approach for biomarker selection and the
integration of repeated measures experiments from two platforms. \emph{BMC Bioinformatics} \bold{13}:325.

Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K. (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA. \emph{Metabolomics}, \bold{6}(1), 119-128.
}

\author{Raivo Kolde <rkolde@gmail.com> (pheatmap), Benoit Liquet, Kim-Anh Le Cao, Benoit Gautier (for mixOmics).}

\seealso{\code{\link{multilevel}}, \code{\link{pheatmap}} and http:www.mixOmics.org for more details..}

\examples{
## First example: one-factor analysis with sPLS-DA
# -------------------
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation

design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)
vac18.splsda.multilevel <- multilevel(X, ncomp = 3, design = design,
                                         method = "splsda", keepX = c(30, 137, 123))


# set up colors for pheatmap
col.samp <- c("lightgreen", "red", "lightblue", "darkorange",
              "purple", "maroon", "blue", "chocolate", "turquoise",
              "tomato1", "pink2", "aquamarine")
col.stimu <- c("darkblue", "purple", "green4","red3")
col.stimu <- col.stimu[as.numeric(Y)]
col.stimu <- unique(col.stimu)

pheatmap.multilevel(vac18.splsda.multilevel, 
                    # colors:
                    col_sample = col.samp, 
                    col_stimulation = col.stimu,
                    #labels:
                    label_annotation = c("Subject", "Stimulus"),
                    # scaling:
                    scale = 'row',
                    # distances and clutering
                    clustering_distance_rows = "euclidean", 
                    clustering_distance_cols = "euclidean", 
                    clustering_method = "complete",
                    #  show col/row names and font
                    show_colnames = FALSE,
                    show_rownames = FALSE, 
                    fontsize = 8, 
                    fontsize_row = 3,
                    fontsize_col = 2,
                    border = FALSE, 
                    width = 10)

## Second example: two-factor analysis with sPLS-DA
# --------------------
\dontrun{
data(vac18.simulated) # on the simulated data this time

X <- vac18.simulated$genes
design <- data.frame(sample = vac18.simulated$sample,
                     stimul = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

vac18.splsda2.multilevel <- multilevel(X, ncomp = 2, design = design,
                            keepX = c(200, 200), method = 'splsda')

# set up colors for each level of pheatmap 
col.sample <- c("lightgreen", "red","lightblue","darkorange","purple","maroon") # 6 samples
col.time <- c("pink","lightblue1") # two time points
col.stimu <- c('green', 'black', 'red', 'blue') # 4 stimulations
# set up labels for the 2 levels in design matrix
label.stimu <- unique(design[, 2])
label.time <- unique(design$time)

pheatmap.multilevel(vac18.splsda2.multilevel,
                                # colors:
                                col_sample=col.sample, 
                                col_stimulation=col.stimu, 
                                col_time=col.time,
                                #labels for each level
                                label_color_stimulation=label.stimu,
                                label_color_time=label.time, 
                                #clustering method
                                clustering_method="ward",
                                #show col/row names and font size
                                show_colnames = FALSE,
                                show_rownames = TRUE,
                                fontsize_row=2)
}
}
\keyword{regression}
\keyword{multivariate}

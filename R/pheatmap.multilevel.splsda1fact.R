# Copyright (C) 2015 
# Benoit Liquet, Universite de Bordeaux, France
# Kim-Anh Le Cao, University of Queensland, Brisbane, Australia
# Benoit Gautier, University of Queensland, Brisbane, Australia
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

pheatmap.multilevel.splsda1fact <- function(result, cluster = NULL, color = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100), 
                                            col_sample = NULL, col_stimulation = NULL, label_annotation = NULL, 
                                            breaks = NA, border_color = "grey60", cellwidth = NA, cellheight = NA, 
                                            scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, 
                                            clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
                                            clustering_method = "complete", treeheight_row = ifelse(cluster_rows, 50, 0), treeheight_col = ifelse(cluster_cols, 50, 0), 
                                            legend = TRUE, annotation = NA, annotation_colors = NA, annotation_legend = TRUE, 
                                            show_rownames = TRUE, show_colnames = TRUE, fontsize = 10, 
                                            fontsize_row = fontsize, fontsize_col = fontsize, filename = NA, 
                                            width = NA, height = NA, order_sample = NULL, tab.prob.gene = NULL, ...) { 
  
  if (result$ncomp == 1) {
    name.probe <- names(result$loadings$X[unique(which(result$loadings$X[, 1] != 0)), 1:result$ncomp])
  } else {
    name.probe <- rownames(result$loadings$X[unique(which(result$loadings$X != 0, arr.ind = TRUE)[, 1]), 1:result$ncomp])
  }
  
  ## note: not sure what the cluster argument is for since it is set to NULL
  if (is.null(order_sample)) 
    order_sample <- 1:dim(result$Xw)[1]
  mat <- result$Xw[order_sample, setdiff(name.probe, cluster)]
  rownames(mat) <- order_sample
  geneX <- colnames(mat)
  
  if (!(is.null(result$tab.prob.gene))) 
    geneX <- result$tab.prob.gene[match(geneX, result$tab.prob.gene[, 1]), 2] 
  
  matt <- t(mat)
  rownames(matt) <- geneX
  
  # extract information contained in the design matrix
  sample = as.character(result$design[, 1])
  cond = factor(result$design[, 2])

  name.sample <- unique(sample)
  nsujet <- length(name.sample)
  
  annotation <- data.frame(Sample = sample, Stimulation = cond)
  rownames(annotation) <- 1:(dim(annotation)[1])
  
  if (is.null(col_sample)) 
    col_sample <- colors()[sample(1:400, nsujet)]
  
  Sample <- col_sample[1:nsujet]
  names(Sample) <- c("1", 2:nsujet)
  
  if (is.null(col_stimulation)) 
    color_stimulation <- colors()[sample(1:400, nlevels(cond))]
  
  Stimulation <- col_stimulation
  names(Stimulation) <- levels(cond)
  annotation_colors <- list(Sample = Sample, Stimulation = Stimulation)
  
  if (!(is.null(label_annotation))) 
    names(annotation_colors) <- names(annotation) <- label_annotation
  
  # note: message supressed for Ward method 
  #'The "ward" method has been renamed to "ward.D"; note new "ward.D2"'
  #'# as maintainers from the phetmap package have not done the update
  suppressMessages(
  pheatmap(matt, color = color, breaks = breaks, border_color = border_color, 
           cellwidth = cellwidth, cellheight = cellheight, scale = scale, 
           cluster_rows = cluster_rows, cluster_cols = cluster_cols, 
           clustering_distance_rows = clustering_distance_rows, 
           clustering_distance_cols = clustering_distance_cols, 
           clustering_method = clustering_method, treeheight_row = treeheight_row, 
           treeheight_col = treeheight_col, legend = legend, annotation = annotation, 
           annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
           show_rownames = show_rownames, show_colnames = show_colnames, 
           fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, 
           filename = filename, width = width, height = height, 
           ...)
  )
}
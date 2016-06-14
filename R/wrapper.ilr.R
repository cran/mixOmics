#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2015
# last modified: 25-02-2016
#
# Copyright (C) 2015
#
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
#############################################################################################################



# logratio.transfo
logratio.transfo = function(X,
logratio = "none" # one of ('none','CLR','ILR')
)
{
    
    if (logratio == 'ILR')
    {
        if (any(class(X) != 'ilr'))
        {   # data are ilr transformed, then the data loose 1 variable, but we'll use V to reconstruct the matrix
            X = ilr.transfo(X)
        }
    }else if (logratio == 'CLR') {
        X = clr.transfo(X)
    }
    #if logratio = "none", do nothing
    
    return(X)
}


# 1 - ilr transform of the data, isoLMR function from robCompositions package, with changes
# -----------------

# KA changed the function to add a min value when many zeroes in data (prob with log and division by 0 otherwise)
ilr.transfo = function(x, fast = TRUE, min.value = min(x[which(x != 0)])*0.01)
{
    
    # ilr transformation
    x.ilr = matrix(NA, nrow = nrow(x), ncol = ncol(x)-1)
    D = ncol(x)
    # KA added: a little something to avoid 0 values
    if (fast)
    {
        for (i in 1 : ncol(x.ilr))
        {
            x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(((apply(as.matrix(x[, (i+1) : D, drop = FALSE]), 1, prod) + min.value)^(1 / (D-i))) / (x[,i]+ min.value))
            #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop = FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
        }
    } else {
        for (i in 1 : ncol(x.ilr))
        {
            x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(apply(as.matrix(x[, (i+1):D]), 1, function(x){exp(log(x))})/(x[, i]+ min.value) + min.value)
            #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(apply(as.matrix(x[,(i+1):D]), 1, function(x){exp(log(x))})/(x[,i]))
        }
    }
    class(x.ilr) = 'ilr'
    return(as.matrix(x.ilr))
}




# 2 - back transformation from ilr to clr space
# -------------------
clr.backtransfo = function(x)
{
    # construct orthonormal basis
    V = matrix(0, nrow = ncol(x), ncol = ncol(x)-1)
    for( i in 1:ncol(V) )
    {
        V[1:i, i] = 1/i
        V[i+1, i] = (-1)
        V[, i] = V[, i] * sqrt(i/(i+1))
    }
    rownames(V) = colnames(x)
    return(V)
    
}


# CLR transformation
clr.transfo = function(x)
{
    # KA added
    min.value = min(x[which(x != 0)])*0.01
    
    
    #if (dim(x)[2] < 2) stop("data must be of dimension greater equal 2")
    if (dim(x)[2] == 1)
    {
        res = list(x.clr = x, gm = rep(1, dim(x)[1]))
    } else{
        geometricmean = function (x) {
            #       if (any(na.omit(x == 0)))
            #         0
            #       else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
            #     }
            # KA changed to
            exp(mean(log(x + min.value)))
        }
        gm = apply(x, 1, geometricmean)
        # KA changed
        x.clr = log((x + min.value) / (gm))
        res = x.clr #list(x.clr = x.clr, gm = gm)
    }
    class(res) = "clr"
    return(res)
}


# 2 - Filter samples with counts that are too low and denoising
# -----------------

# KA changed the function to add a nan value when many zeroes in data (prob with log and division by 0 otherwise)
#data: data.raw (count data.frame )
#indiv is a data.frame with information on the sample
## sample is the colname use as ID for example RSID
mixMC.filter = function(data, indiv , sample , affiliation = NULL, taxonomy = NULL)
{
    
    if (is.matrix(data))
    data = as.data.frame(data)
    
    if (!is.data.frame(data))
    stop("'data' should be a data.frame.", call. = FALSE)
    
    if (is.matrix(indiv) || is.vector(indiv))
    indiv = as.data.frame(indiv)
    
    if (!is.data.frame(indiv))
    stop("'indiv' should be a data.frame.", call. = FALSE)
    
    #   if ((!is.vector(sample)) || length(sample)>1)
    #     stop("'sample' should be a character vector of length = 1.", call. = FALSE)
    #
    if (!sample %in% colnames(indiv) || (!is.vector(sample)) || length(sample)>1)
    stop("'sample' should be in:",paste0(as.vector(colnames(indiv)),collapse = ';'), call. = FALSE)
    
    # to determine a sensible cutoff
    remove.patient = which(apply(data, 1, sum) < 10)
    #match to sample ID
    remove.RSID = indiv[remove.patient,sample ]
    # remove all matching (unique) patients to keep a repeated measures balanced design
    remove.final = which(indiv[,sample] %in%  remove.RSID)
    indiv.diverse = indiv[-c(remove.final), ]
    data.diverse.unfiltered = data[-c(remove.final), ]
    
    # KA added alternative denoising method using as proposed by:
    # http://enterotype.embl.de/enterotypes.html (Between Class Analysis)
    # KA changed: colSums as the data are transposed in our case (OTU in columns)
    # KA amended output
    
    keep.otu = which(colSums(data.diverse.unfiltered)*100/(sum(colSums(data.diverse.unfiltered))) > 0.01)
    data.filter.diverse = data.diverse.unfiltered[, keep.otu]
    
    result = list()
    
    result$keep.otu = keep.otu
    result$data.filter.diverse = data.filter.diverse
    result$indiv.diverse = indiv.diverse
    # update affiliation for diverse; check OK
    if (!is.null(affiliation)){
        affiliation.diverse = affiliation[keep.otu]
        result$affiliation.diverse = affiliation.diverse
    }
    # storing taxo info as matrix; check OK
    if (!is.null(taxonomy)){
        taxonomy.diverse = taxonomy[keep.otu, ]
        result$taxonomy.diverse = taxonomy.diverse
    }
    # -----------------------
    # return data
    # ----------------------
    return(invisible(result))
    
}

# 3 - Using TSS (Total-Sum-Scaling) on the filtered raw counts for the diverse OTU set.
# -----------------

# KA changed the function to add a nan value when many zeroes in data (prob with log and division by 0 otherwise)
#data: data.filter (filter count data.frame )
#log_normalisation is the choice between log TSS normalisation or not
normalisation.TSS = function(data, log_normalisation = FALSE)
{
    
    if (is.matrix(data))
    data = as.data.frame(data)
    
    if (!is.data.frame(data))
    stop("'data' should be a data.frame.", call. = FALSE)
    
    if (!is.logical(log_normalisation)) {
        stop("'log_normalisation' should be either a logical value.",
        call. = FALSE)
    }
    
    
    TSS.divide = function(x){
        x/sum(x)
    }
    
    if (!isTRUE(log_normalisation))
    {
        #     TSS normalises count data by dividing each OTU read count by the total number of reads in each sample.
        #     It converts counts to scaled ratios.
        #     The first 10 values of the function are checked so that they return same values as the output data matrix.
        result = t(apply(data, 1, TSS.divide))
        
    }
    else
    {
        # Apply the same function to the data set dividing the log of the value by the sum of all log values.
        result = t(apply(log(data+1), 1, TSS.divide))
    }
    
    
    return(invisible(result))
    
}

# # Function to create the Annotation file for a graphiAn input
#
# taxo.output = function(contrib, Y, name.var, taxonomy){
#   # contrib is a result from plotCOntrib$contrib
#   # Y indicates the sample class
#   # name.var is the OTU name matching with contrib
#   # taxonomy is the input taxonomy file (so far a matrix for the different classification levels, we only consider the family level here)
#
#   # taxonomy info
#   taxonomy = as.data.frame(taxonomy)
#   taxo.input = cbind(taxonomy[rownames(contrib),], OTU = sub('.', '_', rownames(contrib),fixed = TRUE))
#   # indicate where NA values are, or blanks
#   na.index = c(which(is.na(taxo.input$Family)), which(as.character(taxo.input$Family) ==  ''))
#
#   # taxonomy at the OTU levels
#   taxo.level = name.var[rownames(contrib)]
#   # node size
#   node.size = vector(length = nrow(contrib))
#   for(k in 1:nrow(contrib)){
#     # based on median: extract the max median from a contrib output
#     node.size[k] = abs(as.numeric(max(contrib[k, 1:nlevels(Y)]))) #[contrib[k, 'GroupContrib']])
#     # based on loading
#     #node.size[k] = as.numeric(abs(contrib[k, 'importance']))
#   }
#   ##node.size = round((node.size)*100, 2)
#   node.size = round(node.size*50)
#
#
#   #--node color: green is pos, red is neg
#   node.color = c(rep('#388ECC', nrow(contrib)))
#   node.color[which(contrib$importance < = 0)] = '#F0E442'
#
#   #-- edge color (no need)
#   edge.color = c(rep(NA, nrow(contrib)))
#
#   #--background color
#   back.color = contrib$Contrib
#
#   #--body site
#   GroupContrib = contrib$GroupContrib
#
#   # summarise the results
#   annot.family = cbind(taxo.level, node.size, node.color, edge.color, back.color, GroupContrib)
#
#   # remove na values: when we dont have the family information
#   if (length(na.index) > 0 ) {
#     taxo.input = taxo.input[-c(na.index), ]
#     annot.family = annot.family[-c(na.index), ]
#   }
#   return(list(taxo = taxo.input, annot = annot.family, na.index = na.index))
# }
#
# #  ================================================== 
# # outputs for cladogram and call graphlan
# #  ================================================== 
#
# # first, name.var: remove the OTU_97 text for example
# # name.var2 = sub('OTU_97.', 'OTU_97_', rownames(taxonomy)) #alternatively could be rownames(taxonomy)
# # names(name.var2) = rownames(taxonomy)
#
#
# graphlan.mixMC = function(data , taxonomy , Y, choice.ncomp = 2,path.output = NULL,col.sample.unique = NULL,path.graphlan = NULL){
#
#   # data
#   # Y indicates the sample class
#   # taxonomy is the input taxonomy file (so far a matrix for the diffrent clasification levels, we only consider the family level here)
#   name.var2 = sub('.', '_', rownames(taxonomy),fixed = TRUE) #alternatively could be rownames(taxonomy)
#   names(name.var2) = rownames(taxonomy)
#    # print(name.var2)
#   # calculate contribution
#   mat.taxo = mat.annot = NULL # initialise
#   # then fill in tree input per component
#   for(k in 1: choice.ncomp){
#     contrib = plotContrib(data, method = 'median', comp = k, plot = F)
#     res = taxo.output(contrib = contrib$contrib, Y = Y, name.var = name.var2, taxonomy = taxonomy)
#     mat.taxo = rbind(mat.taxo, res$taxo)
#     mat.annot = rbind(mat.annot, res$annot)
#   }
#   mat.taxo2 = matrix(nrow = dim(mat.taxo)[1],ncol = dim(mat.taxo)[2])
#
#   if (!is.null(col.sample.unique) && length(col.sample.unique) != nlevels(Y))
#     stop("'col.sample.unique' should have a length ",nlevels(Y), call. = FALSE)
#
#   # suppress double or triple points
#   mat.taxo.final = matrix(ncol = 1,nrow = nrow(mat.taxo))
#
#   for(i in 1:nrow(mat.taxo))
#   {
#     mat.taxo.level = NULL
#     for(j in 1:ncol(mat.taxo[i,]))
#     {
#       if ( !is.na(mat.taxo[i,j][[1]][1]))
#       {
#         if ((!as.character(mat.taxo[i,j][[1]][1]) == "") && (!as.character(mat.taxo[i,j][[1]][1]) == 'NA'))
#         {
#
#           if (j == 1)
#           mat.taxo.level = mat.taxo[i,j]
#         else
#           mat.taxo.level = paste(mat.taxo.level,mat.taxo[i,j],sep = '.')
#         mat.taxo2[i,j] = as.character(mat.taxo[i,j])
#         }
#         else
#         {
#           mat.taxo2[i,j] = 'NA'
#         }
#       }
#       else
#       {
#         mat.taxo2[i,j] = 'NA'
#       }
#     }
#     mat.taxo.final[i,] = mat.taxo.level
#   }
#   # input ready for graphLAn
#   write.table(mat.taxo.final, sep = '.', row.names = F, col.names = F, na = '', quote = F, file = paste(path.output, 'treeInput.txt', sep = ''))
#
#   # ------------------
#   # #prepare background annotation
#   # ------------------
#
#   if (!is.null(col.sample.unique))
#   {
#     for (i in 1:nlevels(Y))
#     {
#       mat.annot[mat.annot[,6] == levels(Y)[i],5] = col.sample.unique[i]
#     }
#   }
#
#   annotation.otu = c(as.character(mat.taxo2[,2]),as.character(mat.taxo2[,3]),as.character(mat.taxo2[,5]))
#   annotation.otu = unique(annotation.otu)
#   red.annotation.otu = c()
#   let = c(letters,LETTERS)
#   ncont = 1
#   list.mat.taxo = c(as.character(mat.taxo2[,2]),as.character(mat.taxo2[,3]),as.character(mat.taxo2[,5]))
#   # print()
#   tb = table(factor(list.mat.taxo,levels = annotation.otu))
#
#   for(i in 1:length(as.vector(tb)))
#   {
#     if ((as.vector(tb)[i]*2)<(nchar(names(tb)[i])))
#     {
#       red.annotation.otu = c(red.annotation.otu,paste(let[ncont],names(tb)[i],sep = ':'))
#       ncont = ncont+1
#     }
#     else
#     {
#       red.annotation.otu = c(red.annotation.otu,names(tb)[i])
#     }
#   }
#
#   n.nrow = 9+3*nrow(mat.annot)+3*length(annotation.otu)
#   mat.annot.final = matrix(nrow = n.nrow,ncol = 3)
#   mat.annot.final[,3] = ''
#
#   mat.annot.final[1,1] = 'clade_separation'
#   mat.annot.final[1,2] = 0.5
#   mat.annot.final[2,1] = 'branch_thickness'
#   mat.annot.final[2,2] = 1.5
#   mat.annot.final[3,1] = 'branch_bracket_depth'
#   mat.annot.final[3,2] = 0.8
#   mat.annot.final[4,1] = 'branch_bracket_width'
#   mat.annot.final[4,2] = 0.25
#   mat.annot.final[5,1] = 'clade_marker_size'
#   mat.annot.final[5,2] = 40
#   mat.annot.final[6,1] = 'clade_marker_edge_color'
#   mat.annot.final[6,2] = '#555555'
#   mat.annot.final[7,1] = 'clade_marker_edge_width'
#   mat.annot.final[7,2] = 1.2
#   mat.annot.final[8,1] = 'annotation_background_alpha'
#   mat.annot.final[8,2] = 0.6
#   mat.annot.final[9,1] = 'annotation_legend_font_size'
#   mat.annot.final[9,2] = 9
#
#   mat.annot.final[1:nrow(mat.annot)+9,2] = 'clade_marker_color'
#   mat.annot.final[1:nrow(mat.annot)+9,1] = as.character(mat.annot[,1])
#   mat.annot.final[1:nrow(mat.annot)+9,3] = as.character(mat.annot[,3])
#   mat.annot.final[1:nrow(mat.annot)+(nrow(mat.annot)+9),2] = 'clade_marker_size'
#   mat.annot.final[1:nrow(mat.annot)+(nrow(mat.annot)+9),1] = as.character(mat.annot[,1])
#   mat.annot.final[1:nrow(mat.annot)+(nrow(mat.annot)+9),3] = mat.annot[,2]
#   mat.annot.final[1:nrow(mat.annot)+(nrow(mat.annot)*2+9),2] = 'annotation_background_color'
#   mat.annot.final[1:nrow(mat.annot)+(nrow(mat.annot)*2+9),1] = as.character(mat.annot[,1])
#   mat.annot.final[1:nrow(mat.annot)+(nrow(mat.annot)*2+9),3] = as.character(mat.annot[,5])
#   mat.annot.final[1:length(annotation.otu)+(nrow(mat.annot)*3+9),2] = 'annotation'
#   mat.annot.final[1:length(annotation.otu)+(nrow(mat.annot)*3+9),1] = as.character(annotation.otu)
#   mat.annot.final[1:length(annotation.otu)+(nrow(mat.annot)*3+9),3] = as.character(red.annotation.otu)
#   mat.annot.final[1:length(annotation.otu)+(nrow(mat.annot)*3+length(annotation.otu)+9),2] = 'annotation_background_color'
#   mat.annot.final[1:length(annotation.otu)+(nrow(mat.annot)*3+length(annotation.otu)+9),1] = annotation.otu
#   mat.annot.final[1:length(annotation.otu)+(nrow(mat.annot)*3+length(annotation.otu)+9),3] = '#e2e2e2'
#   mat.annot.final[1:length(annotation.otu)+(nrow(mat.annot)*3+2*length(annotation.otu)+9),2] = 'annotation_font_size'
#   mat.annot.final[1:length(annotation.otu)+(nrow(mat.annot)*3+2*length(annotation.otu)+9),1] = annotation.otu
#   mat.annot.final[1:length(annotation.otu)+(nrow(mat.annot)*3+2*length(annotation.otu)+9),3] = 11
#
#
#   write.table(mat.annot.final, sep = '\t', row.names = F, col.names = F, quote = F, file = paste(path.output, 'annotation.txt', sep = ''))
#
#
#   # --------------
#   # call graphLan with bash code
#   # -------------
#
#
#   # was checking whether those were individual OTU selected on each component
#   #table(mat.annot[, "taxo.level"])
#
#
#   
#   x = paste(path.graphlan,"graphlan_annotate.py --annot ",path.output,"annotation.txt ",path.output,"treeInput.txt ",path.output,"treeInput.xml",sep = '')
#   cat(x)
#   
#   system(x)
#   
#   x = paste(path.graphlan,"graphlan.py ",path.output,"treeInput.xml ",path.output,"Cladogram.png --dpi 300 --size 7 --pad 0.75",sep = '')
#   cat(x)
#   
#   system(x)
#   
#   
# }



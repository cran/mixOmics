#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 24-05-2015
# last modified: 27-05-2016
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

# ========================================================================================================
# internal_predict.DA: prediction module for Discriminant Analysis. Used in 'predict.mint.block.pls.R'
# ========================================================================================================



#out
#object$X
#class(objet)
#dist
#variatesX
#ncomp
#Y.hat
#newdata
#t.pred
#Y


internal_predict.DA = function(object, out, q, dist)
{
    
    if (length(grep("plsda",class(object)))==0) # a DA analysis (mint).(block).(s)plsda
    stop("'Object' is not from a Discriminant Analysis", call.=FALSE)
    
    out.DA = list()
    J = length(object$X) #at this stage we have a list of blocks
    p = lapply(object$X, ncol)
    t.pred = out$variates
    Y.hat = out$predict
    newdata = out$newdata #actually concat.newdata
    variatesX = object$variates[-(J + 1)];
    ncomp = object$ncomp
    
    Y = object$Y
    Y.prim = unmap(object$Y)
    G = cls = list()
    for (i in 1 : J)
    {
        G[[i]] = sapply(1:q, function(x) {apply(as.matrix(variatesX[[i]][Y.prim[, x] == 1,,drop=FALSE]), 2, mean)})
        if (ncomp[i] == 1)
        G[[i]] = t(t(G[[i]]))
        else
        G[[i]] = t(G[[i]])
        colnames(G[[i]]) = paste("dim", c(1:ncomp[i]), sep = " ")
        
    }
    names(G)=names(object$X)
    
    ### Start: Maximum distance
    if (any(dist == "all") || any(dist == "max.dist"))
    {
        cls$max.dist = lapply(1:J, function(x){matrix(sapply(1:ncomp[x], ### List level
            function(y){apply(Y.hat[[x]][, , y, drop = FALSE], 1,  ### component level
                function(z){
                    paste(levels(Y)[which(z == max(z))], collapse = "/")
                }) ### matrix level
            }), nrow = nrow(newdata[[x]]), ncol = ncomp[x])
        })
        cls$max.dist = lapply(1:J, function(x){colnames(cls$max.dist[[x]]) = paste(rep("comp", ncomp[x]), 1 : ncomp[[x]], sep = " ");
            rownames(cls$max.dist[[x]]) = rownames(newdata[[x]]); return(cls$max.dist[[x]])})
        names(cls$max.dist)=names(object$X)
    }
    
    
    ### Start: Centroids distance
    if (any(dist == "all") || any(dist == "centroids.dist"))
    {
        cl = list()
        centroids.fun = function(x, G, h, i) {
            q = nrow(G[[i]])
            x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
            
            if (h > 1) {
                d = apply((x - G[[i]][, 1:h])^2, 1, sum)
            }
            else {
                d = (x - G[[i]][, 1])^2
            }
            cl.id = paste(levels(Y)[which(d == min(d))], collapse = "/")
        }
        
        for (i in 1 : J)
        {
            cl[[i]] = matrix(nrow = nrow(newdata[[1]]), ncol = ncomp[i])
            
            for (h in 1 : ncomp[[i]])
            {
                cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 1, function(x) {centroids.fun(x = x, G = G, h = h, i = i)})
                cl[[i]][, h] = cl.id
            }
        }
        
        cls$centroids.dist = lapply(1:J, function(x){colnames(cl[[x]]) = paste(rep("comp", ncomp[x]), 1 : ncomp[[x]], sep = " ");
            rownames(cl[[x]]) = rownames(newdata[[x]]); return(cl[[x]])})
        names(cls$centroids.dist)=names(object$X)
    }### End: Centroids distance
    
    
    ### Start: Mahalanobis distance
    if (any(dist == "all") || any(dist == "mahalanobis.dist"))
    {
        cl = list()
        Sr.fun = function(x, G, Yprim, h, i) {
            q = nrow(G[[i]])
            Xe = Yprim %*% G[[i]][, 1:h]
            #Xr = object$variates$X[, 1:h] - Xe
            Xr = variatesX[[i]][, 1:h] - Xe
            Sr = t(Xr) %*% Xr/nrow(Yprim)
            Sr.inv = solve(Sr)
            x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
            if (h > 1) {
                mat = (x - G[[i]][, 1:h]) %*% Sr.inv %*% t(x - G[[i]][, 1:h])
                d = apply(mat^2, 1, sum)
            } else {
                d = drop(Sr.inv) * (x - G[[i]][, 1])^2
            }
            cl.id = paste(levels(Y)[which(d == min(d))], collapse = "/")
        }
        
        for (i in 1 : J){
            cl[[i]] = matrix(nrow = nrow(newdata[[1]]), ncol = ncomp[i])
            
            for (h in 1:ncomp[[i]]) {
                cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 1, Sr.fun, G = G, Yprim = Y.prim, h = h, i = i)
                cl[[i]][, h] = cl.id
            }
        }
        
        cls$mahalanobis.dist = lapply(1:J, function(x){colnames(cl[[x]]) = paste(rep("comp", ncomp[x]), 1 : ncomp[[x]], sep = " ");
            rownames(cl[[x]]) = rownames(newdata[[x]]);return(cl[[x]])})
        names(cls$mahalanobis.dist)=names(object$X)
    } ### End: Mahalanobis distance
    
    out.DA$class = cls
    
    if (length(grep("block",class(object)))!=0 & J>1) # a block
    {
        out.DA$centroids = G
    }else{ #not a block
        out.DA$centroids = G[[1]]
        out.DA$class = lapply(out.DA$class,function(x){x[[1]]})
    }
    if (any(dist == "all"))
    dist = "all"
    
    out.DA$dist = dist
    ### End if discriminant analysis is performed
    
    # at this stage, we have the classification of each sample for each dataset of object$X
    # now we need to combine the classification by vote (majority wins), only when more than one block, otherwise 'vote' is classic classification
    if (length(object$X)>1)
    {
        for (ijk in 1:length(out.DA$class))# loop on the dist
        {
            # create a temporary array to make computation on the lists easier
            temp=array(, c(nrow(newdata[[1]]), min(ncomp), J))
            for(i in 1:J)
            {
                temp[, , i] = out.DA$class[[ijk]][[i]][, 1:min(ncomp)]
                
            }
            # look at the majority vote for all dataset of object$X (with table), if more than a unique max, we put NA
            table.temp = apply(temp,c(1,2), function(x){a=table(x); if (length(which(a==max(a)))==1) {b=names(which.max(a))}else{b=NA}; b})
            colnames(table.temp) = colnames(out.DA$class[[ijk]][[i]])[1:min(ncomp)]
            rownames(table.temp) = rownames(out.DA$class[[ijk]][[i]])
            out.DA$vote[[ijk]] = table.temp
        }
        names(out.DA$vote) = names(out.DA$class)
    }else{
        out.DA$vote = out.DA$class
    }
    out.DA
    
}






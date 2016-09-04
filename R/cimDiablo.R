#############################################################################################################
# Authors:
#   Amrit Singh, the University of British Columbia, Vancouver, Canada
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2015
# last modified: 19-08-2016
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



################################################
#
## 2) cim_diablo
#
################################################

# This function is a small wrapper of cim. For more customisation, please use cim

cimDiablo = function(object,
ncomp=1,
margins = c(2, 15),
legend.position="topright",
transpose = FALSE,
row.names = TRUE,
col.names = TRUE,
size.legend=1.5)
{
    
    # check input object
    if (!any(class(object) == "block.splsda"))
    stop("cimDiablo is only available for 'block.splsda' objects")

    if (length(object$X) <= 1)
    stop("This function is only available when there are more than 3 blocks") # so 2 blocks in X + the outcome Y

    if(ncomp > min(object$ncomp))
    stop("'ncomp' needs to be higher than object$ncomp")

    X = object$X
    Y = object$Y

    #need to reorder variates and loadings to put 'Y' in last
    indY = object$indY
    object$variates = c(object$variates[-indY], object$variates[indY])
    object$loadings = c(object$loadings[-indY], object$loadings[indY])
    
    #reducing loadings for ncomp
    object$loadings = lapply(object$loadings, function(x){x[, 1:ncomp, drop=FALSE]})
    
    keepA = lapply(object$loadings, function(i) apply(abs(i), 1, sum) > 0)
    XDatList = mapply(function(x, y){
        x[, y]
    }, x=X, y=keepA[-length(keepA)], SIMPLIFY=FALSE)
    XDat = do.call(cbind, XDatList)
    XDat[which(XDat > 2)] = 2
    XDat[which(XDat < -2)] = -2
    
    dark = brewer.pal(n = 12, name = 'Paired')[seq(2, 12, by = 2)]
    VarLabels = factor(rep(names(X), lapply(keepA[-length(keepA)], sum)), levels = names(X))#[order(names(X))])
    
    ## Plot heatmap
    opar = par()[! names(par()) %in% c("cin", "cra", "csi", "cxy", "din", "page")]
    par(mfrow=c(1,1))
    cim(XDat,transpose= transpose,
    row.names = row.names, col.names = col.names,
    col.sideColors = dark[as.numeric(VarLabels)],
    row.sideColors = color.mixo(as.numeric(Y)), margins = margins)
    
    legend(legend.position, c("Rows", c(levels(Y)[order(levels(Y))], "",
    "Columns", names(X))), col = c(1, color.mixo(1:nlevels(Y)), 1,
    1, dark[1:nlevels(VarLabels)][match(levels(VarLabels), names(X))]),
    pch = c(NA, rep(19, nlevels(Y)), NA, NA, rep(19, nlevels(VarLabels))), bty="n", cex = size.legend,
    text.font = c(2, rep(1, nlevels(Y)), NA, 2, rep(1, nlevels(VarLabels))))
    
    par(opar)
    
    return(invisible(XDat))
}


#############################################################################################################
# Author :
#   Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
#
# created:  27-07-2017
# last modified: 27-07-2017
#
# Copyright (C) 2017
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

t.test.process = function(mat.error.rate, alpha = 0.01)
{
    # mat.error.rate has nrep rows and ncomp columns
    # we test successively whether adding a component improves the results
    
    max = ncol(mat.error.rate) #number max of components included
    pval = NULL
    opt = 1 #initialise the first optimal number of genes
    for(j in 2:max)
    {
        pval[j] = NA
        temp = try(t.test(mat.error.rate[,opt],mat.error.rate[,j],alternative="greater")$p.value, silent=T) #t.test of "is adding X comp improves the overall results"
        if(any(class(temp) == "try-error") || is.na(temp)) # temp can be NaN when error.keepX is constant
        {
            ncomp_opt = NULL
            break
        } else {
            pval[j] = temp #had to be two steps (temp then pval =temp) otherwise the class is lost
        }
        if( (pval[j]< (alpha)))
        opt=j #if the p-value is lower than 0.05, the optimal number of comp is updated
    }
    if(all(class(temp) != "try-error") & !is.na(temp))
    ncomp_opt= opt
    
    return(ncomp_opt)
}


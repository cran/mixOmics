#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 15-04-2015
# last modified: 12-04-2016
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
# Calculate the explained variance of one dataset based on its variates
# ========================================================================================================


explained_variance = function(data, variates, ncomp)
{
    #check input data
    check = Check.entry.single(data, ncomp)
    data = check$X
    ncomp = check$ncomp
    
    isna = is.na(data)
    if (sum(isna > 0))
    {
        warning("NA values put to zero, results will differ from PCA methods used with NIPALS")
        data[isna] = 0
    }
    nor2x = sum((data)^2) # total variance in the data
	exp.varX = NULL
	for (h in 1:ncomp)
	{
		temp = variates[, h] / drop(t(variates[, h]) %*% (variates[, h]))
	    exp_var_new = as.numeric(t(variates[, h]) %*% data %*% t(data) %*% temp )/nor2x
	    exp.varX = append(exp.varX, exp_var_new)
	
	}
    names(exp.varX) = paste("comp", 1:ncomp)
    
    # result: vector of length ncomp with the explained variance per component
    exp.varX
}


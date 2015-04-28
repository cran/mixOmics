# Copyright (C) 2009 
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
# Leigh Coonan, Queensland Faculty for Advanced Bioinformatics, Australia
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

plot.pca <- function(x, ...){
    per.var = (x$sdev)/sum(x$sdev)
    per.var.vec=as.vector(per.var[1: x$ncomp])
    barplot(per.var.vec, names.arg = seq(1, x$ncomp,by=1), xlab="Principal Components", ylab="Proportion of Explained Variance")
}

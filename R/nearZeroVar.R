# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
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


nearZeroVar <- 
function (x, freqCut = 95/5, uniqueCut = 10) 
{
    if (is.vector(x)) 
        x = matrix(x, ncol = 1)
    freqRatio = apply(x, 2, function(data) {
        t = table(data[!is.na(data)])
        if (length(t) <= 1) {
            return(0)
        }
        w = which.max(t)
        return(max(t, na.rm = TRUE)/max(t[-w], na.rm = TRUE))
    })
    lunique = apply(x, 2, function(data) length(unique(data[!is.na(data)])))
    percentUnique = 100 * lunique/apply(x, 2, length)
    zeroVar = (lunique == 1) | apply(x, 2, function(data) all(is.na(data)))
	
    out = list()
	out$Position = which((freqRatio > freqCut & percentUnique <= uniqueCut) | zeroVar)
	names(out$Position) = NULL
    out$Metrics = data.frame(freqRatio = freqRatio, percentUnique = percentUnique)
    out$Metrics = out$Metrics[out$Position, ]
    out
}

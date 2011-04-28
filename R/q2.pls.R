# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
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


q2.pls <-
function(X, Y, ncomp, mode, M, fold, max.iter, tol)
{
	p = ncol(X)
	q = ncol(Y)
	n = nrow(X)
     
    X.temp = scale(X)
    Y.temp = scale(Y)
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    PRESS = Q2 = matrix(nrow = ncomp, ncol = q)
     
    #-- boucle sur h --#
    for (h in 1:ncomp) {
        
        #-- initialisation --#
        u = Y.temp[, 1] 
        a.old = 0
        b.old = 0
        iter = 1
         
        repeat {
            #-- compute loading vectors and variates associated to X --#
            a = crossprod(X.temp, u) / drop(crossprod(u))
            a = a / drop(sqrt(crossprod(a)))
            t = X.temp %*% a / drop(crossprod(a))
            
            #-- compute loading vectors and variates associated to Y --#		
            b = crossprod(Y.temp, t) / drop(crossprod(t))
            u = Y.temp %*% b / drop(crossprod(b))
				
            if ((crossprod(a - a.old) < tol) || (iter == max.iter)) break
             
            a.old = a
            iter = iter + 1
        }
		 
        RSS[h + 1, ] = colSums((Y.temp - t %*% t(b))^2) 
		
        #-- compute PRESS --#
        press = matrix(nrow = n, ncol = q)
		
        for (i in 1:M) {
            omit = fold[[i]]
            X.train = as.matrix(X.temp[-omit, ])
            Y.train = as.matrix(Y.temp[-omit, ])
            X.test = matrix(X.temp[omit, ], nrow = length(omit))
            Y.test = matrix(Y.temp[omit, ], nrow = length(omit))
            u.cv = Y.train[, 1] 
            a.old.cv = 0
            iter.cv = 1
         
            repeat {
                a.cv = crossprod(X.train, u.cv) / drop(crossprod(u.cv))
                a.cv = a.cv / drop(sqrt(crossprod(a.cv)))
                t.cv = X.train %*% a.cv
            		
                b.cv = crossprod(Y.train, t.cv) / drop(crossprod(t.cv))
                u.cv = Y.train %*% b.cv / drop(crossprod(b.cv))
				
                if ((crossprod(a.cv - a.old.cv) < tol) || (iter.cv == max.iter)) break
             
                a.old.cv = a.cv
                iter.cv = iter.cv + 1
            }
		    	
            Y.hat.cv = (X.test %*% a.cv) %*% t(b.cv)
            press[omit, ] = (Y.test - Y.hat.cv)^2	
        }
		 
        #-- deflation des matrices --#
        c = crossprod(X.temp, t) / drop(crossprod(t))		
        X.temp = X.temp - t %*% t(c)   
         
        #-- mode canonique --#
        if (mode == "canonical") {
            e = crossprod(Y.temp, u) / drop(crossprod(u))
            Y.temp = Y.temp - u %*% t(e)
        }
         
        #-- mode classic --#
        if(mode == "classic") Y.temp = Y.temp - t %*% t(b)                 
         
        #-- mode regression --#
        if(mode == "regression") {
            d = crossprod(Y.temp, t) / drop(crossprod(t))
            Y.temp = Y.temp - t %*% t(d)
        }
		
        #-- mode invariant --#
        if (mode == "invariant") Y.temp = Y	
          
        PRESS[h, ] = colSums(press)
        Q2[h, ] = 1 - PRESS[h, ]/RSS[h, ]
     
    } #-- fin boucle sur h --#
     
    #-- valeurs sortantes --#
    if (q > 1) {
        Q2.total = 1 - rowSums(PRESS/rowSums(RSS[-ncomp, ]))
        Q2 = cbind(Q2, Q2.total)
    }
	 
    return(invisible(Q2))
}

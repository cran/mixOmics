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


q2.spls <-
function(X, Y, ncomp, mode, keepX, keepY, M, fold, max.iter, tol)
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
        nx = p - keepX[h]
        ny = q - keepY[h]
         
        #-- svd de M = t(X)*Y --#         
        XY = crossprod(X.temp, Y.temp)
        svd.XY = svd(XY, nu = 1, nv = 1)
        a.old = svd.XY$u
        b.old = svd.XY$v
         
        #-- latent variables --#
        t = X.temp %*% a.old / drop(crossprod(a.old))
        t = t / drop(sqrt(crossprod(t)))
         
        u = Y.temp %*% b.old / drop(crossprod(b.old))
        u = u / drop(sqrt(crossprod(u)))
         
        iter = 1
         
        repeat {
            a = t(X.temp) %*% u
            b = t(Y.temp) %*% t
             
            if (nx != 0) { 
                a = ifelse(abs(a) > abs(a[order(abs(a))][nx]), 
                    (abs(a) - abs(a[order(abs(a))][nx])) * sign(a), 0)
            }
            a = a / drop(sqrt(crossprod(a)))
		     
            if (ny != 0) {
                b = ifelse(abs(b) > abs(b[order(abs(b))][ny]),
                    (abs(b) - abs(b[order(abs(b))][ny])) * sign(b), 0)
            }
            b = b / drop(sqrt(crossprod(b)))
			 
            t = X.temp %*% a / drop(crossprod(a))
            t = t / drop(sqrt(crossprod(t)))
               
            u = Y.temp %*% b / drop(crossprod(b))
            u = u / drop(sqrt(crossprod(u)))
             
            if ((crossprod(a - a.old)) < tol || (iter == max.iter)) break
             
            a.old = a
            b.old = b
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
                    
            XY = crossprod(X.train, Y.train)
            svd.XY = svd(XY, nu = 1, nv = 1)
            a.old.cv = svd.XY$u
            b.old.cv = svd.XY$v
             
            #-- latent variables --#
            t.cv = X.train %*% a.old.cv / drop(crossprod(a.old.cv))
            t.cv = t.cv / drop(sqrt(crossprod(t.cv)))
             
            u.cv = Y.train %*% b.old.cv / drop(crossprod(b.old.cv))
            u.cv = u.cv / drop(sqrt(crossprod(u.cv)))
             
            iter.cv = 1
             
            repeat {
                a.cv = t(X.train) %*% u.cv
                b.cv = t(Y.train) %*% t.cv
                 
                if (nx != 0) { 
                    a.cv = ifelse(abs(a.cv) > abs(a.cv[order(abs(a.cv))][nx]), 
                           (abs(a.cv) - abs(a.cv[order(abs(a.cv))][nx])) * sign(a.cv), 0)
                }
                a.cv = a.cv / drop(sqrt(crossprod(a.cv)))
		         
                if (ny != 0) {
                    b.cv = ifelse(abs(b.cv) > abs(b.cv[order(abs(b.cv))][ny]),
                           (abs(b.cv) - abs(b.cv[order(abs(b.cv))][ny])) * sign(b.cv), 0)
                }
                b.cv = b.cv / drop(sqrt(crossprod(b.cv)))
			     
                t.cv = X.train %*% a.cv / drop(crossprod(a.cv))
                t.cv = t.cv / drop(sqrt(crossprod(t.cv)))
                 
                u.cv = Y.train %*% b.cv / drop(crossprod(b.cv))
                u.cv = u.cv / drop(sqrt(crossprod(u.cv)))
                 
                if ((crossprod(a.cv - a.old.cv)) < tol || (iter.cv == max.iter)) break
                 
                a.old.cv = a.cv
                b.old.cv = b.cv
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
         
        #-- mode regression --#
        if(mode == "regression") {
            d = crossprod(Y.temp, t) / drop(crossprod(t))
            Y.temp = Y.temp - t %*% t(d)
        }
          
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

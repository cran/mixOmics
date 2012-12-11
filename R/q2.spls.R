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
function (X, Y, ncomp, mode, M, fold, max.iter, tol, keepX, keepY) 
{
    p = ncol(X)
    q = ncol(Y)
    n = nrow(X)

    X.temp = scale(X)
    Y.temp = scale(Y)
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    PRESS = Q2 = matrix(nrow = ncomp, ncol = q)
    keepX = !keepX
    keepY = !keepY
	
    #-- boucle sur h --#
    for (h in 1:ncomp) {
        #-- svd de Mat = t(X)*Y --#
        X.aux = X.temp        
        Y.aux = Y.temp         
         
        Mat = crossprod(X.aux, Y.aux)
        svd.M = svd(Mat, nu = 1, nv = 1)
        a.old = svd.M$u
        b.old = svd.M$v
         
        #-- latent variables --#
        t = X.temp %*% a.old / drop(crossprod(a.old))
        t = t / drop(sqrt(crossprod(t)))
         
        u = Y.temp %*% b.old / drop(crossprod(b.old))
        u = u / drop(sqrt(crossprod(u)))
         
        iter = 1
         
        #-- boucle jusqu'à convergence de a et de b --#
        repeat {
            a = t(X.temp) %*% u
            b = t(Y.temp) %*% t
             
            a[keepX[, h]] = 0
            a = a / drop(crossprod(u))
            a = a / drop(sqrt(crossprod(a)))
		     
            b[keepY[, h]] = 0
            b = b / drop(crossprod(t))
			
            t = X.temp %*% a / drop(crossprod(a))
            t = t / drop(sqrt(crossprod(t)))
             
            u = Y.temp %*% b / drop(crossprod(b))
            u = u / drop(sqrt(crossprod(u)))
           
            if ((crossprod(a - a.old)) < tol || (iter == max.iter)) break
             
            a.old = a
            b.old = b
            iter = iter + 1
        }

        #-- added for sparse --#		
        Y.sel = Y.temp
        Y.sel[, keepY[h, ]] = 0		
        RSS[h + 1, ] = colSums((Y.sel - t %*% t(b))^2)

        #-- compute PRESS --#
        press = matrix(nrow = n, ncol = q)
		
        for (i in 1:M) {
            omit = fold[[i]]
            X.train = as.matrix(X.temp[-omit, ])
            Y.train = as.matrix(Y.temp[-omit, ])
            X.test = matrix(X.temp[omit, ], nrow = length(omit))
            Y.test = matrix(Y.temp[omit, ], nrow = length(omit))
            #u.cv = Y.train[, 1]
            #a.old.cv = 0
            #iter.cv = 1
			
            #-- svd de Mat = t(X)*Y --#         
            Mat = crossprod(X.train, Y.train)
            svd.M = svd(Mat, nu = 1, nv = 1)
            a.old.cv = svd.M$u
            b.old.cv = svd.M$v
         
            #-- latent variables --#
            t.cv = X.train %*% a.old.cv / drop(crossprod(a.old.cv))
            t.cv = t.cv / drop(sqrt(crossprod(t.cv)))
         
            u.cv = Y.train %*% b.old.cv / drop(crossprod(b.old.cv))
            u.cv = u.cv / drop(sqrt(crossprod(u.cv)))
         
            iter.cv = 1
			
            repeat {
                a.cv = crossprod(X.train, u.cv)/drop(crossprod(u.cv))
                a.cv = a.cv/drop(sqrt(crossprod(a.cv)))
                t.cv = X.train %*% a.cv
                b.cv = crossprod(Y.train, t.cv)/drop(crossprod(t.cv))
                u.cv = Y.train %*% b.cv/drop(crossprod(b.cv))

                a.old.cv = a.cv
                iter.cv = iter.cv + 1
                a.cv = t(X.train) %*% u.cv
                b.cv = t(Y.train) %*% t.cv
             
                a.cv[keepX[, h]] = 0
                a.cv = a.cv / drop(crossprod(u.cv))
                a.cv = a.cv / drop(sqrt(crossprod(a.cv)))
		     
                b.cv[keepY[, h]] = 0
                b.cv = b.cv / drop(crossprod(t.cv))
			
                t.cv = X.train %*% a.cv / drop(crossprod(a.cv))
                t.cv = t.cv / drop(sqrt(crossprod(t.cv)))
             
                u.cv = Y.train %*% b.cv / drop(crossprod(b.cv))
                u.cv = u.cv / drop(sqrt(crossprod(u.cv)))
           
                if ((crossprod(a.cv - a.old.cv) < tol) || (iter.cv == max.iter)) 
                    break
             
                a.old.cv = a.cv
                b.old.cv = b.cv
                iter.cv = iter.cv + 1				
            }
			
            Y.hat.cv = (X.test %*% a.cv) %*% t(b.cv)
            press[omit, ] = (Y.test - Y.hat.cv)^2

            Y.hat.cv = (X.test %*% a.cv) %*% t(b.cv)
            Y.sel = Y.test
            Y.sel[, keepY[h, ]] = 0
            press[omit, ] = (Y.sel - Y.hat.cv)^2
        }

        #-- deflation des matrices --#
        c = crossprod(X.temp, t)/drop(crossprod(t))
        X.temp = X.temp - t %*% t(c)

        #-- mode canonique --#
        if (mode == "canonical") {
            e = crossprod(Y.temp, u)/drop(crossprod(u))
            Y.temp = Y.temp - u %*% t(e)
        }

        #-- mode classic --#
        if (mode == "classic") 
            Y.temp = Y.temp - t %*% t(b)

        #-- mode regression --#
        if (mode == "regression") {
            d = crossprod(Y.temp, t)/drop(crossprod(t))
            Y.temp = Y.temp - t %*% t(d)
        }

        #-- mode invariant --#
        if (mode == "invariant") 
            Y.temp = Y

        PRESS[h, ] = colSums(press)
        PRESS[h, keepY[h, ]] = NA
        RSS[h + 1, keepY[h, ]] = NA
        Q2[h, ] = 1 - PRESS[h, ] / RSS[h, ]

    } #-- fin boucle sur h --#

    #-- valeurs sortantes --#
    if (q > 1) {
      Q2.total = 1 - rowSums(PRESS, na.rm = TRUE) / rowSums(RSS[-ncomp, ], na.rm = TRUE)
      Q2 = cbind(Q2, Q2.total)
    }
    
    return(invisible(Q2))
}

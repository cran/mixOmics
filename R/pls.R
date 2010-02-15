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



`pls` <-
function(X, Y, ncomp = 3, mode = c("regression", "canonical", "invariant", "classic"),
max.iter = 500, tol = 1e-06, scaleY = TRUE)
{


#-- validation des arguments --#
if (length(dim(X)) != 2) 
        stop("'X' must be a numeric matrix.")

X = as.matrix(X)
Y = as.matrix(Y)


#if (length(dim(Y)) == 0)
#Y = as.matrix(Y, ncol = 1)
#else
#Y = as.matrix(Y)

if (!is.numeric(X) || !is.numeric(Y)) 
        stop("'X' and/or 'Y' must be a numeric matrix.")

n = nrow(X)
p = ncol(X)
q = ncol(Y)

if ((n != nrow(Y))) 
stop("unequal number of rows in 'X' and 'Y'.")

#if (missing(ncomp))
#    ncomp = mat.rank(X)$rank

if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
stop("invalid number of variates, 'ncomp'")

if(ncomp > p) {
warning("Reset maximum number of variates 'ncomp' to rank(X).")
ncomp = mat.rank(X)$rank
}

ncomp = round(ncomp)
mode = match.arg(mode)

#-- initialisation des matrices --#
X.names = dimnames(X)[[2]]
if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")

if (dim(Y)[2] == 1) Y.names = "Y"
else {
Y.names = dimnames(Y)[[2]]
if (is.null(Y.names)) Y.names = paste("Y", 1:q, sep = "")
}

ind.names = dimnames(X)[[1]]
if (is.null(ind.names)) {
ind.names = dimnames(Y)[[1]]
rownames(X) = ind.names
}
if (is.null(ind.names)) {
ind.names = 1:n
rownames(X) = rownames(Y) = ind.names
}

#-- centrer et réduire les données --#
X = scale(X, center = TRUE, scale = TRUE)
if(scaleY == TRUE) Y = scale(Y, center = TRUE, scale = TRUE) 

X.temp = X.na = X
Y.temp = Y.na = Y
#-- if NA values, set them to zero
ind.X.na = is.na(X.na)
ind.Y.na = is.na(Y.na)
X.temp[ind.X.na] = 0
Y.temp[ind.Y.na] = 0

mat.t = matrix(nrow = n, ncol = ncomp)
mat.u = matrix(nrow = n, ncol = ncomp)
mat.a = matrix(nrow = p, ncol = ncomp)
mat.b = matrix(nrow = q, ncol = ncomp)
mat.c = matrix(nrow = p, ncol = ncomp)
#KA: I add this matrix for regression mode
mat.d = matrix(nrow = q, ncol = ncomp)

#-- boucle sur h --#
for (h in 1:ncomp) {

#-- initialisation --#
u = Y.temp[, 1] 
a.old = 0
b.old = 0
iter = 1

repeat {
#--compute loading vectors associated to X
a = t(X.temp) %*% u / drop(crossprod(u))
a = a / drop(sqrt(crossprod(a)))

#--compute latent variable associated to X
t = X.temp %*% a / drop(crossprod(a))

#--compute loading vectors associated to Y
# PLS2 mode classic alone as in tenenhaus p.128

#if(mode != "classic") {b = t(Y.temp) %*% t / drop(crossprod(t)); b = b / drop(sqrt(crossprod(b)))} else {b = t(Y.temp) %*% t / drop(crossprod(t))}
b = t(Y.temp) %*% t / drop(crossprod(t))

#--compute latent variable associated to Y
# mode classic alone
u = Y.temp %*% b / drop(crossprod(b))

if (crossprod(a - a.old) < tol & crossprod(b - b.old) < tol) break

if (iter == max.iter) {
warning(paste("Maximum number of iterations reached for dimension", h),
call. = FALSE)
break
}

a.old = a
b.old = b
iter = iter + 1
}

#-- deflation des matrices --#
c = t(X.temp) %*% t / drop(crossprod(t))
X.na = X.na - t %*% t(c)   
#X.temp = X.temp - t %*% t(c)   

#-- mode canonique --#
if (mode == "canonical") {
e = t(Y.temp) %*% u / drop(crossprod(u))
Y.na = Y.na - u %*% t(e)
#Y.temp = Y.temp - u %*% t(e)

}

#-- mode classic --#
if(mode == "classic") Y.na = Y.na - t %*% t(b)                 


#-- mode regression --#
if(mode == "regression") {
d = t(Y.temp) %*% t /drop((t(t) %*% t))
Y.na = Y.na - t %*% t(d)                 
}

#-- mode invariant --#
if (mode == "invariant") {
Y.na = Y 
#Y.temp = Y 
}

X.temp = X.na
Y.temp = Y.na
X.temp[ind.X.na] = 0
Y.temp[ind.Y.na] = 0

mat.t[, h] = t
mat.u[, h] = u
mat.a[, h] = a
mat.b[, h] = b
mat.c[, h] = c
##mat.d[, h] = d    #added

} #-- fin boucle sur h --#

#-- valeurs sortantes --#
rownames(mat.a) = rownames(mat.c) = X.names
rownames(mat.b) = Y.names
rownames(mat.t) = rownames(mat.u) = ind.names

comp = paste("comp", 1:ncomp)
colnames(mat.t) = colnames(mat.u) = comp
colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = comp 

cl = match.call()
cl[[1]] = as.name('pls')
result = list(
	call=cl,
	X = X, 
	Y = Y, 
	ncomp = ncomp, 
	mode = mode, 
	mat.c = mat.c,
	variates = list(X = mat.t, Y = mat.u),
	loadings = list(X = mat.a, Y = mat.b), 
	names = list(X = X.names, Y = Y.names, indiv = ind.names)
)

class(result) = "pls"
return(invisible(result))
}


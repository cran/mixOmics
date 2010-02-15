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

`spls` <-
function(
	X, 
	Y, 
	ncomp = 2, 
	mode = c("regression", "canonical"),
	max.iter = 500, 
	tol = 1e-06,
	keepX = c(rep(ncol(X), ncomp)), 
	keepY = c(rep(ncol(Y), ncomp)),
	scaleY = TRUE    
	)
{

mode = match.arg(mode)

##if (missing(keepX) & missing(keepY)) {
##result = pls(X, Y, ncomp = ncomp, max.iter = max.iter, tol = tol, mode = mode)
##}
##else {
#-- validation des arguments --#
if (length(dim(X)) != 2) 
stop("'X' must be a numeric matrix.")

X = as.matrix(X)
Y = as.matrix(Y)

if (!is.numeric(X) || !is.numeric(Y)) 
stop("'X' and/or 'Y' must be a numeric matrix.")

n = nrow(X)
p = ncol(X)
q = ncol(Y)

if ((n != nrow(Y))) 
stop("unequal number of rows in 'X' and 'Y'.")

#if (missing(ncomp))
#ncomp = mat.rank(X)$rank

if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
stop("invalid number of variates, 'ncomp'")

if(ncomp > p) 
stop("the number of variates 'ncomp' must be lower or equal than ", p, ".", 
call. = FALSE)

#-- initialisation des matrices --#
X.names = dimnames(X)[[2]]
if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")

Y.names = dimnames(Y)[[2]]
if (is.null(X.names)) Y.names = paste("Y", 1:q, sep = "")

ind.names = dimnames(X)[[1]]
if (is.null(ind.names)) {
ind.names = dimnames(Y)[[1]]
rownames(X) = ind.names
}

if (is.null(ind.names)) {
ind.names = 1:n
rownames(X) = rownames(Y) = ind.names
}

if (length(keepX) != ncomp) 
stop("length of 'keepX' must be equal to ", ncomp)

if (length(keepY) != ncomp) 
stop("length of 'keepY' must be equal to ", ncomp)

if (any(keepX > p)) 
stop("each component of 'keepX' must be lower or equal than ", p, ".",
call. = FALSE)

if (any(keepY > q)) 
stop("each component of 'keepY' must be lower or equal than ", q, ".",
call. = FALSE)

#-- centrer et réduire les données --#
X = scale(X, center = TRUE, scale = TRUE)
if(scaleY == TRUE) Y = scale(Y, center = TRUE, scale = TRUE) 


X.temp = X.na = X
Y.temp = Y.na = Y
ind.X.na = is.na(X.na)
ind.Y.na = is.na(Y.na)
X.temp[ind.X.na] = 0
Y.temp[ind.Y.na] = 0

mat.t = matrix(nrow = n, ncol = ncomp)
mat.u = matrix(nrow = n, ncol = ncomp)
mat.a = matrix(nrow = p, ncol = ncomp)
mat.b = matrix(nrow = q, ncol = ncomp)
mat.c = matrix(nrow = p, ncol = ncomp)
#KA: added:
mat.d = matrix(nrow = q, ncol = ncomp)

#-- boucle sur h --#
for (h in 1:ncomp) {
nx = p - keepX[h]
ny = q - keepY[h]

#-- svd de M = t(X)*Y --#
M = crossprod(X.temp, Y.temp)
svd.M = svd(M, nu = 1, nv = 1)
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

if (nx != 0) { 
a = ifelse(abs(a) > abs(a[order(abs(a))][nx]), 
(abs(a) - abs(a[order(abs(a))][nx])) * sign(a), 0)
a = a / drop(sqrt(crossprod(a)))
}

if (ny != 0) {
b = ifelse(abs(b) > abs(b[order(abs(b))][ny]),
(abs(b) - abs(b[order(abs(b))][ny])) * sign(b), 0)
b = b / drop(sqrt(crossprod(b)))
}

t = X.temp %*% a / drop(crossprod(a))
t = t / drop(sqrt(crossprod(t)))

u = Y.temp %*% b / drop(crossprod(b))
u = u / drop(sqrt(crossprod(u)))

if (crossprod(a - a.old) < tol & crossprod(b - b.old) < tol) break

if (iter == max.iter) {
warning(paste("Maximum number of iterations reached for the component", h),
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

#-- mode canonique --#
if (mode == "canonical") {
e = t(Y.temp) %*% u / drop(crossprod(u))
Y.na = Y.na - u %*% t(e)
}

#-- mode regression --#
if (mode == "regression") {
d = t(Y.temp) %*% t /drop((t(t) %*% t))
Y.na = Y.na - t %*% t(d)
}

X.temp = X.na
Y.temp = Y.na
X.temp[ind.X.na] = 0
Y.temp[ind.Y.na] = 0

mat.a[, h] = a
mat.b[, h] = b
mat.c[, h] = c
if (mode == "regression") mat.d[, h] = d   #KA: added
mat.t[, h] = t
mat.u[, h] = u  

} #-- fin boucle sur h --#

#-- valeurs sortantes --#
rownames(mat.a) = rownames(mat.c) = X.names
rownames(mat.b) = Y.names
rownames(mat.t) = rownames(mat.u) = ind.names

dim = paste("comp", 1:ncomp)
colnames(mat.t) = colnames(mat.u) = dim
colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = dim 

cl = match.call()
cl[[1]] = as.name('spls')

result = list(
call = cl,
X = X, Y = Y, ncomp = ncomp, mode = mode, 
keepX = keepX,
keepY = keepY,
mat.c = mat.c, mat.t = mat.t, 
variates = list(X = mat.t, Y = mat.u),
loadings = list(X = mat.a, Y = mat.b),
names = list(X = X.names, Y = Y.names, indiv = ind.names))

class(result) = c("spls", "pls") 
return(invisible(result))
##}  #end if keepX
}


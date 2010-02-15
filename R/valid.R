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


# -------------------- RMSEP/q2/PRESS/RSS FOR BOTH PLS AND sPLS------------------

#Note: RMSEP computation should be the same for the three modes when using PLS if same seed (for CV) or if LOOCV

#'valid' = function(X, ...) UseMethod('valid')


`valid` <-
function(
	X, 
	Y, 
	ncomp = 3, 
	mode = c("regression", "invariant", "classic"),
	max.iter = 500, 
	tol = 1e-06,
        criterion = c('rmsep', 'q2'), 
	method = c('pls', 'spls'),
	keepX = if(method =='pls') NULL else c(rep(ncol(X), ncomp)),
	keepY = if(method =='pls') NULL else c(rep(ncol(Y), ncomp)),
	scaleY = TRUE,
	validation = c('loo', 'Mfold'),
	M = if(validation == 'Mfold') 10 else nrow(X)
	){

#do warning for mode + other warnings ?
if (missing(validation)) stop("Choose between cross-validation: 'Mfold' or loo-cv: 'loo'")
if (missing(method)) stop("Choose a method: 'pls' or 'spls' ")
if (missing(criterion)) stop("Choose a validation criterion")
if (missing(mode)) stop("Choose a mode")
if (mode == 'canonical') stop("Only regression, classic or invariant mode are allowed !")
if ((method == 'spls') & (mode == 'invariant')) stop("No invariant mode with sPLS")


if (length(dim(X)) != 2) 
stop("'X' must be a numeric matrix.")

X = as.matrix(X)
Y = as.matrix(Y)

if (!is.numeric(X) || !is.numeric(Y)) 
stop("'X' and/or 'Y' must be a numeric matrix.")


n = nrow(X)
p = ncol(X)
q = ncol(Y)


#regression coeff
B.all <- array(0, dim = c(p, q, ncomp))      #on all data  (for RSS)
B.hat <- array(0, dim = c(p, q, ncomp))      # on test set (for RMSEP)

#Y prediction
Y.all <- array(0, dim = c(n, q, ncomp))  #on all data
Y.hat.press <- array(0, dim = c(n, q, ncomp))  # on test set  to be consistent with tenenhaus
Y.hat.rmsep <- array(0, dim = c(n, q, ncomp))  # on test set  to be consistent with package pls

coef = matrix(0, ncol = p + 1, nrow = q)

# ---criteria
PRESS = matrix(0, nrow = q, ncol = ncomp)
RSS = matrix(0, nrow = q, ncol = ncomp) 
q2 = vector(length=ncomp)
q2V = matrix(0, nrow = q, ncol = ncomp)
q2cum = vector(length=ncomp)
q2Vcum = matrix(0, nrow = q, ncol = ncomp)


#--center and scale data --#
X = scale(X, center = TRUE, scale = TRUE)
if(scaleY == TRUE) Y = scale(Y, center = TRUE, scale = TRUE) 



#compute B.hat for all variables to compute RSS
if(any(criterion =='q2')){

if(method == 'pls') {object = pls(X = X , Y = Y, ncomp = ncomp, mode = mode, max.iter = max.iter, tol = tol)} 
else 
{object = spls(X = X , Y = Y, ncomp = ncomp, mode = mode, 
max.iter = max.iter, tol = tol, keepX = keepX, keepY = keepY)}

a.all = object$loadings$X
b.all = object$loadings$Y
c.all = object$mat.c

for (h in 1:ncomp) {
	W.all = a.all[, 1:h] %*% solve(t(c.all[, 1:h]) %*% a.all[, 1:h]) 
	if(q==1){B.all[,,h] = W.all %*% as.vector(t(b.all[, 1:h]))} else {B.all[,,h] = W.all %*% t(b.all[, 1:h])}
	Y.all[,,h] = X %*% B.all[,,h]

	RSS[, h] = apply((Y - Y.all[,,h])^2, 2, sum)    # on all data
}  #end h
} #end if


#-----M fold cross validation or loo
# define the folds
if (validation=='Mfold') {fold = sample(c(1:M), n, replace=TRUE)} else {fold=c(1:n)}


for(i in 1:M){

if(validation=='Mfold'){
	X.train = X[fold!=i,]
	Y.train = Y[c(fold!=i),]
	X.test = X[fold==i,]
}

if(validation=='loo'){
	X.train = X[-i,]
	Y.train = Y[-i,]
	X.test = X[i,]
}

# -- pls or spls
if(method == 'pls') {object = pls(X = X.train , Y = Y.train, ncomp = ncomp, 
mode = mode, max.iter = max.iter, tol = tol)} 
else 
{object = spls(X = X.train , Y = Y.train, ncomp = ncomp, 
mode = mode, max.iter = max.iter, tol = tol, keepX = keepX, keepY = keepY)}

a = object$loadings$X
b = object$loadings$Y
c = object$mat.c

for(h in 1:ncomp){
	W = a[, 1:h] %*% solve(t(c[, 1:h]) %*% a[, 1:h]) 
	if(q==1){B.hat[,,h] = W %*% as.vector(t(b[, 1:h]))} else {B.hat[,,h] = W %*% t(b[, 1:h])}
	Y.hat.rmsep[fold==i,,h] = X.test %*% B.hat[,,h] + colMeans(as.matrix(Y.train))   #to be consistent with pls package
	Y.hat.press[fold==i,,h] = X.test %*% B.hat[,,h]                       #to be consistent with Tenenhaus
	
	PRESS[, h] = apply((Y - Y.hat.press[,,h])^2, 2, sum)  # on test data

}  #end h

}  # end i


# ----compute RMSEP --------------------
if(any(criterion =='rmsep')){
rmsep.mat = matrix(nrow=q, ncol=ncomp)
for (j in 1:q){
	if(ncomp==1) {rmsep.mat[j,] = sqrt(mean((Y.hat.rmsep[,j,1] - Y[,j])^2))} else {rmsep.mat[j,] = sqrt(colMeans((Y.hat.rmsep[,j,1:ncomp] - Y[,j])^2))}
	}

intercept = sqrt(mean((Y[,1])^2))
rmsep.mat = cbind(c(rep(intercept, q)), rmsep.mat)
colnames(rmsep.mat) = c('intercept', paste('dim', c(1:ncomp), sep=''))
rownames(rmsep.mat) = colnames(Y)
}  #end if


# ----- compute q2 --------------------
if(any(criterion =='q2')){

if(q==1){RSS.0 = c(rep(n - 1, q), RSS[, -ncomp])}else{RSS.0 = cbind(rep(n - 1, q), RSS[, -ncomp])}

# --q2 = 1 - PRESS/RSS
if(q==1){q2 = 1- PRESS/RSS.0}else{q2 = 1- apply(PRESS, 2, sum)/apply(RSS.0, 2, sum)}

# --q2V = 1 - PRESS/RSS   option 1
q2V = 1- PRESS/RSS.0      #if q=1, then q2V = q2


if(q!=1){rownames(q2V) = rownames(q2Vcum) = colnames(Y)}
rownames(PRESS) = rownames(RSS) = colnames(Y)
colnames(PRESS) = colnames(RSS) = paste("comp", 1:ncomp)
if(q==1){colnames(q2) = colnames(q2V) = paste("comp", 1:ncomp)}else{names(q2) = colnames(q2V) = paste("comp", 1:ncomp)}

}  #end if


return(invisible(list(
Y.hat = Y.hat.rmsep, 
fold=fold,
rmsep = if (any(criterion == 'rmsep')) rmsep.mat else NULL, 
Q2 = if (any(criterion == 'q2')) { list(RSS=RSS, PRESS=PRESS, q2 = q2, q2V = q2V)} else NULL 

)))
}




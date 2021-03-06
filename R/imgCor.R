#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
#
# created: 2009
# last modified: 2015
#
# Copyright (C) 2009
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


# -----------------------------------------
# --- small example but see help file
# ----------------------------------------
#library(mixOmics)
#data(liver.toxicity)
#X <- liver.toxicity$gene
#Y <- liver.toxicity$clinic
#imgCor(X=X[,1:10],Y=Y[,1:4])


imgCor <-
function(X,
         Y,  
         type = "combine", 
         X.var.names = TRUE, 
         Y.var.names = TRUE,  	 
         sideColors = TRUE,
         interactive.dev = TRUE,
         title = TRUE,
         color, row.cex, col.cex,symkey, keysize, 
         xlab, ylab, margins, lhei, lwid) 
{
  #-- checking general input arguments ---------------------------------------#
  #---------------------------------------------------------------------------#
  
  #-- check that the user did not enter extra arguments
  arg.call = match.call()
  user.arg = names(arg.call)[-1]
  
  err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())), 
                 error = function(e) e)
  
  if ("simpleError" %in% class(err))
    stop(err[[1]], ".", call. = FALSE)
  
  default.arg = c("color", "row.cex", "col.cex","symkey", "keysize", 
                  "xlab", "ylab", "margins", "lhei", "lwid")
  function.arg = c(names(mget(names(formals()), sys.frame(sys.nframe()))),
                   default.arg)
  not.arg = !(user.arg %in% function.arg)
  
  if (any(not.arg)) {
    unused.arg = user.arg[not.arg]
    not.arg = which(not.arg) + 1
    output = rep("", length(not.arg))
    
    for (i in 1:length(not.arg)) {
      output[i] = paste0(unused.arg[i], " = ", arg.call[[not.arg[i]]])
    }
    
    output = paste0("(",paste(output, collapse = ", "), ").")
    msg = "unused argument "
    if (length(not.arg) > 1) msg ="unused arguments "  
    stop(msg, output, call. = FALSE)
  }
  
  #-- X and Y matrices
  if (is.data.frame(X)) X = as.matrix(X)
  if (is.data.frame(Y)) Y = as.matrix(Y)
  
  if (!is.matrix(X) || !is.matrix(Y)) 
    stop("'X' and/or 'Y' must be a numeric matrix.", call. = FALSE)
  
  if (!is.numeric(X) || !is.numeric(Y)) 
    stop("'X' and/or 'Y' must be a numeric matrix.", call. = FALSE)
  
  p = ncol(X)
  q = ncol(Y)  
  
  #-- X.var.names
  if (!is.logical(X.var.names)) 
    stop("'X.var.names' must be a logical constant (TRUE or FALSE).",
         call. = FALSE)
  
  #-- Y.var.names
  if (!is.logical(Y.var.names)) 
    stop("'Y.var.names' must be a logical constant (TRUE or FALSE).",
         call. = FALSE)
  
  #-- type
  choices = c("combine", "separate")
  type = choices[pmatch(type, choices)]
  
  if (is.na(type)) 
    stop("'type' should be one of 'combine' or 'separate'.", 
         call. = FALSE)
  
  #-- sideColors
  if (is.logical(sideColors)) {
    if(isTRUE(sideColors)) sideColors = c("blue", "red") else sideColors = NULL
  }
  else {
    if (length(sideColors) != 2)
      stop("'sideColors' must be a character vector of length ", 2, ".", 
           call. = FALSE)
  }
  
  isColor = sapply(sideColors, function(x) { 
    tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE) })
  
  if (any(!isColor)) 
    stop("'sideColors' must be a character vector of recognized colors.", 
         call. = FALSE)
  
  #-- interactive.dev
  if (!is.logical(interactive.dev))
    stop("'interactive.dev' must be a logical constant (TRUE or FALSE).",
         call. = FALSE)
  
  #-- title
  if (!is.logical(title))
    stop("'title' must be a logical constant (TRUE or FALSE).",
         call. = FALSE)
  
  #-- end checking --#
  #------------------#
  
  
  #-- calling the cim function for mapping -----------------------------------#
  #---------------------------------------------------------------------------#
  
  #-- if type = "combine" --#
    if (type == "combine") {		
      
      if (isTRUE(X.var.names)) X.var.names = colnames(X)
      else X.var.names = rep("", p)
      
      if (isTRUE(Y.var.names)) Y.var.names = colnames(Y)
      else Y.var.names = rep("", q)
      
      mat.lab = c(X.var.names, Y.var.names)
      
      matcor = cor(cbind(X, Y), use = "pairwise")
      colnames(matcor) = rownames(matcor) = mat.lab
        matcor=matcor[(p + q):1,]
        
        row.sideColors = col.sideColors = NULL
        
        if (!is.null(sideColors)) {
          bg.col = par("bg")
          row.sideColors = c(rep(sideColors[2], q),rep(sideColors[1], p))
          print(length(row.sideColors))
          col.sideColors = rev(row.sideColors)
        }
        
        cim(matcor, cluster = "none", 
            title = if(title) "[X,Y] correlation matrix" else NULL,
            col.sideColors = col.sideColors,
            row.sideColors = row.sideColors,
            row.names = TRUE, col.names = TRUE)
    }
     
    # representation des matrices de correlation de #
    # X, Y et entre X et Y                          #
    #-----------------------------------------------#
    if (type == "separate") {
        Xcor = cor(X, use = "pairwise")
        Ycor = cor(Y, use = "pairwise")
        XYcor = cor(X, Y, use = "pairwise")

        Xcor = Xcor[p:1,]#t(Xcor[p:1, ])
        Ycor = Ycor[q:1,]#t(Ycor[q:1, ])
        XYcor = XYcor[p:1,]
		
        if (interactive.dev) {   
          cim(Xcor, cluster = "none",
              title = if(title) "X correlation matrix" else NULL,
              row.names = X.var.names, col.names = X.var.names)
          
          devAskNewPage(TRUE)
          cim(Ycor, cluster = "none", 
              title = if(title) "Y correlation matrix" else NULL,
              row.names = Y.var.names, col.names = Y.var.names)
          
          cim(XYcor, cluster = "none",
              title = if(title) "XY correlation matrix" else NULL,			
              row.names = X.var.names, col.names = Y.var.names)
        }
		else {
		  GD = getOption("device")
		  if (is.character(GD)) {
		    if (GD == "RStudioGD") GD = FALSE
		  }
		  else 
		    GD = TRUE
		  
		  if (isTRUE(GD)) getOption("device")()
		  cim(Xcor, cluster = "none",
		      title = if(title) "X correlation matrix" else NULL,
		      row.names = X.var.names, col.names = X.var.names)
		  
		  if (isTRUE(GD)) getOption("device")()
		  cim(Ycor, cluster = "none", 
		      title = if(title) "Y correlation matrix" else NULL,
		      row.names = Y.var.names, col.names = Y.var.names)
		  
		  if (isTRUE(GD)) getOption("device")()
		  cim(XYcor, cluster = "none",
		      title = if(title) "XY correlation matrix" else NULL,
		      row.names = X.var.names, col.names = Y.var.names)
		}
    }	
}

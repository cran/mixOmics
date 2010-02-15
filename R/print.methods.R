# ------------- print pls -------------------------------------

'print.pls' <-
function(x, ...){

    mode = paste("'", x$mode, "'", sep = "")
	
    cat("\nCall:\n", deparse(x$call), "\n\n")
	
    cat(" PLS with a", mode, "mode with", x$ncomp, "PLS components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")

    cat(" No variable selection \n\n")
	
    cat(" Available components: \n", 
        "-------------------- \n")
	
	cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}


# ----------------print spls ----------------------------------

'print.spls' <-
function(x, ...){

    mode = paste("'", x$mode, "'", sep = "")
    keepX = paste("[", x$keepX, "]", sep = "")
    keepY = paste("[", x$keepY, "]", sep = "")
	
	cat("\nCall:\n", deparse(x$call), "\n\n")
	
    cat(" sPLS with a", mode, "mode with", x$ncomp, "sPLS components. \n")
    cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")
	
    cat(" Selection of", keepX, "variables on each of the sPLS components on the X data set. \n")
    cat(" Selection of", keepY, "variables on each of the sPLS components on the Y data set. \n\n")

    cat(" Available components: \n", 
        "-------------------- \n")
	
	cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}


# ---------------- print rcc ----------------------------------

'print.rcc' <-
function(x, ...){

    cat("\nCall:\n", deparse(x$call), "\n\n")
  
    cat(" rCCA with regularization parameters", x$lambda[1], "and", x$lambda[2], "for the X and Y data. \n")
    cat(" You entered data X of dimensions :", nrow(x$X), ncol(x$X), "\n")
    cat(" You entered data Y of dimensions :", nrow(x$Y), ncol(x$Y), "\n\n")
	
    cat(" Available components: \n", 
        "-------------------- \n")
	
	cat(" canonical correlations: see object$cor \n")
	cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
}


# ------- print for summary with (s)PLS object or rcc ----------

'print.summary' <-
function(x, ...){

    print.gap = 4
    what = x$what
    digits = x$digits

	# --------------------- output pls/spls ---------------------------
	
    if(x$method == "pls" | x$method == "spls"){

        if (x$method == "pls") {
            cat(" PLS mode:", x$mode)
			cat("\n Number of variates considered:", x$ncomp, "\n")
        }	
        else {
            cat(" sPLS mode:", x$mode)
			cat("\n Number of variates considered:", x$ncomp)
			cat("\n Number of X-variables selected on each of the sPLS components:", x$keepX)
			cat("\n Number of Y-variables selected on each of the sPLS components:", x$keepY, "\n")
        }			
        

        #---------- affichage communauté ----------#
        if (any(what == "all") || any(what == "communalities")) { 
            cat("\n\n Communalities Analysis:\n",
                "----------------------")
				
            cat("\n X-Variables vs their own Variates: see object$CM.X$own \n")
            cat("\n X-Variables vs the opposite Variates: see object$CM.X$opp \n")
            cat("\n Y-Variables vs their own Variates: see object$CM.Y$opp \n")
            cat("\n Y-Variables vs the opposite Variates: see object$CM.Y$opp \n")
        }

        #--------- affichage redondance -----------#
        if (any(what == "all") || any(what == "redundancy")) {
            cat("\n\n Redundancy Analysis:\n",
                "-------------------\n")
				
            cat("\n X-Variables vs their own Variates: see object$Rd.X$own \n")
            cat("\n X-Variables vs the opposite Variates: see object$Rd.X$opp \n")
            cat("\n Y-Variables vs their own Variates: see object$Rd.Y$opp \n")
            cat("\n Y-Variables vs the opposite Variates: see object$Rd.Y$opp \n")
        }

        #---------- tableau VIP ---------#
        if (any(what == "all") || any(what == "VIP")) {
            cat("\n\n", "Variable Importance in the Projection (VIP): see object$VIP \n",
                        "------------------------------------------- \n\n")
        }

    }  #end if pls


    # --------------------------- output rcc ---------------------------
    if(x$method == "rcc" ){
        print.gap = 4
        if (any(what == "all")) {
            cat(" Number of canonical variates considered:", x$ncomp, "\n")
            cat("\n Canonical correlations:",
                "\n ----------------------\n")
            print(round(x$can.cor, digits = digits), print.gap = print.gap)
        }

        #-- affichage communauté --#
        if (any(what == "all") || any(what == "communalities")) { 
            cat("\n\n Canonical Communalities Analysis:\n",
                "--------------------------------")

            cat("\n X-Variables vs their own Canonical Variates: see object$Cm.X$own \n")
            cat("\n X-Variables vs the opposite Canonical Variates: see object$Cm.X$opp \n")
            cat("\n Y-Variables vs their own Canonical Variates: see object$Cm.Y$own \n")
            cat("\n Y-Variables vs the opposite Canonical Variates: see object$Cm.Y$opp \n")
        }

        #--------- affichage redondance -----------#
        if (any(what == "all") || any(what == "redundancy")) {
            cat("\n\n Redundancy Analysis:\n",
                "-------------------\n")
				
            cat("\n X-Variables vs their own Variates: see object$Rd.X$own \n")
            cat("\n X-Variables vs the opposite Variates: see object$Rd.X$opp \n")
            cat("\n Y-Variables vs their own Variates: see object$Rd.Y$opp \n")
            cat("\n Y-Variables vs the opposite Variates: see object$Rd.Y$opp \n")
        }

    }  #end rcc
}







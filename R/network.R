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


network <-
function(...) UseMethod("network")


network.default <-
function(mat, 
         threshold = 0.5, 
         X.names = NULL, 
         Y.names = NULL,
         color.node = c("white", "white"), 
         shape.node = c("circle", "rectangle"), 
         color.edge = c("blue", "red"), 
         lty.edge = c("solid", "solid"), 
         lwd.edge = c(1, 1), 
         show.edge.labels = FALSE,
         show.color.key = TRUE, 
         symkey = TRUE, 
         keysize = 1,
         breaks, 
         interactive = FALSE, 
         alpha = 1,
	 ...) 
{

    # validation des arguments #
	#--------------------------#
	if (length(dim(mat)) != 2) 
        stop("'mat' must be a numeric matrix.")

    mat = as.matrix(mat)

    if (!is.numeric(mat)) 
        stop("'mat' must be a numeric matrix.")

    if (length(color.node) != 2) 
		stop("'color.node' must be a vector of length 2.")

	if (length(shape.node) != 2) 
		stop("'shape.node' must be a vector of length 2.")

	if (length(color.edge) < 2 & (class(color.edge) != "function")) 
		stop("'color.edge' must be a vector of length larger than or equal to 2.")
		
	if (length(lty.edge) != 2) 
		stop("'lty.edge' must be a vector of length .")

	if (length(lwd.edge) != 2) 
		stop("'lwd.edge' must be a vector of length 2.")
		
	p = nrow(mat)
	q = ncol(mat)
    dim = min(p, q)	

	if (is.null(X.names)) X.names = rownames(mat)
	if (is.null(Y.names)) Y.names = colnames(mat)
	
	if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")
	if (is.null(Y.names)) Y.names = paste("Y", 1:q, sep = "")
	
	if (is.null(threshold) || !is.numeric(threshold) || threshold < 0 || threshold > max(abs(mat)))
        stop("Invalid value for 'threshold', it must be a positive numeric value < ", max(abs(mat)))
		
    id = bin.color(t(mat), threshold = threshold, breaks = breaks, col = color.edge, symkey = symkey)
    mat = as.vector(t(mat))
    col.id = as.vector(id$bin)	
	color.edge = id$col[col.id[1]]

	for (i in 2:length(mat)) {
		color.edge = c(color.edge, id$col[col.id[i]])
	}

	# Définition des sommets #
	#------------------------#
    Xn = paste("X", 1:p, sep = "")
    Yn = paste("Y", 1:q, sep = "")
    nodes = data.frame(name = c(Xn, Yn), group = c(rep("x", p), rep("y", q)))
	
	node.X = rep(Xn, each = q)
	node.Y = rep(Yn, p)

	# Définition des arêtes #
	#-----------------------#
	relations = data.frame(from = node.X, to = node.Y, weight = mat)

	# Décide quels sont les arêtes à incluir dans le réseau #
	#-------------------------------------------------------#
	idx = abs(mat) >= threshold
	relations = relations[idx, ]
	color.edge = color.edge[idx]

	# Génère un graphe avec toutes les arêtes signifiantes #
	#------------------------------------------------------#
	gR = graph.data.frame(relations, directed = FALSE, vertices = nodes)
	
	# Attributs des sommets #
	#-----------------------#
	V(gR)$label = c(X.names, Y.names)
	
	V(gR)$label.color = "black"
	
	V(gR)$color = color.node[1]
	V(gR)$color[V(gR)$group == "y"] = color.node[2]

	V(gR)$shape = shape.node[1]
	V(gR)$shape[V(gR)$group == "y"] = shape.node[2]
	
	# Attributs des arêtes #
	#----------------------#
	if (show.edge.labels) E(gR)$label = round(E(gR)$weight, 2)
	
	E(gR)$label.color = "black"
	
	E(gR)$color = color.edge 
	
	E(gR)$lty = lty.edge[1]
	E(gR)$lty[E(gR)$weight < 0] = lty.edge[2]
	
	E(gR)$width = lwd.edge[1]
	E(gR)$width[E(gR)$weight < 0] = lwd.edge[2]
	
	gR = delete.vertices(gR, which(degree(gR) == 0) - 1)
	
	# Attributs pour le plot #
	#------------------------#
	lhei = c(keysize, 4) 
	lwid = c(keysize, 4)    
	lmat = matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)	
	nc = length(id$col)
	x = seq(0, 1, length = nc + 2)
	z.mat = seq(0, 1, length = nc + 1)
	z.mat = matrix(z.mat, ncol = 1)
	
	if ((id$lim[1] < -threshold) & (id$lim[2] < threshold)) {
		xv = c(0, x[nc + 1])
		lv = round(c(id$lim[1], -threshold), 2)
		col = c(id$col, "white")
	}
	
	if ((id$lim[1] > -threshold) & (id$lim[2] > threshold)) {
		xv = c(x[2], 1)
		lv = round(c(threshold, id$lim[2]), 2)
		col = c("white", id$col)
	}
				
	if ((id$lim[1] < -threshold) & (id$lim[2] > threshold)) {
		idn = max(which(id$breaks < 0))	
		idp = min(which(id$breaks > 0))	
		xv = c(0, x[idn + 1], x[idp], 1)	
		lv = round(c(id$lim[1], -threshold, threshold, id$lim[2]), 2)
		col = c(id$col[1:idn], "white", id$col[(idn + 1):nc])
	}
	
	#----------------------------------#
	# Construction du graphe de départ #
	#----------------------------------#
	nn = vcount(gR)
	V(gR)$label.cex = min(2/log(nn), 1)
	E(gR)$label.cex = min(2.25/log(nn), 1)
	cex0 = 2*V(gR)$label.cex
	
	def.par = par(no.readonly = TRUE)
	
	par(pty = "s", mar = c(0, 0, 0, 0))
	plot(1:100, 1:100, type = "n", axes = FALSE, xlab = "", ylab = "")
	cha = V(gR)$label
	cha = paste(" ", cha, " ")
	xh = strwidth(cha, cex = cex0)
	yh = strheight(cha, cex = cex0) * 5/2.75

	V(gR)$size = xh
	V(gR)$size2 = yh
	
	par(def.par)
	
	l = layout.fruchterman.reingold(gR, weights = (1 - abs(E(gR)$weight))^alpha)

	if (isTRUE(!interactive)) {
		if (isTRUE(show.color.key)) {
			layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
			par(mar = c(5, 4, 2, 1), cex = 0.75)
			image(z.mat, col = col, xaxt = "n", yaxt = "n")
			box()
			par(usr = c(0, 1, 0, 1))						
			axis(1, at = xv, labels = lv)
			title("Color key", font.main = 1)
			par(def.par)
			par(new = TRUE)
		}
	
		par(pty = "s", mar = c(0, 0, 0, 0))
		plot(gR, layout = l)
		par(def.par)
	}	
	
	#----------------------#
	# Procédure interactif #
	#----------------------#
	gE.none = FALSE
	if (isTRUE(interactive)) {
        dev.off()
		
		# Barre de contrôle #
		#-------------------#
		min.cut = threshold
		max.cut = max(mat)
		threshold.old = threshold

		getOption("device")("width" = 4, "height" = 2, "xpos" = 250, "ypos" = 0) 
		def.par = par(no.readonly = TRUE)

		cuts = seq(0, 1, length = 21)
		par(mai = c(0.25, 0.15, 0.3, 0.15), bg = gray(0.95))
		layout(matrix(c(0, 1, 0), ncol = 1, nrow = 3), 
		widths = 1, heights = c(0.25, 1, 0.25))

		plot(cuts, type = "n", rep(0, 21), xlab = "", ylab = "",
		xlim = c(-0.10, 1.10), axes = FALSE)
		title("threshold control", cex.main = 1.9, font.main = 1)
		text(0.5, -0.6, "value", cex = 1.5)
		text(0, -0.6, round(min.cut, 2), cex = 1.4)
		text(1, -0.6, round(max.cut, 2), cex = 1.4)
		mtext(min.cut, side = 1, line = -0.4, outer = FALSE, cex = 0.95)

		rect(-0.1, -0.3, -0.02, 0.3, col = "white") 
		rect(1.02, -0.3, 1.1, 0.3, col = "white")
		points(1.06, 0, pch = 3, cex = 2.4)
		lines(c(-0.085, -0.035), c(0, 0))

		for (i in seq(0, 1, length = 21)) lines(c(i, i), c(-0.22, 0.2))

		x = pos = 0
		rect(-0.01, -0.045, x, 0.04, col = "red")
		rect(x, -0.045, 1.01, 0.04, col = "white")

		getOption("device")()

		# Plot du graphe de départ #
		#--------------------------#
		if (isTRUE(show.color.key)) {
			layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
			par(mar = c(5, 4, 2, 1), cex = 0.75)
			image(z.mat, col = col, xaxt = "n", yaxt = "n")
			box()
			par(usr = c(0, 1, 0, 1))						
			axis(1, at = xv, labels = lv)
			title("Color key", font.main = 1)
			par(def.par)	
			par(new = TRUE)
		}
		
		par(pty = "s", mar = c(0, 0, 0, 0))	
		plot(gR, layout = l)
		par(def.par)
		
		nD = dev.cur()
		dev.set(dev.prev())
		DEV.x = grDevices::dev.cur()

		repeat {
			grDevices::dev.set(DEV.x)

			z = .Internal(locator(1, type = "n"))
			x = z[[1]]
			y = z[[2]]
			flag = z[[3]]

			if (flag == 0) break

			if (0 <= x & x <= 1 & -0.22 <= y & y <= 0.22) {
				rect(0, -0.045, x, 0.04, col = "red")
				rect(x, -0.045, 1.01, 0.04, col = "white")
				pos = x
			}

			if (1.02 <= x & x <= 1.1 & -0.3 <= y & y <= 0.3) {
				x = pos + 0.05
				idx = which.min(abs(cuts - x))
				x = cuts[idx]
				pos = x
				rect(0, -0.045, x, 0.04, col = "red")
				rect(x, -0.045, 1.01, 0.04, col = "white")
			}

			if (-0.1 <= x & x <= -0.02 & -0.3 <= y & y <= 0.3) {
				x = pos - 0.05
				idx = which.min(abs(cuts - x))
				x = cuts[idx]
				pos = x
				rect(0, -0.045, x, 0.04, col = "red")
				rect(x, -0.045, 1.01, 0.04, col = "white")
			}

			mtext(round(threshold, 3), side = 1, line = -0.4, cex = 0.9, 
			col = gray(0.95), font = 2)
			threshold = (max.cut - min.cut) * pos + min.cut
			mtext(round(threshold, 3), side = 1, line = -0.4, cex = 0.9)

			grDevices::dev.set(nD) 

			# Plot du nouveau graphe #
			#------------------------#
			if (threshold >= threshold.old) {

				# Décide quels sont les arêtes à supprimer du réseau #
				#----------------------------------------------------#
				supp.edge = E(gR)[abs(E(gR)$weight) < threshold]

				# Génère un graphe avec toutes les arêtes signifiantes #
				#------------------------------------------------------#
				gE = delete.edges(gR, supp.edge)
				gE = delete.vertices(gE, which(degree(gE) == 0) - 1)
				
				# Plot du graphe #
				#----------------#
				nn = vcount(gE)
				V(gE)$label.cex = min(2/log(nn), 1)
				E(gE)$label.cex = min(2.25/log(nn), 1)
				cex0 = 2*V(gE)$label.cex
				
				def.par = par(no.readonly = TRUE)
				
				par(pty = "s", mar = c(0, 0, 0, 0))
				plot(1:100, 1:100, type = "n", xaxt = "n")
				cha = V(gE)$label
				cha = paste(" ", cha, " ")
				xh = strwidth(cha, cex = cex0)
				yh = strheight(cha, cex = cex0) * 5/2.75
	
				V(gE)$size = xh
				V(gE)$size2 = yh	
				
				par(def.par)
	
				l = layout.fruchterman.reingold(gE, weights = (1 - abs(E(gE)$weight))^alpha)

				if (isTRUE(show.color.key)) {
					layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
					par(mar = c(5, 4, 2, 1), cex = 0.75)
					image(z.mat, col = col, xaxt = "n", yaxt = "n")
					box()
					par(usr = c(0, 1, 0, 1))						
					axis(1, at = xv, labels = lv)
					title("Color key", font.main = 1)
					par(def.par)
					par(new = TRUE)
				}
				
				par(pty = "s", mar = c(0, 0, 0, 0))
				plot(gE, layout = l)
				par(def.par)
				
				threshold.old = threshold
			}
			else {
				#-------------------------------------------------------#
				# Décide quels sont les arêtes à incluir dans le réseau #
				#-------------------------------------------------------#
				supp.edge = E(gR)[abs(E(gR)$weight) < threshold]

				# Génère un graphe avec toutes les arêtes signifiantes #
				#------------------------------------------------------#
				gE = delete.edges(gR, supp.edge)
				gE = delete.vertices(gE, which(degree(gE) == 0) - 1)

				# Plot du graphe #
				#----------------#
				nn = vcount(gE)
				V(gE)$label.cex = min(2/log(nn), 1)
				E(gE)$label.cex = min(2.25/log(nn), 1)
				cex0 = 2*V(gE)$label.cex
				
				def.par = par(no.readonly = TRUE)
				
				par(pty = "s", mar = c(0, 0, 0, 0))
				plot(1:100, 1:100, type = "n", xaxt = "n")
				cha = V(gE)$label
				cha = paste(" ", cha, " ")
				xh = strwidth(cha, cex = cex0)
				yh = strheight(cha, cex = cex0) * 5/2.75
	
				V(gE)$size = xh
				V(gE)$size2 = yh

				par(def.par)	
	
				l = layout.fruchterman.reingold(gE, weights = (1 - abs(E(gE)$weight))^alpha)

				if (isTRUE(show.color.key)) {
					layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
					par(mar = c(5, 4, 2, 1), cex = 0.75)
					image(z.mat, col = col, xaxt = "n", yaxt = "n")
					box()
					par(usr = c(0, 1, 0, 1))						
					axis(1, at = xv, labels = lv)
					title("Color key", font.main = 1)
					par(def.par)
					par(new = TRUE)
				}
				
				par(pty = "s", mar = c(0, 0, 0, 0))
				plot(gE, layout = l)
				par(def.par)
				
				threshold.old = threshold
			}

		grDevices::dev.set(DEV.x)
		gE.none = TRUE
		} # fin du bucle
		
		if (gE.none != FALSE) gR = gE
	}
	
	return(invisible(gR))
}	



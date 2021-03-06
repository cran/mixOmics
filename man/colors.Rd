\name{color.jet}
\encoding{UTF-8}
\alias{color.jet}
\alias{color.spectral}
\alias{color.GreenRed}
\alias{color.mixo}

\title{Color Palette for mixOmics}

\description{
The functions create a vector of \code{n} "contiguous" colors (except the \code{color.mixo} which are colors used internally to fit our logo colors). 
}

\usage{
color.jet(n, alpha = 1)
color.spectral(n, alpha = 1)
color.GreenRed(n, alpha = 1)
color.mixo(num.vector)
}

\arguments{
  \item{n}{an integer, the number of colors \eqn{(\geq 1)} to be in the palette.}
  \item{alpha}{a numeric value between 0 and 1 for alpha channel (opacity).}
  \item{num.vector}{for \code{color.mixo} an integer vector specifying which colors to use in the mixOmics palette (there are only 10 colors available.}
}

\details{
The function \code{color.jet(n)} create color scheme, beginning with dark blue, 
ranging through shades of blue, cyan, green, yellow and red, 
and ending with dark red. This colors palette is suitable for displaying ordered (symmetric) data, 
with \code{n} giving the number of colors desired.
}

\value{
For \code{color.jet(n)}, \code{color.spectral(n)}, \code{color.GreenRed(n)} a character vector, \code{cv}, of color names. This can be used either to create a user-defined color palette for subsequent graphics by \code{palette(cv)}, a \code{col=} specification in graphics functions or in \code{par}.

For \code{color.mixo}, a vector of colors matching the mixOmics logo (10 colors max.)
}

\seealso{
\code{\link{colorRamp}}, \code{\link{palette}}, \code{\link{colors}} 
for the vector of built-in "named" colors; 
\code{\link{hsv}}, \code{\link{gray}}, \code{\link{rainbow}}, 
\code{\link{terrain.colors}}, ... to construct colors; and \code{\link{heat.colors}}, 
\code{\link{topo.colors}} for images. 
}

\examples{
# -----------------------
# jet colors
# ----------------------
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
  image(matrix(z, ncol = 1), col = color.jet(n), 
        xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
  box()
  par(usr = c(-1, 1, -1, 1))    
  axis(1, at = c(-1, 0, 1))
}

# -----------------------
# spectral colors
# ----------------------
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
  image(matrix(z, ncol = 1), col = color.spectral(n), 
        xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
  box()
  par(usr = c(-1, 1, -1, 1))    
  axis(1, at = c(-1, 0, 1))
}

# -----------------------
# GreenRed colors
# ----------------------
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
  image(matrix(z, ncol = 1), col = color.GreenRed(n), 
        xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
  box()
  par(usr = c(-1, 1, -1, 1))    
  axis(1, at = c(-1, 0, 1))
}

# # --------------------------------
# mixOmics colors
# # -------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

my.colors = color.mixo(1:5)
my.pch = ifelse(nutrimouse$genotype == 'wt', 16, 17)
#plotIndiv(nutri.res, ind.names = FALSE, group = my.colors, pch = my.pch, cex = 1.5)

}


\keyword{color}



\name{jet.colors}
\encoding{latin1}
\alias{jet.colors}

\title{Jet Colors Palette}

\description{
Create a vector of \code{n} "contiguous" colors. 
}

\usage{
jet.colors(n)
}

\arguments{
  \item{n}{an integer, the number of colors \eqn{(\geq 1)} to be in the palette.}
}

\details{
The function \code{jet.colors(n)} create color scheme, beginning with dark blue, 
ranging through shades of blue, cyan, green, yellow and red, 
and ending with dark red. This colors palette is suitable for displaying ordered (symmetric) data, 
with \code{n} giving the number of colors desired.
}

\value{
A character vector, \code{cv}, of color names. This can be used either to create a 
user-defined color palette for subsequent graphics by \code{palette(cv)}, a 
\code{col=} specification in graphics functions or in \code{par}.
}

\seealso{\code{\link{colorRamp}}, \code{\link{palette}}, \code{\link{colors}} 
for the vector of built-in "named" colors; 
\code{\link{hsv}}, \code{\link{gray}}, \code{\link{rainbow}}, 
\code{\link{terrain.colors}}, \eqn{\ldots} to construct colors; and \code{\link{heat.colors}}, 
\code{\link{topo.colors}} for images. 
}

\examples{
## Jet Color Scales Strips
def.par <- par(no.readonly = TRUE)
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
    image(matrix(z, ncol = 1), col = jet.colors(n), 
          xaxt = "n", yaxt = "n", main = paste("n = ", n))
    box()
    par(usr = c(-1, 1, -1, 1))		
    axis(1, at = c(-1, 0, 1))
}
par(def.par)
}

\keyword{color}
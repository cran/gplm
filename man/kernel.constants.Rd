\name{kernel.constants}
\alias{kernel.constants}
\title{Kernel constants}
\description{
Calculates several constants of a (product) kernel function.
}
\usage{
kernel.constants(kernel = "biweight", d = 1, product = TRUE)
}
\arguments{
  \item{kernel}{text string, see \code{\link{kernel.function}}}
  \item{d}{integer (dimension of the kernel)}
  \item{product}{(if d>1) product or spherical kernel}
}
\details{
  The constants which are calculated are the second moment,
  the square norm and the canonical bandwidth of the kernel
  (only the two latter terms depend on the dimension d).
}
\value{
  List with components:
  \item{m2}{second moment}
  \item{c2}{square norm}
  \item{d0}{canonical bandwidth}
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{kernel.function}}}
\examples{
  kernel.constants()                  ## default (biweight), d=1
  kernel.constants("epanechnikov",1)  ## epanechnikov, d=1
  kernel.constants("epanechnikov",2)  ## product epanechnikov, d=2
}
\keyword{smooth}


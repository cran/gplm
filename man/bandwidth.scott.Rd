\name{bandwidth.scott}
\alias{bandwidth.scott}
\title{Scott's rule of thumb}
\description{
Calculates Scott's rule of thumb bandwidth vector.
}
\usage{
bandwidth.scott(x, kernel = "biweight", product = TRUE)
}
\arguments{
  \item{x}{n x d matrix, data}
  \item{kernel}{text string, see \code{\link{kernel.function}}}
  \item{product}{(if d>1) product or spherical kernel}
}
\details{
  The default
  bandwidth vector is computed by Scott's rule of thumb for the Gaussian
  kernel and adapted to the chosen kernel function.
}
\value{
  d x 1 bandwidth vector used for calculation
}
\references{
  Scott, D.W. (1992).
  \emph{Multivariate Density Estimation: Theory, Practice, and
    Visualization}. New York, Chichester: Wiley.
}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{kernel.function}}, \code{\link{kde}}}
\examples{
  ## two-dimensional data
  n <- 1000
  u <- runif(n)
  thresh <- 0.4
  x1 <- rnorm(n)*(u<thresh) +rnorm(n,mean=3)*(u>=thresh)
  x2 <- rnorm(n)*(u<thresh) +rnorm(n,mean=9)*(u>=thresh)
  bandwidth.scott( cbind(x1,x2) )
}
\keyword{smooth}


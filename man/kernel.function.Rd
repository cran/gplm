\name{kernel.function}
\alias{kernel.function}
\title{Kernel function}
\description{
Calculates several kernel functions (uniform, triangle, epanechnikov,
biweight, triweight, gaussian).
}
\usage{
kernel.function(u, kernel = "biweight", product = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{n x d matrix}
  \item{kernel}{text string}
  \item{product}{(if d>1) product or spherical kernel}
}
\details{
  The kernel parameter is  a text string specifying
  the univariate kernel function which is either the gaussian pdf
  or proportional to \eqn{(1-|u|^p)^q}.
  Possible text strings are "triangle" (p=q=1),
  "uniform" (p=1, q=0), "epanechnikov" (p=2, q=1),
  "biweight" or "quartic" (p=q=2),
  "triweight" (p=2, q=3), "gaussian" or "normal" (gaussian pdf).

  The multivariate kernels are obtained by a
  product of unvariate kernels \eqn{K(u_1)...K(u_d)}
  or by a spherical (radially symmetric) kernel
  proportional to \eqn{K(||u||)}. (The resulting kernel
  is a density, i.e. integrates to 1.)
}
\value{
  n x 1 vector of kernel weights
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
%%\seealso{}
\examples{
  kernel.function(0)                         ## default (biweight)
  kernel.function(0, kernel="epanechnikov")  ## epanechnikov
  kernel.function(0, kernel="gaussian")      ## equals dnorm(0)
}
\keyword{smooth}

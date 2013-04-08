\name{convol}
\alias{convol}
\title{Kernel convolution}
\description{
Calculates the convolution of data with a kernel function.
}
\usage{
convol(x, h = 1, grid = NULL, y = 1, w = 1, p = 2, q = 2,
       product = TRUE, sort = TRUE)
}
\arguments{
  \item{x}{n x d matrix, data}
  \item{h}{scalar or 1 x d, bandwidth(s)}
  \item{grid}{m x d matrix, where to calculate the convolution (default = x)}
  \item{y}{n x c matrix, optional responses}
  \item{w}{scalar or  n x 1 or 1 x m or n x m, optional weights}
  \item{p}{integer or text, see \code{\link{kernel.function}}}
  \item{q}{integer, see \code{\link{kernel.function}}}
  \item{product}{(if d>1) product or spherical kernel}
  \item{sort}{logical, TRUE if data need to be sorted}
}
\details{
  The kernel convolution which is calculated is
  \eqn{\sum_i K_h(x_i - grid_{j})\,y_i\,w_{ij}} for
  \eqn{i=1,...,n} and \eqn{j=1,...,m}. The kernel function is determined
  by the kernel parameters p and q, see
  \code{\link{kernel.function}}. The default kernel is the biweight
  (quartic) kernel function.   Note that the DLL requires the data matrix
  to be sorted by its first column.
}
\value{
  m x c matrix
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{kernel.function}}, \code{\link{kde}}, \code{\link{kreg}}}
\examples{
  n <- 100
  x <- rnorm(n)
  convol(x,h=0.8,grid=-3:3)/n  ## estimates density of x at points -3:3
}
\keyword{smooth}

\name{kde}
\alias{kde}
\title{Kernel density estimation}
\description{
  Calculates a kernel density estimate (univariate or multivariate).
}
\usage{
kde(x, bandwidth = NULL, grid = TRUE, kernel = "biweight",
    product = TRUE, sort = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{n x d matrix, data}
  \item{bandwidth}{scalar or 1 x d, bandwidth(s)}
  \item{grid}{logical or m x d matrix (where to calculate the density)}
  \item{kernel}{text string, see \code{\link{kernel.function}}}
  \item{product}{(if d>1) product or spherical kernel}
  \item{sort}{logical, TRUE if data need to be sorted}
}
\details{
  The kernel density estimator is calculated as
  \eqn{\frac{1}{n} \sum_i K_h(x_i - grid_{j})} for
  \eqn{i=1,...,n} and \eqn{j=1,...,m}. The default
  bandwidth is computed by Scott's rule of thumb
  (adapted to the chosen kernel function).

}
\value{
  List with components:
  \item{x}{m x d matrix, where density has been calculated}
  \item{y}{m x 1 vector, density estimates}
  \item{bandwidth}{bandwidth used for calculation}
  \item{rearrange}{if sort=TRUE, index to rearrange x and y to its original order.}
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{kernel.function}}, \code{\link{convol}}, \code{\link{kreg}}}
\examples{
  n <- 1000
  x <- rnorm(n)
  plot(kde(x), type="l")

  ## mixed normal data
  n <- 1000
  u <- runif(n)
  thresh <- 0.4
  x <- rnorm(n)*(u<thresh) +rnorm(n,mean=3)*(u>=thresh)
  h <- 1
  fh <- kde(x,bandwidth=h)
  plot(kde(x,bandwidth=h),type="l",lwd=2); rug(x)
  lines(kde(x,bandwidth=h*1.2),col="red")
  lines(kde(x,bandwidth=h*1.4),col="orange")
  lines(kde(x,bandwidth=h/1.2),col="blue")
  lines(kde(x,bandwidth=h/1.4),col="cyan")

  ## two-dimensional data
  n <- 1000
  u <- runif(n)
  thresh <- 0.4
  x1 <- rnorm(n)*(u<thresh) +rnorm(n,mean=3)*(u>=thresh)
  x2 <- rnorm(n)*(u<thresh) +rnorm(n,mean=9)*(u>=thresh)

  grid1 <- seq(min(x1),max(x1),length=20)  ## grid for x1
  grid2 <- seq(min(x2),max(x2),length=25)  ## grid for x2

  fh <- kde( cbind(x1,x2), grid=create.grid(list(grid1,grid2)) )
  o <- order(fh$x[,2],fh$x[,1])
  density <- (matrix(fh$y[o],length(grid1),length(grid2)))
  
  par(mfrow=c(2,2))
  plot(kde(x1),type="l",main="x1"); rug(x1)
  plot(kde(x2),type="l",main="x2"); rug(x2)
  persp(grid1,grid2,density,main="KDE",
        theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
  contour(grid1,grid2,density, main="KDE Contours")
  points(x1,x2,col="red",pch=18,cex=0.5)
  par(mfrow=c(1,1))
}
\keyword{smooth}

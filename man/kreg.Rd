\name{kreg}
\alias{kreg}
\title{Kernel regression}
\description{
  Calculates a kernel regression estimate (univariate or multivariate).
}
\usage{
kreg(x, y, bandwidth = NULL, grid = TRUE, kernel = "biweight",
     product = TRUE, sort = TRUE)
}
\arguments{
  \item{x}{n x d matrix, data}
  \item{y}{n x 1 vector, responses}
  \item{bandwidth}{scalar or 1 x d, bandwidth(s)}
  \item{grid}{logical or m x d matrix (where to calculate the regression)}
  \item{kernel}{text string, see \code{\link{kernel.function}}}
  \item{product}{(if d>1) product or spherical kernel}
  \item{sort}{logical, TRUE if data need to be sorted}
}
\details{
  The estimator is calculated by Nadaraya-Watson kernel regression.
  Future extension to local linear (d>1) or polynomial (d=1) estimates
  is planned. The default bandwidth is computed by Scott's rule of thumb
  for kde (adapted to the chosen kernel function).
}
\value{
  List with components:
  \item{x}{m x d matrix, where regression has been calculated}
  \item{y}{m x 1 vector, regression estimates}
  \item{bandwidth}{bandwidth used for calculation}
  \item{df.residual}{approximate degrees of freedom (residuals)}
  \item{rearrange}{if sort=TRUE, index to rearrange x and y to its original order.}
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{kernel.function}}, \code{\link{convol}}, \code{\link{kde}}}
\examples{
  n <- 1000
  x <- rnorm(n)
  m <- sin(x)
  y <- m + rnorm(n)
  plot(x,y,col="gray")
  o <- order(x); lines(x[o],m[o],col="green")
  lines(kreg(x,y),lwd=2)

  ## two-dimensional
  n <- 100
  x <- 6*cbind(runif(n), runif(n))-3
  m <- function(x1,x2){ 4*sin(x1) + x2 }
  y <- m(x[,1],x[,2]) + rnorm(n)
  mh <- kreg(x,y)##,bandwidth=1)

  grid1 <- unique(mh$x[,1])
  grid2 <- unique(mh$x[,2])
  est.m  <- t(matrix(mh$y,length(grid1),length(grid2)))
  orig.m <- outer(grid1,grid2,m)
  par(mfrow=c(1,2))
  persp(grid1,grid2,orig.m,main="Original Function",
        theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
  persp(grid1,grid2,est.m,main="Estimated Function",
	theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
  par(mfrow=c(1,1))
  
  ## now with normal x, note the boundary problem,
  ## which can be somewhat reduced by a gaussian kernel
  n <- 1000
  x <- cbind(rnorm(n), rnorm(n))
  m <- function(x1,x2){ 4*sin(x1) + x2 }
  y <- m(x[,1],x[,2]) + rnorm(n)
  mh <- kreg(x,y)##,p="gaussian")

  grid1 <- unique(mh$x[,1])
  grid2 <- unique(mh$x[,2])
  est.m  <- t(matrix(mh$y,length(grid1),length(grid2)))
  orig.m <- outer(grid1,grid2,m)
  par(mfrow=c(1,2))
  persp(grid1,grid2,orig.m,main="Original Function",
        theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
  persp(grid1,grid2,est.m,main="Estimated Function",
	theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
  par(mfrow=c(1,1))
}
\keyword{smooth}

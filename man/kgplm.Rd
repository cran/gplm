\name{kgplm}
\alias{kgplm}
\title{Generalized partial linear model}
\description{
  Fits a generalized partial linear model (kernel-based)
  using the (generalized) Speckman estimator or backfitting (in the
  generalized case combined with local scoring) for two additive
  component functions. 
}
\usage{
kgplm(x, t, y, h, family, link,
          b.start=NULL, m.start=NULL, grid = NULL, m.grid.start = NULL, 
          offset = 0, method = "speckman", sort = TRUE, weights = 1,
          weights.trim = 1, weights.conv = 1, max.iter = 25, eps.conv = 1e-8,
          kernel = "biweight", kernel.product = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{n x p matrix, data for linear part}
  \item{y}{n x 1 vector, responses}
  \item{t}{n x q matrix, data for nonparametric part}
  \item{h}{scalar or 1 x q, bandwidth(s)}
  \item{family}{text string, family of distributions (e.g.
    "gaussian" or "bernoulli", see details for \code{\link{glm.ll}})}
  \item{link}{text string, link function (depending on family,
    see details for \code{\link{glm.ll}})}
  \item{b.start}{p x 1 vector, start values for linear part}
  \item{m.start}{n x 1 vector, start values for nonparametric part}
  \item{grid}{m x q matrix, where to calculate the nonparametric function (default = t)}
  \item{m.grid.start}{m x 1 vector, start values for nonparametric part on grid}
  \item{offset}{offset}
  \item{method}{"speckman" or "backfit"}
  \item{sort}{logical, TRUE if data need to be sorted}
  \item{weights}{binomial weights}
  \item{weights.trim}{trimming weights for fitting the linear part}
  \item{weights.conv}{weights for convergence criterion}
  \item{max.iter}{maximal number of iterations}
  \item{eps.conv}{convergence criterion}
  \item{kernel}{text string, see \code{\link{kernel.function}}}
  \item{kernel.product}{(if p>1) product or spherical kernel}
  \item{verbose}{print additional convergence information}
}
%%\details{
%%}
\value{
  List with components:
  \item{b}{p x 1 vector, linear coefficients}
  \item{b.cov}{p x p matrix, linear coefficients}
  \item{m}{n x 1 vector, nonparametric function estimate}
  \item{m.grid}{m x 1 vector, nonparametric function estimate on grid}
  \item{it}{number of iterations}
  \item{deviance}{deviance}
  \item{df.residual}{approximate degrees of freedom (residuals)}
  \item{aic}{Akaike's information criterion}
}
\references{
  Mueller, M. (2001)
  Estimation and testing in generalized partial linear models -- A
  comparative study. \emph{Statistics and Computing,} 11:299--309.

  Hastie, T. and Tibshirani, R. (1990)
  \emph{Generalized Additive Models.} London: Chapman and Hall.
}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{kernel.function}}, \code{\link{kreg}}}
\examples{
  ## data 
  n <- 1000; b <- c(1,-1); rho <- 0.7
  m <- function(t){ 1.5*sin(pi*t) }
  x1 <- runif(n,min=-1,max=1); u  <- runif(n,min=-1,max=1)
  t  <- runif(n,min=-1,max=1); x2 <- round(m(rho*t + (1-rho)*u))
  x  <- cbind(x1,x2)
  y  <- x \%*\% b + m(t) + rnorm(n)

  ## partial linear model (PLM)
  gh <- kgplm(x,t,y,h=0.25,family="gaussian",link="identity")
  o <- order(t)
  plot(t[o],m(t[o]),type="l",col="green")
  lines(t[o],gh$m[o]); rug(t)

  ## partial linear probit model (GPLM)
  y <- (y>0)
  gh <- kgplm(x,t,y,h=0.25,family="bernoulli",link="probit")

  o <- order(t)
  plot(t[o],m(t[o]),type="l",col="green")
  lines(t[o],gh$m[o]); rug(t)

  ## data with two-dimensional m-function 
  n <- 1000; b <- c(1,-1); rho <- 0.7
  m <- function(t1,t2){ 1.5*sin(pi*t1)+t2 }
  x1 <- runif(n,min=-1,max=1); u  <- runif(n,min=-1,max=1)
  t1 <- runif(n,min=-1,max=1); t2 <- runif(n,min=-1,max=1)
  x2 <- round( m( rho*t1 + (1-rho)*u , t2 ) )
  x  <- cbind(x1,x2); t  <- cbind(t1,t2)
  y  <- x \%*\% b + m(t1,t2) + rnorm(n)

  ## partial linear model (PLM)
  grid1 <- seq(min(t[,1]),max(t[,1]),length=20)
  grid2 <- seq(min(t[,2]),max(t[,2]),length=25)
  grid  <- create.grid(list(grid1,grid2))

  gh <- kgplm(x,t,y,h=0.5,grid=grid,family="gaussian",link="identity")

  o <- order(grid[,2],grid[,1])
  est.m  <- (matrix(gh$m.grid[o],length(grid1),length(grid2)))
  orig.m <- outer(grid1,grid2,m)
  par(mfrow=c(1,2))
  persp(grid1,grid2,orig.m,main="Original Function",
        theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
  persp(grid1,grid2,est.m,main="Estimated Function",
        theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
  par(mfrow=c(1,1))
}
\keyword{smooth}

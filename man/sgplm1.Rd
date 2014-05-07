\name{sgplm1}
\alias{sgplm1}
\title{Generalized partial linear model}
\description{
  Fits a generalized partial linear model (based on smoothing spline)
  using the (generalized) Speckman estimator or backfitting (in the
  generalized case combined with local scoring) for two additive
  component functions.
  In contrast to \code{\link{kgplm}}, this function can be used
  only for a 1-dimensional nonparametric function. 
}
\usage{
sgplm1(x, t, y, spar, df=4, family, link,
       b.start=NULL, m.start=NULL, grid = NULL, offset = 0, 
       method = "speckman", weights = 1, weights.trim = 1, 
       weights.conv = 1, max.iter = 25, eps.conv = 1e-8,
       verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{n x p matrix, data for linear part}
  \item{y}{n x 1 vector, responses}
  \item{t}{n x 1 matrix, data for nonparametric part}
  \item{spar}{scalar smoothing parameter, as in \code{\link{smooth.spline}}}
  \item{df}{scalar equivalent number of degrees of freedom (trace of
          the smoother matrix), as in \code{\link{smooth.spline}}}
  \item{family}{text string, family of distributions (e.g.
    "gaussian" or "bernoulli", see details for \code{\link{glm.ll}})}
  \item{link}{text string, link function (depending on family,
    see details for \code{\link{glm.ll}})}
  \item{b.start}{p x 1 vector, start values for linear part}
  \item{m.start}{n x 1 vector, start values for nonparametric part}
  \item{grid}{m x q matrix, where to calculate the nonparametric function (default = t)}
  \item{offset}{offset}
  \item{method}{"speckman" or "backfit"}
  \item{weights}{binomial weights}
  \item{weights.trim}{trimming weights for fitting the linear part}
  \item{weights.conv}{weights for convergence criterion}
  \item{max.iter}{maximal number of iterations}
  \item{eps.conv}{convergence criterion}
  \item{verbose}{print additional convergence information}
  \item{...}{further parameters to be passed to \code{\link{smooth.spline}}}
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
\note{This function is mainly implemented for comparison. It is not
  really optimized for performance, however since it is spline-based, it
  should be sufficiently fast. Nevertheless, there might be several
  possibilities to improve for speed, in particular I guess that the
  sorting that \code{\link{smooth.spline}} performs in every iteration
  is slowing down the procedure quite a bit.
}
\seealso{\code{\link{kgplm}}}
\examples{
  ## generate data
  n <- 1000; b <- c(1,-1); rho <- 0.7
  mm <- function(t){ 1.5*sin(pi*t) }
  x1 <- runif(n,min=-1,max=1); u  <- runif(n,min=-1,max=1)
  t  <- runif(n,min=-1,max=1); x2 <- round(mm(rho*t + (1-rho)*u))
  x  <- cbind(x1,x2)
  y  <- x \%*\% b + mm(t) + rnorm(n)

  ## fit partial linear model (PLM)
  k.plm <- kgplm(x,t,y,h=0.35,family="gaussian",link="identity")
  s.plm <- sgplm1(x,t,y,spar=0.95,family="gaussian",link="identity")

  o <- order(t)
  ylim <- range(c(mm(t[o]),k.plm$m,s.plm$m),na.rm=TRUE)
  plot(t[o],mm(t[o]),type="l",ylim=ylim)
  lines(t[o],k.plm$m[o], col="green")
  lines(t[o],s.plm$m[o], col="blue")
  rug(t); title("Kernel PLM vs. Spline PLM")

  ## fit partial linear probit model (GPLM)
  y <- (y>0)
  k.gplm <- kgplm(x,t,y,h=0.35,family="bernoulli",link="probit")
  s.gplm <- sgplm1(x,t,y,spar=0.95,family="bernoulli",link="probit")

  o <- order(t)
  ylim <- range(c(mm(t[o]),k.gplm$m,s.gplm$m),na.rm=TRUE)
  plot(t[o],mm(t[o]),type="l",ylim=ylim)
  lines(t[o],k.gplm$m[o], col="green")
  lines(t[o],s.gplm$m[o], col="blue")
  rug(t); title("Kernel GPLM vs. Spline GPLM (Probit)")
}
\keyword{smooth}

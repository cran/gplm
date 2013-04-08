\name{glm.ll}
\alias{glm.ll}
\title{Log-likelihood for GLM}
\description{
  Calculates the log-likelihood function of a GLM.
  Currently only the gaussian and the bernoulli family
  are implemented.
}
\usage{
glm.ll(mu, y, phi=1, family="gaussian", k=1)
}
\arguments{
  \item{mu}{n x 1, predicted regression function}
  \item{y}{n x 1, responses}
  \item{phi}{scalar, nuisance parameter (sigma^2 for the gaussian and
    inverse gaussian families, nu for the gamma family)}
  \item{family}{text string, family of distributions (e.g.
    "gaussian" or "bernoulli", see details below)}
  \item{k}{integer > 0, parameter for the negative binomial}
}
\details{Implemented are the "gaussian" family
  (with links "identity" and "log"), the "bernoulli" family
  (with links "logit" and "probit"), the "gamma" family
  (with link "inverse"), the "poisson" family (with link "log"),
  the "inverse.gaussian" family (with link "inverse.squared") and
  the "negative.binomial" (with its canonical "log" type link).

  The default value k=1 leads to the geometric distribution (as a
  special case of the negative binomial). 
}
\value{
  log-likelihood value
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{glm.lld}}, \code{\link{glm.link}}}
\examples{
  glm.ll(rep(0.4,2), c(0,1), family="bernoulli")
}
\keyword{smooth}

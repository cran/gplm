\name{glm.lld}
\alias{glm.lld}
\title{Log-likelihood derivatives for GLM}
\description{
  Computes first and second derivatives of the individual 
  log-likelihood with respect to the linear predictor.
  Currently only the gaussian (with identity link) and the bernoulli family
  (with logit and probit links) are implemented.
}
\usage{
glm.lld(eta, y, family="gaussian", link="identity", k=1)
}
\arguments{
  \item{eta}{n x 1, linear predictors}
  \item{y}{n x 1, responses}
  \item{family}{text string, family of distributions (e.g.
    "gaussian" or "bernoulli", see details for \code{\link{glm.ll}})}
  \item{link}{text string, link function (depending on family,
    see details for \code{\link{glm.ll}})}
  \item{k}{integer > 0, parameter for the negative binomial}
}
\details{See details for \code{\link{glm.ll}}.
}
\value{
  List with components:
  \item{ll1}{n x 1, vector of first derivatives}
  \item{ll2}{n x 1, vector of second derivatives}
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{glm.ll}}, \code{\link{glm.link}}}
\examples{
  glm.lld(c(-1,2), c(0,1), family="bernoulli", link="logit")
}
\keyword{smooth}

\name{glm.inverse.link}
\alias{glm.inverse.link}
\title{Link function for GLM}
\description{
  Defines the link function for a GLM.
}
\usage{
glm.inverse.link(mu, family="gaussian", link="identity", k=1)
}
\arguments{
  \item{mu}{n x 1, linear predictors}
  \item{family}{text string, family of distributions (e.g.
    "gaussian" or "bernoulli", see details for \code{\link{glm.ll}})}
  \item{link}{text string, link function (depending on family,
    see details for \code{\link{glm.ll}})}
  \item{k}{integer > 0, parameter for the negative binomial}
}
%%\details{
%%}
\value{
  n x 1, vector eta (predictors)
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{glm.ll}}, \code{\link{glm.lld}}, \code{\link{glm.link}}}
\examples{
  glm.inverse.link(c(0.25,0.5), family="bernoulli", link="logit")
}
\keyword{smooth}

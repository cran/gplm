\name{glm.link}
\alias{glm.link}
\title{(Inverse) Link function for GLM}
\description{
  Defines the inverse link function for a GLM.
  Currently only the gaussian and the bernoulli family
  are implemented.
}
\usage{
glm.link(eta, family="gaussian", link="identity", k=1)
}
\arguments{
  \item{eta}{n x 1, linear predictors}
  \item{family}{text string, family of distributions (e.g.
    "gaussian" or "bernoulli", see details for \code{\link{glm.ll}})}
  \item{link}{text string, link function (depending on family,
    see details for \code{\link{glm.ll}})}
  \item{k}{integer > 0, parameter for the negative binomial}
}
%%\details{
%%}
\value{
  n x 1, vector mu = glm.link( eta )
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{glm.ll}}, \code{\link{glm.lld}}}
\examples{
  glm.link(c(-1,2), family="bernoulli", link="logit")
}
\keyword{smooth}

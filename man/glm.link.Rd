\name{glm.link}
\alias{glm.link}
\title{(Inverse) Link function for GLM}
\description{
  Defines the inverse link function for a GLM.
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
  n x 1, vector mu (responses)
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{glm.ll}}, \code{\link{glm.lld}}, \code{\link{glm.inverse.link}}}
\examples{
  glm.link(c(-1,2), family="bernoulli", link="logit")
}
\keyword{smooth}

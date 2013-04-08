\name{kbackfit}
\alias{kbackfit}
\title{Backfitting for an additive model using kernel regression}
\description{
  Implements kernel-based backfitting in an additive model,
  optional with a partial linear term.
}
\usage{
kbackfit(t, y, h, x = NULL, grid = NULL, weights.conv = 1,
           offset = 0, method = "generic",
           max.iter = 50, eps.conv = 1e-04, m.start = NULL,
           kernel = "biweight")
}
\arguments{
  \item{y}{n x 1 vector, responses}
  \item{t}{n x q matrix, data for nonparametric part}
  \item{h}{scalar or 1 x q, bandwidth(s)}
  \item{x}{optional, n x p matrix, data for linear part}
  \item{grid}{m x q matrix, where to calculate the nonparametric function (default = t)}
  \item{weights.conv}{weights for convergence criterion}
  \item{offset}{offset}
  \item{method}{one of \code{"generic"}, \code{"linit"} or \code{"modified"}}
  \item{max.iter}{maximal number of iterations}
  \item{eps.conv}{convergence criterion}
  \item{m.start}{n x q matrix, start values for m}
  \item{kernel}{text string, see \code{\link{kernel.function}}}
}
%%\details{
%%}
\value{
  List with components:
  \item{c}{constant}
  \item{b}{p x 1 vector, linear coefficients}
  \item{m}{n x q matrix, nonparametric marginal function estimates}
  \item{m.grid}{m x q matrix, nonparametric marginal function estimates
    on grid}
  \item{rss}{residual sum of squares}
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{kernel.function}}, \code{\link{kreg}}}
%%\examples{}
\keyword{smooth}

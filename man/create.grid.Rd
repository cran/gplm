\name{create.grid}
\alias{create.grid}
\title{Create a grid for kernel estimation}
\description{
  Helps to define a grid for kernel denity or regression estimates
  (univariate or multivariate).
}
\usage{
create.grid(grid.list, sort=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid.list}{list of 1-dimensional vectors containing the grid
    values for each dimension}
  \item{sort}{sort the vectors (can be set to FALSE if vectors are
    already sorted in ascending order)}
}
\details{
  This function allows easily to define grids for the "gplm" package.
  If the data are d-dimensional and the grid vector lengths are n1,
  ... nd, then the output is a (n1*...*nd) x d matrix with each row
  corresponding to one d-dimensional data point at which the function
  estimate is to be calculated.
}
\value{
  m x d grid matrix
}
%%\references{}
\author{Marlene Mueller}
%%\note{}
\seealso{\code{\link{expand.grid}}, \code{\link{kde}}, \code{\link{kreg}}}
\examples{
  v1 <- 1:5
  v2 <- 3:1
  grid <- create.grid(list(v1,v2))

  x <- matrix(rnorm(60),30,2)
  v1 <- seq(min(x[,1]),max(x[,1]),length=10)
  v2 <- seq(min(x[,2]),max(x[,2]),length=5)
  grid <- create.grid(list(v1,v2))
}
\keyword{smooth}

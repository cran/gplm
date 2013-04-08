"kernel.function" <- function(u,kernel="biweight",product=TRUE){

  p <- 2; q <- 2  ## biweight kernel 
  if(kernel=="triangle"){                      p <- q <- 1 }    
  if(kernel=="uniform"){                       p <- 1; q <- 0 } 
  if(kernel=="epanechnikov"){                  p <- 2; q <- 1 } 
  if(kernel=="biweight" || kernel=="quartic"){ p <- 2; q <- 2 } 
  if(kernel=="triweight"){                     p <- 2; q <- 3 } 
  if(kernel=="gaussian"|| kernel=="normal"){   p <- 0; q <- 0 }
  
  u <- as.matrix(u)
  d <- ncol(u)
  if (p > 0){
    if (product){
      x <- 1-sqrt(u*u)^p
      c <- p*gamma(1/2)*gamma(1+q+1/p)/( 2*pi^(1/2)*gamma(1+q)*gamma(1/p) )
      k <- (c^d) * apply(x^q,1,prod) * apply(x>=0,1,prod)
    }else{
      x <- 1-sqrt(rowSums(u*u))^p
      c <- p*gamma(d/2)*gamma(1+q+d/p)/( 2*pi^(d/2)*gamma(1+q)*gamma(d/p) )
      k <- c * x^q * (x>=0)
    }
  }else{
    k <- apply(dnorm(u),1,prod)
  }
  return(k)
}

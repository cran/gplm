"kernel.constants" <- function(kernel="biweight",d=1,product=TRUE){

  p <- 2; q <- 2  ## biweight kernel 
  if(kernel=="triangle"){                      p <- q <- 1 }    
  if(kernel=="uniform"){                       p <- 1; q <- 0 } 
  if(kernel=="epanechnikov"){                  p <- 2; q <- 1 } 
  if(kernel=="biweight" || kernel=="quartic"){ p <- 2; q <- 2 } 
  if(kernel=="triweight"){                     p <- 2; q <- 3 } 
  if(kernel=="gaussian"|| kernel=="normal"){   p <- 0; q <- 0 }

  m2 <- c2 <- d0 <- NA
  v <- 2*pi^(d/2)/( d*gamma(d/2))   ## volume of the d-dimensional unit sphere 
  
  if (product || p==0){
    if (p > 0){
      m2 <-  p*gamma(1+q+1/p)*gamma(1+3/p) / ( 3*gamma(1/p)*gamma(1+3/p+q) )
      c2 <-  p^2 *(gamma(1+1/p+q))^2 *gamma(1+1/p)*gamma(1+2*q)/
        ( 2*gamma(1+1/p+2*q)*(gamma(1/p)*gamma(1+q))^2 )
      c2 <- c2^d
    }else{
      m2 <- 1
      c2 <- 1/((2*sqrt(pi))^2)
    }
    d0 <- (c2/(m2^2))^(1/(d+4))
  }else{
    if (p==1 && q==1){
      m2 <- c(1/6, 0.15, 2/15, 0.12, 0.1072)[d]
      c2 <- (2+2*d)/((2+d)*v)
    }
    if (p==2 && q>=0){
      m2 <- 1/(2*q+d+2)
      if (q==0){ c2 <- 1/v }
      if (q==1){ c2 <- (4+2*d)/((4+d)*v)}
      if (q==2){ c2 <- c(5/7, 9/(5*pi), 35/(22*pi), 24/(5*pi^2))[d]}
      if (q==3){ c2 <- c(350/429, 16/(7*pi), 315/(143*pi), 50/(7*pi^2))[d]}
    }
  }
    
  return(list(m2=m2,c2=c2,d0=d0))
}

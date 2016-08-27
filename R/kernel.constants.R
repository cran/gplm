"kernel.constants" <- function(kernel="biweight",d=1,product=TRUE){

  if (kernel=="triangular"){ kernel <- "triangle" }
  if (kernel=="rectangle" || kernel=="rectangular"){ kernel <- "uniform" }
  if (kernel=="quartic"){ kernel <- "biweight" }
  if (kernel=="normal"){  kernel <- "gaussian" }

  kernel.names <- c("triangle","uniform","epanechnikov","biweight",
                    "triweight","gaussian")

  m2 <- n2 <- d0 <- NA
  volume.d <- pi^(d/2)/gamma(d/2+1)  ## volume of d-dim. unit sphere

  mm2 <- c(1/6, 1/3, 0.2, 1/7, 1/9, 1)
  nn2 <- c(2/3, 0.5, 0.6, 5/7, 350/429, 1/(2*sqrt(pi)))
                              ## cf. Wand & Jones, p. 176
  pp <- c(1,2,2,2,2,0)
  qq <- c(1,0,1,2,3,NA)
  names(mm2) <- names(nn2) <- names(pp) <- names(qq) <- kernel.names

  p <- pp[kernel]
  q <- qq[kernel]

  if (product || p==0){
    m2 <-  mm2[kernel]
    n2 <-  nn2[kernel]
    d0 <- (n2^d/(m2^2))^(1/(d+4))
  }else{
    if (p==1){
      m2 <- (d+1)/((d+2)*(d+3))
      n2 <- 2*(d+1)/((d+2)*volume.d)
    }
    if (p==2){
      m2 <- 1/(2*q+d+2)
      if (q==0){ fac <- 1 }               ## uniform
      if (q==1){ fac <- 2*(d+2)/(d+4) }   ## epanechnikov
      if (q==2){ fac <- 6*(d+2)*(d+4)/((d+6)*(d+8)) }                ## biweight
      if (q==3){ fac <- 20*(d+2)*(d+4)*(d+6)/((d+8)*(d+10)*(d+12)) } ## triweight
      n2 <- fac/volume.d
    }
    d0 <- (n2/(m2^2))^(1/(d+4))
  }

  return(list(m2=m2,n2=n2,d0=d0))
}

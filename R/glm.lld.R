"glm.lld" <- function(eta, y, family="gaussian", link="identity", k=1){

  if (family=="bernoulli"){
    if (link=="logit"){
      e <- exp(eta); g <- e/(1+e)
      ll1  <-  y - g; ll2  <-  -g/(1+e) 
    }
    if (link=="probit"){
      g    <-  pnorm(eta); gg <- g*pnorm(-eta)
      gg8  <-  pnorm(8)*(1-pnorm(8)) 
      gg   <-  gg*(gg>gg8) + gg8*(gg<= gg8)  
      g1   <-  dnorm(eta)
      ll1  <-  (y - g)*g1/gg
      g11  <-  g1*g1
      ll2  <-  -eta*g1/gg - (1-2*g)*g11/(gg*gg)
      ll2  <-  ll2*(y-g); ll2 <- ll2 - g11/gg
    }
  }

  if (family=="gaussian"){
    if (link=="identity"){
      ll1 <- (y-eta)
      ll2 <- rep(-1,length(eta))
    }
    if (link=="log"){
      g <- exp(eta)
      ll1 <- (y-g)*g
      ll2 <- (y-2*g)*g
    }
  }

  if (family=="poisson"){
    if (link=="log"){
      g <- exp(eta)
      ll1 <- y - g
      ll2 <- -g
    }
  }

  if (family=="gamma"){
    if (link=="inverse"){
      g <- 1/eta
      ll1 <- - y + g
      ll2 <- - g*g
    }
  }

  if (family=="inverse.gaussian"){
    if (link=="inverse.squared"){
      g <- 1/sqrt(eta)
      ll1 <- 0.5*(-y + g)
      ll2 <- -0.25 *g/eta
    }
  }

  if (family=="negative.binomial"){
    if (link=="log"){
      e <- exp(eta); ee <- 1-e
      ll1 <- y - k*e/ee
      ll2 <- -k/(ee*ee)
    }
  }
  
  l <- data.frame(ll1,ll2)
  return(l)
}


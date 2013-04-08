"glm.link" <- function(eta, family="gaussian", link="identity", k=1){
  if (family=="bernoulli"){
    if (link=="logit"){
      e <- exp(-eta); mu<- 1/(1+e)
    }
    if (link=="probit"){
      mu <- pnorm(eta)
    }
  }
  if (family=="gaussian"){
    if (link=="identity"){
      mu <- eta
    }
  }
  if (link=="log"){  ## gaussian+poisson
    mu <- exp(eta)
  }
  if (family=="gamma"){
    if (link=="reciprocal"){
      mu <- 1/eta
    }
  }
  if (family=="inverse.gaussian"){
    if (link=="quadratic.reciprocal"){
      mu <- 1/sqrt(eta)
    }
  }
  if (family=="negative.binomial"){
    if (link=="log"){
      e <- exp(eta)
      mu <- e/(k*(1-e))
    }
  }
  
  return(mu)
}


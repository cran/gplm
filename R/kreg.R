"kreg" <- function(x,y,bandwidth=NULL,grid=TRUE,kernel="biweight",
                   product=TRUE,sort=TRUE){

  p <- 2; q <- 2  ## biweight kernel 
  if(kernel=="triangle"){                      p <- q <- 1 }    
  if(kernel=="uniform"){                       p <- 1; q <- 0 } 
  if(kernel=="epanechnikov"){                  p <- 2; q <- 1 } 
  if(kernel=="biweight" || kernel=="quartic"){ p <- 2; q <- 2 } 
  if(kernel=="triweight"){                     p <- 2; q <- 3 } 
  if(kernel=="gaussian"|| kernel=="normal"){   p <- 0; q <- 0 }
  
  x <- as.matrix(x)  ## nxd
  y <- as.matrix(y)  ## nxc
  n <- nrow(x)
  d <- ncol(x)
  c <- ncol(y)
  
  if (sort){
    or <- order(x[,1])
    ro <- order((1:n)[or])
    x  <- x[or,,drop=FALSE]
    y  <- y[or,,drop=FALSE]
  }
  
  if (is.logical(grid)){  ## grid TRUE, but not given
    if (grid){
      ng <- 400
      n.grid <- max(c(ng^(1/d),5))
      grids <- vector("list", d)
      for (j in d:1){
        grids[[d-j+1]] <- seq(min(x[,j]),max(x[,j]),length=n.grid)
      }
      grid <- as.matrix(expand.grid(grids))
      grid <- grid[,d:1,drop=FALSE]
      m <- nrow(grid)
      ro.grid <- NULL
    }else{                ## grid FALSE
      grid <- x
      m <- n
      ro.grid <- ro
    }
  }else{                  ## grid given
    grid <- as.matrix(grid)
    if (sort){
      m <- nrow(grid)
      or.grid <- order(grid[,1])
      ro.grid <- order((1:m)[or.grid])
      grid  <- grid[or.grid,,drop=FALSE]
    }
  }
  
  ##  if (is.null(grid)){
  ##    grid <- x
  ##    m <- n
  ##    org <- or
  ##    rog <- ro
  ##  }else{
  ##    grid <- as.matrix(grid)     ## mxd
  ##    m <- nrow(grid)
  ##    if (sort){
  ##      org <- order(grid[,1])
  ##      rog <- order((1:m)[org])
  ##      grid  <- grid[org,,drop=FALSE]
  ##    }
  ##  }
  
  if (missing(bandwidth)) {
    if (p > 0){  ## non-gaussian
      fac <- kernel.constants(p,q,d)$d0 / (1/(2*sqrt(pi))^d) ^(1/(d+4))
    }else{       ## gaussian
      fac <- 1
    }
    s <- rbind( sd(x), diff( apply(x,2,quantile,c(0.25,0.75)) )/1.349 )
    s <- apply(s,2,min)
    bandwidth <- fac *s* n^(-1/(d+4)) ## Scott's ROT
  }
  
  tmp <- convol(x,bandwidth,grid=grid,y=cbind(rep(1,n),y),
                p=p,q=q,product=product,sort=FALSE)
  mh <- as.matrix((tmp[,2:ncol(tmp)])/(tmp[,1]))
  
  kconst <- kernel.function(t(rep(0,q)), kernel=kernel, product=product)
  df <- n - kconst*sum(1/tmp[,1])/(n*prod(bandwidth))
  
  if (sort){
    return(list(x=grid,y=mh,bandwidth=bandwidth,df=df,rearrange=ro.grid))
  }else{
    return(list(x=grid,y=mh,bandwidth=bandwidth,df=df))
  }
}



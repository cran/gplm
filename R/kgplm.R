"kgplm" <- function(x,t,y,h,family,link,
                    b.start=NULL,m.start=NULL,grid=NULL,m.grid.start=NULL,
                    offset=0,method="speckman",sort=TRUE,weights=1,
                    weights.trim=1,weights.conv=1,max.iter=25,eps.conv=1e-8,
                    kernel="biweight",kernel.product=TRUE,verbose=FALSE){

  kernel.p <- 2; kernel.q <- 2  ## biweight kernel 
  if(kernel=="triangle"){                      kernel.p <- kernel.q <- 1 }    
  if(kernel=="uniform"){                       kernel.p <- 1; kernel.q <- 0 } 
  if(kernel=="epanechnikov"){                  kernel.p <- 2; kernel.q <- 1 } 
  if(kernel=="biweight" || kernel=="quartic"){ kernel.p <- 2; kernel.q <- 2 } 
  if(kernel=="triweight"){                     kernel.p <- 2; kernel.q <- 3 } 
  if(kernel=="gaussian"|| kernel=="normal"){   kernel.p <- 0; kernel.q <- 0 }

  x <- as.matrix(x)
  t <- as.matrix(t)
  y <- as.matrix(y)

  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(t)

  if (length(h)==1){  h  <- rep(h,q) }##matrix(h,1,q) }

  if (length(weights)==1){ weights <- rep(weights,n) }
  if (length(weights.trim)==1){ weights.trim <- rep(weights.trim,n) }
  if (length(weights.conv)==1){ weights.conv <- rep(weights.conv,n) }
  if (length(offset)==1){ offset<-rep(offset,n) }

  n.grid <- 0
  m.grid <- NULL
  havegrid<- !is.null(grid)
  if (havegrid) {
    grid <- as.matrix(grid)
    n.grid <- nrow(grid)
  }

  if (sort) { 
    or <- order(t[,1]) 
    ro <- order((1:n)[or])
    t <- as.matrix(t[or,])

    x <- as.matrix(x[or,])
    y <- as.matrix(y[or,])
    weights <- as.matrix(weights[or])
    weights.trim <- as.matrix(weights.trim[or])
    weights.conv  <- as.matrix(weights.conv[or])
    offset <- as.matrix(offset[or])
    if (!missing(m.start)) { m.start <- as.matrix(m.start[or]) }

    if (havegrid){
      or.grid <- order(grid[,1]) 
      ro.grid <- order((1:n.grid)[or.grid])
      grid <- as.matrix(grid[or.grid,])
      if (!missing(m.grid.start)) { m.grid.start <- as.matrix(m.grid.start[or,]) }
    }
  }

  one <- rep(1,n)

  if (missing(b.start)) { b.start<- rep(0,p) }
  if (missing(m.start)) { m.start<- rep(0,n) }

  xnew <- matrix(rep(0,n*p),nrow=n)
  t<-t(t(t)/h)
  if (havegrid){
    grid<-t(t(grid)/h)
    if (missing(m.grid.start)) { m.grid.start <- rep(0,n.grid) }
  }

  it <- 0
  stop.crit <- FALSE
  xb <- x %*% b.start + offset
  b <- b.start
  dev.start <- Inf

  while ( (stop.crit==FALSE)&(it< max.iter)){
    it <- it +1
    if (verbose){
      print( paste("kgplm: iteration no.",as.character(it)) )
    }
    
    ll <- glm.lld(xb+m.start,y,family=family,link=link)
    zm <- m.start-ll$ll1/ll$ll2
    z  <- xb+zm
    wnew <- as.vector(weights* ll$ll2)
    
    tmp <- convol(t,y=as.matrix(cbind(one,zm,z,x))*as.vector(wnew),
                  p=kernel.p,q=kernel.q,product=kernel.product,sort=FALSE)
    denom <- tmp[,1]
    m <- tmp[,2]/denom
    ##print(data.frame(tmp[,1:2]/h,z,xb,m.start,y,ll)[1:10,])
    xnew <- as.matrix( x-tmp[,4:ncol(tmp)]/denom )
    znew <- as.matrix( z-tmp[,3]/denom )
    
    if (havegrid){
      tmp <- convol(t,grid=grid,y=cbind(one,zm)*as.vector(wnew),
                    p=kernel.p,q=kernel.q,product=kernel.product,sort=FALSE)
      ##print(dim(tmp))
      m.grid <- tmp[,2]/tmp[,1]
    }
    
    if (method=="backfitting"){ 
      wnew<- x*as.vector(weights.trim*weights*ll$ll2)
    }
    else{
      wnew<- xnew*as.vector(weights.trim*weights*ll$ll2)
    }
    
    B  <- t(wnew) %*% xnew
    ##print(B)
    bv <- chol2inv(chol(-B))
    b  <- solve(B,t(wnew) %*% znew)

    eta <- xb+m
    mu  <- glm.link(eta, family=family, link=link)
    dev <- 2*sum(weights *glm.ll(y,y,family=family)-weights*glm.ll(mu,y,family=family))
    chgd <- abs((dev-dev.start)/dev.start)
    db   <- b-b.start
    chgb <- sqrt( sum(db*db)/sum(b.start*b.start) )
    dm   <- m-m.start
    chgm <- sqrt( sum(dm*dm)/sum(m.start*m.start) )

    if (it==1){
      stop.crit<- FALSE
    }
    else{
      stop.crit <- ( ((chgb<eps.conv)&(chgm<eps.conv)) | (chgd<eps.conv) )
    }
    if (verbose){
      print(paste("  Deviance: ",as.character(dev)))
      print(paste("  Change in b: ",as.character(chgb)))
      print(paste("  Change in m: ",as.character(chgm)))
      print(paste("  Change in deviance: ",as.character(chgd)))
    }
    b.start <- b
    m.start <- m
    dev.start <- dev
    xb<- x %*% b.start+offset
    if (havegrid==1) { m.grid.start <- m.grid }
  }
  
  tmp  <- convol(t,y=(xnew*as.vector(weights* ll$ll2)))
  xnew <- xnew-tmp/denom
  ##denom <- denom*n*prod(h)
  tmp  <- t(wnew) %*% xnew
  kconst <- kernel.function(t(rep(0,q)), kernel=kernel, product=kernel.product)
  df   <- n + sum(diag(bv %*% tmp)) - kconst*sum(as.vector(weights* ll$ll2)/denom)
  aic <- 2*(n-df)+dev
  ##print("kgplm")
  ##print(paste("df",df))
  ##print(paste("aic",aic))
  ##print(sum(diag(bv %*% tmp)))
  ##print(kconst*sum(as.vector(weights* ll$ll2)/denom))

  m<-m[ro]
  return(list(b=b,bv=bv,m=m,m.grid=m.grid,it=it,df.residual=df,deviance=dev,aic=aic))
}


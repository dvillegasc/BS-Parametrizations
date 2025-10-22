dBS4 <- function(x, mu=1, sigma=0.5, log=FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS4 to BS (original)
  new_mu    <- (sigma/mu)
  new_sigma <- 1 / sqrt(mu*sigma)
  
  res <- dBS(x=x, mu=new_mu, sigma=new_sigma, log=log)
  return(res)
}


pBS4 <- function(q, mu=1, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS4 to BS (original)
  new_mu    <- (sigma/mu)
  new_sigma <- 1 / sqrt(mu*sigma)
  
  cdf <- pBS(q=q, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  
  return(cdf)
}

qBS4 <- function(p, mu=1, sigma=0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0 | sigma>=1)) 
    stop(paste("sigma must be positive", "\n", ""))
  
  # Changing from BS4 to BS (original)
  new_mu    <- (sigma/mu)
  new_sigma <- 1 / sqrt(mu*sigma)
  
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  
  q <- qBS(p=p, mu=new_mu, sigma=new_sigma, lower.tail=lower.tail, log.p=log.p)
  return(q)
}


rBS4 <- function(n, mu=1, sigma=0.5){
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
      
      # Changing from BS4 to BS (original)
      new_mu    <- (sigma/mu)
      new_sigma <- 1 / sqrt(mu*sigma)
      
      r <- rBS(n=n, mu=new_mu, sigma=new_sigma)
      r
}


hBS4 <- function(x, mu, sigma){
  if (any(x < 0)) 
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 )) 
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
      
      h <- dBS4(x, mu, sigma) / pBS4(x, mu, sigma, lower.tail=FALSE)
      h
}


BS10 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "dBS10", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "dBS10", substitute(sigma.link),
                      c("log", "logit", "probit", "own"))
  structure(
    list(family = c("dBS10", "Birnbaum-Saunders - tenth parameterization"),
         parameters = list(mu=TRUE, sigma=TRUE),
         nopar = 2,
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         
         # First derivatives
         dldm = function(y, mu, sigma) {
           a0 <- mu / sqrt(2 * sigma)
           b0 <- (mu^2) / 2
           da_dm <- 1 / sqrt(2 * sigma)
           db_dm <- mu
           
           term1 <- (-1 / a0) * da_dm
           term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_dm
           term3 <- (1 / (y + b0)) * db_dm
           term4 <- (-1 / (2 * b0)) * db_dm
           term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           result <- term1 + term2 + term3 + term4 + term5
           return(result)
         },
         
         dldd = function(y, mu, sigma) { 
           a0 <- mu / sqrt(2 * sigma)
           b0 <- (mu^2) / 2
           da_ds <- -mu / (2 * sigma * sqrt(2 * sigma))
           db_ds <- 0 # Beta no depende de sigma
           
           term1 <- (-1 / a0) * da_ds
           term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds

           result <- term1 + term2 
           return(result)
         },
         
         # Second derivatives
         
         d2ldm2 = function(y, mu, sigma) {
           a0 <- mu / sqrt(2 * sigma)
           b0 <- (mu^2) / 2
           da_dm <- 1 / sqrt(2 * sigma)
           db_dm <- mu
           
           # dldm
           t1 <- (-1 / a0) * da_dm
           t2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_dm
           t3 <- (1 / (y + b0)) * db_dm
           t4 <- (-1 / (2 * b0)) * db_dm
           t5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           dldm <- t1 + t2 + t3 + t4 + t5
           return(-dldm * dldm) 
         },
         
         d2ldd2 = function(y, mu, sigma) {
           a0 <- mu / sqrt(2 * sigma)
           b0 <- (mu^2) / 2
           
           da_ds <- -mu / (2 * sigma * sqrt(2 * sigma))
           db_ds <- 0 # <--
           
           # dldd
           t1 <- (-1 / a0) * da_ds
           t2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           
           dldd <- t1 + t2
           return(-dldd * dldd)
         },
         
         d2ldmdd = function(y, mu, sigma) {
           a0 <- mu / sqrt(2 * sigma)
           b0 <- (mu^2) / 2
           
           da_dm <- 1 / sqrt(2 * sigma)
           db_dm <- mu
           da_ds <- -mu / (2 * sigma * sqrt(2 * sigma))
           db_ds <- 0 # <--
           
           # dldm
           m1 <- (-1 / a0) * da_dm
           m2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_dm
           m3 <- (1 / (y + b0)) * db_dm
           m4 <- (-1 / (2 * b0)) * db_dm
           m5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           dldm <- m1 + m2 + m3 + m4 + m5
           
           # dldd
           d1 <- (-1 / a0) * da_ds
           d2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           dldd <- d1 + d2
           
           return(-dldm * dldd)
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS10(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS10", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(mean(y), length(y))}),
         sigma.initial = expression({sigma <- rep(2, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

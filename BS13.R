BS13 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "dBS13", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "dBS13", substitute(sigma.link),
                      c("log", "logit", "probit", "own"))
  structure(
    list(family = c("dBS13", "Birnbaum-Saunders - thirteenth parameterization"),
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
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           db_dm <- 1 / sigma
           da_dm <- 0 # <--
           
           term3 <- (1 / (y + b0)) * db_dm
           term4 <- (-1 / (2 * b0)) * db_dm
           term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           result <- term3 + term4 + term5
           return(result)
         },
         
         dldd = function(y, mu, sigma) { 
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           db_ds <- -b0 / sigma
           da_ds <- -1 / (2 * sigma * sqrt(sigma))
           
           term1 <- (-1 / a0) * da_ds
           term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           term3 <- (1 / (y + b0)) * db_ds
           term4 <- (-1 / (2 * b0)) * db_ds
           term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           
           result <- term1 + term2 + term3 + term4 + term5
           return(result)
         },
         
         # Second derivatives
         
         d2ldm2 = function(y, mu, sigma) {
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           db_dm <- 1 / sigma
           
           t3 <- (1 / (y + b0)) * db_dm
           t4 <- (-1 / (2 * b0)) * db_dm
           t5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           dldm <- t3 + t4 + t5
           
           return(-dldm * dldm) 
         },
         
         d2ldd2 = function(y, mu, sigma) {
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           
           db_ds <- -b0 / sigma
           da_ds <- -1 / (2 * sigma * sqrt(sigma))
           
           t1 <- (-1 / a0) * da_ds
           t2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           t3 <- (1 / (y + b0)) * db_ds
           t4 <- (-1 / (2 * b0)) * db_ds
           t5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           dldd <- t1 + t2 + t3 + t4 + t5
           
           return(-dldd * dldd)
         },
         
         d2ldmdd = function(y, mu, sigma) {
           b0 <- mu / sigma
           a0 <- 1 / sqrt(sigma)
           
           db_dm <- 1 / sigma
           da_ds <- -1 / (2 * sigma * sqrt(sigma))
           db_ds <- -b0 / sigma
           
           # dldm
           m3 <- (1 / (y + b0)) * db_dm
           m4 <- (-1 / (2 * b0)) * db_dm
           m5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           dldm <- m3 + m4 + m5
           
           # dldd
           d1 <- (-1 / a0) * da_ds
           d2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           d3 <- (1 / (y + b0)) * db_ds
           d4 <- (-1 / (2 * b0)) * db_ds
           d5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           dldd <- d1 + d2 + d3 + d4 + d5
           
           return(-dldm * dldd)
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS13(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS13", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(mean(y), length(y))}),
         sigma.initial = expression({sigma <- rep(2, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

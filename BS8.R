BS8 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "dBS8", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "dBS8", substitute(sigma.link),
                      c("log", "logit", "probit", "own"))
  structure(
    list(family = c("dBS8", "Birnbaum-Saunders - five parameterization"),
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
           a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
           b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
           db_dm <- b0 / (2 * mu)
           
           term1 <- (1 / (y + b0)) * db_dm
           term2 <- -1 / (2 * b0) * db_dm
           term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           result <- term1 + term2 + term3
           return(result)
         },
         
         dldd = function(y, mu, sigma) { 
           a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
           b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
           
           da_ds <- 1 / sqrt(5 * (sigma - 1))
           db_ds <- -b0 * (2 * sigma - 1) / (2 * sigma * (sigma - 1))
           
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
           a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
           b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
           db_dm <- b0 / (2 * mu)
           
           term1 <- (1 / (y + b0)) * db_dm
           term2 <- -1 / (2 * b0) * db_dm
           term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           
           dldm <- term1 + term2 + term3
           
           return(-dldm * dldm) 
         },
         
         d2ldd2 = function(y, mu, sigma) {
           a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
           b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
           
           da_ds <- 1 / sqrt(5 * (sigma - 1))
           db_ds <- -b0 * (2 * sigma - 1) / (2 * sigma * (sigma - 1))
           
           term1 <- (-1 / a0) * da_ds
           term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           term3 <- (1 / (y + b0)) * db_ds
           term4 <- (-1 / (2 * b0)) * db_ds
           term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           
           dldd <- term1 + term2 + term3 + term4 + term5
           
           return(-dldd * dldd)
         },
         
         d2ldmdd = function(y, mu, sigma) {
           a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
           b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
           
           db_dm <- b0 / (2 * mu)
           da_ds <- 1 / sqrt(5 * (sigma - 1))
           db_ds <- -b0 * (2 * sigma - 1) / (2 * sigma * (sigma - 1))
           
           # dldm
           m1 <- (1 / (y + b0)) * db_dm
           m2 <- -1 / (2 * b0) * db_dm
           m3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
           dldm <- m1 + m2 + m3
           
           # dldd
           d1 <- (-1 / a0) * da_ds
           d2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
           d3 <- (1 / (y + b0)) * db_ds
           d4 <- (-1 / (2 * b0)) * db_ds
           d5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
           dldd <- d1 + d2 + d3 + d4 + d5
           
           return(-dldm * dldd)
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS8(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS8", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(mean(y), length(y))}),
         sigma.initial = expression({sigma <- rep(2, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 1),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

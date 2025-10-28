BS4 <- function(mu.link = "log", sigma.link = "logit"){
  mstats <- checklink("mu.link", "BS4", substitute(mu.link),
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "BS4", substitute(sigma.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  structure(
    list(family = c("BS4", "Birnbaum-Saunders - four parameterization"),
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
         dldm = function(y, sigma, mu) {
           
           t <- y / mu
           term1 <- -1 / mu
           term2 <- (t + 3) / (2 * mu * (t + 1))
           term3 <- -(1 / (2 * sigma^2 * mu)) * (1/t - t)
           
           result <- term1 + term2 + term3
           return(result)
         },
         
         dldd = function(y, sigma, mu) {
           
           t <- y / mu
           B <- t + (1/t) - 2
           term1 <- -1 / sigma
           term2 <- (1 / (sigma^3)) * B
           
           result <- term1 + term2
           return(result)
         },
         
         # Second derivatives
         
         d2ldm2 = function(y, sigma, mu) {
           
           t <- y / mu
           term1 <- -1 / mu
           term2 <- (t + 3) / (2 * mu * (t + 1))
           term3 <- -(1 / (2 * sigma^2 * mu)) * (1/t - t)
           result <- term1 + term2 + term3
           
           return(-result * result)
         },
         
         d2ldd2 = function(y, sigma, mu) {
           
           t <- y / mu
           B <- t + (1/t) - 2
           term1 <- -1 / sigma
           term2 <- (1 / (sigma^3)) * B
           result <- term1 + term2
           
           return(-result * result)
         },
         
         d2ldmdd = function(y, sigma, mu) {
           t <- y / mu
           
           term1_m <- -1 / mu
           term2_m <- (t + 3) / (2 * mu * (t + 1))
           term3_m <- -(1 / (2 * sigma^2 * mu)) * (1/t - t)
           dldm <- term1_m + term2_m + term3_m
           
           B <- t + (1/t) - 2
           term1_d <- -1 / sigma
           term2_d <- (1 / (sigma^3)) * B
           dldd <- term1_d + term2_d
           
           d2ldmdd <- -dldm * dldd
           return(d2ldmdd)
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS4(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS4", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(mean(y), length(y))}),
         sigma.initial = expression({sigma <- rep(0.5, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
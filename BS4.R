BS4 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "BS4", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "BS4", substitute(sigma.link),
                      c("log", "logit", "probit", "cloglog", "cauchit", "own"))
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
         dldm = function(y, mu, sigma) {
           result <- (y / (sigma + mu * y)) + sigma - (mu * y)
           return(result)
         },
         
         dldd = function(y, sigma, mu) {
           result <- (1 / (sigma + mu * y)) + mu - (sigma / y)
           return(result)
         },
         
         # Second derivatives
         
         d2ldm2 = function(y, sigma, mu) {
           result <- (y / (sigma + mu * y)) + sigma - (mu * y)
           return(-result * result)
         },
         
         d2ldd2 = function(y, sigma, mu) {
           result <- (1 / (sigma + mu * y)) + mu - (sigma / y)
           return(-result * result)
         },
         
         d2ldmdd = function(y, sigma, mu) {
           
           dldm <- (y / (sigma + mu * y)) + sigma - (mu * y)

           dldd <- (1 / (sigma + mu * y)) + mu - (sigma / y)
           
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

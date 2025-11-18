BS5 <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "BS5", substitute(mu.link),
                      c("log", "inverse", "identity", "own"))
  dstats <- checklink("sigma.link", "BS5", substitute(sigma.link),
                      c("log", "logit", "probit", "cloglog", "cauchit", "own"))
  structure(
    list(family = c("BS5", "Birnbaum-Saunders - four parameterization"),
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
           
           term1 <- sigma / (y * (sigma + 1) + sigma * mu)
           term2 <- -1 / (2 * mu)
           term3 <- (y * (sigma + 1)) / (4 * mu^2)
           term4 <- -(sigma^2) / (4 * y * (sigma + 1))
           
           result <- (term1 + term2 + term3 + term4)
           return(result)
         },
         
         dldd = function(y, mu, sigma) { 
           
           term1 <- 0.5 + (1 / (2 * (sigma + 1)))
           term2 <- mu / ((y * (sigma + 1) + sigma * mu) * (sigma + 1))
           term3 <- -y / (4 * mu)
           term4 <- -(mu * (sigma^2 + 2 * sigma)) / (4 * y * (sigma + 1)^2)
           
           result <- (term1 + term2 + term3 + term4)
           return(result)
         },
         
         # Second derivatives
         
         d2ldm2 = function(y, mu, sigma) {
           
           term1 <- sigma / (y * (sigma + 1) + sigma * mu)
           term2 <- -1 / (2 * mu)
           term3 <- (y * (sigma + 1)) / (4 * mu^2)
           term4 <- -(sigma^2) / (4 * y * (sigma + 1))
           dldm <- term1 + term2 + term3 + term4
           
           return(-dldm * dldm) 
         },
         
         d2ldd2 = function(y, mu, sigma) {
           term1 <- 0.5 + (1 / (2 * (sigma + 1)))
           term2 <- mu / ((y * (sigma + 1) + sigma * mu) * (sigma + 1))
           term3 <- -y / (4 * mu)
           term4 <- -(mu * (sigma^2 + 2 * sigma)) / (4 * y * (sigma + 1)^2)
           dldd <- term1 + term2 + term3 + term4
           
           return(-dldd * dldd)
         },
         
         d2ldmdd = function(y, mu, sigma) {

           # dldm
           dm1 <- sigma / (y * (sigma + 1) + sigma * mu)
           dm2 <- -1 / (2 * mu)
           dm3 <- (y * (sigma + 1)) / (4 * mu^2)
           dm4 <- -(sigma^2) / (4 * y * (sigma + 1))
           dldm <- dm1 + dm2 + dm3 + dm4
           
           # dldd
           dd1 <- 0.5 + (1 / (2 * (sigma + 1)))
           dd2 <- mu / ((y * (sigma + 1) + sigma * mu) * (sigma + 1))
           dd3 <- -y / (4 * mu)
           dd4 <- -(mu * (sigma^2 + 2 * sigma)) / (4 * y * (sigma + 1)^2)
           dldd <- dd1 + dd2 + dd3 + dd4
           
           return(-dldm * dldd)
         },
         
         
         G.dev.incr = function(y,mu,sigma,...) -2*dBS5(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS5", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial    = expression({mu    <- rep(mean(y), length(y))}),
         sigma.initial = expression({sigma <- rep(0.5, length(y)) }),
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

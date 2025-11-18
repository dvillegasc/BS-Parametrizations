#Verificacion de la dBS4

integrate(dBS4, lower=0, upper=99, mu=1.7, sigma=2.3) #1 with absolute error < 2.8e-06

#Verificacion de las derivadas


#dBS4 <- function(x, mu, sigma, log = FALSE) {
  # Validaciones
  #if (any(x <= 0)) stop("x must be positive")
  #if (any(mu <= 0)) stop("mu (mu_A) must be positive")
  #if (any(sigma <= 0)) stop("sigma (lambda_A) must be positive")
  
  #sqrt_x <- sqrt(x)
  #term_A <- (sigma / (x * sqrt_x)) + (mu / sqrt_x)
  #term_B <- (sigma * sqrt_x) - (mu / sqrt_x)
  #K <- 1 / (2 * sqrt(2 * pi)) # Constante
  
  #if (log) {
    #log_pdf <- log(K) + log(term_A) - 0.5 * (term_B^2)
    #return(log_pdf)
  #} else {
   # pdf <- K * term_A * exp(-0.5 * (term_B^2))
    #return(pdf)
  #}
#}


#Derivadas manuales

library(gamlss)


dldm_manual <- function(y, mu, sigma) {
  result <- (y / (sigma + mu * y)) + sigma - (mu * y)
  return(result)
}

dldd_manual <- function(y, sigma, mu) {
  result <- (1 / (sigma + mu * y)) + mu - (sigma / y)
  return(result)
}


# --- Second derivates ---

d2ldm2_manual <- function(y, sigma, mu) {
  deriv1 <- dldm_manual(y, mu, sigma)
  return(-deriv1 * deriv1)
}

d2ldd2_manual <- function(y, sigma, mu) {
  deriv1 <- dldd_manual(y, mu, sigma)
  return(-deriv1 * deriv1)
}

d2ldmdd_manual <- function(y, sigma, mu) {
  dldm <- dldm_manual(y, mu, sigma)
  dldd <- dldd_manual(y, mu, sigma)
  return(-dldm * dldd)
}

#Derivadas computacionales

    dldm_compu <- function(y, mu, sigma) {
      
      dm <- gamlss::numeric.deriv(
        expr = dBS4(y, mu, sigma, log = TRUE), 
        theta = "mu",                          
        delta = 1e-04)
      
      # Extrae el gradiente
      dldm <- as.vector(attr(dm, "gradient"))
      return(dldm)
    }
    
    
    dldd_compu <- function(y, mu, sigma) {
      
      ds <- gamlss::numeric.deriv(
        expr = dBS4(y, mu, sigma, log = TRUE), 
        theta = "sigma",                       
        delta = 1e-04)
      
      dldd <- as.vector(attr(ds, "gradient"))
      return(dldd)
    }
    
    #Segundas derivadas compu
    d2ldm2_compu <- function(y, mu, sigma) {
      
      dm <- gamlss::numeric.deriv(
        expr = dBS4(y, mu, sigma, log = TRUE), 
        theta = "mu",
        deltha= 1e-04)
      
      d2ldm2 <- as.vector(attr(dm, "hessian"))[1, 1]
      return(d2ldm2)
    }
    
    
    d2ldd2_compu <- function(y, mu, sigma) {
      
      ds <- gamlss::numeric.deriv(
        expr = dBS4(y, mu, sigma, log = TRUE), 
        theta = "sigma")
      
      d2ldd2 <- as.vector(attr(ds, "hessian"))[1, 1]
      return(d2ldd2)
    }
    
    
    
    
    
    # PRUEBA
    
    y_test     <- c(1, 2, 5, 15)
    mu_test    <- 0.7
    sigma_test <- 0.75
    
    cat("--- Verificación de dldm (derivada de mu) ---\n")
    manual_mu <- dldm_manual(y = y_test, mu = mu_test, sigma = sigma_test)
    compu_mu  <- dldm_compu(y = y_test, mu = mu_test, sigma = sigma_test)
    
    print(data.frame(y = y_test, manual = manual_mu, computacional = compu_mu))
    
    
    cat("\n--- Verificación de dldd (derivada de sigma) ---\n")
    manual_sigma <- dldd_manual(y = y_test, mu = mu_test, sigma = sigma_test)
    compu_sigma  <- dldd_compu(y = y_test, mu = mu_test, sigma = sigma_test)
    
    print(data.frame(y = y_test, manual = manual_sigma, computacional = compu_sigma))
    

    
#-------------------------- Validation Familia ----------------------------------
    
    
    
n <- 1000
    
# True parameters are:
true_mu <- 4
true_si <- 5
    
y <- rBS4(n=n, mu=true_mu, sigma=true_si)
    
library(gamlss)
mod <- gamlss(y ~ 1, family=BS4,
              n.cyc = 100)
    
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)
exp(confint(mod)) # Intervalos de confianza reales
    
hist(y, freq = FALSE,
     xlim = c(0, 3), ylim = c(0, 10),
     main = "Histogram of Data-I",
     xlab = "y", ylab = "f(y)",
     col = "burlywood1",
     border = "darkgoldenrod4")
    
    
curve(dBS4(x, mu = 4.5, sigma= 5), from = 0.0000001, to = 10,
      add= TRUE,
      ylim = c(0, 1.7),
      col = "navy",       
      lwd = 2,            
      las = 1,
      type= "l",
      ylab = "f(x)",      
      xlab = "x")          
    
# Curvas adicionales
curve(dBS4(x, mu = 2, sigma= 5), add = TRUE, col = "deepskyblue", type= "l", lwd = 2)
    
curve(dBS4(x, mu = 3, sigma= 5), add = TRUE, col = "gold", type= "l", lwd = 2)
    
curve(dBS4(x, mu = 3.5, sigma= 5), add = TRUE, col = "red", type= "l", lwd = 2)
    
curve(dBS4(x, mu = 4, sigma= 5), add = TRUE, col = "#F28E2B", type= "l", lwd = 2)

curve(dBS4(x, mu = 4.5, sigma= 5), add = TRUE, col = "#F96F9B", type= "l", lwd = 2)
    

legend("topright",
       col = c("navy", "deepskyblue", "gold", "red", "#F28E2B","#F96F9B"),
       lty = 1,
       bty="n",
       cex = 0.9,       
       legend = c("μA = 1.5","μA = 2", "μA = 3", "μA = 3.5", "μA = 4", "μA = 4.5"))
    
    
    
    
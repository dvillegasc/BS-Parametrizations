#Verificacion de la dBS5

integrate(dBS5, lower=0, upper=99, mu=1.7, sigma=2.3) #1 with absolute error < 5.4e-05

#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)


dldm_manual <- function(y, mu, sigma) {
  
  term1 <- sigma / (y * (sigma + 1) + sigma * mu)
  term2 <- -1 / (2 * mu)
  term3 <- (y * (sigma + 1)) / (4 * mu^2)
  term4 <- -(sigma^2) / (4 * y * (sigma + 1))
  
  result <- (term1 + term2 + term3 + term4)
  return(result)
}

dldd_manual <- function(y, mu, sigma) { 
  
  term1 <- 0.5 + (1 / (2 * (sigma + 1)))
  term2 <- mu / ((y * (sigma + 1) + sigma * mu) * (sigma + 1))
  term3 <- -y / (4 * mu)
  term4 <- -(mu * (sigma^2 + 2 * sigma)) / (4 * y * (sigma + 1)^2)
  
  result <- (term1 + term2 + term3 + term4)
  return(result)
}


# --- Second derivates ---

d2ldm2_manual <- function(y, mu, sigma) {
  
  term1 <- sigma / (y * (sigma + 1) + sigma * mu)
  term2 <- -1 / (2 * mu)
  term3 <- (y * (sigma + 1)) / (4 * mu^2)
  term4 <- -(sigma^2) / (4 * y * (sigma + 1))
  dldm <- term1 + term2 + term3 + term4
  
  return(-dldm * dldm) 
}

d2ldd2_manual <- function(y, mu, sigma) {
  term1 <- 0.5 + (1 / (2 * (sigma + 1)))
  term2 <- mu / ((y * (sigma + 1) + sigma * mu) * (sigma + 1))
  term3 <- -y / (4 * mu)
  term4 <- -(mu * (sigma^2 + 2 * sigma)) / (4 * y * (sigma + 1)^2)
  dldd <- term1 + term2 + term3 + term4
  
  return(-dldd * dldd)
}

d2ldmdd_manual <- function(y, mu, sigma) {
  
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
}

#Derivadas computacionales

dldm_compu <- function(y, mu, sigma) {
  
  dm <- gamlss::numeric.deriv(
    expr = dBS5(y, mu, sigma, log = TRUE), 
    theta = "mu",                          
    delta = 1e-04)
  
  # Extrae el gradiente
  dldm <- as.vector(attr(dm, "gradient"))
  return(dldm)
}


dldd_compu <- function(y, mu, sigma) {
  
  ds <- gamlss::numeric.deriv(
    expr = dBS5(y, mu, sigma, log = TRUE), 
    theta = "sigma",                       
    delta = 1e-04)
  
  dldd <- as.vector(attr(ds, "gradient"))
  return(dldd)
}

#Segundas derivadas compu
d2ldm2_compu <- function(y, mu, sigma) {
  
  dm <- gamlss::numeric.deriv(
    expr = dBS5(y, mu, sigma, log = TRUE), 
    theta = "mu",
    deltha= 1e-04)
  
  d2ldm2 <- as.vector(attr(dm, "hessian"))[1, 1]
  return(d2ldm2)
}


d2ldd2_compu <- function(y, mu, sigma) {
  
  ds <- gamlss::numeric.deriv(
    expr = dBS5(y, mu, sigma, log = TRUE), 
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
true_mu <- 1
true_si <- 5

y <- rBS5(n=n, mu=true_mu, sigma=true_si)

library(gamlss)
mod <- gamlss(y ~ 1, family=BS5,
              n.cyc = 100)

exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)

#------------------------ Grafica 1 ------------------------------------

curve(dBS5(x, mu = 1, sigma= 2), from = 0.0000001, to = 2,
      #add= TRUE,
      ylim = c(0, 3),
      col = "deepskyblue",       
      lwd = 2,            
      las = 1,
      type= "l",
      ylab = "f(x)",      
      xlab = "x")          

curve(dBS5(x, mu = 1, sigma= 5), add = TRUE, col = "gold", type= "l", lwd = 2)

curve(dBS5(x, mu = 1, sigma= 10), add = TRUE, col = "red", type= "l", lwd = 2)

curve(dBS5(x, mu = 1, sigma= 25), add = TRUE, col = "#F28E2B", type= "l", lwd = 2)

curve(dBS5(x, mu = 1, sigma= 50), add = TRUE, col = "#F96F9B", type= "l", lwd = 2)

curve(dBS5(x, mu = 1, sigma= 100), add = TRUE, col = "navy", type= "l", lwd = 2)


legend("topright",
       col = c("deepskyblue", "gold", "red", "#F28E2B","#F96F9B", "navy"),
       lty = 1,
       bty="n",
       cex = 0.9,       
       legend = c("δ = 2","δ = 5", "δ = 10", "δ = 25", "δ = 50", "δ = 100"))
  

#------------------------ Grafica 2 ------------------------------------



curve(dBS5(x, mu = 1, sigma= 5), from = 0.0000001, to = 4,
      #add= TRUE,
      ylim = c(0, 1),
      col = "deepskyblue",       
      lwd = 2,            
      las = 1,
      type= "l",
      ylab = "f(x)",      
      xlab = "x")          

curve(dBS5(x, mu = 1.5, sigma= 5), add = TRUE, col = "gold", type= "l", lwd = 2)

curve(dBS5(x, mu = 2, sigma= 5), add = TRUE, col = "red", type= "l", lwd = 2)

curve(dBS5(x, mu = 2.5, sigma= 5), add = TRUE, col = "#F28E2B", type= "l", lwd = 2)

curve(dBS5(x, mu = 3, sigma= 5), add = TRUE, col = "#F96F9B", type= "l", lwd = 2)

curve(dBS5(x, mu = 3.5, sigma= 5), add = TRUE, col = "navy", type= "l", lwd = 2)


legend("topright",
       col = c("deepskyblue", "gold", "red", "#F28E2B","#F96F9B", "navy"),
       lty = 1,
       bty="n",
       cex = 0.9,       
       legend = c("μ = 1","μ = 1.5", "μ = 2", "μ = 2.5", "μ = 3", "μ = 3.5"))

#-------------------------------- Grafica 3 --------------------------------

varBS5 <- function(mu, sigma) {
  numerador <- (mu^2) * (2 * sigma + 5)
  denominador <- (sigma + 1)^2
  return(numerador / denominador)
}


mu <- 2
sigma <- seq(0, 20, length.out = 200) 

var_values <- varBS5(mu = mu, sigma = sigma)


plot(sigma, var_values, 
     type = "l",           
     lwd = 2,              
     ylim = c(0, 20),      
     xlim = c(0, 20),      
     xlab = expression(delta),  # Símbolo griego delta
     ylab = "Var[T]",      
     main = "",            
     las = 1)             


legend(x= 8, y= 11,
       lty = 1,
       bty="n",
       cex = 0.9,       
       legend = expression(mu == 2))


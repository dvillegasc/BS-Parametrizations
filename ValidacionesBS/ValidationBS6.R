#Verificacion de la dBS6

integrate(dBS6, lower=0, upper=99, mu=1.7, sigma=2.3) #1 with absolute error < 0.00011

#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)


dldm_manual = function(y, mu, sigma) {
  b <- (2 * mu) / (2 + sigma^2)
  db_dm <- 2 / (2 + sigma^2)
  
  term1 <- (1 / (y + b)) * db_dm
  term2 <- -1 / (2 * mu)
  term3 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_dm
  
  result <- term1 + term2 + term3
  return(result)
}

dldd_manual = function(y, mu, sigma) { 
  b <- (2 * mu) / (2 + sigma^2)
  db_ds <- -(4 * mu * sigma) / ((2 + sigma^2)^2)
  
  term1 <- -1 / sigma
  term2 <- (1 / sigma^3) * ((y / b) + (b / y) - 2)
  term3 <- (1 / (y + b)) * db_ds
  term4 <- (-1 / (2 * b)) * db_ds
  term5 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_ds
  
  result <- term1 + term2 + term3 + term4 + term5
  return(result)
}


# --- Second derivates ---

d2ldm2_manual <- function(y, mu, sigma) {
  b <- (2 * mu) / (2 + sigma^2)
  db_dm <- 2 / (2 + sigma^2)
  
  term1 <- (1 / (y + b)) * db_dm
  term2 <- -1 / (2 * mu)
  term3 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_dm
  
  dldm <- term1 + term2 + term3
  
  return(-dldm * dldm) 
}

d2ldd2_manual = function(y, mu, sigma) {
  b <- (2 * mu) / (2 + sigma^2)
  db_ds <- -(4 * mu * sigma) / ((2 + sigma^2)^2)
  
  term1 <- -1 / sigma
  term2 <- (1 / sigma^3) * ((y / b) + (b / y) - 2)
  term3 <- (1 / (y + b)) * db_ds
  term4 <- (-1 / (2 * b)) * db_ds
  term5 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_ds
  
  dldd <- term1 + term2 + term3 + term4 + term5
  
  return(-dldd * dldd)
}

d2ldmdd_manual = function(y, mu, sigma) {
  # Auxiliares
  b <- (2 * mu) / (2 + sigma^2)
  db_dm <- 2 / (2 + sigma^2)
  db_ds <- -(4 * mu * sigma) / ((2 + sigma^2)^2)
  
  # dldm
  m1 <- (1 / (y + b)) * db_dm
  m2 <- -1 / (2 * mu)
  m3 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_dm
  dldm <- m1 + m2 + m3
  
  # dldd
  d1 <- -1 / sigma
  d2 <- (1 / sigma^3) * ((y / b) + (b / y) - 2)
  d3 <- (1 / (y + b)) * db_ds
  d4 <- (-1 / (2 * b)) * db_ds
  d5 <- (1 / (2 * sigma^2)) * ((y / b^2) - (1 / y)) * db_ds
  dldd <- d1 + d2 + d3 + d4 + d5
  
  return(-dldm * dldd)
}


#Derivadas computacionales

dldm_compu <- function(y, mu, sigma) {
  
  dm <- gamlss::numeric.deriv(
    expr = dBS6(y, mu, sigma, log = TRUE), 
    theta = "mu",                          
    delta = 1e-04)
  
  # Extrae el gradiente
  dldm <- as.vector(attr(dm, "gradient"))
  return(dldm)
}


dldd_compu <- function(y, mu, sigma) {
  
  ds <- gamlss::numeric.deriv(
    expr = dBS6(y, mu, sigma, log = TRUE), 
    theta = "sigma",                       
    delta = 1e-04)
  
  dldd <- as.vector(attr(ds, "gradient"))
  return(dldd)
}

#Segundas derivadas compu
d2ldm2_compu <- function(y, mu, sigma) {
  
  dm <- gamlss::numeric.deriv(
    expr = dBS6(y, mu, sigma, log = TRUE), 
    theta = "mu",
    deltha= 1e-04)
  
  d2ldm2 <- as.vector(attr(dm, "hessian"))[1, 1]
  return(d2ldm2)
}


d2ldd2_compu <- function(y, mu, sigma) {
  
  ds <- gamlss::numeric.deriv(
    expr = dBS6(y, mu, sigma, log = TRUE), 
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

y <- rBS6(n=n, mu=true_mu, sigma=true_si)

library(gamlss)
mod <- gamlss(y ~ 1, family=BS6,
              n.cyc = 100)

exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)

#------------------------ Grafica 1 ------------------------------------

curve(dBS6(x, mu = 2, sigma= 0.1), from = 0.0000001, to = 3,
      #add= TRUE,
      ylim = c(0, 2),
      col = "deepskyblue",       
      lwd = 2,            
      las = 1,
      type= "l",
      ylab = "f(x)",      
      xlab = "x")          

curve(dBS6(x, mu = 2, sigma= 0.3), add = TRUE, col = "gold", type= "l", lwd = 2)

curve(dBS6(x, mu = 2, sigma= 0.5), add = TRUE, col = "red", type= "l", lwd = 2)

curve(dBS6(x, mu = 2, sigma= 0.75), add = TRUE, col = "#F28E2B", type= "l", lwd = 2)

curve(dBS6(x, mu = 2, sigma= 1), add = TRUE, col = "#F96F9B", type= "l", lwd = 2)

curve(dBS6(x, mu = 2, sigma= 1.5), add = TRUE, col = "navy", type= "l", lwd = 2)


legend("topright",
       col = c("deepskyblue", "gold", "red", "#F28E2B","#F96F9B", "navy"),
       lty = 1,
       bty="n",
       cex = 0.9,       
       legend = c("α = 0.1","α = 0.3", "α = 0.5", "α = 0.75", "α = 1", "α = 1.5"))


#------------------------ Grafica 2 ------------------------------------



curve(dBS6(x, mu = 0.75, sigma= 0.1), from = 0.0000001, to = 4,
      #add= TRUE,
      ylim = c(0, 1),
      col = "deepskyblue",       
      lwd = 2,            
      las = 1,
      type= "l",
      ylab = "f(x)",      
      xlab = "x")          

curve(dBS6(x, mu = 1, sigma= 0.1), add = TRUE, col = "gold", type= "l", lwd = 2)

curve(dBS6(x, mu = 1.5, sigma= 0.1), add = TRUE, col = "red", type= "l", lwd = 2)

curve(dBS6(x, mu = 2, sigma= 0.1), add = TRUE, col = "#F28E2B", type= "l", lwd = 2)

curve(dBS6(x, mu = 2.5, sigma= 0.1), add = TRUE, col = "#F96F9B", type= "l", lwd = 2)

curve(dBS6(x, mu = 3, sigma= 0.1), add = TRUE, col = "navy", type= "l", lwd = 2)


legend("topright",
       col = c("deepskyblue", "gold", "red", "#F28E2B","#F96F9B", "navy"),
       lty = 1,
       bty="n",
       cex = 0.9,       
       legend = c("μ = 0.1", "μ = 1","μ = 1.5", "μ = 2", "μ = 2.5", "μ = 3"))

#-------------------------------- Grafica 3 --------------------------------

varBS6 <- function(mu, sigma) {
  numerador <- ((mu*sigma)^2 * (4+5*sigma^2))
  denominador <- ((2 + sigma^2)^2)
  return(numerador / denominador)
}


mu <- 2
sigma <- seq(0, 25, length.out = 200) 

var_values <- varBS6(mu = mu, sigma = sigma)


plot(sigma, var_values, 
     type = "l",           
     lwd = 2,              
     ylim = c(0, 20),      
     xlim = c(0, 20),      
     xlab = expression(alpha),  # Símbolo griego delta
     ylab = "Var[T]",      
     main = "",            
     las = 1)             


legend(x= 8, y= 11,
       lty = 1,
       bty="n",
       cex = 0.9,       
       legend = expression(mu == 2))


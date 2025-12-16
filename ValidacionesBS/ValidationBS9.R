#Verificacion de la dBS9

integrate(dBS9, lower=0, upper=999, mu=10, sigma=1.5) #1 with absolute error < 0.00011

#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)


dldm_manual <- function(y, mu, sigma) {
  a0 <- sqrt(2 * (sigma - 1))
  b0 <- mu / sigma
  db_dm <- 1 / sigma
  
  term1 <- (1 / (y + b0)) * db_dm
  term2 <- -1 / (2 * b0) * db_dm
  term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  return(term1 + term2 + term3)
}

dldd_manual <- function(y, mu, sigma) { 
  a0 <- sqrt(2 * (sigma - 1))
  b0 <- mu / sigma
  
  da_ds <- 1 / a0
  db_ds <- -b0 / sigma
  
  term1 <- (-1 / a0) * da_ds
  term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  term3 <- (1 / (y + b0)) * db_ds
  term4 <- (-1 / (2 * b0)) * db_ds
  term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
  
  return(term1 + term2 + term3 + term4 + term5)
}


# --- Computacionales ---

dldm_compu <- function(y, mu, sigma) {
  dm <- gamlss::numeric.deriv(expr = dBS9(y, mu, sigma, log = TRUE), theta = "mu", delta = 1e-04)
  dldm <- as.vector(attr(dm, "gradient"))
  return(dldm)
}

dldd_compu <- function(y, mu, sigma) {
  ds <- gamlss::numeric.deriv(expr = dBS9(y, mu, sigma, log = TRUE), theta = "sigma", delta = 1e-04)
  dldd <- as.vector(attr(ds, "gradient"))
  return(dldd)
}





# PRUEBA

y_test     <- c(1, 2, 5, 15)
mu_test    <- 10 # Varianza
sigma_test <- 1.5 # Gamma (>1)

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

y <- rBS9(n=n, mu=true_mu, sigma=true_si)

library(gamlss)
mod <- gamlss(y ~ 1, family=BS9,
              n.cyc = 100)

exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)

#------------------------ Grafica 1 ------------------------------------

curve(dBS9(x, mu = 10, sigma= 1.05), from = 0.0000001, to = 25,
      #add= TRUE,
      ylim = c(0, 0.25),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS9(x, mu = 10, sigma= 1.1), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS9(x, mu = 10, sigma= 1.3), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS9(x, mu = 10, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS9(x, mu = 10, sigma= 1.7), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS9(x, mu = 10, sigma= 1.9), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("γ = 1.05","γ = 1.1", "γ = 1.3", "γ = 1.5", "γ = 1.7", "γ = 1.9"))


#------------------------ Grafica 2 ------------------------------------


curve(dBS9(x, mu = 5, sigma= 1.5), from = 0.0000001, to = 4,
      #add= TRUE,
      ylim = c(0, 0.85),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS9(x, mu = 10, sigma= 1.5), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS9(x, mu = 15, sigma= 1.5), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS9(x, mu = 20, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS9(x, mu = 25, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS9(x, mu = 30, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("σ² = 5", "σ² = 10","σ² = 15", "σ² = 20", "σ² = 25", "σ² = 30"))

#-------------------------------- Grafica 3 --------------------------------

varBS9 <- function(mu, sigma) {
  numerador <- ((mu^2) * (sigma - 1) * (5*sigma - 3))
  denominador <- sigma^2
  return(numerador / denominador)
}


mu <- 2
sigma <- seq(1.05, 10, length.out = 100) 

var_values <- varBS9(mu = mu, sigma = sigma)


plot(mu, var_values, 
     type = "l",           
     lwd = 2,              
     ylim = c(0, 20),      
     xlim = c(1, 250),      
     xlab = expression(mu[A]),  # Símbolo griego delta
     ylab = "Var[T]",      
     main = "",            
     las = 1)             


legend(x= 2.2, y= 1.7,
       lty = 1,
       bty="n",
       cex = 0.9,       
       legend = expression(lambda[A] == 2))

















#------------------------ Grafica 1 ------------------------------------

curve(dBS9(x, mu = 2, sigma= 1.1), from = 0.0000001, to = 3,
      ylim = c(0, 0.9),
      col = "deepskyblue",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(x)",      
      xlab = "x",
      main = "")          

curve(dBS9(x, mu = 2, sigma= 1.3), add = TRUE, col = "gold", type= "l", lwd = 2)
curve(dBS9(x, mu = 2, sigma= 1.5), add = TRUE, col = "red", type= "l", lwd = 2)
curve(dBS9(x, mu = 2, sigma= 1.7), add = TRUE, col = "#F28E2B", type= "l", lwd = 2)
curve(dBS9(x, mu = 2, sigma= 1.9), add = TRUE, col = "#F96F9B", type= "l", lwd = 2)
curve(dBS9(x, mu = 2, sigma= 2.1), add = TRUE, col = "navy", type= "l", lwd = 2)

legend("topright",
       col = c("deepskyblue", "gold", "red", "#F28E2B","#F96F9B", "navy"),
       lty = 1, bty="n", cex = 0.9,        
       legend = c("φ = 1.1","φ = 1.3", "φ = 1.5", "φ = 1.7", "φ = 1.9", "φ = 2.1"))


#------------------------ Grafica 2 ------------------------------------

curve(dBS9(x, mu = 0.75, sigma= 1.5), from = 0.0000001, to = 2,
      ylim = c(0, 1.35),
      col = "deepskyblue",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(x)",      
      xlab = "x",
      main = "(b) Variando Media (Phi=2)")          

curve(dBS9(x, mu = 1, sigma= 1.5), add = TRUE, col = "gold", type= "l", lwd = 2)
curve(dBS9(x, mu = 1.25, sigma= 1.5), add = TRUE, col = "red", type= "l", lwd = 2)
curve(dBS9(x, mu = 1.5, sigma= 1.5), add = TRUE, col = "#F28E2B", type= "l", lwd = 2)
curve(dBS9(x, mu = 1.75, sigma= 1.5), add = TRUE, col = "#F96F9B", type= "l", lwd = 2)
curve(dBS9(x, mu = 2, sigma= 1.5), add = TRUE, col = "navy", type= "l", lwd = 2)

legend("topright",
       col = c("deepskyblue", "gold", "red", "#F28E2B","#F96F9B", "navy"),
       lty = 1, bty="n", cex = 0.9,        
       legend = c("μ = 0.75", "μ = 1","μ = 1.25", "μ = 1.5", "μ = 1.75", "μ = 2"))

#-------------------------------- Grafica 3 --------------------------------

varBS9 <- function(mu, sigma) {
  numerador <- ((mu^2) * (sigma - 1) * (5*sigma - 3))
  denominador <- sigma^2
  return(numerador / denominador)
}

mu <- 2
sigma <- seq(1.01, 250, length.out = 500) 

var_values <- varBS9(mu = mu, sigma = sigma)

# 3. Graficar
plot(sigma, var_values, 
     type = "l",            
     lwd = 2,              
     ylim = c(0, 21),       # Límite en Y un poco más arriba de 20
     xlim = c(0, 250),      # Límite en X igual a la imagen
     xlab = expression(phi), 
     ylab = "Var[T]",      
     main = "",            
     las = 1)

# Linea de la asintota teórica: 5 * mu^2
# Para mu=2, el limite es 5*(2^2) = 20
abline(h = 5 * mu^2, col="red", lty=2)

text(125, 10, labels = expression(mu == 2), cex = 1.2)

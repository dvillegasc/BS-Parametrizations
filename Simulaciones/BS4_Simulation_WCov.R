# To perform the simulation -----------------------------------------------
if (!dir.exists("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls")) {
  dir.create("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls")
}


library("parSim")

parSim(
  ### SIMULATION CONDITIONS
  n = c(200, 400, 600, 800, 1000),
  mu = c(2, 3.5, 4.5),
  sigma = c(0.4, 5),
  
  reps = 1000,                                # repetitions
  write = TRUE,                               # Writing to a file
  name = "C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls/BS4_Sim_WCov",  # Name of the file
  nCores = 6,                                 # Number of cores to use
  progressbar = TRUE,               # Progress bar
  export = c("dBS", "pBS", "qBS", "rBS", "hBS", "dBS4", "pBS4", "qBS4", "rBS4", "hBS4", "BS4"),
  
  expression = {
    library(gamlss)
    library(gamlss2)
    
    # True parameter values
    y <- rBS4(n=n, mu, sigma)
    
    f   <- y ~ 1
    mod <- try(suppressMessages(gamlss2(f, family=BS4)), silent = TRUE)
    
    if (class(mod)[1] == "try-error") {
      mu_hat    <- NA
      sigma_hat <- NA
    }
    else {
      mu_hat <- exp(coef(mod, what="mu")["(Intercept)"])
      sigma_hat <- exp(coef(mod, what="sigma")["(Intercept)"])
    }
    
    # Results list:
    Results <- list(
         mu_hat = mu_hat,
      sigma_hat = sigma_hat
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------

archivos <- list.files(pattern = "^BS4_Sim_WCov.*\\.txt$", 
                       path="C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls",
                       full.names = TRUE)

archivos

lista_datos <- lapply(archivos, read.table, header = TRUE, 
                      sep = "", stringsAsFactors = FALSE)
datos <- do.call(rbind, lista_datos)

datos$case <- with(datos, ifelse(mu==2 & sigma==0.4, 1, 
                                 ifelse(mu==2 & sigma==5, 2,
                                        ifelse(mu==3.5 & sigma==0.4, 3,
                                               ifelse(mu==3.5 & sigma==5, 4, 
                                                      ifelse(mu==4.5 & sigma==0.4, 5, 6))))))
datos$case <- as.factor(datos$case)

# To analize the results --------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

trim <- 0.10 # percentage of values to be trimmed

dat <- datos %>% group_by(n, mu, case) %>% 
  summarise(nobs = n(),
            mean_mu = mean(mu_hat, trim=trim, na.rm=TRUE),
            mean_si = mean(sigma_hat, trim=trim, na.rm=TRUE),
            mse_mu = mean((mu_hat - mu)^2, trim=trim, na.rm=TRUE),
            mse_si = mean((sigma_hat - sigma)^2, trim=trim, na.rm=TRUE),
            bias_mu = mean(mu_hat-mu, trim=trim, na.rm=TRUE),
            bias_si = mean(sigma_hat-sigma, trim=trim, na.rm=TRUE),
  )

dat


# Plots

if (!dir.exists("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Figs")) {
  dir.create("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Figs")
}


p1 <- ggplot(dat, aes(x=n, y=bias_mu, colour=case)) +
  geom_line() + 
  ylab(expression(paste("Bias for ", mu))) +
  ylim(min(dat$bias_mu), 0.1)

p2 <- ggplot(dat, aes(x=n, y=bias_si, colour=case)) +
  geom_line() + 
  ylab(expression(paste("Bias for ", sigma)))+
  ylim(min(dat$bias_mu), 0.1)

p1_final <- p1 + theme_bw(base_size = 13)
p2_final <- p2 + theme_bw(base_size = 13)


ggsave(filename="C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Figs/bias_BS5_Sim_WCov.pdf", width=12, height=6,
       plot=p1_final+p2_final)


p3 <- ggplot(dat, aes(x=n, y=mse_mu, colour=case)) +
  geom_line() + 
  ylab(expression(paste("MSE for ", mu)))

p4 <- ggplot(dat, aes(x=n, y=mse_si, colour=case)) +
  geom_line() + 
  ylab(expression(paste("MSE for ", sigma)))

p3_final <- p3 + theme_bw(base_size = 13)
p4_final <- p4 + theme_bw(base_size = 13)

ggsave(filename="C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Figs/mse_BS5_Sim_WCov.pdf", width=12, height=6,
       plot=p3_final+p4_final)


# Tables

trim <- 0.10 # percentage of values to be trimmed

dat <- datos %>% group_by(n, mu) %>% 
  summarise(nobs = n(),
            mean_mu = mean(mu_hat, trim=trim, na.rm=TRUE),
            ab_mu = mean(abs(mu_hat-mu), trim=trim, na.rm=TRUE),
            mse_mu = mean((mu_hat - mu)^2, trim=trim, na.rm=TRUE)
  )

dat


dat |> filter(mu == 1, sigma == 1) |> 
  select(mean_mu, mean_si, ab_mu, ab_si, mse_mu, mse_si) -> a
a[, -1]

library(xtable)
xtable(a[, -1])



dat_summary <- datos %>% 
  group_by(n, mu, case) %>%
  summarise(
    mean_mu = mean(mu_hat, trim=trim, na.rm=TRUE),
    mse_mu = mean((mu_hat - mu)^2, trim=trim, na.rm=TRUE),
    bias_mu = mean(mu_hat - mu, trim=trim, na.rm=TRUE),
    .groups = 'drop'
  )

tabla_comparativa <- dat_summary %>%
  filter(near(mu, 0.25)) %>%
  select(n, "θ̂" = mean_mu, "Bias" = bias_mu, "MSE" = mse_mu)

print(tabla_comparativa)
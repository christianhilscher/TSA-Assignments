rm(list=ls())
library(dict)
library(ggplot2)
library(reshape)
############################################
# Defining all functions I need later
############################################

# Function for shifting the variables back - i.e. get e_{t-1} from e_t
shift <- function(x, lag) {
  n <- length(x)
  xnew <- rep(NA, n)
  if (lag < 0) {
    xnew[1:(n-abs(lag))] <- x[(abs(lag)+1):n]
  } else if (lag > 0) {
    xnew[(lag+1):n] <- x[1:(n-lag)]
  } else {
    xnew <- x
  }
  return(xnew)
}

# Function for getting x (one sample)
get_x <- function(length, a, beta){
  # Need more epsilons since I have two lags
  epsilon <- rnorm(length+2)
  # Calculating e
  e_tmp <- epsilon - (5/6)*shift(epsilon,1) + (1/6)*shift(epsilon,2)
  e <- e_tmp[3:length(e_tmp)]
  
  #Calculating u
  u_initial <- e[1]
  u <- u_initial
  for (i in seq(2,ts_length)){
    u_t <- a*u[i-1] + e[i]
    u <- c(u, u_t)
  }
  
  # Calculating x
  x <- beta_0 + u
  return(x)
}

get_coef_se <- function(x){
  x_t1 <- shift(x,1)
  x_t1 <- x_t1[2:length(x_t1)]
  
  y <-  x[2:length(x)]
  ols_frame <- data.frame(y, x_t1)
  modl <- lm(y ~ x_t1, data=ols_frame)
  sum <- summary(modl)
  
  d <- dict()
  #Adding the coefficient as well as SE
  d[['coef_and_se']] <- sum$coefficients[2,1:2]
  #Now including varaince of errors; needed for strategy 2 and 3
  d[['gamma_e']] <- sum$sigma
  d[['e']] <- sum$residuals
  return(d)
}

k <- function(z){
  c <- (6*pi*z)/5
  value <- 3/(c**2) * (sin(c)/(c)-cos(c))
  return(value)
}

covs <- function(x, h){
  x_th <- shift(x,h)[(h+1):length(shift(x,h))]
  x_t <- x[(h+1):length(x)]
  
  
  return(cov(x_t, x_th))
}

get_omega <- function(e, strategy){
  
  gamma_0 <- var(e)
  T <- length(e)
  sum_T <- length(e) - 2
  
  if (strategy == "strategy2"){
    q <- floor(T**(4/5))
  } else {
    q <- floor(T**(1/3))
  }
  
  sum_overh <- 0
  for (h in seq(sum_T)){
    z <- (h/(q+1))
    value <-  k(z)*covs(e,h)*2
    sum_overh <- sum_overh + value
  }
  
  omega_hat <- (T/(T-2))*(gamma_0 + sum_overh)
  return(omega_hat)
}

strategy1 <- function(a, a1_hat, a1_hat_se){
  t_ratio <- (a1_hat-a)/a1_hat_se
  return(t_ratio)
}

strategies_lrvar <- function(a, a1_hat, a1_hat_se, gamma_e, e, strategy){
  omega_hat <- get_omega(e, strategy)
  
  correction <- sqrt(gamma_e/omega_hat)
  secondterm <- ((omega_hat-gamma_e)/(2*sqrt(omega_hat))) * ((1000*a1_hat_se)/(gamma_e))
  
  t_ratio <- correction * (a1_hat-a)/a1_hat_se - secondterm
  return(t_ratio)
}

get_tratios <- function(x_mat, a, sample_size){
  t_ratios1 <- vector()
  t_ratios2 <- vector()
  t_ratios3 <- vector()
  
  for (sample in seq(sample_size)){
    
    # Frist getting the coefficients, gamma_e and the errors themselves
    d <- get_coef_se(x=x_mat[sample,])
    
    t1 <- strategy1(a=1, 
                    a1_hat=d[['coef_and_se']][1], 
                    a1_hat_se=d[['coef_and_se']][2])
    t_ratios1 <- cbind(t_ratios1, t1)
    
    t2 <- strategies_lrvar(a=1, 
                           a1_hat = d[['coef_and_se']][1], 
                           a1_hat_se = d[['coef_and_se']][2],
                           gamma_e = d[['gamma_e']],
                           e = d[['e']],
                           strategy="strategy2")
    t_ratios2 <- cbind(t_ratios2, t2)
    
    t3 <- strategies_lrvar(a=1, 
                           a1_hat = d[['coef_and_se']][1], 
                           a1_hat_se = d[['coef_and_se']][2],
                           gamma_e = d[['gamma_e']],
                           e = d[['e']],
                           strategy="strategy3")
    t_ratios3 <- cbind(t_ratios3, t3)
  }
  
  out_mat <- matrix(c(t_ratios1, 
                      t_ratios2, 
                      t_ratios3), nrow=3, 
                    ncol=1000, byrow = TRUE)
  return(out_mat)
}

calculate_P <- function(ratio_mat, critical_val){
  samples <- dim(ratio_mat)[1]
  p_list <- vector()
  
  for (sample in seq(samples)){
    values <- ratio_mat[sample,]
    p <- sum((values<(critical_val))/length(values))
    p_list <- cbind(p_list, p)
  }
  return(p_list)
}

############################################
# Setting the parameters and calulating
############################################


ts_length <- 1000
sample_size <- 1000
a <- 1
beta_0 <- 1

x_mat <- matrix(NA, sample_size, ts_length)

# Filling up the matrix; each row is sample
for (sample in seq(sample_size)){
  x_mat[sample,] <- get_x(ts_length, a, beta_0)
}

ratios <- get_tratios(x_mat, 1, sample_size)
calculate_P(ratios, (-2.86))




 plotting <- function(e){
   
   T <- length(e) + 1
   sum_T <- length(e) - 1
   
   horizon <- seq(1, sum_T)
   df = data.frame('horizon' = horizon)
   
   for (i in cbind('strategy2', 'strategy3')){
     if (i == "strategy2"){
       q <- floor(T**(4/5))
     } else {
       q <- floor(T**(1/3))
     }
     
     cov_list <- vector()
     kernel_list <- vector()
     sum_overh <- 0
     for (h in seq(sum_T)){
       z <- (h/(q+1))
       value <-  k(z)
       kernel_list <- rbind(kernel_list, value)
       
       cov <- covs(e, h-1)
       cov_list <- rbind(cov_list, cov)
     }
     df[i] <- kernel_list
   }
   df['cov'] <- cov_list
   return(df)
 }
   
 results <- get_coef_se(x_mat[109,])
 res <- results[['e']]
 values <- seq(1, 998)
 
 df1 <- plotting(res)
 
p1 = ggplot(data=df1, aes(x=horizon)) + 
   geom_line(aes(y=strategy2), color="#0072B2") + 
   geom_line(aes(y=strategy3), color="#D55E00") + 
   geom_point(aes(y=cov), alpha=0.2, colour="black", size=0.2) +
   xlab('Horizon') + 
   ylab('Values') + 
   ggtitle('Kernel values and covariances')

print(p1)

random



df1$cov
covs(res,1)



results[['gamma_e']]**2

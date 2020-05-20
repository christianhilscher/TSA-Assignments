rm(list=ls())
library(ggplot2)
library(tidyr)
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

spectral_dens <- function(q, lambda, Mod=0){
  i = complex(real = 0, imaginary = 1)
  
  total_l <- 0
  for (l in seq((-q), q)){
    value <- exp(i*lambda*(-l))
    total_l <- total_l + value
  }
  
  total_j <- 0
  for (j in seq((-q), q)){
    value <- exp(i*lambda*j)
    total_j <- total_j + value
  }
  
  T_c <- 1/(2*q + 1)**2 * (total_j*total_l)
  spectral_d <- T_c * 1/(2*pi)
  
  if (Mod==0){
    return(spectral_d)
  } else if (Mod==1){
    return(Mod(spectral_d))
  }
  else{
    stop('Please provide either 0 or 1 as argument for Mod')
  }
}


f_lambda <- function(lambdas, q_list){
  df <- data.frame(lambdas)
  for (q in q_list){
    
    sp_densitites <- mapply(spectral_dens, q=q, lambda=lambdas)
    num_sp_densitites <- Mod(sp_densitites)
    
    name <- "density_q"
    fullname <- paste(name, toString(q), sep="")
    
    df[[fullname]] <- num_sp_densitites
  }
  return(df)
}

VR <- function(lambdas, q_list){
  df <- data.frame(lambdas)
  for (q in q_list){
    s_lambda0 <- vector()
    for (lambda in lambdas){
      s_current <- 2* integrate(spectral_dens, 0, lambda, q=q, Mod=1)$value
      s_lambda0 <- rbind(s_lambda0, s_current)
    }
    
    s_pi <- 2* integrate(spectral_dens, 0, pi, q=q, Mod=1)$value
    
    name <- "VR_q"
    fullname <- paste(name, toString(q), sep="")
    
    results <- s_lambda0/s_pi
    df[[fullname]] <- results
  }
  
  return(df)
}

prepare_plot <- function(df){
  df_plot <- pivot_longer(df, -lambdas, names_to="q", values_to="values")
  return(df_plot)
}
############################################

points <- 500
lambdas <- seq(from=0, 
               to=pi, 
               length.out = points)
q_list = c(1, 3, 10)


#Defining color palette
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df_f <- f_lambda(lambdas, q_list)
df_VR <- VR(lambdas, q_list)

df_f_plot <- prepare_plot(df_f)
df_VR_plot <- prepare_plot(df_VR)

p1=ggplot(data=df_f_plot, aes(x=lambdas, y=values, colour=q)) + 
  geom_line(alpha=0.8, size=0.7) + 
  scale_colour_manual(values = cbPalette) +
  xlab('Lambda') + 
  ylab(expression(f[x](lambda))) + 
  ggtitle("Values of the spectal density function over the interval [0,pi]")

p2=ggplot(data=df_VR_plot, aes(x=lambdas, y=values, colour=q)) + 
  geom_line(alpha=0.8, size=0.7) + 
  scale_colour_manual(values = cbPalette) +
  xlab('Lambda') + 
  ylab(expression(VR(lambda[0],q))) +
  ggtitle("Values of the VR as a function of lambda between [0,pi]")

print(p1)
print(p2)
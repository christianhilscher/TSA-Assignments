rm(list=ls())
library(dict)
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

generate_x <- function(tslength, errors, b){
  
  e_shifted <- shift(errors,1)
  x <- (errors+b*e_shifted)[2:1001]
  return(x)
}

DGP <- function(tslength, samplesize, b, models){
  
  storage <- dict()
  for (model in models){
    storage[[model]] <- vector()
    
    for (i in seq(1, samplesize)){
      if (model=="i"){
        errors <- rnorm(1001)
      }else if (model=="ii"){
        errors <- rt(1001, 5)
      }else if(model=="iii"){
        errors <- runif(1001)
      }else{
        stop('Please provide one of the models: i, ii, iii')
      }
      
      storage[[model]] <- rbind(storage[[model]], 
                                generate_x(tslength, errors, b))
    }
  }
  return(storage)
}

delta <- function(b, lambdaj){
  val <- (1 + b**2 + 2*b*cos(lambdaj))
  return(val)
}

Ix <- function(lambdaj, x){
  
  bigT <- length(x)
  little_t_seq <- seq(bigT)-1
  
  prefactor <- 1/(2*pi*bigT)
  
  
  sequence_sin <- cos(lambdaj * little_t_seq)
  sequence_cos <- cos(lambdaj * little_t_seq)
  
  val <- prefactor*(sum(x*sequence_cos)**2 + sum(x*sequence_sin)**2)
  return(val)
}

sigma <- function(lambda_list, b, x){
  
  bigM <- length(lambda_list)
  
  sumoverIx <- sapply(lambda_list, Ix, x=x)
  sumoverdelta <- delta(b=b, lambda_list)
  
  val <- (1/bigM) * sum((sumoverIx/sumoverdelta))
  return(val)
}

objf <- function(b, lambda_list, x){
  
  firstterm <- delta(b, lambda_list)
  secondterm <- sigma(lambda_list, b, x)
  
  thirdterm_Ix <- sapply(lambda_list, Ix, x)
  
  tmp <- -log(firstterm) - log(secondterm) - (thirdterm_Ix/firstterm)*(1/secondterm)
  val <- sum(tmp)
  
  return(val)
}
#######################################################
TS_length <- 1000
Sample_size <- 1000
b <- 0.25
models <- cbind('i', 'ii', 'iii')
sigma_list <- cbind(1, 5/3, 1/12)
#Grid to choose the optimal b from
b_list <- seq(-1,1, by=0.05)

# Initialising final storage space
df_final <- data.frame()

# Generating the data
data <- DGP(TS_length, Sample_size, b, models)
M <- floor((TS_length-1)/2)
lambdas <- (2*pi/TS_length) * seq(1:M)


# Looping through the models
for (model in models){
  pos <- which(models==model)
  
  data_used <- data[[model]]
  sample_size <- dim(data_used)[1]
  ts_length <- dim(data_used)[2]
  
  df <- data.frame()
  
  # Looping through the samples
  for (s in seq(sample_size)){
    
    res <- sapply(b_list, objf, lambda_list=lambdas, x=data_used[s,])
    
    b_opt <- b_list[which.max(res)]
    sigma_opt <- sigma(lambdas, b_opt, x=data_used[s,]) * 2 * pi
    
    df <- rbind(df, cbind(b_opt, sigma_opt))
  }
  colnames(df) <- c("b", "sigma")
  
  bias_b <- (mean(df$b) - b)*ts_length
  bias_sigma <- (mean(df$sigma) - sigma_list[pos])*ts_length
  
  # Getting the MSE
  sigma_b <- var(df$b)
  sigma_sigma <- var(df$sigma)
  
  df_final <- rbind(df_final, cbind(bias_b, bias_sigma, sigma_b, sigma_sigma))
  print(model)
  
}

df_i <- df_final[1,]
df_ii <- df_final[2,]
df_iii <- df_final[3,]

# Getting the MSE for b and sigma
calc_MSE <- function(df){
  val_b <- ((df$bias_b/1000)**2 + df$sigma_b)*1000
  val_sigma <- ((df$bias_sigma/1000)**2 + df$sigma_sigma)*1000
  
  return(c(val_b, val_sigma))
}

calc_MSE(df_i)


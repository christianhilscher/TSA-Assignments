rm(list=ls())
library(ggplot2)
library(tidyr)
library(dict)
library(tidyverse)

############################################
# Defining all functions I need later
############################################

f <- function(x, z){
  value <- x**(z-1)*exp(-x)
  return(value)
}

intf <- function(z){
  if (z != 0){
    value <- integrate(f, lower=0,upper=Inf, z=z)$value
  }else{
    value <- integrate(f, lower=1e-3,upper=Inf, z=z)$value
  }
  
  return(value)
}

psi <- function(j, d){
  
  if (d<0){
    upper <- j**(-d-1)
    lower <- intf(abs(d))
  }else{
    upper <- j**(d-1)
    lower <- intf(d)
  }
  
  value <- upper/lower
  return(value)
}


get_data <- function(overall_length, d_list){
  y <- matrix(nrow=overall_length, ncol=length(d_list))
  for (i in seq_along(d_list)){
    res <- rep(NA, length(overall_length))
    
    for (t in seq_len(overall_length)){
      val <- 0
      for (j in seq_len(t-1)){
        tmp <- psi(j,d_list[i]) * et[t-j]
        val <- val + tmp
      }
      res[t] <- sum(val)
    }
    
    y[,i] <- res
  }
  return(y)
}

get_gamma0 <- function(data, d_list, d){
  
  pos <- which(d_list==d)
  sigma2 <- var(data[,pos])
  upper <- intf(1-2*d)
  lower <- intf(1-d)**2
  
  final <- sigma2*(upper/lower)
  return(final)
}

get_corrs <- function(data, laglength, d_list){
  
  sto = dict()
  for (i in seq_along(d_list)){
    lagframe <- data.frame()
    
    gamma0 <- get_gamma0(data, d_list, d_list[i])
    lagframe <- cbind(0, gamma0)
    
    
    for (h in seq(laglength)){
      pre <- (h-1+d_list[i])/(h-d_list[i])
      gammah <- pre * lagframe[h,2]
      
      lagframe <- rbind(lagframe, cbind(h, gammah))
    }
    
    lagframe[,2] <- lagframe[,2]/lagframe[1,2]
    colnames(lagframe) <- c("h", 'real')
    
    estimated <- acf(data[,i], lag.max=laglength, 
                     plot=FALSE, demean = FALSE)$acf
    final <- cbind(lagframe, estimated)
    sto[[toString(i)]] <- data.frame(final)
  }
  
  return(sto)
}

prepplot <- function(data){
  df_plot <- pivot_longer(data, -h, names_to = 'type', values_to = 'values')
  return(df_plot)
}
############################################
############################################
############################################

M <- 2000
TS <- 1000
overall_length <- M+TS
et <- rnorm(M+TS)

d_list <- c(-0.45, -0.25, 0, 0.25, 0.45)

y <- get_data(overall_length, d_list)
yshort <- y[M:(TS+M), ]

yshort <- load(file="y_short.Rda")
# Dictionary with 5 dataframes for each d
a <- get_corrs(yshort, 25, d_list)


d_list <- c(-0.45, -0.25, 0, 0.25, 0.45)
yshort <- read.csv2("yshort.csv")
a <- get_corrs(yshort, 25, d_list)


cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
for (key in a$keys()){
  df <- prepplot(a[[key]])
  spec <- as.numeric(key)
  
  d_spec <- d_list[spec]
  name <- ("Autocorrelations with d=")
  
  p <- ggplot(data=df, aes(x=h, y=values, fill=type)) + 
    geom_col(position="dodge") + 
    scale_colour_manual(values = cbPalette) + 
    ggtitle(paste(name, d_spec, sep=""))
  plot(p)
}


df_y <- data.frame(yshort)
colnames(df_y) <- c("X", "d=-0.45", "d=-0.25", "d=0", "d=0.25", "d=0.45")
df_y_plot <- pivot_longer(df_y, -X, names_to = 'type', values_to = 'values')

for (key in a$keys()){
  spec <- as.numeric(key)
  
  d_spec <- d_list[spec]
  namehighl <- ("d=")
  name <- ("Time series of y with d=")
  
  highl <- paste(namehighl, d_spec, sep="")
  
  z <- mutate(df_y_plot, highlight=ifelse(type==highl, highl, "Other"))
  p2 <- ggplot(data=z, aes(x=X, y=values, group=type, color=highlight)) +  
    geom_line() + 
    scale_color_manual(values = c("#009E73", "lightgrey")) +
    scale_size_manual(values=c(1.5,0.02)) +
    theme(legend.position="none") +
    ggtitle(paste(name, d_spec, sep=""))
  plot(p2)
}

d_spec <- d_list[1]
namehighl <- ("d=")
highl <- paste(namehighl, d_spec, sep="")
highl

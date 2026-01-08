################################################
## RV estimation under microstructure effects ##
################################################

# This script illustrates the impact of market microstructure noise on the latent 
# price process, specifically focusing on bid-ask bounces. For this purpose we simulate a 
# latent price path following a GBM. This path is then distorted by adding artificial 
# bid-ask bounces. We follow by generating Volatility Signature Plots to demonstrate how to find an
# sampling frequency.



rm(list=ls())

# Load packages
library(tidyverse)

################################################################################


# We start by defining the function to simulate the GBM price paths
set.seed(2026)


FUN_GBM <-  function(N) {
  
  mu <- 0.05 # drift component (constant)
  sigma <- 0.2 # volatility (constant)
  
  
  # create iid innovations
  innovations <- rnorm(N, mean = 0, sd = 1)
  
  # cumulate innovations 
  cum_innovations <- cumsum(innovations)
  
  # Define Time Grid
  
  j <- 1:N
  time_vector <- j / N
  
  
  # Compute exponent of geometric brownian motion
  expo <- (mu - 0.5*sigma^2) * time_vector + (sigma/sqrt(N))*cum_innovations
  
  # Compute Price
  P_0 <- 25
  
  P_j_N <- P_0*exp(expo)
  
  # Create dataframe
  data <-data.frame(
    Time = c( time_vector),
    P_j_N = c( P_j_N),
    N_steps = as.integer(N ))
  
  return(data)
  
}

################################################################################
# Simulate GBM with 390 observations corresponding to minute prices on a trading day

simulate_day <- function(n = 390, P0 = 25, mu = 0.05, sigma = 0.2, tick = 1) {
  
  innovations <- rnorm(n)
  cum_innovations <- cumsum(innovations)
  time_vector <- 1:n / n
  expo <- (mu - 0.5*sigma^2)*time_vector + (sigma/sqrt(n))*cum_innovations
  P_j_n <- P0 * exp(expo)
  
  # Bid/Ask
  Bid <- floor(P_j_n / tick) * tick
  Ask <- ceiling(P_j_n) * tick
  
  # Observed Price
  indicator <- sample(c(0,1), n, replace = TRUE)
  P_obs <- ifelse(indicator == 1, Bid, Ask)
  
  data.frame(t = 1:n, P_j_n = P_j_n, P_obs = P_obs)
}

df<- simulate_day(n = 390, P0 = 25, mu = 0.05, sigma = 0.2, tick = 1)



ggplot(df, aes(x = t)) +
  geom_line(aes(y = P_j_n, color = "True Price"), linewidth = 1) +
  geom_line(aes(y = P_obs, color = "Observed Price"), linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = c(
    "True Price" = "midnightblue",
    "Observed Price" = "darkorange"
  ), name = NULL) +
  labs(
    x = "Intraday Observations",
    y = "Price",
    color = "Legend",
    title = "Simulated GBM and observed Bid-Ask Prices"
  ) +
  theme_minimal(base_size = 14)


# create plot for paper
p1<-ggplot(df, aes(x = t)) +
  geom_line(aes(y = P_j_n, color = "True Price"), linewidth = 1) +
  geom_line(aes(y = P_obs, color = "Observed Price"), linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = c(
    "True Price" = "midnightblue",
    "Observed Price" = "darkorange"
  ), name = NULL) +
  labs(
    x = "Intraday Observations",
    y = "Price",
    color = "Legend",
  ) +
  theme_minimal(base_size = 14)+
  theme(
    strip.text = element_text(face = "bold", color = "black", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.text.y = element_text(face = "bold", color = "black", size = 11),
    axis.title.y = element_text(face = "bold", color = "black", size = 13),
    axis.title.x = element_text(face = "bold", color = "black", size = 13),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.position = "none"
  )

ggsave(
  filename = "BID-ASK-Plot.png",  
  plot = p1,                       
  width = 9, height = 3,          
  dpi = 600                       
)


################################################################################
# We now create the volatility signature plot 

# Define function to compute Realize Variances
compute_RV <- function(price_series, delta) {
  log_ret <- diff(log(price_series[seq(1, length(price_series), by = delta)]))
  sum(log_ret^2)
}

#We simulate 500 trading days with 390 observations each day and compute the Realized
#Variances for each day using different sampling frequencies

days <- 500       # number of day
n <- 390          # min per day
freq <- c(1,5,10,15,30,60,120)# Sampling frequencies
RV_mat <- matrix(0, nrow = days, ncol = length(freq))

for (d in 1:days) {
  day_data <- simulate_day(n = n, tick = 1 )
  for (f in 1:length(freq)) {
    RV_mat[d,f] <- compute_RV(day_data$P_obs, freq[f])
  }
}

# Take average 
RV_avg <- colMeans(RV_mat)
signature_plot_data <- data.frame(freq = freq, RV = RV_avg)

ggplot(signature_plot_data, aes(x = freq, y = RV)) +
  geom_line(linetype = "dashed",color = "darkgrey", size = 1.5) +
  geom_point(color = "darkred", size = 3) +
  labs(x = "Sampling Frequency in minutes", 
  y = "Average Realized Variance", 
  title = "Volatility Signature Plot") +
  theme_minimal(base_size = 14)



p2<-ggplot(signature_plot_data, aes(x = freq, y = RV)) +
  geom_line(linetype = "dashed",color = "darkgrey", size = 1.5) +
  geom_point(color = "darkred", size = 3) +
  labs(x = "Sampling Frequency", 
  y = "Average RV")+
  theme_minimal(base_size = 14)+
  theme(
    strip.text = element_text(face = "bold", color = "black", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.text.y = element_text(face = "bold", color = "black", size = 11),
    axis.title.y = element_text(face = "bold", color = "black", size = 13),
    axis.title.x = element_text(face = "bold", color = "black", size = 13),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
  )

ggsave(
  filename = "Vola-Signature-Plot.png",  
  plot = p2,                      
  width = 9, height = 3,          
  dpi = 600 )                      


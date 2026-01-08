###########################################
# Simulation of Geometric Brownian Motion##
###########################################

#In this script we present the simulation conducted in Section 2.3 of the seminar paper, 
#where a Geometric Brownian Motion is simulated to illustrate how a continuous-time process 
#can be approximated in discrete time. By increasing the number of discrete steps $N$ on the interval $[0,1]$, 
#the simulated price path becomes smoother and more closely resembles a continuous process.  
#This demonstrates how continuous-time models can be made practically usable for numerical analysis and simulation-based studies.

################################################################################

rm(list=ls()) 

library(tidyverse)    # to organize data

set.seed(2025) # set seed to reproduce results

# Create Function simulating GBM price path
FUN_GBM <-  function(N) {
  
  mu <- 0.05 # drift component (constant)
  sigma <- 0.3 # volatility (constant)
  
  
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
  P_j_N <- exp(expo)
  
  # Create dataframe
  data <-data.frame(
    Time = c( time_vector),
    P_j_N = c( P_j_N),
    N_steps = as.integer(N ))
  
  return(data)
  
}


# set sample sizes
sample_size = c(10,100,1000,100000)
results_list <- list() # create list to store results

# loop over sample sizes to create price paths
for(n in sample_size){
  
  # call function 
  GBM <- FUN_GBM(n)
  
  # store result
  results_list[[as.character(n)]] <-GBM
  
}

# Safe results and prepare for plotting
GBM_data <- bind_rows(results_list)
GBM_data <- GBM_data %>%
  mutate(
    Label = paste0("N = ", N_steps)
  )


# Create Plot for price path
ggplot(GBM_data, aes(x = Time, y = P_j_N)) +
  geom_point(color = "darkblue", alpha = 0.8, size = 1) +
  facet_wrap(~ Label, ncol = 2, scales = "free_y") +
  labs(
    title = "Prices following a Geometric Brownian Motion",
    x = "Time",
    y = expression(P[j/N])
  ) +
  theme_minimal(base_size = 12)



# Create Plot for returns (for illustration only, not included in the paper)
GBM_returns <- GBM_data %>% 
  filter(N_steps == 1000) %>%  
  mutate(return = log(P_j_N) - lag(log(P_j_N))) 


ggplot(GBM_returns, aes(x = Time, y = return)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "GBM Log Returns",
    x = "Time",
    y = "Return"
  )


# Create Plot used in Paper
GBM_data <- GBM_data %>%
  mutate(PointSize = case_when(
    N_steps <= 100 ~ 2,
    N_steps <= 1000 ~ 0.7,
    N_steps == 100000 ~ 0.3   
    
  ))

p <- ggplot(GBM_data, aes(x = Time, y = P_j_N)) +
  geom_point(aes(size = PointSize), color = "darkblue", alpha = 0.7) +  
  scale_size_identity() +  
  facet_wrap(~ Label, ncol = 2, scales = "fixed") +
  labs(
    title = NULL,
    x = NULL,
    y =  expression(bold(P[t[j]]))
  ) +
  theme_light(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", color = "black", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.text.y = element_text(face = "bold", color = "black", size = 11),
    axis.title.y = element_text(face = "bold", color = "black", size = 17),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

ggsave(
  filename = "GBM-SIM Plot.png",  
  plot = p,                      
  width = 9, height = 5,          
  dpi = 600                       
)

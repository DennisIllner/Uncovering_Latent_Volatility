##############################################################
## Consistency of Realized Variance using GBM price path ##
##############################################################


#This section demonstrates the consistency of the Realized Variance ($RV$) estimator 
#within a Geometric Brownian Motion framework. We visualize the consistency property 
#by simulating 10 trading days. Each day represents a standard trading session from 09:30 to 16:00 (390 minutes).
#To approximate a continuous-time process, we use high-frequency discretization at the one-second level, 
#resulting in 23,400 observations per trading day. The objective is to compare the daily $RV$ estimates 
#across various sampling frequencies ranging from ultra-high-frequency (1 second) to low-frequency (60 minutes) 
#against the true, latent integrated variance


################################################################################
# We start by generating the price path

rm(list=ls())

library(tidyverse)


set.seed(2025)


mu <- 0.05 # drift component (constant)
sigma <- 0.2 # volatility (constant)


# Observations
price_obs_daily <- 60 * 390   #60 Seconds times 390 min each day
n_days <- 10
N <- price_obs_daily* n_days

# create iid innovations
innovations <- rnorm(N, mean = 0, sd = 1)

# cumulate innovations 
cum_innovations <- cumsum(innovations)


# Define Time Grid
j <- 1:N
time_vector <- (j / N) * n_days


# Compute Step Size bc. we do not have grid from [0,1]
dt <- n_days / N   # duration divided by number of obs

# Compute exponent of geometric brownian motion
expo <- (mu - 0.5*sigma^2) * time_vector + (sigma*sqrt(dt))*cum_innovations

# Compute Price
P_j_N <- exp(expo)

# Create dataframe
data <-data.frame(
  Time = c( time_vector),
  P_j_N = c( P_j_N),
  N_steps = as.integer(N ))


# Create Plot (not included in paper, shown here for illustration)
ggplot(data, aes(x = Time, y = P_j_N)) +
  geom_point(color = "darkblue", alpha = 0.8, size = 1) +
  labs(
    title = "Geometric Brownian Motion",
    x = "Time",
    y = expression(bold(P[j/N]))
  ) +
  
  theme_light(base_size = 12) +
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


# Now we gonna compute realized variances 

# Add Intraday times second frequencies for illsutration
intraday_times <- format(
  seq(
    from = as.POSIXct("2020-01-01 09:30:00"),
    to   = as.POSIXct("2020-01-01 15:59:59"),
    by   = "1 sec"
  ),
  "%H:%M:%S"  
)


data <- data %>%
  mutate(
    Day = rep(1:n_days, each = price_obs_daily),
    Time = rep(intraday_times, times = n_days))


# Create function to compute RV
FUN_RV <- function(k, dataset) {
  
  data %>%
    group_by(Day) %>%
    filter(((row_number() - 1) %% k) == 0) %>%  # filter according to sampling frequency
    mutate(log_ret = log(P_j_N) - lag(log(P_j_N))) %>% # compute log return
    summarise(RV = sum(log_ret^2, na.rm = TRUE)) %>%  # compute realized variance
    ungroup() %>% 
    mutate(k = k)
}




# Set sampling frequencies
k_value <- c(1, 60, 300, 900, 1800, 3600) # (1 second, 1min, 5min, 15min, 30min, 60min) 

results_list <- list() # initialize empty list to store reuslts
i <- 0 # initialize index counter to store results



for(k in k_value){
  
  i = i+1 # increase index counter
  rv <- FUN_RV(k, data) #call function
  results_list[[i]] <- rv # store results
}

data_final <- dplyr::bind_rows(results_list)



# Prepare Data for nice plotting
data_plot <- data_final %>%
  
  mutate(Day = factor(Day),  # to get categorical axis
         True_IV = sigma^2, # Assign true IV for reference
         # Assign Frequency and sort
         Frequency = factor(
           case_when(
             k == 1    ~ "1 sec",
             k == 60   ~ "1 min",
             k == 300  ~ "5 min",
             k == 900  ~ "15 min",
             k == 1800 ~ "30 min",
             k == 3600 ~ "60 min"),
           levels = c("60 min","30 min", "15 min", "5 min", "1 min", "1 sec"),
           ordered = TRUE
         )
  ) 


# creat plot
ggplot(data_plot, aes(x = Day, y = RV)) +
  
  # True integrated variance IV
  geom_point(aes(y = True_IV), 
             shape = 15,      # 15 = Quadrat
             color = "darkblue", 
             size = 4) +
  
  # reaized variances
  geom_point(color = "darkorange", size = 3, alpha = 1) +
  
  # create grid
  facet_wrap(~ Frequency, nrow = 5) + 
  
  # label
  labs(
    title = "Realized Variance and IV using different sampling frequencies",
    x = "day",
    y = "RV and IV",
    color = "",   
    shape = ""    
  ) +
  theme_bw(base_size = 14)


# Plot used for Paper
p<-ggplot(data_plot, aes(x = Day)) +
  
  # True IV 
  geom_point(aes(y = True_IV, shape = "True IV", color = "True IV"),
             size = 4, alpha = 1) +
  
  # realized variances
  geom_point(aes(y = RV, shape = "RV", color = "RV"),
             size = 3, alpha = 1) +
  
  # create grid
  facet_wrap(~ Frequency, nrow = 5) +
  
  # labels
  labs(
    #title = "Realized Variance and IV using different sampling frequencies",
    x = "Trading Day",
    y = " Daily Volatility",
    color = "",
    shape = ""
  ) +
  
  # Scale manually
  scale_color_manual(values = c("True IV" = "darkblue", "RV" = "darkorange")) +
  scale_shape_manual(values = c("True IV" = 15, "RV" = 16)) +
  
  theme_light(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", color = "black", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.text.y = element_text(face = "bold", color = "black", size = 11),
    axis.title.y = element_text(face = "bold", color = "black", size = 11),
    axis.title.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.position = "none"
  )

ggsave(
  filename = "Consistency-GBM-Plot.png",  
  plot = p,                       
  width = 9, height = 5,          
  dpi = 600                      
)


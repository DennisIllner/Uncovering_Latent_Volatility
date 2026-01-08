##################################################################
## Asymptotic Properties of RV using GARCH(1,1) diffusion-limit ##
##################################################################

#This script investigates the convergence and distribution of the realized 
#variance estimator. The following code starts by estimating a  GARCH(1,1)
#model to calibrate the parameters  for the GARCH(1,1) diffusion limit.

#We then follow by generating  sample paths using the Euler Maruyama discretization with 
#and examine the asymptotic properties of the realized variance estimator focusing on consistency
# and the asymptotic normality.
################################################################################


# clean R environment
rm(list = ls())



# Load packages
library(tidyverse)    
library(gridExtra)    
library(rugarch)      

################################################################################
# We start by loading and preparing the data for estimation
# Note that we used high-frequency data consisting of 1-minute prices


#Read high-frequency 1-minute prices - IBM
data <- read.csv("IBM.txt", header = FALSE)
colnames(data) <- c("Date", "Time", "Open", "High", "Low", "Close","Volume")
data$Date <- as.Date(data$Date, format = "%m/%d/%Y") # convert date format

# Check length of dataset
max_date <- max(data$Date, na.rm = TRUE)
min_date <- min(data$Date, na.rm = TRUE)


#Reduce TIME SPAN
data = data %>% filter(Date >="2006-01-01")


adj_close_price <- data %>% group_by(Date) %>%
  slice_tail(n = 1) %>%  # take last close price observed each day
  ungroup()

#daily returns
returns <- adj_close_price %>% 
  mutate(daily_returns = Close/lag(Close) - 1,
         daily_log_returns = log(Close) - log(lag(Close)),
         daily_squared_return = daily_log_returns^2) %>% 
  select(Date, daily_returns, daily_log_returns, daily_squared_return) %>% 
  drop_na(daily_log_returns)

################################################################################
# We define in-sample for estimation of GARCH(1,1)
estimation_start <- as.Date("2015-10-01")
estimation_end   <- as.Date("2022-09-30")

# Filtern in-sample for estimation
estimation_data <- returns %>%
  filter(Date >= estimation_start & Date <= estimation_end)

################################################################################
# Now we estimate the GARCH(1,1) model

# Define GARCH(1,1)
garch_1_1_spec = ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm")

garch_1_1_fit <- ugarchfit(spec = garch_1_1_spec, data = estimation_data$daily_log_returns)
garch_1_1_fit

# Get GARCH parameters
garch_coefs <- coef(garch_1_1_fit)
omega_garch <- garch_coefs["omega"]
alpha_garch <- garch_coefs["alpha1"]
beta_garch  <- garch_coefs["beta1"]

###############################################################################
# Define the mapping function to get the parameters of the GARCH(1,1) diffusion-limit


get_continuous_params <- function(omega, alpha, beta, m) {
  
  # Theta 
  theta_cont <- -m * log(alpha + beta)
  
  # Omega
  omega_cont <- m * omega / (1 - alpha - beta)
  
  # Lambda 
  numerator <- 2 * alpha * log(alpha + beta)^2 * (1-beta*(alpha+beta))
  
  term1 <- (1 - (alpha + beta)^2) * (1-beta)^2
  
  term2 <-  (-alpha*(1-beta*(alpha+beta)))
  
  term3 <- (6 * log(alpha+beta) + 2*log(alpha+beta)^2 + 4*(1-alpha-beta))
  
  denominator <- term1 + term2*term3
  
  lambda_cont = numerator / denominator
  
  
  return(list(theta=theta_cont, omega=omega_cont, lambda=lambda_cont))
}

# Call function
cont_params <- get_continuous_params(omega_garch, alpha_garch, beta_garch, m=1)

################################################################################
# Define function to simlate data of the calibrated GARCH(1,1) diffusion-limit 
# using the Euler-Discretization

simulate_fun <- function(cont_params, sim_days = 2500,
                         steps_per_day = 390, seed = 2025) {
  
  theta  <- cont_params$theta
  omega  <- cont_params$omega
  lambda <- cont_params$lambda
  
  dt <- 1 / steps_per_day
  
  sim_spot_var       <- numeric(sim_days * steps_per_day)
  sim_integrated_var <- numeric(sim_days)
  sim_prices <- numeric(sim_days * steps_per_day)
  
  
  
  set.seed(seed)
  
  sigma2 <- omega          # mean
  log_p  <- log(100)       # start at 100
  
  for (i in 1:sim_days) {
    
    w_price <- rnorm(steps_per_day)
    w_vola  <- rnorm(steps_per_day)
    
    intraday_sigma2 <- numeric(steps_per_day)
    
    
    for (t in 1:steps_per_day) {
      
      # price (log)
      log_p <- log_p + sqrt(sigma2) * sqrt(dt) * w_price[t]
      
      # volatility process
      sigma2 <- theta * omega *dt + sigma2*(1-theta*dt + (2*lambda*theta*dt)^0.5 * w_vola[t])
      
      # store vola
      intraday_sigma2[t] <- sigma2
      
      # store minute and vola
      idx <- (i - 1) * steps_per_day + t
      sim_prices[idx] <- log_p
      sim_spot_var[idx] <- sigma2
      
      
    }
    
    sim_integrated_var[i] <- sum(intraday_sigma2) * dt
    
    
  }
  
  list(
    integrated_var = sim_integrated_var,
    prices         = sim_prices,
    spot_var       = sim_spot_var
  )
  
}

################################################################################
# After defining the function, we simulate 2500 trading days from the calibrated GARCH(1,1) diffusion limit,
# with 390 observations per day corresponding to minute prices within a trading day

sim_days <- 2500
steps_per_day <- 390
seed = 1996
sim_data <- simulate_fun(cont_params, sim_days,
                         steps_per_day , seed )


#Time Grid
time_days <- (0:(sim_days * steps_per_day -1)) / steps_per_day

################################################################################
# Create Plot for Spot VOlatility


# Create Dataframe with Spot Volatility
df_spot <- data.frame(
  time = time_days,
  spot_var = sim_data$spot_var
)

# Create plot
ggplot(df_spot, aes(x = time, y = sqrt(spot_var)*100)) +
  geom_line(linewidth = 0.3) +
  labs(
    x = "Time (days)",
    y = "Spot volatility (%)",
    title = "Spot Volatility "
  ) +
  coord_cartesian(xlim = c(0, 2500)) +
  theme_light(base_size = 12) 
  
# Create Plot for Paper
p1<-ggplot(df_spot, aes(x = time, y = sqrt(spot_var)*100)) +
  geom_line(linewidth = 0.3) +
  labs(
    #x = "Time (days)",
    y = "Spot volatility (%)",
    #title = "Spot Volatility "
  ) +
  coord_cartesian(xlim = c(0, 2500)) +
  theme_light(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3),
    axis.text = element_text(face = "bold", color = "black", size = 20),
    axis.title.y = element_text(face = "bold", color = "black", size = 20),
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = "Spot-Vola-GARCH-DIFF.png",
  plot = p1,
  width = 9, height = 6,
  dpi = 600
)

################################################################################
# Create Plot for Integrated VOlatility

ggplot(
  data.frame(day = 1:2500,
             vol = 100 * sqrt(sim_data$integrated_var)),
  aes(x = day, y = vol)
) +
  geom_line() +
  labs(
    x = "Day",
    y = "Daily Integrated volatility (%)",
    title = "Daily Integrated volatility"
  ) +
  coord_cartesian(xlim = c(1, 2500)) +
  theme_light(base_size = 12)


# Create Plot for Paper
p2<-ggplot(
  data.frame(day = 1:2500,
             vol = 100 * sqrt(sim_data$integrated_var)),
  aes(x = day, y = vol)
) +
  geom_line() +
  labs(
    x = "Day",
    y = "Daily Integrated volatility (%)",
    title = "Daily Integrated Volatility"
  ) +
  coord_cartesian(xlim = c(1, 2500)) +
  theme_light(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3),
    axis.text = element_text(face = "bold", color = "black", size = 20),
    axis.title.y = element_text(face = "bold", color = "black", size = 20),
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = "IV-GARCH-DIFF.png",
  plot = p2,
  width = 9, height = 6,
  dpi = 600
)

################################################################################
# Create Plot for log price process

df_min_prices <- data.frame(
  time = time_days,
  min_prices = sim_data$prices
)

ggplot(df_min_prices, aes(x = time, y = min_prices)) +
  geom_line(linewidth = 0.3) +
  labs(
    x = "Time (days)",
    y = "log price",
    title = "log price process"
  ) +
  coord_cartesian(xlim = c(0, 2500)) +
  theme_light(base_size = 12)

p3<-ggplot(df_min_prices, aes(x = time, y = min_prices)) +
  geom_line(linewidth = 0.3) +
  labs(
    x = "Time (days)",
    y = "log price",
    title = "log price process"
  ) +
  coord_cartesian(xlim = c(0, 2500)) +
  theme_light(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3),
    axis.text = element_text(face = "bold", color = "black", size = 20),
    axis.title.y = element_text(face = "bold", color = "black", size = 20),
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = "log-Price-GARCH-DIFF.png",
  plot = p3,
  width = 9, height = 6,
  dpi = 600
)


################################################################################
#Create daily returns

# daily returns
opens  <- sim_data$prices[seq(1, sim_days*steps_per_day, by=steps_per_day)]
closes <- sim_data$prices[seq(steps_per_day, sim_days*steps_per_day, by=steps_per_day)]
daily_returns <- closes - opens

# Create plot for daily returns (plot not used in paper but shown for illustration)
ggplot(
  data.frame(day = 1:2500,
             daily_ret = daily_returns),
  aes(x = day, y = daily_ret)
) +
  geom_line(linewidth = 0.3) +
  labs(
    x = "Day",
    y = "Daily returns",
    title = "Daily returns"
  ) +
  coord_cartesian(xlim = c(1, 2500)) +
  theme_light(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3),
    axis.text = element_text(face = "bold", color = "black", size = 11),
    axis.title = element_text(face = "bold", color = "black", size = 13),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_blank(),
    legend.position = "none"
  )


################################################################################
# Create Plot for daily squared returns

# comüute r^2
daily_squared_returns <- daily_returns^2

# Create Plot
ggplot(
  data.frame(day = 1:2500,
             daily_ret2 = sqrt(daily_squared_returns)*100),
  aes(x = day, y = daily_ret2)
) +
  geom_line(linewidth = 0.3) +
  labs(
    x = "Day",
    y = "Daily squared returns (%)",
    title = "Daily squared returns"
  ) +
  coord_cartesian(xlim = c(1, 2500)) +
  theme_light(base_size = 12)

# Create plot used for paper
p4<-ggplot(
  data.frame(day = 1:2500,
             daily_ret2 = sqrt(daily_squared_returns)*100),
  aes(x = day, y = daily_ret2)
) +
  geom_line(linewidth = 0.3) +
  labs(
    x = "Day",
    y = "Daily squared returns (%)",
    title = "Daily squared returns"
  ) +
  coord_cartesian(xlim = c(1, 2500)) +
  #theme_minimal()
  theme_light(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3),
    axis.text = element_text(face = "bold", color = "black", size = 20),
    axis.title.y = element_text(face = "bold", color = "black", size = 20),
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = "squared-return-GARCH-DIFF.png",
  plot = p4,
  width = 9, height = 6,
  dpi = 600
)


################################################################################
# We now turn to create to compute RV for increasing sampling frequencies to
# show convergence of RV to IV


# Create function to compute RV for different sampling frequencies
realized_var <- function(prices, steps_per_day, n) {
  rv <- numeric(sim_days)
  
  for (i in 1:sim_days) {
    
    idx <- ((i-1)*steps_per_day + 1):(i*steps_per_day)
    
    p_day <- prices[idx]
    
    # use every n-th price
    p_sampled <- p_day[seq(1, length(p_day), by = n)]
    
    # compute intraday return
    r_n <- diff(p_sampled) 
    
    # compute RV
    rv[i] <- sum(r_n^2, na.rm = TRUE)
  }
  rv
}



n_vec <- c(1, 5,10, 15, 30, 60)  # initialize sampling frequencies

min_prices <- sim_data$prices

rv_list <- list() # initialize list to store RV

for (n in n_vec) {
  rv_list[[paste0("RV_", n)]] <- realized_var(min_prices, steps_per_day, n)
}


# store IV and
df_iv <- data.frame(
  Day = 1:sim_days,
  True_IV = 100 * sqrt(sim_data$integrated_var)
)

# change long format
rv_df <- data.frame(Day = 1:sim_days) %>%
  bind_cols(as.data.frame(rv_list)) %>%
  rename_with(~ paste0(gsub("RV_", "", .), "-min"), starts_with("RV_")) %>%
  pivot_longer(cols = -Day, names_to = "Frequency", values_to = "RV_Value")

# scale RV to plot them in %
rv_df$RV_Value <- 100 * sqrt(rv_df$RV_Value)

# Set order for grid
rv_df$Frequency <- factor(rv_df$Frequency, 
                          levels = c("60-min", "30-min", "15-min", "10-min", "5-min", "1-min"))

# Create Plot
ggplot(rv_df, aes(x = Day)) +
  geom_line(data = df_iv, aes(y = True_IV, color = "True IV"), size = 0.5, alpha = 1.7) +
  geom_line(aes(y = RV_Value, color = "RV"), size = 0.2, alpha = 0.75) +
  facet_wrap(~ Frequency, ncol = 2, scales = "fixed") +
  scale_color_manual(values = c("True IV" = "darkblue", "RV" = "darkorange")) +
  labs(x = "Days", y = " Daily Volatility (%)",title = "RV and IV series", color = "") +
  theme_light(base_size = 12)

# 4.Create Plot used for Paper
p5 <-ggplot(rv_df, aes(x = Day)) +
  geom_line(data = df_iv, aes(y = True_IV, color = "True IV"), size = 0.5, alpha = 1.7) +
  geom_line(aes(y = RV_Value, color = "RV"), size = 0.2, alpha = 0.75) +
  facet_wrap(~ Frequency, ncol = 2, scales = "fixed") +
  scale_color_manual(values = c("True IV" = "darkblue", "RV" = "darkorange")) +
  labs(x = "", y = " Daily Volatility (%)", color = "") +
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

ggsave(
  filename = "Consistency-Panel.png",
  plot = p5,
  width = 9, height = 6,
  dpi = 600
)

################################################################################
# We now simulate 2500 trading days from the calibrated GARCH(1,1) diffusion limit,
# with 23400 observations per day corresponding to second-level prices within a trading day

sim_days <- 2500
steps_per_day <- 23400

sim_data <- simulate_fun(cont_params, sim_days,
                         steps_per_day , seed = 1996)


time_days <- (0:(sim_days * steps_per_day -1)) / steps_per_day


# Price Process
df_sec_prices <- data.frame(
  time = time_days,
  sec_prices = sim_data$prices
)

# Initialize sampling frequencies
n_vec <- c(1, 60, 300, 900, 1800, 3600)
rv_list <- list()

for (n in n_vec) {
 
  rv_raw <- realized_var(df_sec_prices$sec_prices, steps_per_day, n)
  
  # create label
  label <- ifelse(n < 60, paste0(n, "-sec"), paste0(n/60, "-min"))
  rv_list[[label]] <- rv_raw
}

# change to long format
rv_df <- rv_list %>% 
  as_tibble() %>% 
  mutate(Day = 1:sim_days) %>% 
  pivot_longer(cols = -Day, names_to = "Frequency", values_to = "RV_Value")

# change to volatility measure
rv_df <- rv_df %>% 
  mutate(RV_Value = 100 * sqrt(RV_Value))

# IV Dataframe
df_iv <- data.frame(
  Day = 1:sim_days,
  True_IV = 100 * sqrt(sim_data$integrated_var)
)

# Join RV and IV to compute samplin error
density_df <- rv_df %>%
  left_join(df_iv, by = "Day") %>%
  mutate(Error = RV_Value - True_IV)



# Set order for grid
density_df$Frequency <- factor(density_df$Frequency, 
                               levels = c("60-min", "30-min", "15-min", "5-min", "1-min", "1-sec"))



ggplot(density_df, aes(x = Error)) +
  geom_density(fill = "lightblue", alpha = 0.5, color = "darkblue", linewidth = 0.75) +
  #geom_density( color = "darkblue", linewidth = 1) +
  # ncol = 1 erzwingt, dass die Plots untereinander stehen!
  facet_wrap(~ Frequency, ncol = 2, scales = "free_y") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  #coord_cartesian(xlim = c(-1, 1)) + 
  labs(
    x = " (RV - IV)",
    y = "Density",
    title = "Degeneration of measurement error",
  ) +
  theme_light(base_size = 12)

# Create plot paper
p6<- ggplot(density_df, aes(x = Error)) +
  geom_density(fill = "lightblue", alpha = 0.5, color = "darkblue", linewidth = 0.75) +
  #geom_density( color = "darkblue", linewidth = 1) +
  # ncol = 1 erzwingt, dass die Plots untereinander stehen!
  facet_wrap(~ Frequency, ncol = 2, scales = "free_y") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  #coord_cartesian(xlim = c(-1, 1)) + 
  labs(
    x = " (RV - IV)",
    y = "Density",
    #title = "Degeneration of measurement error",
    
  ) +
  theme_light(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", color = "black", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.text.y = element_text(face = "bold", color = "black", size = 11),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold", color = "black", size = 11),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.position = "none"
  )

ggsave(
  filename = "Error-Distribution-Panel.png",
  plot = p6,
  width = 9, height = 5,
  dpi = 600
)


################################################################################
# Now we examine the asymptotic normality property of the RV estimator

# set delta
dt_sim <- 1 / steps_per_day

# Compute Integrated Quarticitiy by squaring and summing spot variance
iq_daily <- colSums(matrix(sim_data$spot_var^2, nrow = steps_per_day)) * dt_sim

# Store IQ inf dataframe
df_iq <- data.frame(
  Day = 1:sim_days,
  True_IQ = iq_daily
)

# Create DF for standardized measurement error using 1 second sampling frequency
df_clt <- rv_df %>%
  filter(Frequency == "1-sec") %>%
  mutate(RV_raw = (RV_Value / 100)^2) %>%
  left_join(df_iv, by = "Day") %>%
  mutate(IV_raw = (True_IV / 100)^2) %>%
  left_join(df_iq, by = "Day") %>%
  mutate(
    Standardized_Error = (RV_raw - IV_raw) / sqrt(2 * dt_sim * True_IQ)
  )

#Create Plot
ggplot(df_clt, aes(x = Standardized_Error)) +
  # density of standardized error
  geom_density(fill = "lightblue", alpha = 0.5, color = "black") +
  
  # standard normal
  stat_function(fun = dnorm, 
                args = list(mean = 0, sd = 1), 
                color = "darkorange", linewidth = 1.2) +
  
  # Restrict x-axis to the relevant range
  coord_cartesian(xlim = c(-4, 4)) +
  
  labs(
    x = expression(paste("Standardized Error: ", frac(RV - IV, sqrt(2 * Delta * IQ)))),
    title = "Asymptotic Normality of standardized errors") +
  theme_light(base_size = 12) 




p7<-ggplot(df_clt, aes(x = Standardized_Error)) +
  # density of standardized error
  geom_density(fill = "lightblue", alpha = 0.5, color = "black") +
  
  # standard normal
  stat_function(fun = dnorm, 
                args = list(mean = 0, sd = 1), 
                color = "darkorange", linewidth = 1.2) +
  
  # Restrict x-axis to the relevant range
  coord_cartesian(xlim = c(-4, 4)) +
  
  labs(
    x = "Stanardized sampling error") +
  theme_light(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", color = "black", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.text.y = element_text(face = "bold", color = "black", size = 11),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold", color = "black", size = 11),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.position = "none")

ggsave(
  filename = "Scaled-Error-Distribution-Panel.png",
  plot = p7,
  width = 9, height = 4,
  dpi = 600
)







library(dplyr)
library(readr)
library(ggplot2)

library(depmixS4)
library(matrixcalc)

# Function to add multiple plots to one
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Set seed for reproduceable results
set.seed(1)

# Pre and post cutoff
cutoff <- as.Date("2017-01-01") 
post_cutoff <- as.Date("2018-01-12")

# Load data from csv (data manually collected from coinmarketcap.com)
df_data <- read_csv("bitcoin_hist.csv") %>%
  mutate(Date = as.Date(Date, "%b %d, %Y")) %>%
  arrange(Date) %>%
  mutate(logret = log(lead(Close)/Close))

# Filter post data
df_post_data <- df_data %>%
  filter(Date > post_cutoff,
         !is.na(logret)) %>%
  dplyr::select(Date, post_close = Close)

# Filter pre data
df_data <- df_data %>%
  filter(Date >= cutoff, Date <= post_cutoff,
         !is.na(logret))

##
# HMM
hmm_model <- depmix(
  logret ~ 1,
  df_data,
  nstates = 3
)
hmm_model_fitted <- fit(hmm_model)

df_data <- df_data %>%
  bind_cols(posterior(hmm_model_fitted)) %>%
  mutate(Close_end=lead(Close),
         logret_end = lead(logret))

# Multiplot
p1 <- ggplot(df_data) +
  geom_segment(aes(x=Date, xend=lead(Date), y=Close, yend=Close_end, colour=factor(state))) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  labs(title = "Bitcoin Closing Prices, Returns, and State", y = "Close (USD)")
p2 <- ggplot(df_data) +
  geom_segment(aes(x=Date, xend=lead(Date), y=logret, yend=logret_end, colour=factor(state))) +
  labs(y = "Log Returns") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank())
p3 <- ggplot(df_data, aes(x=Date, y=state)) +
  geom_line() +
  labs(y = "HMM State") +
  theme_classic()
multiplot(plotlist = list(p1, p2, p3), cols = 1)

# Shapiro-Wilk Test
# normality tests of logret show not significantly different from normal distribution
shapiro.test((df_data %>% filter(state == 1))$logret)
shapiro.test((df_data %>% filter(state == 2))$logret)
shapiro.test((df_data %>% filter(state == 3))$logret)

# Plot a density of logreturns, they look kind of normal
ggplot(df_data, aes(x=logret, colour=factor(state))) +
  geom_density() +
  stat_function(fun = dnorm,
                args = list(
                  mean = mean((df_data %>% filter(state == 1))$logret, na.rm=TRUE),
                  sd =   sd((df_data %>% filter(state == 1))$logret, na.rm=TRUE)),
                colour = "black") +
  stat_function(fun = dnorm,
                args = list(
                  mean = mean((df_data %>% filter(state == 2))$logret, na.rm=TRUE),
                  sd =   sd((df_data %>% filter(state == 2))$logret, na.rm=TRUE)),
                colour = "black") +
  stat_function(fun = dnorm,
                args = list(
                  mean = mean((df_data %>% filter(state == 3))$logret, na.rm=TRUE),
                  sd =   sd((df_data %>% filter(state == 3))$logret, na.rm=TRUE)),
                colour = "black")


##
# Predict forward
# Get transition matrix from hmm model
trans_matrix <- matrix(c(
  hmm_model_fitted@transition[[1]]@parameters$coefficients,
  hmm_model_fitted@transition[[2]]@parameters$coefficients,
  hmm_model_fitted@transition[[3]]@parameters$coefficients
  ), ncol = 3, byrow = TRUE
)

# Get response coefficients from hmm model
hmm_means <- c(
  hmm_model_fitted@response[[1]][[1]]@parameters$coefficients[1],
  hmm_model_fitted@response[[2]][[1]]@parameters$coefficients[1],
  hmm_model_fitted@response[[3]][[1]]@parameters$coefficients[1]
)
hmm_sds <- c(
  hmm_model_fitted@response[[1]][[1]]@parameters$sd[1],
  hmm_model_fitted@response[[2]][[1]]@parameters$sd[1],
  hmm_model_fitted@response[[3]][[1]]@parameters$sd[1]
)

# Function to simulate the HMM model forward
pred_state_seq <- function(data, trans, length)
{
  # Get initial state
  df_pred <- NULL
  curr_state <- last(data$state)
  curr_date  <- last(data$Date)
  curr_close <- last(data$Close)
  
  for (i in 1:length)
  {
    # Get next state
    probs <- cumsum(trans[curr_state, ])
    ran <- runif(1, 0, 1)
    if (ran < probs[1]) {curr_state <- 1}
    else if (ran < probs[2]) {curr_state <- 2}
    else {curr_state <- 3}
    
    # Get realization
    realized_logret <- rnorm(1, hmm_means[curr_state], hmm_sds[curr_state])
    
    # Get Closing price
    curr_close <- curr_close * exp(realized_logret)
    
    # Increment date
    curr_date <- curr_date + 1
    
    df_pred <- df_pred %>%
      bind_rows(
        data.frame(
          Date = curr_date,
          state = curr_state,
          pred_close = curr_close
        )
      )
  }
  
  return(df_pred)
}

# Do 10000 30-day forward simulations to get confidence intervals
df_pred <- NULL
for (i in 1:10000)
{
  tmp <- pred_state_seq(df_data, trans_matrix, 30)
  df_pred <<- df_pred %>%
    bind_rows(
      bind_cols(
        data.frame(trial = rep(i, 30)),
        tmp
      )
    )
  
  if (i %% 100 == 0) print(i)
}

# Group trials, calculate intervals, join with actual data
full_data <- df_pred %>%
  group_by(Date) %>%
  summarize(lower_5 = quantile(pred_close, 0.025),
            lower_10 = quantile(pred_close, 0.1),
            median = quantile(pred_close, 0.5),
            upper_10 = quantile(pred_close, 0.9),
            upper_5 = quantile(pred_close, 0.975),
            mean = mean(pred_close)) %>%
  full_join(df_data, by=c("Date")) %>%
  full_join(df_post_data, by=c("Date"))

# Final chart, BTC price with prediction
ggplot(full_data, aes(x=Date)) +
  geom_line(aes(y=Close), colour = "black") +
  geom_ribbon(aes(ymin = lower_5, ymax = upper_5), alpha = 0.3) +
  geom_ribbon(aes(ymin = lower_10, ymax = upper_10), alpha = 0.3) +
  geom_line(aes(y=mean), colour = "blue") +
  geom_line(aes(y=post_close), colour = "red") +
  labs(title = "Bitcoin Daily Close Price and HMM Prediction",
       subtitle = "HMM is a 3 state model on log value of daily returns",
       x = "Date", y = "Close Price (USD)") +
  theme_bw()

ggsave(filename = "btc_hmm.png", width = 10, height = 5)

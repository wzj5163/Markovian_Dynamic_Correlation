rm(list = ls())
setwd("~/Documents/Code/Markov_Dynamic_Correlation/0_6_0_CSI300_2013to2023_TtestDiscretization")

# Load necessary libraries
library(dplyr)

# Function to simulate future states using P1
simulate_with_P1 <- function(P, start_state, N) {
  states <- c("-1", "0", "1")
  current_state <- which(states == start_state)
  future_states <- character(N)
  
  for (i in 1:N) {
    probs <- P[current_state, ]
    future_states[i] <- sample(states, size = 1, prob = probs)
    current_state <- which(states == future_states[i])
  }
  
  return(future_states)
}

# Function to simulate future states using P2
simulate_with_P2 <- function(P, start_states, N) {
  states <- c("-1", "0", "1")
  current_states <- start_states
  future_states <- character(N)
  
  for (i in 1:N) {
    sequence_key <- paste(current_states, collapse = ",")
    probs <- P[sequence_key, ]
    next_state <- sample(states, size = 1, prob = probs)
    future_states[i] <- next_state
    current_states <- c(tail(current_states, 1), next_state)
  }
  
  return(future_states)
}

# Assume P1 and P2 are already loaded or calculated
P1 <- read.csv(paste0("Transitions/P1_ttest_alpha0d05.csv"),row.names = 1)
colnames(P1) = c("-1", "0", "1")
P2 <- read.csv(paste0("Transitions/P2_ttest_alpha0d05.csv"),row.names = 1)
colnames(P2) = c("-1", "0", "1")
# Simulate future states
simulation_times = 1000
N <- 10  # Number of periods to simulate

# Initialize lists to store simulation outcomes
future_states_P1 <- vector("list", simulation_times)
future_states_P2 <- vector("list", simulation_times)

# Initial probability distributions for starting states
P1_population <- colSums(P1) / sum(colSums(P1))
P2_population <- colSums(P2) / sum(colSums(P2))

# Fix the random seed for reproducibility
set.seed(123)

# Perform simulations
for (i in 1:simulation_times){
  start_state_P1 <- sample(c("-1", "0", "1"), size = 1, prob = P1_state0)  # Starting state for P1
  start_states_P2 <- sample(c("-1", "0", "1"), size = 2, prob = P1_state0)   # Starting states for P2
  
  future_states_P1[[i]] <- simulate_with_P1(P1, start_state_P1, N)
  future_states_P2[[i]] <- simulate_with_P2(P2, start_states_P2, N)
  }

# Comparison
# For a basic comparison, you might want to analyze the distribution of final states
final_states_P1 <- sapply(future_states_P1, function(x) tail(x, 1))
final_states_P2 <- sapply(future_states_P2, function(x) tail(x, 1))

# Calculate the frequency of final states for each simulation set
final_states_freq_P1 <- table(final_states_P1)
final_states_freq_P2 <- table(final_states_P2)

# Print the frequencies for comparison
print("Frequency of Final States for P1:")
print(final_states_freq_P1)
print("Frequency of Final States for P2:")
print(final_states_freq_P2)

rm(list = ls())
setwd("~/Documents/Code/Markov_Dynamic_Correlation/0_6_0_CSI300_2013to2023_TtestDiscretization")
dir.create("Transitions/",showWarnings = FALSE)
library(xts)
library(PortfolioAnalytics)
library(Spillover)
library(lubridate)
library(writexl)
library(huge)
library(readxl)
source('MDC_functions.R')

discretized_ecor = read.csv("Correlations/Discrezation/t_test/ecor.csv")

corr_to_transition_probability <- function(corr_matrix, lag=1) {
  P_count = matrix(0, ncol = 3, nrow = 3^lag)
  colnames(P_count) = c("-1", "0", "1")
  # Iterate through the matrix to count transitions
  for (i1 in 1:(nrow(corr_matrix)-lag)) {
    for (j in 1:(ncol(corr_matrix))) {
      S1 = corr_matrix[i1, j]
      current_state = corr_matrix[i1+lag, j]
      if (!is.na(S1) & !is.na(current_state)){
        S1_index = which(c("-1","0","1") == as.character(S1))
        current_index = which(c("-1","0","1") == as.character(current_state))
        # Increment the count for this state transition
        P_count[S1_index, current_index] = P_count[S1_index, current_index]+1
      }
    }
  }  
  # Calculate transition probabilities
  P <- sweep(P_count, 1, rowSums(P_count), "/")
  return(P)
}

P = corr_to_transition_probability(discretized_ecor, lag=1)
rownames(P) = c("-1", "0", "1")
write.csv(P, paste0("Transitions/P1_ttest_alpha0d05.csv"),row.names = TRUE)


corr_to_transition_probability <- function(corr_matrix, lag = 1) {
  # Number of states for each lag position
  num_states <- 3
  total_states <- num_states^lag
  
  # Initialize the transition count matrix
  P_count <- array(0, dim = c(total_states, num_states))
  dimnames(P_count) <- list(rep("", total_states), c("-1", "0", "1"))
  
  # Generate all possible states combinations for the given lag
  states_combinations <- expand.grid(replicate(lag, c(-1, 0, 1), simplify = FALSE))
  rownames(P_count) <- apply(states_combinations, 1, paste, collapse = ",")
  
  # Iterate through the matrix to count transitions
  for (i in 1:(nrow(corr_matrix) - lag)) {
    for (j in 1:ncol(corr_matrix)) {
      # Extract the sequence of states leading up to the current state
      sequence <- corr_matrix[i:(i + lag - 1), j]
      if (all(!is.na(sequence))) {
        current_state <- corr_matrix[i + lag, j]
        # Convert states to a string key
        sequence_key <- paste(sequence, collapse = ",")
        if (!is.na(current_state)) {
          # Find the index for the sequence and the current state
          sequence_index <- which(rownames(P_count) == sequence_key)
          current_index <- which(c("-1", "0", "1") == as.character(current_state))
          # Increment the count for this state transition
          P_count[sequence_index, current_index] <- P_count[sequence_index, current_index] + 1
        }
      }
    }
  }
  
  # Calculate transition probabilities
  P <- sweep(P_count, 1, rowSums(P_count), "/")
  P[is.na(P)] <- 0
  return(P)
}

# Example usage
P1 <- corr_to_transition_probability(as.matrix(discretized_ecor), lag = 1)
write.csv(P1, paste0("Transitions/P1_ttest_alpha0d05.csv"),row.names = TRUE)
P2 <- corr_to_transition_probability(as.matrix(discretized_ecor), lag = 2)
write.csv(P2, paste0("Transitions/P2_ttest_alpha0d05.csv"),row.names = TRUE)
P3 <- corr_to_transition_probability(as.matrix(discretized_ecor), lag = 3)
write.csv(P3, paste0("Transitions/P3_ttest_alpha0d05.csv"),row.names = TRUE)
P4 <- corr_to_transition_probability(as.matrix(discretized_ecor), lag = 4)
write.csv(P4, paste0("Transitions/P4_ttest_alpha0d05.csv"),row.names = TRUE)
P5 <- corr_to_transition_probability(as.matrix(discretized_ecor), lag = 5)
write.csv(P5, paste0("Transitions/P5_ttest_alpha0d05.csv"),row.names = TRUE)
P6 <- corr_to_transition_probability(as.matrix(discretized_ecor), lag = 6)
write.csv(P6, paste0("Transitions/P6_ttest_alpha0d05.csv"),row.names = TRUE)


# Initialize L
L=1
# Initially run the function with L=1
P <- corr_to_transition_probability(as.matrix(discretized_ecor), lag = L)
# Loop to increment L until condition is met or a max limit is reached to prevent infinite loops
max_lag = 12
while(!all(P %in% c(0,1)) && L<=max_lag){
  L=L+1
  P <- corr_to_transition_probability(as.matrix(discretized_ecor), lag = L)
}

if (L > max_lag_limit) {
  cat("Reached max lag limit without fulfilling condition. Last L:", L, "\n")
} else {
  cat("Condition met at L:", L, "\n")
}



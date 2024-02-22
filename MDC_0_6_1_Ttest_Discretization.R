rm(list = ls())
setwd("~/Documents/Code/Markov_Dynamic_Correlation/0_6_0_CSI300_2013to2023_TtestDiscretization")
dir.create("Correlations/",showWarnings = FALSE)
dir.create("Correlations/Discrezation/",showWarnings = FALSE)
dir.create("Correlations/Discrezation/t_test/",showWarnings = FALSE)
library(xts)
library(PortfolioAnalytics)
library(lubridate)
library(writexl)
# library(huge)
library(readxl)
# source('MDC_functions.R')


count_matrix = read.csv(paste0("Correlations/Corr_count_matrix.csv"))
count_matrix = as.matrix(count_matrix)
ecor_matrix = read.csv(paste0("Correlations/Corr_matrix.csv"))
ecor_matrix = as.matrix(ecor_matrix)

t_stat_matrix <- ecor_matrix * sqrt(count_matrix - 2) / sqrt(1 - ecor_matrix^2)

# Calculate degrees of freedom
df <- count_matrix - 2

# Calculate the absolute t-statistics
abs_t_stat_matrix <- abs(t_stat_matrix)

# Calculate the cumulative probability for the lower tail
lower_tail_prob <- pt(-abs_t_stat_matrix, df)

# Calculate two-tailed p-value
p_value_matrix <- 2 * lower_tail_prob

alpha = 0.05

discretized_ecor = matrix(0, nrow = nrow(ecor_matrix), ncol = ncol(ecor_matrix))

positive_indices <- which(ecor_matrix > 0 & p_value_matrix <= alpha)
negative_indices <- which(ecor_matrix < 0 & p_value_matrix <= alpha)
na_indices <- which(is.na(p_value_matrix) )

discretized_ecor[positive_indices] = 1
discretized_ecor[negative_indices] = -1
discretized_ecor[na_indices] = NA

# Flatten the matrix to a vector for frequency calculation
discretized_vec <- as.vector(discretized_ecor)
# Replace NA with a placeholder
discretized_vec[is.na(discretized_vec)] <- "NA"

# Factor and count all values including the placeholder for NA
value_counts <- table(factor(discretized_vec, levels = c("-1", "0", "1", "NA")))
# Total count NAs for proportion calculation
total_count <- sum(value_counts)
# Calculate proportions
proportions <- value_counts / total_count
# Display proportions
print(proportions)

# Factor and count all values excluding the placeholder for NA
nonNA_value_counts <- table(factor(discretized_vec, levels = c("-1", "0", "1")))
total_nonNA_count <- sum(nonNA_value_counts)
proportions <- nonNA_value_counts / total_nonNA_count
print(proportions)

write.csv(discretized_ecor, "Correlations/Discrezation/t_test/ecor.csv")


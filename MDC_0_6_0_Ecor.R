rm(list = ls())
setwd("~/Documents/Code/Markov_Dynamic_Correlation/0_6_0_CSI300_2013to2023_TtestDiscretization")
version = "0_5_0/"
dir.create("Correlations",showWarnings = FALSE)
dir.create("Excess_returns",showWarnings = FALSE)
dir.create("Transitions",showWarnings = FALSE)
library(xts)
library(PortfolioAnalytics)
library(Spillover)
library(lubridate)
library(writexl)

price<-read.csv("CSI300_2000to20231219_BackwardAdjust.csv")
price_xts<-xts(zoo(price[,-1],order.by = as.Date(price$DateTime,format='%Y/%m/%d')))
return_xts<- Return.calculate(price_xts, method="log")[-1,]
return_xts<- return_xts['2012-12/2023']
return_xts<- return_xts[,!is.na(return_xts[1,])]

#### Calculate empirical correlation matrices of these stocks. ####
# I want to focus on each month of their returns, so the empirical correlation matrices would be non-overlapping. 
# For each month, if a stock has a 5-day consecutive returns equal to 0, we believe this stock is suspended for trading, 
# so we set the returns of this stock in this month to NAs and this stock has NA correlation values to other stocks.
# Notice the empirical correlation matrices are calculated by excess returns.
# Finally, I expect to get 132 non-overlapping empirical correlation matrices, and save them in separate csv files.

# Create a sequence of dates
dates <- seq(ymd("2012-12-01"), ymd("2023-11-01"), by = "1 month")
dates <- format(dates, "%Y-%m")
# Convert to a list
paired_dates_list <- as.list(dates)
ecor_matrix = c()
count_matrix = c()
for (i in 1:length(paired_dates_list)) {
  return_xts_1month <- return_xts[paired_dates_list[[i]]]
  # If a stock has a 5-day consecutive returns equal to 0, we define this stock is suspended for trading
  suspended_stocks <- apply(return_xts_1month, 2, function(x) any( rle(x==0)$lengths >= 5 & rle(x==0)$values))
  
  # Check if there are any columns with all zeroes
  if (any(suspended_stocks)) {
    cat(paired_dates_list[[i]],"Suspended stocks:\n")
    print(names(return_xts_1month)[suspended_stocks])
  }
  
  # Set returns of suspended stocks to NA
  return_xts_1month[,suspended_stocks] <- NA
  
  # Calculate market return, as average return of all stocks
  market_return = rowMeans(return_xts_1month, na.rm = TRUE)
  
  # Calculate excess return of each stock
  excess_return_1month = sweep(return_xts_1month, 1, market_return, "-")
  excess_return_1month[,suspended_stocks] = NA
  df_excess_return_1month = data.frame(excess_return_1month)
  colnames(df_excess_return_1month) = colnames(return_xts)
  write_xlsx(df_excess_return_1month, paste0("Excess_returns/Return_",paired_dates_list[[i]],".xlsx"))
  
  
  
  # Calculate empirical correlation matrix by excess return and set the suspended stocks' correlations to NA
  ecor_1month = cor(excess_return_1month, use = "pairwise.complete.obs")
  ecor_1month[suspended_stocks,] = NA
  ecor_1month[,suspended_stocks] = NA
  df_ecor_1month = data.frame(ecor_1month)
  colnames(df_ecor_1month)=colnames(return_xts)
  # write_xlsx(df_ecor_1month, paste0("Correlations/Corr_",paired_dates_list[[i]],".xlsx"))
  write.csv(df_ecor_1month, paste0("Correlations/Corr_",paired_dates_list[[i]],".csv"),row.names = FALSE)
  
  # Identify indices of the upper triangular off-diagonal part
  ecor_row = ecor_1month[upper.tri(ecor_1month)]
  ecor_matrix = rbind(ecor_matrix, ecor_row)
  
  # Initialize the matrix for storing the number of observations
  obs_count_matrix <- matrix(NA, ncol = ncol(ecor_1month), nrow = ncol(ecor_1month))
  for (i in 1:(ncol(ecor_1month)-1)) {
    for (j in (i+1):ncol(ecor_1month)) {
      # Calculate the number of non-NA observations for each pair
      valid_obs <- sum(!is.na(excess_return_1month[,i]) & !is.na(excess_return_1month[,j]))
      obs_count_matrix[i,j] <- valid_obs
      obs_count_matrix[j,i] <- valid_obs  # Making the matrix symmetric
    }
  }
  # Set rows and columns for suspended stocks to NA
  obs_count_matrix[suspended_stocks,] <- NA
  obs_count_matrix[,suspended_stocks] <- NA
  # Extract the upper triangular part of the matrix as a vector
  obs_count_upper_tri <- obs_count_matrix[upper.tri(obs_count_matrix)]
  # To store this information
  count_matrix <- rbind(count_matrix, obs_count_upper_tri)  # Assuming ecor_matrix is prepared for this operation
  
}
ecor_matrix
df_ecor_matrix = data.frame(ecor_matrix)
write.csv(df_ecor_matrix, paste0("Correlations/Corr_matrix.csv"),row.names = FALSE)
count_matrix
df_count_matrix = data.frame(count_matrix)
write.csv(df_count_matrix, paste0("Correlations/Corr_count_matrix.csv"),row.names = FALSE)

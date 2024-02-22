rm(list = ls())
setwd("~/Documents/Code/Markov_Dynamic_Correlation/0_5_0_CSI300_2013to2023_MeanReverting")
dir.create("Correlations/",showWarnings = FALSE)
dir.create("Correlations/Discrezation/",showWarnings = FALSE)
dir.create("Correlations/Discrezation/shrinkage/",showWarnings = FALSE)
dir.create("Correlations/Discrezation/truncation/",showWarnings = FALSE)
dir.create("Correlations/huge",showWarnings = FALSE)
dir.create("Correlations/huge/shrinkage",showWarnings = FALSE)
dir.create("Correlations/huge/truncation",showWarnings = FALSE)
library(xts)
library(PortfolioAnalytics)
library(lubridate)
library(writexl)
library(huge)
library(readxl)
source('MDC_functions.R')

# price<-read.csv("CSI300_2000to20231219_BackwardAdjust.csv")
# price_xts<-xts(zoo(price[,-1],order.by = as.Date(price$DateTime,format='%Y/%m/%d')))
# return_xts<- Return.calculate(price_xts, method="log")[-1,]
# return_xts<- return_xts['2012-12/2023']
# return_xts<- return_xts[,!is.na(return_xts[1,])]

# Create a sequence of dates
dates <- seq(ymd("2012-12-01"), ymd("2023-11-01"), by = "1 month")
dates <- format(dates, "%Y-%m")
# Convert to a list
paired_dates_list <- as.list(dates)
lambda_shrinkage_t = list()
lambda_truncation_t = list()

for (i in 1:length(paired_dates_list)) {
  excess_return_1month = read_excel(paste0("Excess_returns/Return_",paired_dates_list[[i]],".xlsx"))
  
  # Nonparanormal 
  stocks.shrinkage = huge.npn(excess_return_1month)
  stocks.truncation = huge.npn(excess_return_1month, npn.func = "truncation")
  # stocks.skeptic = huge.npn(stocks, npn.func = "skeptic")
  out.shrinkage=huge(coredata(stocks.shrinkage))
  out.truncation=huge(coredata(stocks.truncation))
  # out.skeptic=huge(coredata(stocks.skeptic))
  #  optimal lambda of out
  out.shrinkage.select=huge.select(out.shrinkage)
  out.truncation.select=huge.select(out.truncation)
  # out.skeptic.select=huge.select(out.skeptic) # not available
  
  lambda_shrinkage_t[[i]]=out.shrinkage.select$opt.lambda
  lambda_truncation_t[[i]]=out.truncation.select$opt.lambda
  
  graph_shrinkage_t=as.matrix(out.shrinkage.select$refit)
  graph_truncation_t=as.matrix(out.truncation.select$refit)
  df = data.frame(graph_shrinkage_t)
  colnames(df)=colnames(excess_return_1month)
  write_xlsx(df, paste0("Correlations/huge/shrinkage/mb_RIC_",paired_dates_list[[i]],".xlsx"))
  df = data.frame(graph_truncation_t)
  colnames(df)=colnames(excess_return_1month)
  write_xlsx(df, paste0("Correlations/huge/truncation/mb_RIC_",paired_dates_list[[i]],".xlsx"))
  
  # read FEVD_norm, FEVD_nonnorm and empirical correlation matrix
  # Emp_corr = read_excel(paste0("Correlations/Corr_",paired_dates_list[[i]],".xlsx"))
  Emp_corr = read.csv(paste0("Correlations/Corr_",paired_dates_list[[i]],".csv"))
  
  # Loop over upper triangular off-diagonal
  Empiric_corr_shrinkage_t = discretize_corr(as.matrix(Emp_corr),graph_shrinkage_t)
  Empiric_corr_truncation_t = discretize_corr(as.matrix(Emp_corr),graph_truncation_t)
  Npn_corr_shrinkage_t = discretize_corr(as.matrix(cor(stocks.shrinkage)),graph_shrinkage_t)
  Npn_corr_truncation_t = discretize_corr(as.matrix(cor(stocks.truncation)),graph_truncation_t)
  
  df3 = data.frame(Empiric_corr_shrinkage_t)
  colnames(df3)=colnames(excess_return_1month)
  write_xlsx(df3, paste0("Correlations/Discrezation/shrinkage/Corr_mb_RIC_",paired_dates_list[[i]],".xlsx"))
  df3 = data.frame(Empiric_corr_truncation_t)
  colnames(df3)=colnames(excess_return_1month)
  write_xlsx(df3, paste0("Correlations/Discrezation/truncation/Corr_mb_RIC_",paired_dates_list[[i]],".xlsx"))
  df4 = data.frame(Npn_corr_shrinkage_t)
  colnames(df4)=colnames(excess_return_1month)
  write_xlsx(df4, paste0("Correlations/Discrezation/shrinkage/Npn_Corr_mb_RIC_",paired_dates_list[[i]],".xlsx"))
  df4 = data.frame(Npn_corr_truncation_t)
  colnames(df4)=colnames(excess_return_1month)
  write_xlsx(df4, paste0("Correlations/Discrezation/truncation/Npn_Corr_mb_RIC_",paired_dates_list[[i]],".xlsx"))
  
}
Lambda_shrinkage_t<-as.matrix(cbind(unlist(lambda_shrinkage_t)))
Lambda_truncation_t<-as.matrix(cbind(unlist(lambda_truncation_t)))
df3 = data.frame(Lambda_shrinkage_t)
df3 = cbind(dates,df3)
colnames(df3)=c("Date","Lambda")
write_xlsx(df3, paste0("Correlations/huge/shrinkage/Lambda_mb_RIC_",paired_dates_list[[i]],".xlsx"))
df3 = data.frame(Lambda_truncation_t)
df3 = cbind(dates,df3)
colnames(df3)=c("Date","Lambda")
write_xlsx(df3, paste0("Correlations/huge/truncation/Lambda_mb_RIC_",paired_dates_list[[i]],".xlsx"))

Fisher_Z(0)

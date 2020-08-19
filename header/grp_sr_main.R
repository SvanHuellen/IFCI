# You need to set the data path for your own directory
primary_directory = 'C:/Users/Ho Tak Yui/Documents/FCI Project/R code_andy'
library_directory = 'C:/Users/Ho Tak Yui/Documents/FCI Project/R code_andy/implementation'
data_directory    = 'C:/Users/Ho Tak Yui/Documents/FCI Project/R code_andy/data'
output_directory  = primary_directory # type whatever path you want here
#######

# input data files
dependent_filename = 'ch_d03lm1.txt'
indicator_filename = 'grouped_all.txt'
grp_filename       = 'sr_group_test.txt'
#######

# Initialisation of global variables
fbgn                = "1993-01-01"
fend                = "2000-12-01"
forecasting_horizon = 12
start_lag           = 1
max_lag             = 6
update_interval     = forecasting_horizon # prefer to be equal
n_update            = 5
p_value             = .1
lag_length          = max_lag - start_lag + 1
######

# if you want to save outputs, set it to TRUE
save_output = TRUE
# file names
fci_filename             = 'srfci.csv'
unconca_filename         = 'sr_unconca.csv'
weight_filename          = 'sr_weight.csv'
filtered_weight_filename = 'sr_filtered_weight.csv'
mean_filename            = 'cum_mean_diff.csv'
######

# script starts: the below are just function calls according to the above set-up
# set path
setwd(primary_directory)
source(paste(library_directory, 'common.R', sep='/'))
source(paste(library_directory, 'SRFCI.R', sep='/'))
source(paste(library_directory, 'grp_sr.R', sep='/'))

# function 'packages' checks for installed packages before running install.packages(); if not working: use install.packages('xts')
packages_check('xts')
library(xts)

# read data
dependent  = data_read(paste(data_directory, dependent_filename, sep='/'))
data       = data_read(paste(data_directory, indicator_filename, sep='/'))
grp        = grp_read(paste(data_directory, grp_filename, sep='/'))
# size of data
n_pc = dim(data)[2]

# generate lag variables
g_lag = data_lag_expansion(data, lags = max_lag, start_lag = start_lag)

# generate a sample matrix for looping
fbgn   = as.Date(fbgn)
fend   = as.Date(fend)
sample = gen_samp(fbgn, fend, forecasting_horizon, max_lag, update_interval, n_update)

# generate weight and factor time series
regr_res       = grp_leading_filter(data_lag = g_lag, dependent, grp, sample, lag_length, p_value)
factor         = regr_res$factor_ts
weight         = regr_res$weight_ts
remain         = regr_res$ind_count_ts
filtered_table = weight_table(weight)

# conca
conca_res = sr_conca(factor, fend, update_interval, forecasting_horizon, n_update)
srfci     = conca_res$conca
unconca   = conca_res$unconca
mean      = conca_res$cum_diff 

# save outputs
if(save_output){
  write.table(filtered_table, file = paste(output_directory, filtered_weight_filename, sep='/'), sep=',', col.names=FALSE, row.names=FALSE, na='')
  write.zoo(srfci, file = paste(output_directory, fci_filename, sep='/'), sep=',', col.names=TRUE, na='')
  write.zoo(unconca, file = paste(output_directory, unconca_filename, sep='/'), sep=',', col.names=TRUE, na='')
  write.table(grp_weight(weight,grp), file = paste(output_directory, weight_filename, sep='/'), sep=',', col.names=NA, row.names=TRUE, na='')
  write.zoo(mean, file = paste(output_directory, mean_filename, sep='/'), sep=',', col.names=TRUE, na='')
}
# You need to set the data path for your own directory
primary_directory = 'C:/Users/Ho Tak Yui/Documents/FCI Project/R code_andy'
library_directory = 'C:/Users/Ho Tak Yui/Documents/FCI Project/R code_andy/implementation'
data_directory    = 'C:/Users/Ho Tak Yui/Documents/FCI Project/R code_andy/data'
output_directory  = primary_directory # type whatever path you want here
#######

# input data files
dependent_filename = 'lcpi.txt'
indicator_filename = 'lr_g_ori_fs.txt'
#######

# Initialisation of global variables
fbgn                = "1980-01-01"
fend                = "2000-12-01"
forecasting_horizon = 12
maxLag              = 6
ncomp               = 1 # some functionality does not work for n>1
update_interval     = 12 # prefer to be equal
n_update            = 16
thresold            = 0.05 # the level for the weight screening
######

# if you want to save outputs, set it to TRUE
save_output = TRUE
# file names
lrfci_filename  = 'lrfci.csv'
weight_filename = 'lr_weight.csv'
filtered_weight_filename = 'lr_filtered_weight.csv'
sig_table_file_name = 'lr_sig_table.csv'
######

# script starts: the below are just function calls according to the above set-up
# set path
setwd(primary_directory)
source(paste(library_directory, 'common.R', sep='/'))
source(paste(library_directory, 'LRFCI.R', sep='/'))

# function 'packages' checks for installed packages before running install.packages(); if not working: use install.packages('xts')
packages_check('xts')
library(xts)
# use pls package
packages_check('pls')
library(pls)

# read data
dependent  = data_read(paste(data_directory, dependent_filename, sep='/'))
data       = data_read(paste(data_directory, indicator_filename, sep='/'))

# generate lag variables
g_lag = data_lag_expansion(data, start_lag=0, lags = maxLag)

# moving window
fbgn   = as.Date(fbgn)
fend   = as.Date(fend)
pls_ts = mov_win_pls(data=g_lag, dependent, fbgn, fend, forecasting_horizon, update_interval, n_update, maxLag, ncomp)

# get weight
weight = flatten_wgt(pls_ts, attr='weight')
# get filtered weight
w_sparse       = flatten_wgt(pls_ts, attr='w_sparse')
filtered_table = weight_table(w_sparse)
# get factor
lrfci = conca(pls_ts, fend, update_interval, forecasting_horizon, attr = 'scaled_factor')

# weight screening
sig_table = weight_screening(weight, thresold = thresold)

# save outputs
if(save_output){
  write.table(weight, file = paste(output_directory, weight_filename, sep='/'), sep=',', col.names=TRUE, row.names=TRUE, na='')
  write.table(filtered_table, file = paste(output_directory, filtered_weight_filename, sep='/'), sep=',', col.names=FALSE, row.names=FALSE, na='')
  write.zoo(lrfci, file = paste(output_directory, lrfci_filename, sep='/'), sep=',', col.names=TRUE, na='')
  write.table(sig_table, file = paste(output_directory, sig_table_file_name, sep='/'), sep=',', col.names=FALSE, row.names=TRUE, na='')
}

# True/False weight table for significance (in any lag)
# weight selection: relative difference
# You need to set the data path for your own directory
primary_directory = 'C:/Users/Ho Tak Yui/Documents/FCI Project/R code_andy'
library_directory = 'C:/Users/Ho Tak Yui/Documents/FCI Project/R code_andy/implementation'
data_directory    = 'C:/Users/Ho Tak Yui/Documents/FCI Project/R code_andy/data'
output_directory  = primary_directory # type whatever path you want here
#######

# input data files
data_filename = 'ip.txt'
#######

# Initialisation of global variables
steps               = 6
n_update            = 2
update_interval     = 12
forecast_length     = 12
fbgy                = "1980"
freq                = 'm' #do you really have another freq?
fbgp                = '1'
fendy               = "2000"
fendp               = '12'
eq_shift            = FALSE
dyn_forecast        = TRUE
# variables
target              = 'd12lip'
control             = 'd12lgdp'
ints                = 'd4ird' #short term interest rate
intl                = 'irl' #long term interest rate
ec_def              = 'ec'
# model specification
var_default = c(target, control, ints, ec_def) # variables included in the default model
lag_default = list(c(0,1), c(0,1), c(2), c(12)) # each column c() is corresponding to a variable; target lag 0 has to be included
var_non_fci = c(target, control) # variables included in the non fci part of the model
lag_non_fci = list(c(0,1), c(0,1))
# fci specification
var_fci = c('ec_07fci', 'srfci')
lag_fci = list(c(12), c(0))
######

# if you want to save outputs, set it to TRUE
save_output = TRUE
# output: all available outputs are ('pval', 'cpval', 'mse_fci', 'mse_default', 'cmse_fci', 'cmse_default')
outputs = c('pval', 'cpval', 'mse_fci', 'mse_default', 'cmse_fci', 'cmse_default')
# file names: in the same order as outputs above
filenames = c('pval.csv', 'cpval.csv', 'mse_fci.csv', 'mse_default.csv', 'cmse_fci.csv', 'cmse_default.csv')
######

###### the program starts ###### 
# set path
setwd(primary_directory)
source(paste(library_directory, 'common.R', sep='/'))
source(paste(library_directory, 'forecast.R', sep='/'))
# pacakge
packages_check('xts')
library(xts)

# read data
data  = data_read(paste(data_directory, data_filename, sep='/'))

####### set up ######
if(freq=='m'){
  fbgn = as.Date(paste(fbgy, fbgp, '01', sep='-'))
  fend = as.Date(paste(fendy, fendp, '01', sep='-'))
} else{
  print('freq not implemented')
}

var_list = list('default'=var_default, 'non_fci'=var_non_fci, 'fci'=var_fci)
lag_list = list('default'=lag_default, 'non_fci'=lag_non_fci, 'fci'=lag_fci)
remove(var_default, var_non_fci, var_fci, lag_default, lag_non_fci, lag_fci)

res=encompass(data, target, var_list, lag_list, fbgn, fend, steps, n_update, update_interval, eq_shift, dyn_forecast)

# save outputs
if(save_output){
  for(i in 1:length(outputs)){
    write.zoo(res[[outputs[i]]], file = paste(output_directory, filenames[i], sep='/'), sep=',', col.names=FALSE, row.names=FALSE, na='')
  }
}
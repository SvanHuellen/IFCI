packages_check=function(x){
  # Install function for packages    
  x=as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

data_read = function(filename) {
  # This is a function that read the data file from .txt into R and format it into xts class
  data = read.delim(filename)
  # transform data from charater into numeric ones
  date = as.Date(data[, 1], "%Y-%m-%d")
  data = sapply(data[, -1], function(x) as.numeric(x))
  return(xts(as.data.frame(data), order.by = date))
}

standardize_data = function(data, insample_interval, outsample_interval = insample_interval) {
  # enable out-of-sample standardization
  mean = apply( data[insample_interval], 2, mean, na.rm = TRUE )
  sd = apply( data[insample_interval], 2, sd, na.rm =TRUE )
  return(scale(data[outsample_interval], center=mean, scale=sd))
}

sample_interval = function(data, start_period = time(first(data)), end_period = time(last(data)), start_increment = 0, end_increment = 0) {
  # create a sample time interval from data or start_period/end_period or both
  start_period = seq(start_period, by = paste(start_increment, 'months'), length = 2)[2]
  end_period = seq(end_period, by = paste(end_increment, 'months'), length = 2)[2]
  paste(start_period, '/', end_period, sep = '')
}

data_lag_expansion = function(data, lags, start_lag = 0, init_sample_interval = sample_interval(data, time(first(data)), time(last(data)), 0, 0),
                                   cut_sample_interval = init_sample_interval){
  # create lags from data
  if (lags == 0) { return(data[cut_sample_interval]) } 
  data = data[init_sample_interval]
  out  = list()
  for(i in 1:ncol(data)){
    # out[[(i - 1)*(lags + 1) + 1]] = data[, i]
    for (j in start_lag:lags){
      # out[[(i - 1) * (lags + 1) + 1 + j]] = lag(data[, i], j)
      lag_var = lag(data[, i], j)
      colnames(lag_var) = paste(colnames(data[ ,i]), j, sep='.')
      out[[(i - 1) * (lags - start_lag + 1) + j - start_lag + 1]] = lag_var
    }
  }
  date = time(data)
  return( xts(as.data.frame(out), order.by = date)[cut_sample_interval] )
}

weight_table = function(weight, table_fmt='multi'){
  # reformat the weight matrix to a table of filtered weight
  # single: single index ; multi: multilevel index
  # single is more reusable
  table1 = as.data.frame(weight)
  name   = unlist(lapply(rownames(table1), function(x) sub('\\..*', '', x)))
  lag    = rownames(table1)
  table1 = cbind(name, lag, table1)
  if(table_fmt == 'single'){
    return(table1)
  }
  if(table_fmt == 'multi'){
    table2 = data.frame(table1[FALSE, -1], stringsAsFactors=FALSE)
    table2[,'lag'] = as.character(table2[,'lag'])
    for(i in unique(name)){
      table2[nrow(table2)+1, ]  = c(i, rep(NA, length=ncol(table2)-1))
      table2[nrow(table2)+1, ]  = colnames(table2)
      subtable                  = table1[table1['name']==i, -1]
      subtable[is.na(subtable)] = 0
      table2                    = rbind(table2, subtable)
      table2[nrow(table2)+1, ]  = NA
    }
    return(table2)
  }
  return(weight)
}
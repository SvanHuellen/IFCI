dtoo = function(data, start_period = start(data), end_period = end(data)){
  # Date-TO-Observation
  idx = paste(start_period, end_period, sep='/')
  return(dim(data[idx])[1])
}

otod = function(data, start_period = start(data), n=dtoo(data)){
  # Observation-TO-Date
  idx = paste(start_period, '', sep='/')
  return(tail(index(first(data[idx], n)),1))
}

df_regr = function(data, var_array, lag_list){
  # create the data matrix for regression
  if (length(var_array)!=length(lag_list)){
    stop('variable size not equal to lag size!')
  }
  else {
    output = lag(data[, var_array[1]], lag_list[[1]])
    name   = paste(var_array[1], lag_list[[1]], sep='.')
    for (i in 2:length(lag_list)){
      output = merge(output, lag(data[, var_array[i]], lag_list[[i]]))
      name   = c(name, paste(var_array[i], lag_list[[i]], sep='.'))
    }
  }
  name = gsub('.0', '', name, fixed=TRUE)
  colnames(output) = name
  return(output)
}

pred_error = function(data, target, steps, fbgn, fend, update_interval, forecast_length, dyn_forecast){
# calculate forecast error in experiment and control
  df_pred    = data
  fit_samp   = sample_interval(data=NULL, start_period = fbgn, end_period = fend, start_increment = 0, end_increment = update_interval)
  eqn        = paste(target, '~.', sep='')
  regr       = lm(eqn, data=df_pred[fit_samp])
  # outputs
  pred_value = list()
  pred_error = list()
  mse        = c()
  for(p in 1:steps){
    pred_samp        = sample_interval(data=NULL, start_period = fend, end_period = fend, start_increment = update_interval+p, end_increment = update_interval+forecast_length)
    pred_value[[p]]  = predict(regr, df_pred[pred_samp])
    # update for dynamic prediction
    if(dyn_forecast){
      df_pred[names(pred_value[[p]]), target] = pred_value[[p]] 
      # need to update lag too
      pred_column = names(df_pred[,grep(target, colnames(df_pred), fixed=TRUE)])
      pred_lags   = as.numeric(sub(".*\\.", "", pred_column))
      for( l in 1:length(pred_lags) ){
        if( !is.na(pred_lags[l]) ){
          df_pred[, pred_column[l]] = lag(df_pred[, target], pred_lags[l])
        }
      }
    }
    # error calculation
    pred_error[[p]] = pred_value[[p]] - data[names(pred_value[[p]]), target]
    mse             = c(mse, sqrt(mean(pred_error[[p]]^2)))
  }
  return(list('mse' = mse, 'pred' = pred_value, 'error' = pred_error))
}

mdm_test = function(err1, err2){
  # Modified Diebold-Marianos statistic 
  dt    = list()
  dm    = c()
  stat  = c()
  pval  = c()
  for(p in 1:length(err1)){
    nobs = length(err1[[p]])
    # distance
    dt[[p]]    = ( err1[[p]] - err2[[p]] ) * err1[[p]]
    # statistics
    dm   = c( dm, mean( dt[[p]] ) / sd( dt[[p]] ) * sqrt(nobs) )
    stat = c( stat, sqrt(1/nobs) * sqrt( nobs+1-2*p+1/nobs*p*(p-1) ) * dm[p] )
    pval = c( pval, pt(stat[p], df = nobs-1, lower.tail = FALSE))
  }
  return(list('dm'=dm, 'stat'=stat, 'pval'=pval))
}

encompass=function(data, target, var_list, lag_list, fbgn, fend, steps, n_update, update_interval, eq_shift, dyn_forecast){
  if(eq_shift==1){
    stop('Equation shift not implemented yet!')
  }
  err_fci      = list()
  err_default  = list()
  mse_fci      = c()
  mse_default  = c()
  pval         = c()
  cpval        = c()
  cmse_fci     = c()
  cmse_default = c()
  for(u in 0:n_update){
    # rename fci variables
    fci_name   = paste(var_list[['fci']], u, sep='')
    var_fci_eq = c(var_list[['non_fci']], fci_name)
    lag_fci_eq = c(lag_list[['non_fci']], lag_list[['fci']])
    df_fci     = df_regr(data = data, var_array = var_fci_eq, lag_list = lag_fci_eq)
    df_default = df_regr(data = data, var_array = var_list[['default']], lag_list = lag_list[['default']])
    # forecasting error
    pred_fci     = pred_error(data=df_fci, target=target, steps=steps, fbgn=fbgn, fend=fend, update_interval = u*update_interval, forecast_length = forecast_length, dyn_forecast = dyn_forecast)
    pred_default = pred_error(data=df_default, target=target, steps=steps, fbgn=fbgn, fend=fend, update_interval = u*update_interval, forecast_length = forecast_length, dyn_forecast = dyn_forecast)
    # output mse
    mse_fci     = cbind(mse_fci, pred_fci[['mse']])
    mse_default = cbind(mse_default, pred_default[['mse']])
    # mdm test: output p-values
    mdm  = mdm_test( err1 = pred_fci[['error']], err2 = pred_default[['error']] )
    pval = cbind(pval, mdm[['pval']])
    #mdminv is not used
    #mdminv = mdm_test( err1 = pred_default[['error']], err2 = pred_fci[['error']] )
    # accumulate errors through update
    for(p in 1:steps){
      tmp_err_fci     = as.numeric(pred_fci[['error']][[p]])
      tmp_err_default = as.numeric(pred_default[['error']][[p]])
      if( p > length(err_fci) ){
        # initialization: this is not user unfriendly
        err_fci[[p]]     = tmp_err_fci
        err_default[[p]] = tmp_err_default
      } else {
        err_fci[[p]]     = c(err_fci[[p]], tmp_err_fci)
        err_default[[p]] = c(err_default[[p]], tmp_err_default)
      }
    }
    # cumulative mdm p-values
    cmdm  = mdm_test( err1 = err_fci, err2 = err_default )
    cpval = cbind(cpval, cmdm[['pval']])
    # cumulative mse
    cmse_fci = c(cmse_fci, sqrt(mean(unlist(err_fci)^2)))
    cmse_default = c(cmse_default, sqrt(mean(unlist(err_default)^2)))
  }
  #output: mdm/cmdm: p-values; mse/cmse
  return(list('pval'=pval, 'cpval'=cpval, 
              'mse_fci'=mse_fci, 'mse_default'=mse_default,
              'cmse_fci'=cmse_fci, 'cmse_default'=cmse_default))
}
# Eviews sampling 'smpl' may lead to missing entries in lags data for the first i^th lags
# This script properly fills those data and results in more indicators become significant

gen_samp = function(fbgn, fend, forecasting_horizon, max_lag, update_interval, n_update){
  # generate the sample pairs for moving window
  sample = matrix(nrow = 1 + n_update, ncol = 2)
  for (i in 0:n_update){
    #an update to exclude max_lag and forecasting_horizon
    #sample[ i+1, 1 ] = sample_interval(data=NULL, fbgn, fend, end_increment = update_interval*i + forecasting_horizon)
    #sample[ i+1, 2 ] = sample_interval(data=NULL, fbgn, fend, start_increment = max_lag, end_increment = forecasting_horizon + update_interval*i + max_lag)
    sample[ i+1, 1 ] = sample_interval(data=NULL, fbgn, fend, end_increment = update_interval*i)
    sample[ i+1, 2 ] = sample_interval(data=NULL, fbgn, fend, start_increment = max_lag, end_increment = forecasting_horizon + update_interval*i)
  }
  return(sample)
}

### leading_filter functions
regr = function(data_lag, dependent, samp, p_value = .1){
  # data = g_lag
  # produce ols regression coefficients and select the most significant three
  g_multi  = standardize_data(data_lag, insample_interval = samp, outsample_interval = samp)
  #X        = as.matrix(g_multi[( max_lag - start_lag + 2 ):nrow(g_multi),]) #is this necessary in Eviews?
  X        = as.matrix(g_multi[samp]) 
  y        = dependent[rownames(X)]
  X_name   = unique(lapply(colnames(g_multi), function(x) sub('\\..*', '', x)))
  # initiate outputs
  weight           = matrix(data=0, nrow=ncol(g_multi))
  rownames(weight) = colnames(g_multi)
  const            = 0
  ind_count        = 0
  for(i in 1:length(X_name)){
    # first regression: on all lags of a single independent variable
    d1     = data.frame(y, X[ , grep(X_name[[i]], colnames(X)) ])
    regr1  = lm(d1)
    anova1 = summary(regr1)$coefficients # columns are ('Esitmate', 'Std. Error', 't value', 'Pr(>|t|)')
    # extract 3 with smallest p-value (including intercept)
    X1_sig = rownames(anova1[ order(anova1[,4])[1:3], ]) 
    X1_sig = X1_sig[! X1_sig %in% c('(Intercept)')] # intercept is not needed for regression in lm
    # second regression: on the 3 selected variables
    d2     = data.frame(y, X[ ,  X1_sig ])
    regr2  = lm(d2)
    anova2 = summary(regr2)$coefficients
    #print(paste(X_name[i],summary(regr2)$adj.r.squared,sep=':'))
    # filtered by p-values
    X2_sig = rownames(anova2[which( anova2[,4] < p_value ),])
    if ('(Intercept)' %in% X2_sig){
      const  = const + anova2[c('(Intercept)'),'Estimate']
    }
    coeff = X2_sig[! X2_sig %in% c('(Intercept)')]
    if (length(coeff)>0){
      ind_count = ind_count + 1
    }
    # save significant coeff as weight
    weight[coeff,] = anova2[coeff, 'Estimate']
  }
  return(list('weight'=weight, 'const'=const, 'ind_count'=ind_count))
}

leading_filter = function(data_lag, dependent, sample, p_value = .1){
  # a moving window implementation to produce weight and factor
  w_ts           = matrix( nrow = dim(data_lag)[2], ncol = dim(sample)[1] )
  rownames(w_ts) = colnames(data_lag)
  colnames(w_ts) = paste('weight', 1:ncol(w_ts)-1, sep='')
  factor_ts      = list()
  const_ts       = c()
  ind_count_ts   = c()
  for (i in 1:dim(sample)[1]){
    regr_res       = regr( data_lag = g_lag, dependent, sample[i, 1], p_value )
    w_ts[, i]      = regr_res$weight
    g_multi        = standardize_data( data_lag, insample_interval = sample[i, 1], outsample_interval = sample[i, 2])
    factor         = as.matrix(g_multi[sample[i, 2]]) %*% w_ts[,i]
    factor         = factor/regr_res$ind_count + regr_res$const/regr_res$ind_count
    factor         = as.xts(factor)
    # standardize and rescale to target?
    #factor         = standardize_data(factor)*sd(dependent[sample[i, 2]])+mean(dependent[sample[i, 2]])
    factor_ts[[i]] = factor
    const_ts       = c(const_ts, regr_res$const)
    ind_count_ts   = c(ind_count_ts, regr_res$ind_count)
  }
  # flatten the list into 2d
  nrow                = max(sapply(factor_ts, length))
  factor_ts           = do.call(cbind,lapply(factor_ts, 'length<-', nrow))
  colnames(factor_ts) = paste('fci', 1:ncol(factor_ts)-1, sep='')
  return(list('factor_ts'=factor_ts, 'weight_ts'=w_ts, 'const_ts'=const_ts, 'ind_count_ts'=ind_count_ts))
}

### conca functions
sr_conca = function(factor, fend, update_interval, forecasting_horizon, n_update = dim(factor)[2]-1){
  # conca function for Short Run Model
  unconca  = as.xts(factor)
  conca    = as.xts(factor)
  for(i in 0:n_update){
    unconca[, i+1] = factor[, i+1]
    conca[, i+1]   = unconca[, 1]
  }
  mean_diff = c() #temp_up
  cum_diff  = c() #cum_up
  # an intercept adjustment
  for(i in 1:n_update){
    samp = sample_interval(data=NULL, fend, fend, start_increment = (i-1)*update_interval, end_increment = i*update_interval)
    mean_diff = c(mean_diff, mean(unconca[samp, i+1] - unconca[samp, i], na.rm=TRUE))
  }
  cum_diff = cumsum(mean_diff)
  # hard code for debugging
  #cum_diff = rep(0, length(cum_diff))
  ####
  for(i in 1:n_update){
    for(j in 0:(i-1)){
      samp = sample_interval(data=NULL, fend, fend, start_increment = j*update_interval+1, end_increment = (j+1)*update_interval)
      conca[samp, i+1] = unconca[samp, j+2] - cum_diff[j+1]
      samp = sample_interval(data=NULL, fend, fend, start_increment = i*update_interval+1, end_increment = i*update_interval + forecasting_horizon)
      conca[samp, i+1] = unconca[samp, j+2] - cum_diff[j+1]
    }
  }
  return(list('conca'=conca, 'unconca'=unconca, 'cum_diff'=cum_diff))
}
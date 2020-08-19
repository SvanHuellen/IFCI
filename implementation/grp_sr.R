
grp_read = function(filename) {
  # read the grouping txt file, sort it, del all na columns
  grp = read.delim(filename, stringsAsFactors = FALSE)
  grp = grp[order(grp[, 'group']),]
  grp = grp[, colSums(is.na(grp)) != nrow(grp)]
  return(grp)
}

grp_regr = function(data_lag, dependent, grp, samp, lag_length, p_value = .1){
  # produce ols regression coefficients and select the most significant three
  g_multi  = standardize_data(data_lag, insample_interval = samp, outsample_interval = samp)
  X        = as.matrix(g_multi[samp]) 
  y        = dependent[rownames(X)]
  X_name   = unique(gsub('\\..*', '', colnames(g_multi)))
  # initiate outputs
  weight           = matrix(data=0, nrow=ncol(g_multi))
  rownames(weight) = colnames(g_multi)
  grp_name         = unique(grp[,'group'])
  const            = 0
  regr_count       = 0
  for(i in 1:length(grp_name)){
    member_name = unname(unlist(grp[ grp[, 'group']==grp_name[i], ]['indicators']))
    nmember     = length(member_name)
    d1_list     = list()
    for(j in 1:nmember){
      exact_name   = paste(member_name[j],'.',sep='') #regr expression
      d1_list[[j]] = data.frame(X[ , grepl(exact_name, colnames(X), fixed=TRUE) ])
    }
    d1 = do.call(cbind, d1_list)
    # first regression: on all lags of a single independent variable
    d1     = data.frame(y, d1)
    regr1  = lm(d1)
    anova1 = summary(regr1)$coefficients # columns are ('Esitmate', 'Std. Error', 't value', 'Pr(>|t|)')
    # extract lag_length*nmember/2 with smallest p-value (including intercept)
    # n_sel  = lag_length*nmember/2
    nlag   = 3
    if ('(Intercept)' %in% rownames(anova1[ order(anova1[,4])[1:(2*nmember+1)], ])){
      nlag = 2
    }
    X1_sig = c()
    for(j in 1:nmember){
      exact_name   = paste(member_name[j],'.',sep='')
      anova_member = anova1[ grepl( exact_name, rownames(anova1), fixed=TRUE ), ]
      X1_sig       = c(X1_sig, rownames(anova_member[ order(anova_member[,4])[1:nlag], ]))
    }
    # second regression: on the selected variables
    d2     = data.frame(y, X[ ,  X1_sig ])
    regr2  = lm(d2)
    anova2 = summary(regr2)$coefficients
    # filtered by p-values
    X2_sig = rownames(anova2[which( anova2[,4] < p_value ),])
    if ('(Intercept)' %in% X2_sig){
      const  = const + anova2[c('(Intercept)'),'Estimate']
    }
    coeff = X2_sig[! X2_sig %in% c('(Intercept)')]
    if (length(coeff)>0){
      regr_count = regr_count + 1
    }
    # save significant coeff as weight
    weight[coeff,] = anova2[coeff, 'Estimate']
  }
  return(list('weight'=weight, 'const'=const, 'regr_count'=regr_count))
}

grp_leading_filter = function(data_lag, dependent, grp, sample, lag_length, p_value = .1){
  # a moving window implementation to produce weight and factor
  w_ts           = matrix( nrow = dim(data_lag)[2], ncol = dim(sample)[1] )
  rownames(w_ts) = colnames(data_lag)
  colnames(w_ts) = paste('weight', 1:ncol(w_ts)-1, sep='')
  factor_ts      = list()
  const_ts       = c()
  regr_count_ts   = c()
  for (i in 1:dim(sample)[1]){
    # use grp_regr
    regr_res       = grp_regr( data_lag = g_lag, dependent, grp, sample[i, 1], lag_length, p_value )
    #######
    w_ts[, i]      = regr_res$weight
    g_multi        = standardize_data( data_lag, insample_interval = sample[i, 1], outsample_interval = sample[i, 2])
    factor         = as.matrix(g_multi[sample[i, 2]]) %*% w_ts[,i]
    factor         = factor/regr_res$regr_count + regr_res$const/regr_res$regr_count
    factor         = as.xts(factor)
    # standardize and rescale to target?
    #factor         = standardize_data(factor)*sd(dependent[sample[i, 2]])+mean(dependent[sample[i, 2]])
    factor_ts[[i]] = factor
    const_ts       = c(const_ts, regr_res$const)
    regr_count_ts  = c(regr_count_ts, regr_res$regr_count)
  }
  # flatten the list into 2d
  nrow                = max(sapply(factor_ts, length))
  factor_ts           = do.call(cbind,lapply(factor_ts, 'length<-', nrow))
  colnames(factor_ts) = paste('fci', 1:ncol(factor_ts)-1, sep='')
  return(list('factor_ts'=factor_ts, 'weight_ts'=w_ts, 'const_ts'=const_ts, 'regr_count_ts'=regr_count_ts))
}

grp_weight = function(weight, grp){
  grp_lag = matrix( nrow = nrow(weight), ncol = 1 )
  rownames(grp_lag) = rownames(weight)
  colnames(grp_lag) = 'group'
  for(i in 1:nrow(grp)){
    grp_lag[grepl(grp[i,'indicators'], rownames(grp_lag)),'group']=grp[i, 'group']
  }
  return(cbind(weight, grp_lag))
}
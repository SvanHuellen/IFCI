# core pls algo
pls1 = function(X, y, l=1){
  # copy from wikipedia
  # initialization
  X = as.matrix(X)
  y = as.matrix(y)
  X[is.na(X)] = 0
  y[is.na(y)] = 0
  # create iterative variables
  w = list()
  p = list()
  q = list()
  # iteration for k^th component
  w[[1]] = t(X)%*%y/norm(t(X)%*%y, '2')
  for ( k in 1:l ) {
    t = X%*%w[[k]]
    tt = as.numeric( t(t)%*%t )
    t = t/tt
    p[[k]] = t(X)%*%t
    q[[k]] = as.numeric( t(y)%*%t )
    if (q[[k]] == 0){
      l = k
      break
    }
    if (k<l){
      X = X - tt*t%*%t(p[[k]])
      w[[k+1]] = t(X)%*%y
    }
  }
  # output
  w  = do.call(cbind, w)
  p  = do.call(cbind, p)
  q  = do.call(cbind, q)
  b  = w %*% solve( t(p) %*% w ) %*% t(q)
  b0 = q[1] - t(p[,1]) %*% b
  return(rbind(b0,b))
}

# core pls algo
matpls1 = function(X,y,ncomp=1){
  # reproduce from eviews
  # initialization
  X = as.matrix(X)
  y = as.matrix(y)
  X[is.na(X)] = 0
  y[is.na(y)] = 0
  # create iterative variables
  w = list()
  p = list()
  t = list()
  # iteration for k^th component
  w[[1]] = t(y) %*%X / norm(t(y)%*%X, '2')
  for ( k in 1:ncomp ) {
    w[[k]] = t(y) %*%X / norm(t(y)%*%X, '2')
    t[[k]] = X %*% t( w[[k]] )
    tt = as.numeric( t(t[[k]]) %*% t[[k]] )
    p[[k]] = t(t[[k]]) %*% X / tt
    if (k<ncomp){
      # b: yloadings
      b = as.numeric( t(y) %*% t[[k]] / tt )
      X = X - t[[k]] %*% p[[k]]
      y = y - b * t[[k]]
    }
  }
  # output
  t = t(do.call('rbind', t))
  w = t(do.call('rbind', w))
  p = t(do.call('rbind', p))
  return(list('t'=t,'w'=w,'p'=p))
}

modifiedW = function(weight, loading, ncomp, n=dim(weight)[1]){
  # the raw weight is corresponding to the deflated independent variables (weight)
  # the modified weight is corresponding to the original independent variables (loading)
  param  = list()
  target = weight
  param[[1]] = diag(n)
  for(i in 1:ncomp){
    param[[i+1]] = param[[i]] - param[[i]] %*% weight[,i] %*% t(loading[,i])
    target[,i] = param[[i]] %*% weight[,i]
  }
  return(target)}

sparseW = function(weight, ncomp, maxLag){
  # pick the lag with largest loading
  maxLag   = maxLag + 1
  vec_temp = numeric(maxLag)
  for(c in 1:ncomp){
    for(i in 1:dim(weight)[1]){
      if(i%%(maxLag) == 0){
        for(j in 1:maxLag){
          vec_temp[j] = weight[i-maxLag+j, c]
        }
        temp = weight[i-maxLag+1, c]
        if(temp<0){
          ctrl = min(vec_temp)
        }
        if(temp>0){
          ctrl = max(vec_temp)
        }
        for(j in 1:maxLag)
          if(weight[i-maxLag+j,c]!=ctrl){
            weight[i-maxLag+j,c]=0
          }
      }
    }
  }
  #weight = weight/norm(weight, '2')
  return(weight)}

insig_stats = function(weight, ncomp){
  # label the insignificant independent variables
  filteredwgt       = list()
  insig_indicators  = list()
  insig_indicators1 = list()
  # filteredwgt[[1]]       = weight[weight[,1]!=0, 1]
  # insig_indicators[[1]]  = filteredwgt[[1]][(filteredwgt[[1]]<0.05)& (filteredwgt[[1]]>-0.05)]
  if(ncomp > 0){
    for(c in 1:ncomp){
      filteredwgt[[c]]       = weight[weight[,c]!=0, c]
      insig_indicators[[c]]  = filteredwgt[[c]][(filteredwgt[[c]]<0.05)& (filteredwgt[[c]]>-0.05)]
      insig_indicators1[[c]] = filteredwgt[[c]][(filteredwgt[[c]] < quantile(filteredwgt[[c]],0.6)) & 
                                                  (filteredwgt[[c]] > quantile(filteredwgt[[c]],0.4))]
    }
  }
  return(list('w_filtered'=filteredwgt, 'insig_indicators'=insig_indicators, 
              'insig_indicators1'=insig_indicators1))}

all_matpls1 = function(data, dependent, in_samp, out_samp=in_samp, ncomp=1, method='eviews'){
  # an aggregate function corresponding to sub_all_matpls1
  # standardize variables
  #y = as.matrix(dependent[in_samp]) # use this line if no standardization
  y = standardize_data(data = dependent, insample_interval=in_samp, outsample_interval=in_samp)
  X = standardize_data(data = data, insample_interval=in_samp, outsample_interval=in_samp)
  y = as.matrix(y)
  X = as.matrix(X)
  if(method=='plsr'){
    packages_check('pls')
    library(pls)
    # plsr package
    plsreg1 = plsr(y~X, ncomp=ncomp)
    p = unclass(loadings(plsreg1))
    w = unclass(loading.weights(plsreg1))
    t = unclass(scores(plsreg1))
  }
  else if(method=='eviews'){
    # plsr from EViews
    pls1 = matpls1(X,y, ncomp = ncomp)
    p = pls1$p
    w = pls1$w
    t = pls1$t  
  }
  else{
    print('method not implemented.')
    return()
  }
  
  # calculate modified weights
  w_modified = modifiedW(w, p, ncomp)
  
  # generate fixed-weight factors
  X2 = standardize_data(data = data, insample_interval=in_samp, outsample_interval=out_samp)
  #factor = X2 %*% w_modified # not used
  
  # generate sparse weight and factors
  w_sparse = sparseW(w_modified, ncomp, maxLag)
  sparse_factor = xts(x = X2 %*% w_sparse, order.by = index(X2))
  sparse_factor = standardize_data(data = sparse_factor, insample_interval=in_samp, outsample_interval=out_samp)
  # scaling
  m1 = mean(dependent[in_samp], na.rm=TRUE)
  m2 = sd(dependent[in_samp], na.rm=TRUE)
  scaled_factor = m1+m2*sparse_factor
  # insig stats
  insig_stats       = insig_stats(w_sparse, ncomp)
  w_filtered        = insig_stats$w_filtered
  insig_indicators  = insig_stats$insig_indicators
  insig_indicators1 = insig_stats$insig_indicators1
  
  # output
  return(list('loadings'=p, 'weight'=w, 'scores'=t, 'w_modified'=w_modified, 'w_sparse'=w_sparse,
              'sparse_factor'=sparse_factor, 'scaled_factor' = scaled_factor, 'w_filtered'=w_filtered,
              'insig_indicators'=insig_indicators, 'insig_indicators1'=insig_indicators1))
}

mov_win_pls = function(data, dependent, fbgn, fend, forecasting_horizon, update_interval, n_update, maxLag, ncomp=1){
  # a moving window implementation of all_matpls1
  pls = list()
  for(i in 1:(n_update+1)){
    in_samp  = sample_interval(data, fbgn, fend, start_increment=maxLag, end_increment = (i-1)*update_interval)
    out_samp = sample_interval(data, fbgn, fend, start_increment=maxLag, end_increment = (i-1)*update_interval+forecasting_horizon)
    pls[[i]] = all_matpls1(data, dependent, in_samp, out_samp, ncomp) 
  }
  return(pls)
}

conca = function(ls, fend, update_interval, forecasting_horizon, attr = 'scaled_factor'){
  # concatenate the factor and align a single output the list of pls results
  output  = ls[[1]][[attr]]
  factor0 = output
  if(length(ls)>1){
    for(i in 2:length(ls)){
      factor = factor0
      for(j in 1:(i-1)){
        samp         = sample_interval(data=NULL, fend, fend, (j-1)*update_interval+1, j*update_interval+forecasting_horizon)
        factor       = c(factor, ls[[j+1]][[attr]][samp]) # need an inner append
        factor       = factor[!duplicated(index(factor))]
        factor[samp] = ls[[j+1]][[attr]][samp]
      }
      # warnings when merge xts
      output           = merge(x = output, y = factor, suffixes = c(i-1, i), by = 'row.names', all = TRUE)
      rownames(output) = output$row.names
      output$row.names = NULL
    }
  }
  colnames(output) = paste(attr, 1:length(ls)-1, sep='')
  return(output)
}

flatten_wgt = function(ls, attr='weight'){
  # flatten a list
  output  = ls[[1]][[attr]]
  if(length(ls)>1){
    for(i in 2:length(ls)){
      output           = merge(x = output, y = ls[[i]][[attr]], suffixes = c(i-1,i), by = 'row.names', all = TRUE)
      rownames(output) = output$Row.names
      output$Row.names = NULL
    }
  }
  colnames(output) = paste(attr, 1:length(ls)-1, sep='')
  return(output)
}

weight_screening = function(weight, thresold=0.01){
  # filtered the variables by the raw weight matrix: 1 for significant
  X_name           = unique(gsub('\\..*', '', rownames(weight)))
  output           = matrix(nrow = length(X_name), ncol = length(weight))
  rownames(output) = X_name
  for(i in 1:length(weight)){
    for(j in 1:length(X_name)){
      group_row = grep(X_name[j], rownames(weight))
      for(k in group_row){
        output[X_name[j], i] = 0
        if(abs(weight[k, i])<=thresold){
          output[X_name[j], i] = 1
          next
        }
      }
    }
  }
  return(output)
}

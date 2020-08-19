rm(list=ls())

library(TSclust)
library(readxl)
library(NbClust)
library(reshape2) # install.packages("reshape2")
library(ggplot2)

market_list = c('money') #'forex','future',,'bank','bond','stock'
sample_period_list = c('1992-2007',"2003-2017",
                       '1992-1996_1999-2007',
                       '2003-2006_2010-2017') #
# 

for (m in 1:length(market_list)){
  market = market_list[m]
  
  setwd('/Users/shine/我的坚果云/macro_index秦老师/1-201809-11/grouping201809') ##where the input data
  filename = paste0(market, "_group_result.txt")
  con <- file(filename, "w")
  
  for (s in 1:length(sample_period_list)){
    sample_period = sample_period_list[s]
    
    cat(market, sample_period,file=con,
        sep="\n",append=TRUE)
    
    print (market)
    print (sample_period)
    
    setwd('/Users/shine/我的坚果云/macro_index秦老师/1-201809-11/grouping201809') ##where the input data
    df = read_excel(paste0(sample_period, market,".xlsx"))
    
    cat(paste('before dropping number of rows',toString(nrow(df))),
        file=con,
        sep="\n",append=TRUE)
    
    print ('before dropping')
    print (nrow(df))
    
    if (max(colMeans(is.na(df))) > 0.5){
      df = df[, -which(colMeans(is.na(df)) > 0.5)]
    }
    #delete columns that have too many nans

    df = df[complete.cases(df), ]
    
    rnames = df[,1]
    df = df[,-1]
    
    print ('after dropping')
    print (ncol(df))
    print (nrow(df))

    cat(paste('after dropping number of columns',toString(ncol(df))),
        paste('number of rows',toString(nrow(df))),
        file=con,
        sep="\n",append=TRUE)
      
    df = as.data.frame(scale(df))
    
    dmat = diss(df,'COR')
    
    indexname = c("kl", "ch", "hartigan", "cindex", 
                  "db", "silhouette", "duda", "pseudot2", 
                  "ratkowsky", "ball", "ptbiserial", "gap", 
                  "mcclain", "gamma", "gplus", "tau", "dunn", 
                  "sdindex",  "sdbw")
    
    k_list = c()
    for (i in 1:length(indexname)){
      nb <- NbClust(t(df), diss = dmat, distance = NULL, min.nc = 2,
                    max.nc = round(ncol(df)/2), method = "centroid", index=indexname[i])
      k_list[i] = as.numeric(nb$Best.nc[1])
      #     print (indexname[i])
    }
    
    freq = table(k_list)
    print (freq)
    k = min(as.numeric(names(freq[freq==max(freq)])))
    print (k)
    
    hc.dpred = hclust(dmat, method='centroid')
    
    cluster = data.frame(cutree(hc.dpred, k=k))
    for (i in 1:k){
      print (paste('group',toString(i)))
      print (rownames(subset(cluster, cluster[,1]==i)))
      
      cat(paste('Group',toString(i)),
          rownames(subset(cluster, cluster[,1]==i)),
          file=con,
          sep="\n",append=TRUE)
    }
    
    #plot by group
    for (i in 1:k){
      df1 = subset(df, select=rownames(subset(cluster, cluster[,1]==i)))
      dd = data.frame(cbind(df1, rowMeans(df1)))
      dd$index = 1:nrow(dd)
      meltR = melt(dd, id = 'index')
      ggplot(meltR, aes(x = index, y = value, group = variable, colour = variable)) + 
        geom_line() + 
        geom_point(aes(shape=variable)) +
        ggtitle(paste(market,'cluster',toString(i),sample_period))+
        theme(legend.position="bottom") +
        guides(col = guide_legend(nrow = 8))
      setwd('/Users/shine/Desktop/png')
      ggsave(paste(market,'cluster',toString(i),sample_period,'.png'))
    }
  }
}


rm(list=ls())

library(readxl)
library(writexl)

market_list = c('money','forex','future','bank','bond','stock') #,

### average to form grouped indicators
setwd('/Users/shine/我的坚果云/macro_index/1-201809-11/grouping201809')
df = read_excel("groups_results_summary.xlsx")

for (m in 1:length(market_list)){
  marketname = market_list[m]
  
  setwd('/Users/shine/我的坚果云/macro_index/1-201809-11/grouping201809')
  dd = read_excel(paste0(marketname,"_market_group.xlsx"))
  
  rownames(dd) = dd$X__1
  
  df1 = subset(df, market==marketname)
  
  grouped = data.frame(matrix(ncol = max(df1$group_version1), nrow = nrow(dd)))
  rownames(grouped) = dd$X__1
  max_group = max(df1$group_version1)
  for (i in 1:max_group){ 
    #if indicator group=0, it naturelly drops out as we start loops from 1
    indicator_list = df1$indicators[df1$group_version1==i]
    # don't know why here need to convert the dataframe, otherwise the data type is not right...
    # stuck here for a long time...

    if (length(unique(indicator_list))==1){
      grouped[,i] = subset(dd, select=indicator_list)
    }else{
      a = subset(dd, select=indicator_list)
      a[] = lapply(a, as.character)
      a[] = lapply(a, as.numeric)
      grouped[,i] = rowMeans(a)
    }
  }
  colnames(grouped) = sapply(seq(1,max_group), 
                             function(i) paste0(marketname,'_g',toString(i)))
  grouped$date = row.names(grouped)
  write_xlsx(grouped, paste0(marketname,'_grouped.xlsx'))
  print (nrow(grouped))
}

# could manully combine the multiple markets grouped series into one xlsx file



## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ZINQL)

data(data_example)

head(metadata)

head(otu.tab[,1:5])

## -----------------------------------------------------------------------------
y.Blautia = otu.tab[,'Blautia']
re = ZINQL_fit(y=y.Blautia, meta=metadata, formula=y~Abx+age+(1|study),
          C='Abx', seed=2024, taus=c(0.1, 0.25, 0.5, 0.75, 0.9),
          method = 'Both')
re

## -----------------------------------------------------------------------------
y.Blautia = otu.tab[,'Blautia']
y.Blautia[1:28] = 0
re = ZINQL_fit(y=y.Blautia, meta=metadata, formula=y~Abx+age+(1|study),
          C='Abx', seed=2024, taus=c(0.1, 0.25, 0.5, 0.75, 0.9),
          method = 'Both')
re

## -----------------------------------------------------------------------------
y.Blautia = otu.tab[,'Blautia']

#MinP
re = ZINQL_fit(y=y.Blautia, meta=metadata, formula=y~Abx+age+(1|study),
          C='Abx', seed=2024, taus=c(0.1, 0.25, 0.5, 0.75, 0.9),
          method = 'MinP')
re

#truncted Cauchy
re = ZINQL_fit(y=y.Blautia, meta=metadata, formula=y~Abx+age+(1|study),
          C='Abx', seed=2024, taus=c(0.1, 0.25, 0.5, 0.75, 0.9),
          method = 'Cauchy')
re

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
#In this segment we store important results from ZINQ-L
model.use = c()
set.seed(2024)
## quantile test 
taus = c(0.1, 0.25, 0.5, 0.75, 0.9)
d.result = data.frame(`ZINQL_MinP`=1, `ZINQL_Cauchy`=1)
model.use = c()
for(t in 1:ncol(otu.tab)){
  print(paste('taxa', t, colnames(otu.tab)[t]))
  yt = otu.tab[, t]
  result.full = suppressMessages(try(
    ZINQL_fit(y=yt,formula=y~Abx+age+(1|study),C='Abx',
                          taus=taus, seed=2024, 
                          meta=metadata, method='Both')
  ))
  if ('try-error' %in% class(result.full)){
    result = rep(NA, 2)
    model.use = c(model.use,NA)
  } else{
    result = result.full$`Final_P_value`
    model.use = c(model.use,result.full$`model`)
  }
  d.result = rbind(d.result,result)
}

d.result = d.result[2:nrow(d.result),]
d.result$model.use = model.use

#we should filter out tests with model.use as None
d.result = d.result[which(!(d.result$model.use=='None')),]

#conduct FDR
d.result$ZINQL_MinP_FDR = p.adjust(d.result$ZINQL_MinP, method = 'fdr')
d.result$ZINQL_Cauchy_FDR = p.adjust(d.result$ZINQL_Cauchy, method = 'fdr')

## -----------------------------------------------------------------------------
head(d.result)


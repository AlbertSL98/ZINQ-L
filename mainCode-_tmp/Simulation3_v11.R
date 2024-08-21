library(tidyverse)
library(Matrix)
library(phyloseq)
library(GUniFrac)
library(quantreg)
library(parallel)
library(lmerTest)

#
#sample
#Rscript /dcs05/hongkai/data/shuai/ZINQL/simulation3_4_1/code/Mysim3_ktx_clean_2024_4_1.R -package_dir '/dcs05/hongkai/data/shuai/ZINQL/code/util3.R' -ni 5 -m 100 -s 1 -seed 1 -outdir /users/sli1/tmp/ZINQL_test 
#ni [5, 10, 15, 20]
#m [50, 100, 200]
#
# Get command-line arguments
#/dcs05/hongkai/data/shuai/ZINQL/simulation3_4_1/code/util2.R
#----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
print(args)
#
# Function to parse arguments
parse_arguments <- function(args) {
  options <- list()
  i <- 1
  while (i <= length(args)) {
    if (substring(args[i], 1, 1) == "-") {
      option <- substring(args[i], 2)
      if (i < length(args) && substring(args[i + 1], 1, 1) != "-") {
        value <- args[i + 1]
        options[[option]] <- value
        i <- i + 1
      } else {
        options[[option]] <- TRUE
      }
    }
    i <- i + 1
  }
  return(options)
}
# Parse command-line arguments
parsed_args <- parse_arguments(args)
#----------------------------------------------------------------

#this is where we have RunZhe's package
tmp = parsed_args[['package_dir']]
package_dir = ifelse( is.null(tmp), '/dcs05/hongkai/data/shuai/ZINQL/code/util5_ignore_NA-no_gee-logit_no_allone.R' ,tmp )
print( paste('Importing source code',package_dir) )
source(package_dir)

#number of longti
#this is where we have RunZhe's package
tmp = parsed_args[['ni']]
ni = as.numeric( ifelse( is.null(tmp), 5 ,tmp ) )
print( paste('ni',ni) )

#number of samples generated from data, your final row number will be m * ni
#this is where we have RunZhe's package
tmp = parsed_args[['m']]
m = as.numeric( ifelse( is.null(tmp), 50 ,tmp ) )
print( paste('m',m) )

#whether dense have signal or rare have signal
#1 means dense have signal
#this is where we have RunZhe's package
tmp = parsed_args[['s']]
scene = as.numeric( ifelse( is.null(tmp), 1 ,tmp ) )
print( paste('scene',scene) )


#outdir
tmp = parsed_args[['outdir']]
outdir = ifelse( is.null(tmp), '/dcs05/hongkai/data/shuai/ZINQL/simulation3_4_1/test' ,tmp )
print( paste('outdir', outdir) )

#seed
tmp = parsed_args[['seed']]
seed = as.numeric( ifelse( is.null(tmp), 2024 ,tmp ) )
print( paste('seed', seed) )

#datadir
tmp = parsed_args[['datadir']]
datadir = ifelse( is.null(tmp), '/dcs05/hongkai/data/shuai/ZINQL/data' ,tmp ) 
print( paste('datadir', datadir) )


if( !(dir.exists(outdir)) ) {
  print(paste('mkdir',outdir))
  dir.create(outdir)
}

#datadir = '/dcs05/hongkai/data/shuai/ZINQL/data'
## load taxa category data
taxa_cat = readRDS( file.path(datadir, 'ktx_dense_group.rds') )
#taxa_cat = readRDS("H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\ktx_dense_group.rds")
#
dense.taxa = taxa_cat %>% filter(dense==1) %>% pull(taxa)
rare.taxa = taxa_cat %>% filter(dense==0) %>% pull(taxa)


#load meta data
#meta = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\ktx_meta.rds')
meta = readRDS( file.path(datadir, 'ktx_meta.rds') )

#process meta
meta = meta %>%
  mutate(AGE = scale(AGE)[,1]) 

  
## rarify
#otu.tab.rff = readRDS("H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\ktx_otu_rff.rds")
otu.tab.rff = readRDS(  file.path(datadir, 'ktx_otu_rff.rds')  )

## step 1: logit reg
finished = 1
if(finished==0){
  logit_coef = lapply(1:ncol(otu.tab.rff), function(j){
    df = data.frame(
      y = otu.tab.rff[, j],
      Abx = meta$Abx_noAbx, 
      age = meta$AGE) %>%
      mutate(ybin = as.integer(y > 0))
    logit.fit = glm(ybin ~ ., 
                    data =df %>% dplyr::select(-y), 
                    family = binomial(link = "logit"))
    logit.coef = coef(logit.fit)
  })
  names(logit_coef) = colnames(otu.tab.rff)
  saveRDS(logit_coef, "H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\logit_coef.rds")
} else{
  #logit_coef = readRDS("H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\logit_coef.rds")
  logit_coef = readRDS(  file.path(datadir, 'logit_coef.rds')  )
}




## step 2: quantile reg
finished = 1
if (finished == 0){
  set.seed(2024)
  taus = seq(0.01, 0.99, by = 0.01)
  quantile_coef = lapply(1:ncol(otu.tab.rff), function(j){
    #print(j)
    df_nz = data.frame(
      y = otu.tab.rff[, j],
      Abx = meta$Abx_noAbx, 
      age = meta$AGE) %>%
      filter(y > 0)
    
    #!!!!!!!!!!!!!!!!!!!!!!!otherwise will Error in rq.fit.br(x, y, tau = tau, ...) : Singular design matrix
    df_nz$y = dither(df_nz$y, type = "right", value = 1)
    quantile.fit = rq(y ~ ., data = df_nz, tau = taus)
    quantile.coef = coef(quantile.fit); coef.name = rownames(quantile.coef)
    quantile.coef = lapply(1:nrow(quantile.coef), function(i) approxfun(taus,quantile.coef[i,]))
    names(quantile.coef) = coef.name
    quantile.coef
  })
  names(quantile_coef) = colnames(otu.tab.rff)
  saveRDS(quantile_coef, "H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\quantile_coef.rds" )
} else{
  #quantile_coef = readRDS("H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\quantile_coef.rds")
  quantile_coef = readRDS(  file.path(datadir, 'quantile_coef.rds')  )
}


zerorate = apply(otu.tab.rff, 2, function(x) mean(x == 0))

#scene = 1
taxa.list = taxa_cat$taxa
if(scene==1){
  #dense taxa has signal
  diff_taxa_idx = taxa.list[which(taxa_cat$dense==1)]
  null_taxa_idx = taxa.list[which(taxa_cat$dense==0)]
} else{
  diff_taxa_idx = taxa.list[which(taxa_cat$dense==0)]
  null_taxa_idx = taxa.list[which(taxa_cat$dense==1)]
}

df = data.frame(Abx = meta$Abx_noAbx, 
                age = meta$AGE)

set.seed(seed)
taus = c(0.1, 0.25, 0.5, 0.75, 0.9)

## step 1: sampling
resample_idx = sample(nrow(otu.tab.rff), size = m, replace = T)

## step 2: generate covariates and outcome
# simply repeat the resampled 'person', manually expand their age, and add random error to their PA and DS
df_longi = lapply(resample_idx, function(id){
  Abx = rep(df$Abx[id],each = ni)
  age = seq(from = df$age[id], to = df$age[id] + 0.1 * (ni-1), by = 0.1)
  data.frame(1, Abx, age)
}) %>% do.call(rbind, .)
#scale age
df_longi$age = scale(df_longi$age)[,1]
study = factor(rep(c(1:m),each=ni))


## step3:generate null
# each person's error is the same
hi = rnorm(m, 0, 1)
h = rep(hi, each = ni)
ybin_null = lapply(null_taxa_idx, function(j){
  #print(j)
  logit_coef_j = logit_coef[[j]]
  logit_coef_j['Abx'] = 0
  eta = as.matrix(df_longi) %*% logit_coef_j + h
  eta_j = log((1-zerorate[j])/zerorate[j])
  eta = scale(eta) + eta_j
  probs = exp(eta)/(1+ exp(eta))
  probs[is.nan(probs)] = 1 # if exp(eta) = Inf
  ybin = sapply(probs, function(prob) rbinom(1, 1, prob = prob))
}) %>% do.call(cbind, .)
colnames(ybin_null) = null_taxa_idx

### generate y
u = runif(m*ni, 0.01, 0.99)
hi = rnorm(m, 0, 1)
h = rep(hi, each = ni)

yij_null = lapply(null_taxa_idx, function(j){
  #print(j)
  func_j = quantile_coef[[j]]
  quantile_coef_j = lapply(func_j, function(f) f(u)) %>%
    do.call(cbind, .)
  quantile_coef_j[,'Abx'] = 0
  #this is equivalent to a matrix multi
  yij = round(rowSums(as.matrix(df_longi) * quantile_coef_j) + h)
}) %>% do.call(cbind, .)
yij_null[yij_null < 0] = 0
colnames(yij_null) = null_taxa_idx
yij_null = yij_null * ybin_null

## step 4: generate alt
hi = rnorm(m, 0, 1)
h = rep(hi, each = ni)
ybin_alt = lapply(diff_taxa_idx, function(j){
  #print(j)
  logit_coef_j = logit_coef[[j]]
  eta = as.matrix(df_longi) %*% logit_coef_j + h
  eta_j = log((1-zerorate[j])/zerorate[j])
  eta = scale(eta) + eta_j
  probs = exp(eta)/(1+ exp(eta))
  probs[is.nan(probs)] = 1 # if exp(eta) = Inf
  ybin = sapply(probs, function(prob) rbinom(1, 1, prob = prob))
}) %>% do.call(cbind, .)
colnames(ybin_alt) = colnames(otu.tab.rff)[diff_taxa_idx]

### generate y
u = runif(m*ni, 0.01, 0.99)
hi = rnorm(m, 0, 1)
h = rep(hi, each = ni)

yij_alt = lapply(diff_taxa_idx, function(j){
  func_j = quantile_coef[[j]]
  quantile_coef_j = lapply(func_j, function(f) f(u)) %>%
    do.call(cbind, .)
  yij = round(rowSums(as.matrix(df_longi) * quantile_coef_j) + h)
}) %>% do.call(cbind, .)
yij_alt[yij_alt < 0] = 0
colnames(yij_alt) = diff_taxa_idx
yij_alt = yij_alt * ybin_alt

taxa_reorder = c(null_taxa_idx, diff_taxa_idx)
yij = cbind(yij_null, yij_alt)
ybin = cbind(ybin_null, ybin_alt)


otu.tab = yij
metadata = df_longi
metadata$study = study
rownames(otu.tab) = rownames(metadata) = rownames(ybin) = paste0("subject", rep(1:m, each=ni), "_visit", rep(1:ni, m))

#save the data for downstream test
data.out.dir = file.path(outdir,'data_sim.RData')
save(otu.tab, metadata, ybin, file=data.out.dir)
#load(data.out.dir)

# study = metadata$study
# X = metadata$Abx
# Z = as.matrix(metadata %>% dplyr::select(-Abx, -study))

## ZINQL
taus = c(0.1, 0.25, 0.5, 0.75, 0.9)
#result = lapply(1:ncol(otu.tab), function(t){
pval.ZINQL = list()
for(t in 1:ncol(otu.tab)){
  print(paste('taxa', t))
  yijt = otu.tab[, t]
  result = suppressMessages(try(
    ZINQL_fit(formula=y~Abx+age+(1|study), y=yijt, meta=metadata, C='Abx', taus = taus, seed=2024, method='both', n.positive.cut=5)
  ))
  if ('try-error' %in% class(result)){
    pval.ZINQL[[t]] = rep(NA, 2)
    names(pval.ZINQL[[t]]) = c('ZINQL_MinP', 'ZINQL_Cauchy')
  } else if(result$model=='None'){
    pval.ZINQL[[t]] = rep(NA, 2)
    names(pval.ZINQL[[t]]) = c('ZINQL_MinP', 'ZINQL_Cauchy')
  } else{
    pval.ZINQL[[t]] = result$Final_P_value 
  }
}
pval.ZINQL = do.call(rbind, pval.ZINQL)
rownames(pval.ZINQL) = colnames(otu.tab)


#-------------------------------------------------------------------
## linear regression and linear mixed model
pval.lmm = list()
for(t in 1:ncol(otu.tab)){
  yijt = otu.tab[, t]
  lmfit = lm(yijt ~ Abx+age ,data=metadata)
  pval_lm = summary(lmfit)$coefficients['Abx', 'Pr(>|t|)']
  #
  lmmfit = suppressMessages(lmerTest::lmer(yijt ~ Abx + age + (1 | study),data=metadata))
  pval_lmm = summary(lmmfit)$coefficients['Abx', 'Pr(>|t|)']
  result_bench = c(pval_lm, pval_lmm); names(result_bench) = c('lm','lmm')
  pval.lmm[[t]] = result_bench
}
#turn it from a list to a dataframe
pval.lmm = do.call(rbind, pval.lmm)
rownames(pval.lmm) = colnames(otu.tab)


## existing methods
tmp = suppressWarnings(
  suppressMessages(
    try(
      LinDA_test(otu.tab, metadata)$pval
    )
  )
)
if('try-error' %in% class(tmp)){pval_linda = rep(NA,ncol(otu.tab))}
if(!('try-error' %in% class(tmp))){pval_linda =tmp}

#test
#---------------------------------------------------------------
#do not solve
# index = which(rowSums(otu.tab)!=0)
# otu.tab2 = otu.tab[index,]
# metadata2=metadata[index,]
# LinDA_test(otu.tab2, metadata2)$pval

#---------------------------------------------------------------

outdir_ma = file.path(outdir,'MaOut')
#
tmp = suppressWarnings(
  suppressMessages(
    try(
      MaAsLin2_test(otu.tab, metadata, output = outdir_ma))$pval
    )
)
if('try-error' %in% class(tmp)){pval_ma = rep(NA,ncol(otu.tab))}
if(!('try-error' %in% class(tmp))){pval_ma =tmp}


#tmp = suppressWarnings(
#  suppressMessages(
#    try(
#      LDM_test(otu.tab, metadata, n.perm.max=5000)$pval
#    )
#  )
#)
#tmp = try(
#      LDM_test(otu.tab, metadata, n.perm.max=5000)$pval
#)
#if('try-error' %in% class(tmp)){pval_ldm = rep(NA,ncol(otu.tab))}
#if(!('try-error' %in% class(tmp))){pval_ldm =tmp}

tmp = suppressWarnings(
  suppressMessages(
    try(
      NBZIMM_test(otu.tab, metadata, method = 'nb')$pval
    )
  )
)
if('try-error' %in% class(tmp)){pval_nb = rep(NA,ncol(otu.tab))}
if(!('try-error' %in% class(tmp))){pval_nb =tmp}

tmp = suppressWarnings(
  suppressMessages(
    try(
      NBZIMM_test(otu.tab, metadata, method = 'zinb')$pval
    )
  )
)
if('try-error' %in% class(tmp)){pval_zinb = rep(NA,ncol(otu.tab))}
if(!('try-error' %in% class(tmp))){pval_zinb =tmp}

tmp = suppressWarnings(
  suppressMessages(
    try(
      NBZIMM_test(otu.tab, metadata, method = 'zig_count')$pval
    )
  )
)
if('try-error' %in% class(tmp)){pval_zig_count = rep(NA,ncol(otu.tab))}
if(!('try-error' %in% class(tmp))){pval_zig_count =tmp}


tmp = suppressWarnings(
  suppressMessages(
    try(
      NBZIMM_test(otu.tab, metadata, method = 'zig_prop')$pval
    )
  )
)
if('try-error' %in% class(tmp)){pval_zig_prop = rep(NA,ncol(otu.tab))}
if(!('try-error' %in% class(tmp))){pval_zig_prop =tmp}


#pval_combo = cbind(pval_linda, pval_ma, pval_nb, pval_zinb, pval_zig_count, pval_zig_prop, pval_ldm)
#colnames(pval_combo) = c('LinDA', 'MaAsLin2', 'NB', 'ZINB', 'ZIG_count', 'ZIG_prop', 'LDM')
pval_combo = cbind(pval_linda, pval_ma, pval_nb, pval_zinb, pval_zig_count, pval_zig_prop)
colnames(pval_combo) = c('LinDA', 'MaAsLin2', 'NBMM', 'ZINBMM', 'ZIGMM', 'ZIG_prop')
pval = cbind(pval.ZINQL, pval.lmm,  pval_combo)
rownames(pval) = colnames(otu.tab)

#rows_with_na <- names( which(rowSums(is.na(pval)) > 0) )
#cols_with_na <- names( which(colSums(is.na(pval)) > 0) )
#a = pval[rows_with_na, cols_with_na]

saveRDS( pval, file.path(outdir,'pval.rds') )

# fdr_result = apply(pval, 2, function(x) p.adjust(x, method = 'fdr')) %>%
#   as.data.frame() %>%
#   mutate(taxa = rownames(.)) %>%
#   mutate(truth = case_when(
#     taxa %in% null_taxa_name ~ 'null',
#     taxa %in% diff_taxa_name ~ 'diff'
#   )) %>%
#   pivot_longer(cols = -c('taxa', 'truth'),
#                names_to = 'Test',
#                values_to = 'fdr') %>%
#   mutate(pred = case_when(
#     fdr < 0.05 ~ 'diff',
#     T ~ 'null'
#   ))
# 
# fdr_summ = fdr_result %>%
#   group_by(Test) %>%
#   summarise(
#     FDP = sum((pred == 'diff') & (truth == 'null'))/sum(pred == 'diff'),
#     TPR = sum((pred == 'diff') & (truth == 'diff'))/sum(truth == 'diff')
#   ) %>% as.data.frame()
# 
# #quantile_tests = c(paste0('indep:', taus), paste0('dep:', taus))
# #print(fdr_summ %>% filter(!Test %in% quantile_tests))
# methods = c("dep:cauchy_unequal",
#             "dep:minP",
#             "indep:cauchy_unequal",
#             "indep:minP",
#             "lm", "lmm", "logit_LRT", 
#             "LinDA", "MaAsLin2", #"LDM",
#             "NB", "ZINB", "ZIG_count", "ZIG_prop", 
#             "Setting")
# print(fdr_summ %>% filter(Test %in% methods))
# print("----------------")
# 
#   
# 
# result = mclapply(1:nsim, simu_func, 
#                   model = model, 
#                   error = error,
#                   setting = setting, 
#                   mc.cores = numCores)






suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(phyloseq)
  library(GUniFrac)
  library(quantreg)
  library(parallel)
  library(lmerTest)
})
#
#sample
#Rscript mySim.R -package_dir /dcs05/hongkai/data/shuai/ZINQL/code/util2.R -ni 5 -m 100 -setting alt -error rnorm -seed 1 -outdir /users/sli1/tmp/ZINQL_test
#ni [5, 10, 15, 20]
#m [50, 100, 200]
#null alt
#
# Get command-line arguments
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
package_dir = ifelse( is.null(tmp), 'C:\\Users\\local_user\\Desktop\\Ktx\\Core Code\\util11.R' ,tmp )
#package_dir = ifelse( is.null(tmp), '/dcs05/hongkai/data/shuai/ZINQL/code/util4.R' ,tmp )
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
m = as.numeric( ifelse( is.null(tmp), 100 ,tmp ) )
print( paste('m',m) )

#alt or null
tmp = parsed_args[['setting']]
setting = ifelse( is.null(tmp), 'null' ,tmp )
print( paste('setting',setting) )

#error term: rnorm, rchisq, rt, rcauchy
tmp = parsed_args[['error']]
error = ifelse( is.null(tmp), 'rnorm' ,tmp )
print( paste('error',error) )

#outdir
tmp = parsed_args[['outdir']]
outdir = ifelse( is.null(tmp), 'C:\\Users\\local_user\\Desktop\\Ktx\\Simulation1\\mysim_7_27_util11\\tmp' ,tmp )
print( paste('outdir', outdir) )
#
if( !(dir.exists(outdir)) ) {
  print(paste('mkdir',outdir))
  dir.create(outdir)
}

#seed
tmp = parsed_args[['seed']]
seed = as.numeric( ifelse( is.null(tmp), 2024 ,tmp ) )
print( paste('seed', seed) )

#datadir
tmp = parsed_args[['datadir']]
datadir = ifelse( is.null(tmp), '/dcs05/hongkai/data/shuai/ZINQL/data' ,tmp ) 
print( paste('datadir', datadir) )

#-------------------------------------------------------------------------
#unread part
outtmpdir = file.path(outdir, 'outTmp')
if( !(dir.exists(outtmpdir)) ) {
  print(paste('mkdir',outtmpdir))
  dir.create(outtmpdir)
}

#sig level
level = 0.05 


#datadir
#load(file = "H:\\Runzhe\\ZINQL\\ZINQ_L-main\\Data\\genus.Rdata")
#datadir = '/dcs05/hongkai/data/shuai/ZINQL/data'
#datadir = 'C:\\Users\\local_user\\Desktop\\Ktx\\Simulation3\\mysim\\data'
otu.tab.rff = readRDS(  file.path(datadir, 'ktx_otu_rff.rds')  )
meta = readRDS( file.path(datadir, 'ktx_meta.rds') )

#make sure the order is the same
otu.tab.rff = otu.tab.rff[meta$sample_ID, ]


set.seed(seed)
taus = c(0.1, 0.25, 0.5, 0.75, 0.9)

# random_err = function(n, error = 'rnorm'){
#   if (error == 'rnorm'){
#     return(rnorm(n))
#   } else if (error == 'rchisq'){
#     return(rchisq(n, df = 2))
#   } else if (error == 'rt'){
#     return(rt(n, df = 2))
#   } else if (error == 'rcauchy'){
#     return(rcauchy(n))
#   }
#   return(0)
# }

#for alt null separation
nonAbx_idx = which(meta$Abx_noAbx == 0)
Abx_idx = which(meta$Abx_noAbx == 1)

#four taxas are run at the same time
idx1 = 'Blautia'
idx2 = 'Dorea'
idx3 = 'Enterococcus'
idx4 = 'Anaerofustis'
#idx = c(idx1, idx2, idx3, idx4) # four representative taxa
idx = c(idx1, idx2, idx3, idx4)
#----------------------------------------------------------------------

#---------------------------------------------------------------------
#simulation generate table
if (setting == 'null'){
  yij = lapply(idx, function(t){
    #print(t)
    y = otu.tab.rff[, t]
    
    ## resampling from non-Abx sample only
    
      #print(i)
    resample_tau = runif(m, 0, 1)
    y = quantile(y[nonAbx_idx], probs = resample_tau, type = 3) + 1
    names(y) = NULL    
    ## step 2: generate longitudinal data
  yij = lapply(y, function(yi){
	#print(yi)
	#this is elementwise operation, here y is a one-dimension vector, and resulting yij is also vector
	yij = exp(log(yi) + rnorm(ni))
	yij =  ifelse(yij < 1, 0, round(yij) - 1)  
	}) %>% do.call(c, .)
	otu_tab = yij
    return(otu_tab)
  }) %>% do.call(cbind, .)
  
  colnames(yij) = idx
  study = factor(rep(c(1:m),each=ni))
  Z = matrix(1, nrow = nrow(yij), ncol = 1); colnames(Z) = 'Z';
  #replicate each value ni times
  #rep(rbinom(3, 1, 0.5), each = 5)
  #[1] 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0
  X = rep(rbinom(m, 1, 0.5), each = ni)  
  ybin = (yij > 0) * 1
}

if (setting == 'alt'){
  yij = lapply(idx, function(t){
    y = otu.tab.rff[, t]
    ## resampling from both Abx and non-Abx samples
      
    resample_tau1 = runif(m/2, 0, 1)
    resample_tau2 = runif(m/2, 0, 1)
    y1 = quantile(y[nonAbx_idx], probs = resample_tau1, type = 3) + 1
    y2 = quantile(y[Abx_idx], probs = resample_tau2, type = 3) + 1
    y = c(y1, y2); names(y) = NULL
      
    ## step 2: generate longitudinal data
	yij = lapply(y, function(yi){
		yij = exp(log(yi) + rnorm(ni))
		yij = ifelse(yij < 1, 0, round(yij) - 1)
	  }) %>% do.call(c, .)	  
    otu_tab = yij
    return(otu_tab)
  }) %>% do.call(cbind, .)
  
  colnames(yij) = idx
  study = factor(rep(c(1:m),each=ni))
  Z = matrix(1, nrow = nrow(yij), ncol = 1); colnames(Z) = 'Z'; 
  X = c(rep(0, m/2 * ni), rep(1, m/2 * ni))
  ybin = (yij > 0) * 1
}

otu.tab = yij
metadata = data.frame(X, Z, study)
rownames(otu.tab) = rownames(metadata) = rownames(ybin) = 
  paste0("subject", rep(1:m, each=ni), "_visit", rep(1:ni, m))
#-------------------------------------------------------------------

data.out.dir = file.path(outdir,'data_sim.RData')
save(otu.tab, metadata, ybin, file=data.out.dir)

#-------------------------------------------------------------------
## ZINQL
taus = c(0.1, 0.25, 0.5, 0.75, 0.9)
#result = lapply(1:ncol(otu.tab), function(t){
pval.ZINQL = list()
for(t in 1:ncol(otu.tab)){
  print(paste('taxa', t))
  yijt = otu.tab[, t]
  result = suppressMessages(try(
	  ZINQL_fit(formula=y~X+(1|study), y=yijt, meta=metadata, C='X', taus = taus, seed=2024, method='both', n.positive.cut=5)
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
  lmfit = lm(yijt ~ X ,data=metadata)
  pval_lm = summary(lmfit)$coefficients['X', 'Pr(>|t|)']
  #
  lmmfit = suppressMessages(lmerTest::lmer(yijt ~ X + (1 | study),data=metadata))
  pval_lmm = summary(lmmfit)$coefficients['X', 'Pr(>|t|)']
  result_bench = c(pval_lm, pval_lmm); names(result_bench) = c('lm','lmm')
  pval.lmm[[t]] = result_bench
}
#turn it from a list to a dataframe
pval.lmm = do.call(rbind, pval.lmm)
rownames(pval.lmm) = colnames(otu.tab)


## existing methods
## return NA is it can not be tested
# tmp = suppressWarnings(
#   suppressMessages(
# 	try(
#     #hush(
#       LinDA_test_unadj(otu.tab, metadata)$pval
#     #)
# 	)
#   )
# )

# if('try-error' %in% class(tmp)){pval_linda=rep(NA,4)}
# if(!('try-error' %in% class(tmp))){pval_linda=tmp}

# tmp = suppressWarnings(
#   suppressMessages(
# 	try(
#     #hush(
#       MaAsLin2_test_unadj(otu.tab, metadata, output = outtmpdir=)$pval
#     #)
# 	)
#   )
# )
# 
# if('try-error' %in% class(tmp)){pval_ma=rep(NA,4)}
# if(!('try-error' %in% class(tmp))){pval_ma=tmp}


# pval_ldm = suppressWarnings(
#   suppressMessages(
#     hush(
#       LDM_test(otu.tab, metadata)$pval
#     )
#   )
# )

pval_nb = c()
for(i in idx){
    p.tmp = suppressWarnings(
      suppressMessages(
    	try( 
    	   NBZIMM_indiv_test(otu.tab[,i], metadata, method = 'nb')
      )
      ))
  
  if('try-error' %in% class(p.tmp)){pval_nb=c(pval_nb,NA)}
  if(!('try-error' %in% class(p.tmp))){pval_nb=c(pval_nb,p.tmp)}
}
names(pval_nb)=idx

pval_zinb = c()
for(i in idx){
  p.tmp = suppressWarnings(
	  try(
	    NBZIMM_indiv_test(otu.tab[,i], metadata, method = 'zinb')
      )
  )
  if('try-error' %in% class(p.tmp)){pval_zinb=c(pval_zinb,NA)}
  if(!('try-error' %in% class(p.tmp))){pval_zinb=c(pval_zinb,p.tmp)}
}
names(pval_zinb)=idx

pval_zig_count = c()
for(i in idx){
  p.tmp = suppressWarnings(
    try(
      NBZIMM_indiv_test(otu.tab[,i], metadata, method = 'zig_count')
      )
  )
  if('try-error' %in% class(p.tmp)){pval_zig_count=c(pval_zig_count,NA)}
  if(!('try-error' %in% class(p.tmp))){pval_zig_count=c(pval_zig_count,p.tmp)}
}
names(pval_zig_count)=idx


pval_zig_prop = c()
for(i in idx){
  p.tmp = suppressWarnings(
    try(
      NBZIMM_indiv_test(otu.tab[,i], metadata, method = 'zig_prop')
      )
  )
  if('try-error' %in% class(p.tmp)){pval_zig_prop=c(pval_zig_prop,NA)}
  if(!('try-error' %in% class(p.tmp))){pval_zig_prop=c(pval_zig_prop,p.tmp)}
}
names(pval_zig_prop)=idx


pval_combo = cbind(pval_nb, pval_zinb, pval_zig_count, pval_zig_prop)
colnames(pval_combo) = c('NBMM', 'ZINBMM', 'ZIGMM', 'ZIG_prop')
pval = cbind(pval.ZINQL, pval.lmm, pval_combo)
rownames(pval) = colnames(otu.tab)
print(str(pval))


#dir.create("./Results")
save.image(file.path(outdir,'simulation_result.rData'))
#
pval.out.dir = file.path(outdir,'pval.rds')
saveRDS(pval, file=pval.out.dir)


#
# a = data.frame(c1=c(1,2,3,4),c2=c(4,5,6,9))
# y = c(2,7,10)
# #column wise operation or elementwise operation
# #return a list
# lapply( 1:ncol(a), function(t){ return(max(a[,t])) } )
# this will use do call to turn list to a vector
# lapply( 1:ncol(a), function(t){ return(max(a[,t])) } ) %>% do.call(c,.)
# lapply( 1:ncol(a), function(t){ return(a[,t]^2) } )
# lapply( y, function(t){ return(t^2) } )
# #this return a vector
# sapply( y, function(t){ return(t^2) } )
# sapply( 1:ncol(a), function(t){ return(max(a[,t])) } )
# sapply( 1:ncol(a), function(t){ return(a[,t]^2) } )
# #this will give a matrix, troublesome
# sapply( c(1,2,3), function(i){i+c(10,11,12)} )

#what is do call do.call(c,.) do
#lapply( c(1,2,3), function(i){i+1} )
#turn it from a list to a vector
#lapply( c(1,2,3), function(i){i+1} ) %>% do.call(c,.)
#try add first col 1, second col 2, etc and cbind 
#a = data.frame(c1=c(1,2,3),c2=c(4,5,6))
#lapply( 1:ncol(a), function(k){y=a[,k]+k} ) %>% do.call(cbind,.)
#if no cbind, return a list
#lapply( 1:ncol(a), function(k){y=a[,k]+k} )

#do call is operations on list, like concate from list to vectors or dataframes 
# a = list()
# a[[1]] = c(1,10,2)
# a[[2]] = c(5,9,2)
# do.call(cbind,a)
# do.call(c,a)



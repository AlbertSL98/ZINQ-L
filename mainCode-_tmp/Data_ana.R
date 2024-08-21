library(tidyverse)
library(Matrix)
library(phyloseq)
library(GUniFrac)
library(quantreg)
library(parallel)
library(lmerTest)
library(UpSetR)
  
package_dir = 'H:\\Runzhe\\ZINQL\\Ktx\\Core Code\\util11.R'
print( paste('Importing source code',package_dir) )
source(package_dir)

load('H:\\Runzhe\\ZINQL\\Ktx\\kx_genus_phyloseq.Rdata')
tax_table = kx_genus@tax_table@.Data
otu.tab = t( as.matrix(kx_genus@otu_table@.Data) )

#for GUniFrac::Rarefy, each column should be a taxa


#mapping
taxa.map = tax_table[,'Genus']
names(taxa.map) = rownames(tax_table)

#manually change several names
#these are xxx, weird
taxa.map['OTU_30;size=32723;'] = 'Lachnospiraceae_unspecified'
taxa.map['OTU_714;size=25;'] = 'Clostridiales_unspecified_unspecified'
taxa.map['OTU_382;size=251;'] = 'Tissierellia_unspecified_unspecified_unspecified'
taxa.map['OTU_4;size=160812;'] = 'Erysipelotrichaceae_unspecified'
taxa.map['OTU_108;size=4439;'] = 'Clostridiales Family XIII. Incertae Sedis_unspecified'
#

#select people that have > 10k count
row.sum = rowSums(otu.tab)
otu.tab = otu.tab[which(row.sum > 10000),]

#get genus name for each column
otu.col = colnames(otu.tab)
otu.row = rownames(otu.tab)

#this part has randomness
set.seed(2024)
otu.tab.rff = GUniFrac::Rarefy(otu.tab)$otu.tab.rff
colnames(otu.tab.rff) = taxa.map[colnames(otu.tab.rff)]
row.names(otu.tab.rff) = otu.row

#saveRDS(otu.tab.rff,'H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\otu_tab_rff.rds')
# 
# a = rowSums(otu.tab)
# a = a[which(a<50000)]
# hist(a,breaks=seq(min(a)-1,max(a)-1+1000,by=1000))
# 
# mean(a<=10000)

# dim(tax_table)
# [1] 259   7
# dim(otu.tab.rff)
# [1] 510 259
# so we have 510 samples? 259 taxa

load('H:\\Runzhe\\ZINQL\\Ktx\\metadata_020524.Rdata')

#----------------------------------------------------------------
#Hanbo Dong's code
metadata <- metadata_020524 %>%
  mutate(sample_ID = rownames(metadata_020524))

# Convert to factors
metadata$RACE <- factor(metadata$RACE)
metadata$DM <- factor(metadata$DM)
metadata$PREOP_ABX <- factor(metadata$PREOP_ABX)
metadata$BACTRIM <- factor(metadata$BACTRIM)
metadata$DDRT <- factor(metadata$DDRT)
metadata$GENDER <- factor(metadata$GENDER)

# Relevel them 
metadata$RACE <- relevel(metadata$RACE, ref = "OTHER")
metadata$DM <- relevel(metadata$DM, ref = "NONE")
metadata$PREOP_ABX <- relevel(metadata$PREOP_ABX, ref = "OTHER")
metadata$BACTRIM <- relevel(metadata$BACTRIM, ref = "OTHER")
metadata$DDRT <- relevel(metadata$DDRT, ref = "L")
metadata$GENDER <- relevel(metadata$GENDER, ref = "MALE")

# Define keywords for each antibiotic group
beta_lactams_keywords <- c("piperacillin tazobactam", "amoxicillin", "ampicillin sulbactam", "amoxicillin clavulanate", "ampicillin", "penicillin v potassium", "ceftriaxone", "cefpodoxime", "cefazolin", "cephalexin", "cefepime", "ceftolozane", "ceftolozane/tazobactam", "cefadroxil",  "meropenem", "ertapenem", "aztreonam")
fluoroquinolones_keywords <- c("levofloxacin", "ciprofloxacin")
other_abx_keywords <- c("vancomycin", "metronidazole", "isoniazid", "linezolid", "azithromycin", "trimethoprim", "nitrofurantoin", "doxycycline", "clindamycin", "amikacin", "daptomycin",  "erythromycin", "fosfomycin")

# Create a dataframe 
Abx_Course_longitudinal <- data.frame(
  Post_Op_Day = metadata$Post_Op_Day,
  Beta_lactams = integer(nrow(metadata)),
  Fluoroquinolones = integer(nrow(metadata)),
  Other = integer(nrow(metadata))
)
row.names(Abx_Course_longitudinal) <- row.names(metadata)

# Process each sample in this loop below
for (i in 1:nrow(metadata)) {
  sample_day <- metadata$Post_Op_Day[i]
  # Initialize binary indicators for antibiotic groups, where 0 is non-exposure
  beta_lactams_use <- 0
  fluoroquinolones_use <- 0
  other_abx_use <- 0
  
  # Check each Abx_Type and Abx_Post_Op_Day pair
  for (j in 1:17) {
    abx_type <- tolower(metadata[i, paste0("Abx_Type", j)])
    abx_day <- metadata[i, paste0("Abx_Post_Op_Day", j)]
    
    # Check abx_type and abx_day are not NA and sample taken day is after abx taken day
    if (!is.na(abx_type) && !is.na(abx_day) &&  (sample_day - abx_day) > 0) {
      # If any keyword is detected, this sample is considered as exposed to the corresponding Abx
      if (any(grepl(paste(beta_lactams_keywords, collapse="|"), abx_type))) {
        beta_lactams_use <- 1
      } else if (any(grepl(paste(fluoroquinolones_keywords, collapse="|"), abx_type))) {
        fluoroquinolones_use <- 1
      } else if (any(grepl(paste(other_abx_keywords, collapse="|"), abx_type))) {
        other_abx_use <- 1
      }
    }
  }
  # 1 is Abx Group, 0 is No Abx Group 
  if (beta_lactams_use == 1 && fluoroquinolones_use == 0) {
    Abx_Course_longitudinal$Abx_Category[i] <- 1  # Beta-lactam
  } else if (fluoroquinolones_use == 1 && beta_lactams_use == 0) {
    Abx_Course_longitudinal$Abx_Category[i] <- 1  # Fluoroquinolone
  } else if (beta_lactams_use == 1 && fluoroquinolones_use == 1) {
    Abx_Course_longitudinal$Abx_Category[i] <- 1  # Beta-lactam and Fluoroquinolone
  } else if (other_abx_use == 1) {
    Abx_Course_longitudinal$Abx_Category[i] <- 1  # Other Abx
  } else {
    Abx_Course_longitudinal$Abx_Category[i] <- 0  # No Abx Group
  }
}

metadata <- cbind(metadata, Abx_noAbx = Abx_Course_longitudinal$Abx_Category)

#because we drop some people, we select from meta data
metadata = metadata %>% filter(sample_ID %in% row.names(otu.tab.rff))

#reorder otu.tab.rff so that order is consistent with meta
otu.tab.rff = otu.tab.rff[metadata$sample_ID,]

saveRDS(otu.tab.rff,'H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\otu_tab_rff.rds')
saveRDS(metadata,'H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\meta.rds')


#read prev processed data
otu.tab.rff = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\otu_tab_rff.rds')
metadata = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\meta.rds')

#ZINQL to identify taxas
taxa.full.list = colnames(otu.tab.rff)



#and other adjusting covariates include age, race, Post-Transplant Day(PostOp), tibiotic prophylaxis (PREOP Day), history of diabetes mellitus(DM), preoperative (preop ABX), and pneumocystis jirovelcii prophylaxis (BACTRIM)

#normalized age
# metadata <- metadata %>%
#   mutate(
#     AGE = scale(AGE)[,1],
#     RACE = ifelse(RACE=='AA',1,0),
#     Post_Op_Day = scale(Post_Op_Day)[,1],
#     DM = ifelse(DM=='DM',1,0),
#     PREOP_ABX = ifelse(PREOP_ABX=='CEFAZOLIN',1,0),
#     BACTRIM = ifelse(BACTRIM=='BACTRIM',1,0),
#     GENDER = ifelse(GENDER=='FEMALE',1,0),
#     DDRT = ifelse(GENDER=='D',1,0)
#   )

# cohort_char = metadata %>%
#   group_by(Abx_noAbx) %>%
#   summarise(
#     count = n(),
#     gender = paste(sum(GENDER == 'FEMALE'), '(', round(mean(GENDER == 'FEMALE') * 100, 1), '%)', sep = ''),
#     PREOP_ABX = paste(sum(PREOP_ABX == 'CEFAZOLIN'), '(', round(mean(PREOP_ABX == 'CEFAZOLIN') * 100, 1), '%)', sep = ''),
#     DDRT = paste(sum(DDRT == 'D'), '(', round(mean(DDRT == 'D') * 100, 1), '%)', sep = ''),
#     BACTRIM = paste(sum(BACTRIM == 'BACTRIM'), '(', round(mean(BACTRIM == 'BACTRIM') * 100, 1), '%)', sep = '')
#   )
# 
# cohort_char = t(cohort_char) %>% as.data.frame()
# colnames(cohort_char) = c('No Antibiotics Group', 'Atibiotics Group')
# cohort_char = cohort_char[-1, ]
# cohort_char = cohort_char[, 2:1]
# 
# Abx_group = metadata %>% filter(Abx_noAbx == 1)
# noAbx_group = metadata %>% filter(Abx_noAbx == 0)
# 
# 
# rownames(cohort_char) <- c("Count", "Female Gender",  "Preoperative Antibiotics (PREOP_ABX)","Deceased Donor Transplantation (DDRT)",  "Bactrim Usage (BACTRIM)")
# 
# kable_table <- kable(cohort_char,  caption = "Table 1")
# 
# library(knitr)
# # Print the table
# print(kable_table)


#this is for downstream analysis
metadata <- metadata %>%
  mutate(
    AGE = scale(AGE)[,1],
    PREOP_ABX = ifelse(PREOP_ABX=='CEFAZOLIN',1,0),
    BACTRIM = ifelse(BACTRIM=='BACTRIM',1,0),
    GENDER = ifelse(GENDER=='FEMALE',1,0),
    DDRT = ifelse(DDRT=='D',1,0)
  )


#Abx_noAbx+AGE+PREOP_ABX+GENDER
#ZINQL
study = as.factor(metadata$ID)
X = metadata$Abx_noAbx
Z = as.matrix(metadata %>% dplyr::select(AGE,GENDER))
ybin = ifelse(otu.tab.rff>0,1,0)

#select non-rare taxa, 0.9 as cut off
p.cut = 0.9
zero.rate = c()
for(t in colnames(otu.tab.rff)){
  zero.rate = c(zero.rate,mean(otu.tab.rff[,t]==0))
}
t.select = colnames(otu.tab.rff)[ which(zero.rate<p.cut) ]

otu.tab.rff=otu.tab.rff[,t.select]

#-------------------------------------------------------------------
#debug
#t = 'Coriobacterium'
#y=yijt
#formula=y~Abx_noAbx+AGE+GENDER+(1|ID)
#formula.logistics = NA
#C='Abx_noAbx'
#taus=taus
#seed=2024
#meta=metadata
#method='both'
#n.positive.cut=5
#-------------------------------------------------------------------

#this time, I allow intercept
model.use = c()
set.seed(2024)
## quantile test 
taus = c(0.1, 0.25, 0.5, 0.75, 0.9)
d.result = data.frame(`ZINQL_MinP`=1, `ZINQL_Cauchy`=1, `lm`=1, `lmm`=1)
model.use = c()
for(t in 1:ncol(otu.tab.rff)){
#for(t in 1:5){
#result = lapply(1:5, function(t){
  print(paste('taxa', t, colnames(otu.tab.rff)[t]))
  yijt = otu.tab.rff[, t]
  #i = which(yijt>0)
  
  result.full = suppressMessages(try(
    ZINQL_fit(y=yijt,formula=y~Abx_noAbx+AGE+GENDER+(1|ID),C='Abx_noAbx',taus=taus, seed=2024, meta=metadata, method='both')
  ))
  if ('try-error' %in% class(result.full)){
    result = rep(NA, 2)
    model.use = c(model.use,NA)
  } else{
    result = result.full$`Final_P_value`
    model.use = c(model.use,result.full$`model`)
  }
  
  
  # result = suppressMessages(try(
  #   ZINQL_fit(y=yijt,X='Abx_noAbx', Z=c('AGE','PREOP_ABX','BACTRIM','GENDER','DDRT','intercept'),random='ID',taus=taus, seed=2024, meta=metadata, method='both')
  # ))
  # if ('try-error' %in% class(result)){
  #   result = rep(NA, 2)
  # }
  
  ## linear regression and linear mixed model
  lmfit = lm(yijt ~ X + Z)
  pval_lm = summary(lmfit)$coefficients['X', 'Pr(>|t|)']
  
  lmmfit = suppressMessages(try(lmerTest::lmer(yijt ~ X + Z + (1 | study))))
  if ('try-error' %in% class(lmmfit)){
    pval_lmm = NA
  } else{
    pval_lmm = summary(lmmfit)$coefficients['X', 'Pr(>|t|)']
  }
  result_bench = c(pval_lm, pval_lmm); names(result_bench) = c('lm','lmm')
  result = c(result, result_bench)
  d.result = rbind(d.result,result)
}

d.result = d.result[2:nrow(d.result),]
d.result$model.use = model.use

rownames(d.result) = colnames(otu.tab.rff)
saveRDS(d.result,'H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\DA_7_19_age_gender\\pval_ZINQL_v11_gender_age.rds')

#Linda
set.seed(2024)
## quantile test 
taus = c(0.1, 0.25, 0.5, 0.75, 0.9)
## existing methods
result.linda = linda(
  t(otu.tab.rff), 
  metadata, 
  formula = '~ Abx_noAbx+AGE+GENDER+(1|ID)',
  feature.dat.type = 'count',
  prev.filter = 0, 
  is.winsor = TRUE, 
  outlier.pct = 0.03,
  p.adj.method = "BH", 
  alpha = 0.05)
pval_linda = result.linda[["output"]][["Abx_noAbx"]]$pvalue
pval.linda = data.frame(linda=pval_linda)
rownames(pval.linda) = colnames(otu.tab.rff)
saveRDS(pval.linda,'H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\DA_7_19_age_gender\\pval_Linda_original.rds')


#MasLin2
set.seed(2024)
res.ma = Maaslin2(
  input_data = otu.tab.rff, 
  input_metadata = metadata, 
  output = 'H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\outtmp', 
  #transform = "AST",
  fixed_effects = c('Abx_noAbx','AGE','GENDER'),
  random_effects = c('ID'),
  plot_heatmap = F,
  plot_scatter = F,
  standardize = FALSE)

summary.tab = res.ma$results %>% filter(name == 'Abx_noAbx') 
summary.tab = summary.tab[match(colnames(otu.tab.rff), summary.tab$feature), ]
pval.Mas = data.frame(Maslin = summary.tab %>% pull(pval))
rownames(pval.Mas) = colnames(otu.tab.rff)
saveRDS(pval.Mas,'H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\DA_7_19_age_gender\\pval_MasLin2_original.rds')

#zinb different versions
N = rowSums(otu.tab.rff)
metadata$N = N + 1
#manually create a numeric ID2
# ID.map = 1:length(unique(metadata$ID))
# names(ID.map) = unique(metadata$ID)
# metadata$ID2 = ID.map[ metadata$ID ]

#nb
#!!!!!!!!!!!!!if you have already done rarify, no need to add log(n) since they are all the same
set.seed(2024)
fit.nb = mms(y = otu.tab.rff, 
          #fixed = ~ Abx_noAbx + AGE + RACE + Post_Op_Day + DM + PREOP_ABX + BACTRIM + stats::offset(log(N)),
          fixed = ~ Abx_noAbx + AGE + GENDER,
          #fixed = ~ Abx + age + stats::offset(log(N)),
          random = ~ 1 | ID, 
          #random = ~ 1 | study, 
          data = metadata,
          min.p = 0, 
          method = 'nb')


res = NBZIMM::fixed(fit.nb)$dist
res = res[grepl('Abx_noAbx', rownames(res)), ]
## handle missing taxa
taxa = res$responses
pval.nb = rep(NA, ncol(otu.tab.rff))
pval.nb[match(taxa, colnames(otu.tab.rff))] = res$pvalue



set.seed(2024)
fit.zinb = mms(y = otu.tab.rff, 
               fixed = ~ Abx_noAbx + AGE + GENDER,
            random = ~ 1 | ID, 
            data = metadata,
            min.p = 0, 
            zi_fixed = ~ 1, 
            zi_random = NULL,
            method = 'zinb')

res = NBZIMM::fixed(fit.zinb)$dist
res = res[grepl('Abx_noAbx', rownames(res)), ]
## handle missing taxa
taxa = res$responses
pval.zinb = rep(NA, ncol(otu.tab.rff))
pval.zinb[match(taxa, colnames(otu.tab.rff))] = res$pvalue

set.seed(2024)
otu.tab.log = log2(otu.tab.rff + 1)
fit.zig_count = mms(y = otu.tab.log, 
                    fixed = ~ Abx_noAbx + AGE + GENDER,
                    random = ~ 1 | ID, 
          data = metadata,
          min.p = 0, 
          zi_fixed = ~ 1,
          zi_random = NULL,
          method = 'zig')


res = NBZIMM::fixed(fit.zig_count)$dist
res = res[grepl('Abx_noAbx', rownames(res)), ]
## handle missing taxa
taxa = res$responses
pval.zig_count = rep(NA, ncol(otu.tab.rff))
pval.zig_count[match(taxa, colnames(otu.tab.rff))] = res$pvalue


pval.nbzimm = data.frame(nb=pval.nb,zinb=pval.zinb,zig_count=pval.zig_count)
rownames(pval.nbzimm) = colnames(otu.tab.rff)
saveRDS(pval.nbzimm,'H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\DA_7_19_age_gender\\pval_nbzimm_original.rds')

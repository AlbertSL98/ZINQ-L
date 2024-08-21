library(ggplot2)
library(tidyverse)
library(phyloseq)
library(GUniFrac)

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

set.seed(2024)
row.sum = rowSums(otu.tab)
otu.tab = otu.tab[which(row.sum > 10000),]

#get genus name for each column
otu.col = colnames(otu.tab)
otu.row = rownames(otu.tab)

otu.tab.rff = GUniFrac::Rarefy(otu.tab)$otu.tab.rff
colnames(otu.tab.rff) = taxa.map[colnames(otu.tab.rff)]
row.names(otu.tab.rff) = otu.row
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
#----------------------------------------------------------------


#take max post.op day
meta.max.post.op.day = metadata %>%
  select(ID,Post_Op_Day) %>%
  group_by(ID) %>%
  summarise(Post_Op_Day = max(Post_Op_Day)) %>%
  inner_join(metadata,by=c('ID','Post_Op_Day'))

#some samples are droped in the first rowsum step, so also drop them in meta data
ID.keep = c()
for(i in meta.max.post.op.day$sample_ID){ID.keep = c(ID.keep,ifelse(i %in% rownames(otu.tab.rff),1,0))}

meta.max.post.op.day$keep = ID.keep
meta.max.post.op.day = meta.max.post.op.day %>%
  filter(keep==1)



otu.tab.rff2 = otu.tab.rff[meta.max.post.op.day$sample_ID,]
#
table(meta.max.post.op.day$Abx_noAbx)
#
index1 = which(meta.max.post.op.day$Abx_noAbx==1)
index2 = which(!(meta.max.post.op.day$Abx_noAbx==1))

#Zero rate
zero.prop = c()
na.prop = c()
for(i in 1:ncol(otu.tab.rff2)){
  zero.prop = c( zero.prop, mean(otu.tab.rff2[,i]==0) )
  na.prop = c( na.prop, mean(is.na(otu.tab.rff2[,i])) )
}
#what if we drop taxas >=0.9 zero rate
cut = 0.9
index = which(zero.prop<cut)
otu.tab.rff3 = otu.tab.rff2[,index]
zero.prop2 = zero.prop[index]
taxa.tmp = colnames(otu.tab.rff3)

hist(zero.prop2)
median(zero.prop2)

#cut off is 0.785, median()
saveRDS(otu.tab.rff3, "H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\ktx_otu_rff.rds")

#group taxa by zero rate
#here 1 means rare
dense = ifelse(zero.prop2>=median(zero.prop2),1,0)
taxa.group = data.frame(taxa=taxa.tmp,dense )
saveRDS(taxa.group, "H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\ktx_dense_group.rds")

saveRDS(meta.max.post.op.day, "H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\ktx_meta.rds")

#variable of interest in meta
#Abx_noAbx, DDRT, PREOP_ABX, BACTRIM
#?shall we add post op day

#some interest genus: enterococcus, lactobacillus, eubacterium, anaerostipes, ruminococcus, streptococcus

# taxa.list = colnames(otu.tab.rff2)
# 
# taxa.list =  c('Enterococcus', 'Lactobacillus', 'Eubacterium', 'Anaerostipes', 'Ruminococcus', 'Streptococcus') 
# 
# taxa.list =  c('Blautia', 'Dorea', 'Enterococcus') 
# 

mean(otu.tab.rff2[,'Blautia']==0)
mean(otu.tab.rff2[,'Enterococcus']==0)
mean(otu.tab.rff2[,'Anaerofustis']==0)
mean(otu.tab.rff2[,'Dorea']==0)

meta.max.post.op.day =  readRDS("H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\ktx_meta.rds")
otu.tab.rff2 = readRDS("H:\\Runzhe\\ZINQL\\Ktx\\Simulation3\\mysim\\data\\ktx_otu_rff.rds")
taxa.list = colnames(otu.tab.rff2)

for(taxa in taxa.list){
  print(taxa)
  #
  x = seq(0,1,by=0.02)
  y1 = quantile( otu.tab.rff2[index1,taxa], x)
  y2 = quantile( otu.tab.rff2[index2,taxa], x)
  #
  d.plot = data.frame (x=rep(x,2),y=c(y1,y2), type=c( rep(0,length(x)),rep(1,length(x)) ) )
  #
  outdir = paste( 'H:\\Runzhe\\ZINQL\\Ktx\\Plot', paste('Max_post_op_taxacol_',taxa,'.tiff',sep=''), sep='\\')
  tiff(outdir)
  g = ggplot() +
    geom_line(data=d.plot,aes(x=x,y=y,color=as.factor(type),lty=as.factor(type)),lwd=1) + ggtitle( paste0(taxa,' max post op') ) +  theme(plot.title = element_text(hjust=0.5, size=20) )
  print(g)
  dev.off()
}



custom_theme <- function() {
  list(
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 25, face = "italic"),
      axis.text.x = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 25),
      axis.text.y = element_text(size = 20),
      axis.title.y = element_text(size = 30),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")  # Add axis lines
    ) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    #for first plot, hide legend
    theme( legend.position = "none")
      ,
    scale_color_manual(
      values = c("1" = "#12685D","0" = "#FCBB44"),
      labels = c("1" = "Abx", "0" = "No Abx")
    ) 
    ,
    scale_linetype_manual(
      values = c("1" = "solid", "0" = "dotdash"),
      labels = c("1" = "Abx", "0" = "No Abx")
    ) 
    ,
    labs(color = "Abx Group",linetype='Abx Group')
  )
}

index1 = meta.max.post.op.day %>% filter(Abx_noAbx==0) %>% pull(sample_ID)
index2 = meta.max.post.op.day %>% filter(Abx_noAbx==1) %>% pull(sample_ID)
x = seq(0,1,by=0.01)

#
taxa = 'Blautia'
print(taxa)
#
#y1 is none abx group
y1 = quantile( otu.tab.rff2[index1,taxa], x)
y2 = quantile( otu.tab.rff2[index2,taxa], x)
#
d.plot = data.frame (x=rep(x,2),y=c(y1,y2), type=c( rep(0,length(x)),rep(1,length(x)) ) )
#
g.blautia = ggplot() +
  geom_line(data=d.plot,aes(x=x,y=y,color=as.factor(type),lty=as.factor(type)),lwd=1.5) + 
  ggtitle( paste0(taxa) ) + 
  custom_theme()
  
#
taxa = 'Dorea'
print(taxa)
#
#y1 is none abx group
y1 = quantile( otu.tab.rff2[index1,taxa], x)
y2 = quantile( otu.tab.rff2[index2,taxa], x)
#
d.plot = data.frame (x=rep(x,2),y=c(y1,y2), type=c( rep(0,length(x)),rep(1,length(x)) ) )
#
g.dorea = ggplot() +
  geom_line(data=d.plot,aes(x=x,y=y,color=as.factor(type),lty=as.factor(type)),lwd=1.5) + 
  ggtitle( paste0(taxa) ) + 
  custom_theme()

#
taxa = 'Enterococcus'
print(taxa)
#
#y1 is none abx group
y1 = quantile( otu.tab.rff2[index1,taxa], x)
y2 = quantile( otu.tab.rff2[index2,taxa], x)
#
d.plot = data.frame (x=rep(x,2),y=c(y1,y2), type=c( rep(0,length(x)),rep(1,length(x)) ) )
#
g.enterococcus = ggplot() +
  geom_line(data=d.plot,aes(x=x,y=y,color=as.factor(type),lty=as.factor(type)),lwd=1.5) + 
  ggtitle( paste0(taxa) ) + 
  custom_theme()

#
taxa = 'Anaerofustis'
print(taxa)
#
#y1 is none abx group
y1 = quantile( otu.tab.rff2[index1,taxa], x)
y2 = quantile( otu.tab.rff2[index2,taxa], x)
#
d.plot = data.frame (x=rep(x,2),y=c(y1,y2), type=c( rep(0,length(x)),rep(1,length(x)) ) )
#
g.anaerofustis = ggplot() +
  geom_line(data=d.plot,aes(x=x,y=y,color=as.factor(type),lty=as.factor(type)),lwd=1.5) + 
  ggtitle( paste0(taxa) ) + 
  custom_theme()

library(grid)
library(gridExtra)
#grid.arrange(p1, p2,  ncol = 2,  left = plot.y , bottom = plot.x)
plot.y = textGrob('Empirical Quantile', rot = 90, gp = gpar(fontsize = 30))
plot.x = textGrob(expression('Quantile level'), rot = 0, gp = gpar(fontsize = 30))
#title.tmp = paste('Simulation 1, ',taxa.a,', ', sim.type,sep='' )
#plot.title = textGrob('Quantile Distribution', rot = 0, gp = gpar(fontsize = 25))
plot.title = ''
#
tiff('C:\\Users\\lenovo\\Desktop\\a.tiff',width=2500,height=2500,res=200)
grid.arrange(g.blautia,g.dorea,g.enterococcus,g.anaerofustis,  ncol = 2, top=plot.title, left = plot.y , bottom = plot.x)
dev.off()

#grid.arrange(g.blautia,g.dorea,g.enterococcus,g.anaerofustis,  ncol = 2, left = plot.y , bottom = plot.x)

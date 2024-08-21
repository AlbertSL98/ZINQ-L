library(UpSetR)
library(tidyverse)

pval.ZINQL = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\DA_7_19_age_gender\\pval_ZINQL_v11_gender_age.rds')

# pval.ZINQL = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\pval_ZINQL_with_int_v7.rds')



pval.linda = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\DA_7_19_age_gender\\pval_Linda_original.rds')


pval.Mas = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\DA_7_19_age_gender\\pval_MasLin2_original.rds')
#let us handle missing taxa


pval.nbzimm = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\DA_7_19_age_gender\\pval_nbzimm_original.rds')

#only consider those that use MinP and Cauchy in ZINQL
#fill in NA for those tests that only use logistics
pval.ZINQL = pval.ZINQL %>%
  mutate(ZINQL_MinP=ifelse(model.use=='Both',ZINQL_MinP,NA)) %>%
  mutate(ZINQL_Cauchy=ifelse(model.use=='Both',ZINQL_Cauchy,NA))

pval = cbind(pval.ZINQL,pval.linda,pval.Mas,pval.nbzimm)
p.val.fdr = pval
for(i in 1:ncol(p.val.fdr)){
  p.val.fdr[,i] = p.adjust(p.val.fdr[,i], method = 'fdr')
}

taxa.feature = list()
for(i in colnames(p.val.fdr)){
  d.tmp = p.val.fdr[,i]
  names(d.tmp) = rownames(p.val.fdr)
  taxa.feature[[i]] = names( which(d.tmp<=0.05) )
}


n = c()
for(i in colnames(p.val.fdr)){
  n = c(n,length(taxa.feature[[i]]))
}

d.plot = data.frame(n=n,method=colnames(p.val.fdr))

ggplot( d.plot, ) +
  geom_bar( aes(x = method, y = n), alpha=0.7, stat = "identity", position = "dodge") +
  #labs(title = "Simulation 1, Dorea, Null Hypothesis, sample size=200", x = "Number of repeats", y = "type 1-error") + 
  labs(title = "Number of taxa identified") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.text.x = element_text(angle = 45, hjust = 0.5, size = 15),
    #axis.title.x = element_text(size = 25),
    axis.text.y = element_text(size = 20),
    #axis.title.y = element_text(size = 30),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")  # Add axis lines
  )  +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) + 
  geom_text(aes(x = method, y = n, label = n), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 5)


#if we choose cauchy_minp as the one
#select taxas in this group but not others
g1 = c('ZINQL_MinP', 'ZINQL_Cauchy')
#g1 = c('dep:cauchy_trunc_zero')
g2 = c("lmm","linda", "Maslin", "zig_count")
#g2 = c("lmm","linda", "Maslin", "zig_count",'zinb')
#g2 = c("lmm","linda", "Maslin", "zig_count",'nb','zinb')

#take intersection of methods in group1
taxa1 = taxa.feature[[g1[1]]]
if(length(g1)>1){
  for(i in g1[2:length(g1)]){
    taxa1 = intersect(taxa1, taxa.feature[[i]])
  }
}
taxa2 = c()
for(i in g2){
  taxa2 = c(taxa2,taxa.feature[[i]])
}
#
taxa1 = unique(taxa1)
taxa2 = unique(taxa2)

#
taxa.select = setdiff(taxa1,taxa2)

#let us plot their quantile
d = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\otu_tab_rff.rds')
meta = readRDS('H:\\Runzhe\\ZINQL\\Ktx\\Data_ana\\meta.rds')

#plot their quantile
g.list = list()
#t='Roseburia'
for(t in taxa.select){
  d.plot = d[,t]
  y.abx = d[meta%>%filter(Abx_noAbx==1)%>%pull(sample_ID),t]
  y.noabx = d[meta%>%filter(Abx_noAbx==0)%>%pull(sample_ID),t]
  x.quantile = seq(0,1,by=0.01)
  d.plot = rbind( data.frame(group=rep('Abx',length(x.quantile)),y=quantile(y.abx,x.quantile),x=x.quantile),data.frame(group=rep('Not Abx',length(x.quantile)),y=quantile(y.noabx,x.quantile),x=x.quantile) )
  
  # g.list[[t]] = ggplot() +
  #   geom_line(data=d.plot,aes(x=x,y=y,color=as.factor(group),lty=as.factor(group)),lwd=1.5) +
  #   scale_color_manual(values = c("Abx" = "#12685D", "Not Abx" = "#E64B35") ) +
  #   labs(color = "Group", lty = "Group",
  #        x = 'Quantile',y='Count') +
  #   ggtitle(t) +
  #   theme_minimal() +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = 40),
  #     axis.text.x = element_text(hjust = 0.5, size = 10),
  #     axis.title.x = element_text(size = 10),
  #     axis.text.y = element_text(size = 10),
  #     axis.title.y = element_text(size = 10),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     axis.line = element_line(color = "black")  # Add axis lines
  # ) +
  #   theme(legend.position = 'None')
  
  g.list[[t]] = ggplot() +
    geom_line(data=d.plot,aes(x=x,y=y,color=as.factor(group),lty=as.factor(group)),lwd=1.5) +
    scale_color_manual(values = c("Abx" = "#12685D", "Not Abx" = "#E64B35") ) +
    labs(color = "Group", lty = "Group",
         x = 'Quantile',y='Count') +
    ggtitle(t) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 50, face = "italic"),
      axis.text.x = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 0),
      axis.text.y = element_text(size = 20),
      axis.title.y = element_text(size = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")  # Add axis lines
    ) +
    theme(legend.position = 'None')
}

library(grid)
library(gridExtra)
plot.y = textGrob(expression('Empirical Quantile'), rot = 90, gp = gpar(fontsize = 50))
plot.x = textGrob(expression('Quantile level'), rot = 0, gp = gpar(fontsize = 50))
#plot.title = textGrob(expression('Distribution of Taxa identified by ZINQ-L'), rot = 0, gp = gpar(fontsize = 55))
plot.title = ''
tiff('C:\\Users\\lenovo\\Desktop\\taxa.tiff',res=200,width=8000,height=6000)
grid.arrange(grobs=g.list,  ncol = 3, top=plot.title, left = plot.y , bottom = plot.x)
dev.off()



#heatmap show overlap of them
methods = c('ZINQL_MinP',"lmm","linda", "Maslin", "zig_count")
methods = c('ZINQL_MinP',"lmm","linda", "Maslin", "zig_count",'nb','zinb')
intersect.matrix = matrix(rep(0,length(methods)*length(methods)), ncol=length(methods))
colnames(intersect.matrix) = methods
rownames(intersect.matrix) = methods
#
#?? divide n of intersect by n of union
for(i in methods){
  for(j in methods){
    intersect.matrix[i,j] = length( intersect(taxa.feature[[i]],taxa.feature[[j]]) ) / length( union(taxa.feature[[i]],taxa.feature[[j]]) )
  }
}

for(i in 1:ncol(intersect.matrix)){intersect.matrix[i,i]=NA}

library(pheatmap)
h = pheatmap(intersect.matrix,cluster_rows = F, cluster_cols = F,main='Overlap of taxas identified',na_col='white',treeheight_row=0,treeheight_col=0, fontsize = 12, show_rownames = T, show_colnames = T, cellheight = 30, cellwidth = 30, annotation_names_row = F, annotation_names_col = F)


#methods = c('dep:minP',"lmm","linda", "Maslin", "zig_count",'nb','zinb')
#methods = c("lmm","linda", "Maslin", "zig_count","dep:cauchy_trunc_zero")

#methods = c("lmm","linda", "Maslin", "zig_count","zinb","dep:cauchy_trunc_zero", "dep:minP")

taxa.feature[['ZIGMM']] = taxa.feature[['zig_count']]
taxa.feature[['ZINQ-L Cauchy']] = taxa.feature[['ZINQL_Cauchy']]
taxa.feature[['ZINQ-L MinP']] = taxa.feature[['ZINQL_MinP']]
taxa.feature[['LinDA']] = taxa.feature[['linda']]
taxa.feature[['MaAsLin2']] = taxa.feature[['Maslin']]
methods = c("lmm","LinDA", "MaAsLin2", "ZIGMM","ZINQ-L Cauchy", "ZINQ-L MinP")

me.select = taxa.feature[methods]

#https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
upset( fromList(me.select), 
       order.by = "freq", 
       number.angles = 0, 
       point.size = 3.5, 
       line.size = 2, 
       mainbar.y.label = "Overlap among methods", 
       sets.x.label = "Number of Taxa selected", 
       text.scale = c(2.5, 1.5, 1.5, 2, 2, 2.5), 
       matrix.color = "#70A5D9", 
       main.bar.color = "#BC6356",
       sets.bar.color = "#4F845C", 
       sets=methods)

#text.scale: 1. mainbar.y.label, 2. y-axis (Here is overlap among methods), 3. sets.x.label, 4. x-axis label (here is the number of taxas selected), 5. x-axis labels (here is zinb, linda and so on) 6. number on top of bar


f.minP = taxa.feature[["dep:minP"]]
f.cauchy = taxa.feature[["dep:cauchy_trunc_zero"]]
for(f in f.minP){
  if(!(f %in% f.cauchy)){
    print(f)
  }
}

mean(is.na(pval.ZINQL[,'dep:minP']))

















total_count <- nrow(metadata)
cohort_char = metadata %>%
  group_by(Abx_Category) %>%
  summarise(
    count = n(),
    age = median(AGE),
    gender = paste(sum(GENDER == 'FEMALE'), '(', round(mean(GENDER == 'FEMALE') * 100, 1), '%)', sep = ''),
    race = paste(sum(RACE == 'AA'), '(', round(mean(RACE == 'AA') * 100, 1), '%)', sep = ''),
    DM = paste(sum(DM == 'DM'), '(', round(mean(DM == 'DM') * 100, 1), '%)', sep = ''),
    PREOP_ABX = paste(sum(PREOP_ABX == 'CEFAZOLIN'), '(', round(mean(PREOP_ABX == 'CEFAZOLIN') * 100, 1), '%)', sep = ''),
    DDRT = paste(sum(DDRT == 'D'), '(', round(mean(DDRT == 'D') * 100, 1), '%)', sep = ''),
    PRIORTX = paste(sum(PRIORTX == '1'), '(', round(mean(PRIORTX == '1') * 100, 1), '%)', sep = ''),
    DGF = paste(sum(DGF == 'DGF'), '(', round(mean(DGF == 'DGF') * 100, 1), '%)', sep = ''),
    THYMO = paste(sum(THYMO == 'THYMO'), '(', round(mean(THYMO == 'THYMO') * 100, 1), '%)', sep = ''),
    PREDNISONE = paste(sum(PREDNISONE == 'PRED'), '(', round(mean(PREDNISONE == 'PRED') * 100, 1), '%)', sep = ''),
    BACTRIM = paste(sum(BACTRIM == 'BACTRIM'), '(', round(mean(BACTRIM == 'BACTRIM') * 100, 1), '%)', sep = '')
  )


cohort_char = t(cohort_char) %>% as.data.frame()
# Remove the first row about Abx_category
cohort_char <- cohort_char[-1, ]

# Add column names
colnames(cohort_char) <- c("No Abx", "Other Abx", "Beta-lactam", "FQ", "Beta-lactam and FQ")

rownames(cohort_char) <- c("Count", "Age", "Female Gender", "AA Race", "Diabetes Mellitus (DM)", "Preoperative Antibiotics (PREOP_ABX)", "Bactrim Usage (BACTRIM)", "Deceased Donor Transplantation (DDRT)", "Prior Transplant (PRIORTX)", "Delayed Graft Function (DGF)", "Anti-thymocyte Globulin (THYMO)", "Prednisone Steroid Maintenance Protocol (PREDNISONE)")



# Create a kable
kable_table <- kable(cohort_char, 
                     caption = "Characteristics by 5 Antibiotics Group in Sample Level")

print(kable_table)


#figure out what happens for NA ZINQ-L
#this is n of non-zero taxa
yyy = rownames(pval.ZINQL)[ which(is.na(pval.ZINQL[,'dep:minP'])) ]
for(y in yyy){
  yijt = otu.tab.rff[, y]
  i = which(yijt>0)
  print( paste( y,length(i) ) )
}

t = 'Hydrogenoanaerobacterium'
t = 'Candidatus Phytoplasma'
t='Alloprevotella'
yijt = otu.tab.rff[, t]
i = which(yijt>0)
metadata$GENDER[i]
metadata$AGE[i]

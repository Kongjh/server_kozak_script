# enviroment --------------------------------------------------------------
setwd("~/rstudio/kozak/lucas/")
getwd()
library(psych)
library(corrplot)#载入两个包
library("grid")
library(Hmisc)
library(PerformanceAnalytics)#加载包
library("tidyverse")
# help(package="psych")
# detach("package:tidyverse")
options(digits=5)
options(tibble.width = Inf) # 表示 tibble 总是打印所有列 ，比如 使用 head 等函数的时候
#last_plot()
p <- ggplot() + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), 
        legend.position = "right",legend.text =  element_text(face="bold", size=23),
        legend.title = element_blank(),
        # axis.text.x = element_text(size = 26,face = "bold"),
        # axis.ticks = element_line(),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(.20,"cm"),
        axis.line = element_line(size = 1.5),
        axis.text = element_text(size = 30,face = "bold"),  
        axis.title = element_text(size = 33, face = "bold"))
pp <- ggplot() + theme_classic()
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
a


# <1>  --------------------------------------------------------------------
# https://github.com/markziemann/dee2/blob/master/AccessDEEfromR.md
# library("devtools")
# devtools::install_github("markziemann/dee2/getDEE2")
library("getDEE2")


# <1.1> 不用再跑 --------------------------------------------------------------
mdat <- getDee2Metadata("scerevisiae") # 11369
unique(factor(mdat$QC_summary))
write_csv(mdat,"~/rstudio/kozak/lucas/mdat.csv")
pmdat <- filter(mdat,QC_summary == "PASS")# 3025
write_csv(pmdat,"~/rstudio/kozak/lucas/pmdat.csv")

# <2> ---------------------------------------------------------------------
# DEE2 data is centred around SRA run accessions numbers, these SRR_accessions can be obtained like this:
# > mdat1<-mdat[which(mdat$GSE_accession %in% "GSE33569"),]
# > SRRlist<-as.vector(mdat1$SRR_accession)
# > SRRlist
# [1] "SRR363796" "SRR363797" "SRR363798" "SRR363799"
head(pmdat)
SRRlist <-as.vector(pmdat$SRR_accession)
write_csv(as.data.frame(SRRlist),"~/rstudio/kozak/lucas/pass_srrlist.csv")
x <- getDEE2("scerevisiae",c("SRR7242111"))
x1 <- x$GeneCounts
x2 <- x$TxCounts
x3 <- x$TxInfo
sum(x1$SRR7242111)


# <3> 大数据准备---------------------------------------------------------------------
gexon <- read_csv("~/rstudio/kozak/lucas/gene_exonlength2.csv") %>% select(-"exon_num")
head(gexon)
gscore <- read_csv("~/rstudio/kozak/gene_21nt_aug_score.csv")
head(gscore)
pmdat <- read_csv("~/rstudio/kozak/lucas/pmdat.csv") # passed 的data
head(pmdat)

srrlist <-read_csv("~/rstudio/kozak/lucas/pass_srrlist.csv")
dat <- read_tsv("~/rstudio/kozak/lucas/scerevisiae_se.tsv",col_names = F)
unique(factor(dat$X1)) # 11369
head(dat)
dat2 <- filter(dat, X1 %in% srrlist$SRRlist)
# head(dat2)
# ?spread
sdat <- spread(dat2, key="X1", value="X3", fill = NA, convert = FALSE, drop = FALSE) 
sdat[,1:10]
write_csv(sdat,"~/rstudio/kozak/lucas/sdat.csv")

# <4> 换成循环的思路---------------------------------------------------------------------
sdat <- read_csv("~/rstudio/kozak/lucas/sdat.csv") %>% filter(grepl('^Y',X2)) %>% rename(gene_id = X2)
test <- sdat[,1:10]
head(test)
test <- sdat
num <- seq(1,3025)
rho <- vector(mode="numeric",length=0)
pvalue <- vector(mode="numeric",length=0)
# gene_num <- vector(mode="numeric",length=0)
for (i in num){
  t1 <- select(test,gene_id,i+1)
  s1 <- sum(t1[,2])/1000000
  t1[,2] <- t1[,2]/s1
  # head(t1)
  # head(gexon)
  t2 <- merge(gexon,t1,by="gene_id")
  # head(t2)
  rpkm <- t3[,3]/t3[,2]
  t3 <- cbind(t2,rpkm)
  # head(t3)
  # head(gscore)
  t4 <- merge(t3,gscore,by = "gene_id")
  res <- cor.test(t4$rpkm,t4$score_21nt,method = "s",exact = F)
  res$p.value -> pvalue[i]
  res$estimate -> rho[i]
}
# rho
cname <- colnames(test)
fdat <- cbind(cname[-1],rho,pvalue) %>% as_tibble() %>% rename(SRR_accession = V1)
# write_csv(fdat,"~/rstudio/kozak/lucas/fdat.csv") # 指 final 最后可用的数据
head(fdat)
head(pmdat)

val_dat <- merge(fdat,pmdat,by="SRR_accession")
# write_csv(val_dat,"~/rstudio/kozak/lucas/val_dat.csv")


# 数据分析 --------------------------------------------------------------------
fdat <- read_csv("~/rstudio/kozak/lucas/fdat.csv")

# fdat$rho <- as.numeric(fdat$rho)
num <- seq(1,nrow(fdat))
(p5 <- p + aes(x = num,y = fdat$rho) 
  + geom_point(size = .5)
  + labs(x= "SRR number",y="rho (mRNA level 与 initiation score"))
plot(density(fdat$rho))

# check test（过去式） ---------------------------------------------------------
head(t3)
102.776/3.483
# head(sdat)
b <- sdat[,1:10]
sdat2 <- apply(sdat[,-1],2,function(x){as.numeric(x*1000000/sum(x))})
sdat3 <- cbind(sdat$X2,sdat2)
b3 <- sdat3[,1:10]
head(b3)
# #check
# s <- sum(sdat$ERR1095152)
# as.data.frame(sdat$ERR1095152/s) ->ss

fsdat <- merge(gexon,sdat3,by.x = "gene_id", by.y="V1",all.x = T) 
f <- fsdat[,1:20] %>% as.numeric()
head(f)
fsdat2 <- apply(fsdat[,-1],1,function(x){as.numeric(x)/as.numeric(x[1])})
fsdat3 <- cbind(fsdat$gene_id,fsdat2)
f3 <- fsdat3[,1:20]
head(f3)


# <5> 用 mRNA looped 那个 mutants 来算 score 与 mRNA level---------------------------------------------------------------------
# 加载 data 运行一次 save就行 -----------------------------------------------------
# 注意 这里 几个 factor 没用3个重复的 mrna 比如 Hh 也有一些只有一个样本 比如 F1这类
o4A <- read_csv("~/rstudio/mrna_looped/Factors/4A/mRNA_level.csv") %>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o4A = mRNA_level) 
CAF20d <- read_csv("~/rstudio/mrna_looped/Factors/CAF20d/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(CAF20d = mRNA_level) 
o5Bd <- read_csv("~/rstudio/mrna_looped/Factors/5Bd/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o5Bd = mRNA_level) 
o4G1d <- read_csv("~/rstudio/mrna_looped/Factors/4G1d/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o4G1d = mRNA_level) 
o4G1dchx <- read_csv("~/rstudio/mrna_looped/Factors/4G1dchx/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o4G1dchx = mRNA_level) 
PAB1D <- read_csv("~/rstudio/mrna_looped/Factors/PAB1D/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(PAB1D = mRNA_level) 
PAB1Dchx <- read_csv("~/rstudio/mrna_looped/Factors/PAB1Dchx/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(PAB1Dchx = mRNA_level) 
PAB1d <- read_csv("~/rstudio/mrna_looped/Factors/PAB1d2/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(PAB1d = mRNA_level) 
# PAB1d$PAB1d <- as.numeric(PAB1d$PAB1d)
# PAB1d <- filter(PAB1d, PAB1d > 0 )
o4ED <- read_csv("~/rstudio/mrna_looped/Factors/4ED/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(o4ED = mRNA_level) 
AD <- read_csv("~/rstudio/mrna_looped/Factors/AD/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(AD = mRNA_level) 
HD <- read_csv("~/rstudio/mrna_looped/Factors/HD/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(HD = mRNA_level) 
Hd <- read_csv("~/rstudio/mrna_looped/Factors/Hd2/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(Hd = mRNA_level) 
ho <- read_csv("~/rstudio/mrna_looped/Factors/ho/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(ho = mRNA_level) 
wt <- read_csv("~/rstudio/mrna_looped/Factors/WT/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(wt = mRNA_level) 
Hh <- read_csv("~/rstudio/mrna_looped/Factors/Hh/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_3)/2)  %>% 
  select(gene_id,mRNA_level) %>% rename(Hh = mRNA_level) 
H1 <- read_csv("~/rstudio/mrna_looped/Factors/H1/mRNA_level.csv") %>% select(gene_id,mRNA_level_1) %>% rename(H1 = mRNA_level_1) 
H1s <- read_csv("~/rstudio/mrna_looped/Factors/H1s/mRNA_level.csv")%>% select(gene_id,mRNA_level_1) %>% rename(H1s = mRNA_level_1) 
F1 <- read_csv("~/rstudio/mrna_looped/Factors/F1/mRNA_level.csv") %>% select(gene_id,mRNA_level_1) %>% rename(F1 = mRNA_level_1) 
F1s <- read_csv("~/rstudio/mrna_looped/Factors/F1s/mRNA_level.csv")%>% select(gene_id,mRNA_level_1) %>% rename(F1s = mRNA_level_1) 
mfs <- read_csv("~/rstudio/mrna_looped/Factors/mfs/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(mfs = mRNA_level) 
gal1 <- read_csv("~/rstudio/mrna_looped/Factors/gal1/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(gal1 = mRNA_level) 
m4g1 <- read_csv("~/rstudio/mrna_looped/Factors/m4g1/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(m4g1 = mRNA_level) 
mho <- read_csv("~/rstudio/mrna_looped/Factors/mho/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(mho = mRNA_level) 
KLRK <- read_csv("~/rstudio/mrna_looped/Factors/KLRK/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(KLRK = mRNA_level) 
OE <- read_csv("~/rstudio/mrna_looped/Factors/OE/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(OE = mRNA_level) 
FS <- read_csv("~/rstudio/mrna_looped/Factors/FS/mRNA_level.csv")%>% mutate(mRNA_level = (mRNA_level_1+mRNA_level_2+mRNA_level_3)/3)  %>% 
  select(gene_id,mRNA_level) %>% rename(FS = mRNA_level) 

head(CAF20d)
dfac_mrna <- merge(o4A,CAF20d,by=c("gene_id")) %>% merge(o5Bd,by=c("gene_id")) %>% merge(o4G1d,by=c("gene_id")) %>%
  merge(o4G1dchx,by=c("gene_id"))%>% merge(PAB1D,by=c("gene_id"))%>% merge(PAB1Dchx,by=c("gene_id"))%>% 
  merge(PAB1d,by=c("gene_id"))%>% merge(o4ED,by=c("gene_id"))%>% merge(AD,by=c("gene_id"))%>% 
  merge(HD,by=c("gene_id"))%>% merge(Hd,by=c("gene_id"))%>% merge(ho,by=c("gene_id"))%>% 
  merge(wt,by=c("gene_id"))%>% merge(Hh,by=c("gene_id"))%>% merge(H1,by=c("gene_id"))%>% 
  merge(H1s,by=c("gene_id"))%>% merge(F1,by=c("gene_id"))%>% merge(F1s,by=c("gene_id"))%>% 
  merge(mfs,by=c("gene_id"))%>% merge(gal1,by=c("gene_id"))%>% merge(m4g1,by=c("gene_id"))%>% 
  merge(mho,by=c("gene_id"))%>% merge(KLRK,by=c("gene_id"))%>% merge(OE,by=c("gene_id"))%>% 
  merge(FS,by=c("gene_id"))
head(dfac_mrna)
write_csv(dfac_mrna,"~/rstudio/mrna_looped/Factors/dfac_mrna.csv")


# 开始分析 --------------------------------------------------------------------
gexon <- read_csv("~/rstudio/gene_exonlength.csv") %>% select(-"exon_num")
head(gexon)
gscore <- read_csv("~/rstudio/kozak/gene_21nt_aug_score.csv")
head(gscore)
dfac_mrna <- read_csv("~/rstudio/mrna_looped/Factors/dfac_mrna.csv",col_names = T)
head(dfac_mrna)
test <- dfac_mrna
num <- seq(1,26)
rho <- vector(mode="numeric",length=0)
pvalue <- vector(mode="numeric",length=0)
# gene_num <- vector(mode="numeric",length=0)
for (i in num){
  t1 <- select(test,gene_id,i+1)
  t2 <- merge(gexon,t1,by="gene_id")
  # head(t2)
  # head(gscore)
  t3 <- merge(t2,gscore,by = "gene_id")
  res <- cor.test(t3[,3],t3$score_21nt,method = "s",exact = F)
  res$p.value -> pvalue[i]
  res$estimate -> rho[i]
}
# rho
cname <- colnames(test)
fdat2 <- cbind(cname[-1],rho,pvalue) %>% as_tibble() #%>% rename(SRR_accession = V1)
write_csv(fdat2,"~/rstudio/kozak/lucas/fdat_factors_from_mrnaloop.csv") # 指 final 最后可用的数据
head(fdat2)


# 数据分析 --------------------------------------------------------------------
fdat2 <- read_csv("~/rstudio/kozak/lucas/fdat_factors_from_mrnaloop.csv")

# fdat$rho <- as.numeric(fdat$rho)
num <- seq(1,nrow(fdat2))
(p5 <- p + aes(x = num,y = fdat2$rho) 
  + geom_point(size = .5)
  + labs(x= "factors number",y="rho (mRNA level 与 initiation score"))
plot(density(fdat2$rho))

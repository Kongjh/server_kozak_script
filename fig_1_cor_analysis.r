# enviroment --------------------------------------------------------------
setwd("G:/kozak_new/initiation_score/")
getwd()
library(psych)
library(corrplot)#载入两个包
library("grid")
library(Hmisc)
library(PerformanceAnalytics)#加载包
library("ggseqlogo")
library("tidyverse")
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
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
a

# arabidopsis -------------------------------------------------------------
darab <- read_csv("G:/kozak_new/initiation_score/arabidopsis/final_arabidopsis.csv",col_names = T)
head(darab);nrow(darab) # 21109
as.numeric(cut2(darab$score_21nt,g = 10)) -> darab$bin_score
tapply(log(darab$rpf), INDEX = darab$bin_score,FUN = mean) -> mean_rpf
tapply(log(darab$mrna), INDEX = darab$bin_score,FUN = mean) -> mean_mrna
tapply(log(darab$te), INDEX = darab$bin_score,FUN = mean) -> mean_te

tapply(darab$rpf, INDEX = darab$bin_score,FUN = mean) -> mean_rpf
tapply(darab$mrna, INDEX = darab$bin_score,FUN = mean) -> mean_mrna
tapply(darab$te, INDEX = darab$bin_score,FUN = mean) -> mean_te

bar1 <- tapply(log2(darab$rpf), INDEX = darab$bin_score,FUN = sd)/sqrt(
  tapply(log2(darab$mrna), INDEX = darab$bin_score,FUN = length))
bar2 <- tapply(log2(darab$mrna), INDEX = darab$bin_score,FUN = sd)/sqrt(
  tapply(log2(darab$mrna), INDEX = darab$bin_score,FUN = length))
bar3 <- tapply(log2(darab$te), INDEX = darab$bin_score,FUN = sd)/sqrt(
  tapply(log2(darab$te), INDEX = darab$bin_score,FUN = length))
bar1

(pa1 <- p + aes(x= seq(1,10), y= mean_rpf) + geom_point(size= 2.5) + 
    scale_x_discrete(limits = seq(1,10))
  + geom_errorbar(aes(ymax= mean_rpf+bar1,ymin= mean_rpf-bar1),width =.1,size =.5)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust = 2, 
             label = "rho 0.90075 \n\rp-value <2e-16" , size = 5) + labs(x = "initiation socor bins",y = "RPF")
)

cor.test(darab$score_21nt,darab$rpf,method = "s",exact = F)
cor.test(darab$score_21nt,darab$mrna,method = "s",exact = F)
cor.test(darab$score_21nt,darab$te,method = "s",exact = F)


grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(7,1))) ####将页面分成2*2矩阵
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p_s19_logo, vp = vplayout(1,1))   ###将（1,1)和(1,2)的位置画图c
print(p_s29_logo, vp = vplayout(2,1))   ###将(2,1)的位置画图b
print(p_s37_logo, vp = vplayout(3,1))  ###将（2,2)的位置画图a
# print(p_gene_logo, vp = vplayout(4,1:3))  ###将（2,2)的位置画图a
print(p_gene_logo, vp = vplayout(4,1))
# print(p_ram_logo, vp = vplayout(5,1:3))  ###将（2,2)的位置画图a
print(p_rs28_logo, vp = vplayout(5,1))  ###将（2,2)的位置画图a
print(p_ssu47_logo, vp = vplayout(6,1))
print(p_ssu55_logo, vp = vplayout(7,1))


# fly ---------------------------------------------------------------------
dfly <- read_csv("G:/kozak_new/initiation_score/fly/final_fly.csv",col_names = T)
dfly <- read_delim("G:/kozak_new/initiation_score/fly/3-score-ribo-fly.txt",delim = "\t") %>% filter(rpf >0)
head(dfly);nrow(dfly) # 5964
as.numeric(cut2(dfly$score_21nt,g = 10)) -> dfly$bin_score
tapply(log(dfly$rpf), INDEX = dfly$bin_score,FUN = mean) -> mean_rpf
tapply(log(dfly$mrna), INDEX = dfly$bin_score,FUN = mean) -> mean_mrna
tapply(log(dfly$te), INDEX = dfly$bin_score,FUN = mean) -> mean_te

tapply(dfly$rpf, INDEX = dfly$bin_score,FUN = mean) -> mean_rpf
tapply(dfly$mrna, INDEX = dfly$bin_score,FUN = mean) -> mean_mrna
tapply(dfly$te, INDEX = dfly$bin_score,FUN = mean) -> mean_te

bar1 <- tapply(log2(dfly$rpf), INDEX = dfly$bin_score,FUN = sd)/sqrt(
  tapply(log2(dfly$mrna), INDEX = dfly$bin_score,FUN = length))
bar2 <- tapply(log2(dfly$mrna), INDEX = dfly$bin_score,FUN = sd)/sqrt(
  tapply(log2(dfly$mrna), INDEX = dfly$bin_score,FUN = length))
bar3 <- tapply(log2(dfly$te), INDEX = dfly$bin_score,FUN = sd)/sqrt(
  tapply(log2(dfly$te), INDEX = dfly$bin_score,FUN = length))
bar1

(pa1 <- p + aes(x= seq(1,10), y= mean_rpf) + geom_point(size= 2.5) + 
    scale_x_discrete(limits = seq(1,10))
  + geom_errorbar(aes(ymax= mean_rpf+bar1,ymin= mean_rpf-bar1),width =.1,size =.5)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.2, vjust = 2, 
             label = "rho 0.90075 \n\rp-value <2e-16" , size = 5) + labs(x = "initiation socor bins",y = "RPF")
)

cor.test(dfly$score_21nt,dfly$rpf,method = "s",exact = F)
cor.test(dfly$score_21nt,dfly$mrna,method = "s",exact = F)
cor.test(dfly$score_21nt,dfly$te,method = "s",exact = F)


grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(7,1))) ####将页面分成2*2矩阵
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p_s19_logo, vp = vplayout(1,1))   ###将（1,1)和(1,2)的位置画图c
print(p_s29_logo, vp = vplayout(2,1))   ###将(2,1)的位置画图b
print(p_s37_logo, vp = vplayout(3,1))  ###将（2,2)的位置画图a
# print(p_gene_logo, vp = vplayout(4,1:3))  ###将（2,2)的位置画图a
print(p_gene_logo, vp = vplayout(4,1))
# print(p_ram_logo, vp = vplayout(5,1:3))  ###将（2,2)的位置画图a
print(p_rs28_logo, vp = vplayout(5,1))  ###将（2,2)的位置画图a
print(p_ssu47_logo, vp = vplayout(6,1))
print(p_ssu55_logo, vp = vplayout(7,1))



# mouse -------------------------------------------------------------------
dmouse <- read_csv("G:/kozak_new/initiation_score/mouse/final_mouse.csv",col_names = T)
dmouse <- read_delim("G:/kozak_new/initiation_score/mouse/3-score-ribo-mouse.txt",delim = "\t") %>% filter(rpf >0)
head(dmouse);nrow(dmouse) # 4523
as.numeric(cut2(dmouse$score_21nt,g = 10)) -> dmouse$bin_score
tapply(log2(dmouse$rpf), INDEX = dmouse$bin_score,FUN = mean) -> mean_rpf
tapply(log2(dmouse$mrna), INDEX = dmouse$bin_score,FUN = mean) -> mean_mrna
tapply(log2(dmouse$te), INDEX = dmouse$bin_score,FUN = mean) -> mean_te

tapply(dmouse$rpf, INDEX = dmouse$bin_score,FUN = mean) -> mean_rpf
tapply(dmouse$mrna, INDEX = dmouse$bin_score,FUN = mean) -> mean_mrna
tapply(dmouse$te, INDEX = dmouse$bin_score,FUN = mean) -> mean_te

bar1 <- tapply(log2(dmouse$rpf), INDEX = dmouse$bin_score,FUN = sd)/sqrt(
  tapply(log2(dmouse$mrna), INDEX = dmouse$bin_score,FUN = length))
bar2 <- tapply(log2(dmouse$mrna), INDEX = dmouse$bin_score,FUN = sd)/sqrt(
  tapply(log2(dmouse$mrna), INDEX = dmouse$bin_score,FUN = length))
bar3 <- tapply(log2(dmouse$te), INDEX = dmouse$bin_score,FUN = sd)/sqrt(
  tapply(log2(dmouse$te), INDEX = dmouse$bin_score,FUN = length))
bar1;bar2;bar3

bar1 <- tapply(log2(dmouse$rpf), INDEX = dmouse$bin_score,FUN = sd)
bar2 <- tapply(log2(dmouse$mrna), INDEX = dmouse$bin_score,FUN = sd)
bar3 <- tapply(log2(dmouse$te), INDEX = dmouse$bin_score,FUN = sd)
bar1;bar2;bar3

cor.test(dmouse$score_21nt,dmouse$rpf,method = "s",exact = F) #  rho 0.12155\n\rp-value = 2.4e-16
cor.test(dmouse$score_21nt,dmouse$mrna,method = "s",exact = F) #  rho 0.1013\n\rp-value = 1.9e-05
cor.test(dmouse$score_21nt,dmouse$te,method = "s",exact = F) #  rho 0.063593\n\rp-value = 1.9e-05

(pa1 <- pp + aes(x= seq(1,10), y= mean_rpf) + geom_point(size= 7) + 
    scale_x_discrete(limits = seq(1,10))
  + geom_errorbar(aes(ymax= mean_rpf+bar1,ymin= mean_rpf-bar1),width =.1,size =1.5)
  + ylim(c(9.1,10.2))
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.02, vjust = 1.0, 
             label = "rho= 0.12155\n\rp-value = 2.4e-16" , size = 14) + labs(x = "initiation socor bins",y = "RPF")
)

grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(7,1))) ####将页面分成2*2矩阵
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p_s19_logo, vp = vplayout(1,1))   ###将（1,1)和(1,2)的位置画图c
print(p_s29_logo, vp = vplayout(2,1))   ###将(2,1)的位置画图b
print(p_s37_logo, vp = vplayout(3,1))  ###将（2,2)的位置画图a
# print(p_gene_logo, vp = vplayout(4,1:3))  ###将（2,2)的位置画图a
print(p_gene_logo, vp = vplayout(4,1))
# print(p_ram_logo, vp = vplayout(5,1:3))  ###将（2,2)的位置画图a
print(p_rs28_logo, vp = vplayout(5,1))  ###将（2,2)的位置画图a
print(p_ssu47_logo, vp = vplayout(6,1))
print(p_ssu55_logo, vp = vplayout(7,1))




# yeast -------------------------------------------------------------------
# yeast MS -------------------------------------------------------------------
dy_ms <- read_delim("G:/kozak_new/initiation_score/yeast/ms/2-protein_intensity.txt",delim = "\t") # 784
head(dy_ms)
cor.test(dy_ms$kozak_score,dy_ms$average_ho,method = "s") #0.19083 p-value = 7.7e-08
cor.test(dy_ms$kozak_score,dy_ms$average_up,method = "s") #0.19239 p-value = 6e-08

#下面是用的我的 score，全Y 的
dscore <- read_csv("G:/kozak/iniation_score/gene_21nt_aug_score.csv", col_names = T ) # 5884
head(dscore)
d1 <- select(dy_ms,c("gene","average_ho","average_up")) %>%
  merge(dscore, by.x = "gene",by.y = "gene_id")
head(d1) # 783
cor.test(d1$score_21nt,d1$average_ho,method = "s") # ρ=0.19248 \n\rP=6.1e-08
cor.test(d1$score_21nt,d1$average_up,method = "s") #ρ=0.19415 \n\rP=4.6e-08

as.numeric(cut2(d1$score_21nt,g = 10)) -> d1$bin_score  ## 可以知道是属于哪个bin的
tapply(log2(d1$average_ho), INDEX = d1$bin_score,FUN = mean) -> mean_ms
# (cor1 <- cor.test(Data$ribo,Data$kozak_score,method = "s"))  #0.08748553
bar1 <- tapply(log2(d1$average_ho), INDEX = d1$bin_score, FUN = sd)/sqrt(
  tapply(log2(d1$average_ho), INDEX = d1$bin_score, FUN = length))
bar1 # max = 

(p1b <- p + aes(x= seq(1,10), y= mean_ms) 
  + geom_point(size= 7) 
  + scale_x_discrete(limits = seq(1,10))
  + geom_errorbar(aes(ymax= mean_ms+bar1,ymin= mean_ms-bar1),width =.2,size =1.5)
  + ylim(c(16,19.5))
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.1, vjust =1.0, 
             label = "rho=0.19248 \n\rP=6.1e-08" , size = 14) + labs(x = "initiation socor bins", y = "protein level (MS)")
)


# MS/mrna ----------------------------------------------------------------
head(d1) # 783
dy_ms <- read_delim("G:/kozak_new/initiation_score/yeast/ms/2-protein_intensity.txt",delim = "\t") # 784
dscore <- read_csv("G:/kozak/iniation_score/gene_21nt_aug_score.csv", col_names = T ) # 5884
drpkm <- read_csv("G:/kozak/riboseq/Ho-1-T_rrpkm.csv",col_names = T) #5764
# drpf <- read_csv("G:/kozak/riboseq/Ho-1-28_rrpf.csv",col_names = T) # 5489

d1 <- select(dy_ms,c("gene","average_ho","average_up")) %>% rename(gene_id = gene) %>%
  merge(dscore, by="gene_id") %>% merge(drpkm,by="gene_id") %>% 
  mutate(MSho_rpkm = average_ho/rrpkm,MSup_rpkm = average_up/rrpkm) #%>% merge(drpf,by="gene_id")
head(d1) # 783

cor.test(d1$score_21nt,d1$MSho_rpkm,method = "s") # ρ=0.090879 \n\rP=0.011
cor.test(d1$score_21nt,d1$MSup_rpkm,method = "s") # ρ=0.092365 \n\rP=0.0097

as.numeric(cut2(d1$score_21nt,g = 10)) -> d1$bin_score  ## 可以知道是属于哪个bin的
tapply(log2(d1$MSho_rpkm), INDEX = d1$bin_score,FUN = mean) -> mean_ms_rpkm
# (cor1 <- cor.test(Data$ribo,Data$kozak_score,method = "s"))  #0.08748553
bar2 <- tapply(log2(d1$MSho_rpkm), INDEX = d1$bin_score, FUN = sd)/sqrt(
  tapply(log2(d1$MSho_rpkm), INDEX = d1$bin_score, FUN = length))
bar2 # max = 

(p2b <- p + aes(x= seq(1,10), y= mean_ms_rpkm) 
  + geom_point(size= 7) 
  + scale_x_discrete(limits = seq(1,10))
  + geom_errorbar(aes(ymax= mean_ms_rpkm+bar2,ymin= mean_ms_rpkm-bar2),width =.2,size =1.5)
  + ylim(c(7,12))
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.1, vjust =1.0, 
             label = "rho=0.090879 \n\rP=0.011" , size = 14) + labs(x = "initiation socor bins", y = "MS/mRNA_level")
)
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(5,9))) ####将页面分成2*2矩阵
print(p1b, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p2b, vp = vplayout(2:4,6:8))
dev.off()
# 18 10 inch



# mRNA half life (degradation rate) ---------------------------------------
#下面是用的我的 score，全Y 的
dscore <- read_csv("~/rstudio/kozak/gene_21nt_aug_score.csv", col_names = T ) %>%
  select(-"reads_21nt")# 5884
degradation <- read_csv("~/rstudio/kozak/mRNA half life.csv", col_names = T) %>% 
  rename(gene_id = "Gene ID",polyA = "poly(A)+ Half-life",total = "Total Half-life") %>%
  select(gene_id,polyA,total)# 3890
head(dscore);head(degradation)
dat <- merge(dscore,degradation,by="gene_id",all.y = T) # 3890
head(dat)

cor.test(dat$score_21nt,dat$polyA,method = "s",exact = F) # ρ=0.090879 P=0.51
cor.test(dat$score_21nt,dat$total,method = "s",exact = F) # ρ=0.073787 P=4.1e-06

as.numeric(cut2(dat$score_21nt,g = 10)) -> dat$bin_score  ## 可以知道是属于哪个bin的
tapply(log2(dat$total), INDEX = dat$bin_score,FUN = mean) -> mean_total
# (cor1 <- cor.test(Data$ribo,Data$kozak_score,method = "s"))  #0.08748553
bar2 <- tapply(log2(dat$total), INDEX = dat$bin_score, FUN = sd)/sqrt(
  tapply(log2(dat$total), INDEX = dat$bin_score, FUN = length))
bar2 # max = 

(p4b <- p + aes(x= seq(1,10), y= mean_total) 
  + geom_point(size= 7) 
  + scale_x_discrete(limits = seq(1,10))
  + geom_errorbar(aes(ymax= mean_total+bar2,ymin= mean_total-bar2),width =.2,size =1.5)
  + ylim(c(2,4))
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.1, vjust =1.0, 
             label = "rho=0.073787 \n\rP=4.1e-06" , size = 14) + labs(x = "initiation socor bins", y = "mRNA half-life(log2)")
)
head(dat)
(p5b <- p%+%dat +aes(x=score_21nt,y=log2(total))
  + geom_point(size=.5)
  + annotate(geom = "text", x = -Inf , y = Inf , hjust = -0.1, vjust =1.0, 
             label = "rho=0.073787 \n\rP=4.1e-06" , size = 14) + labs(x = "initiation socor bins", y = "mRNA half-life(log2)")
)

grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(5,9))) ####将页面分成2*2矩阵
print(p4b, vp = vplayout(2:4,2:4))   ###将(2,1)的位置画图b
print(p5b, vp = vplayout(2:4,6:8))
dev.off()
# 18 10 inch






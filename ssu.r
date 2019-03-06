# enviroment --------------------------------------------------------------
setwd("G:/kozak/ssu/")
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
# file chunk --------------------------------------------------------------
# fs5utr <- "G:/kozak/ssu/last_utr_final.csv"
# f5utr <- "G:/common_used/Nagalakshmi_2008_5UTRs_V64.csv"
# fuorf <- "G:/common_used/uORF.csv"
###
dH028 <- read_csv("G:/kozak/riboseq/last_Ho-1-28_final.csv", col_names = T )
dH0T <- read_csv("G:/kozak/riboseq/last_Ho-1-T_final.csv", col_names = T )

dssu <- read_csv("G:/kozak/ssu/last_ssu_final.csv", col_names = T )
d40srs <- read_csv("G:/kozak/ssu/last_rs_final.csv", col_names = T )
d40smrna <- read_csv("G:/kozak/ssu/last_mrna_final.csv", col_names = T )

d21nt <- read_csv("G:/kozak/iniation_score/gene_21nt_aug_score.csv", col_names = T )
dexon <- read_csv("G:/common_used/gene_exonlength.csv", col_names = T )
drrpkm <- read_csv("G:/kozak/riboseq/Ho-1-T_rrpkm.csv", col_names = T )
drrpf <- read_csv("G:/kozak/riboseq/Ho-1-28_rrpf.csv", col_names = T )
dsrpkm <- read_csv("G:/kozak/ssu/ssu_srpkm.csv", col_names = T )
###

###
dH028 <- cbind(dH028, id = rep("rrs", nrow(dH028)))
dH0T <- cbind(dH0T, id = rep("rmrna", nrow(dH0T)))
dssu <- cbind(dssu, id = rep("ssu", nrow(dssu)))
d40srs <- cbind(d40srs, id = rep("srs", nrow(d40srs)))
d40smrna <- cbind(d40smrna, id = rep("smrna", nrow(d40smrna)))

data <- rbind(dH028,dH0T,dssu, d40srs,d40smrna ) #%>% as.tibble()
head(data)

# basic analysis ----------------------------------------------------------
# if(find("%+%")[1]=="package:psych") stop("oops")
#### 看一下 riboseq 数据
####### rrs
(p_rs <- p%+%(filter(data, id == "rrs")) + 
    aes(length, fill = factor(frame)) + geom_bar(position = "dodge") + xlim(c(15,35)) + 
   labs(title = "r rs frame distribution"))

(p_rs <- p%+%(filter(data, id == "rrs")) + xlim(c(15,35)) + 
    aes(length) + geom_bar() +
    labs(title = "r rs length distribution"))

(p_p5_rs <- p%+%(filter(data, id == "rrs", length == 28)) + aes(p5) + geom_bar()  + xlim(c(-25,200))
  + labs(title = "ribo-seq 5' end", x = "5' end position"))
(p_p5_rs <- p%+%(filter(data, id == "rrs", length == 21)) + aes(p5) + geom_bar() +ylim(c(0,2000))  + xlim(c(-25,200))
  + labs(title = "ribo-seq 5' end", x = "5' end position"))
(p_p5_rs <- p%+%(filter(data, id == "rrs")) + aes(pp5) + geom_bar()  
  + labs(title = "ribo-seq 5' end", x = "5' end position"))
### 这里使用 facet 
#### srs
(p_rs <- p%+%(filter(data, id == "srs")) + 
    aes(length, fill = factor(frame)) + geom_bar(position = "dodge") + xlim(c(10,45))+
    labs(title = "s rs length frame distribution"))

(p_rs <- p%+%(filter(data, id == "srs")) + xlim(c(10,45)) + 
    aes(length) + geom_bar() +
    labs(title = "s rs length distribution"))
(p_p5_rs <- p%+%(filter(data, id == "srs", length == 28)) + aes(p5) + geom_bar()  + xlim(c(-25,200))
  + labs(title = "ribo-seq 5' end", x = "5' end position"))
# (p_p5_rs <- p%+%(filter(data, id == "srs", length == 21)) + aes(p5) + geom_bar()  + xlim(c(-25,200))
#   + labs(title = "ribo-seq 5' end", x = "5' end position"))
(p_p5_rs <- p%+%(filter(data, id == "srs", length == 31)) + aes(p5) + geom_bar()  + xlim(c(-25,200))
  + labs(title = "ribo-seq 5' end", x = "5' end position"))
(p_p5_rs <- p%+%(filter(data, id == "srs")) + aes(pp5) + geom_bar()  
  + labs(title = "ribo-seq 5' end", x = "5' end position"))

#### ssu
(p_ <- p%+%(filter(data, id == "ssu")) + aes(length, fill = factor(frame)) 
  + geom_bar(position = "dodge") + xlim(c(10,60)) + labs(title = "ssu length frame distribution"))

(p_rs <- p%+%(filter(data, id == "ssu")) + 
    aes(length) + geom_bar() + xlim(c(10,60))+
    labs(title = "ssu length distribution"))

(p_p5_rs <- p%+%(filter(data, id == "ssu")) + aes(p5) + geom_bar() +xlim(c(-40,10))
  + labs(title = "ssu 5' end", x = "5' end position") )
(p_p5_rs <- p%+%(filter(data, id == "ssu")) + aes(p3) + geom_bar() +xlim(c(0,40))
  + labs(title = "ssu 3' end", x = "3' end position") )

(p_p5_rs <- p%+%(filter(data, id == "ssu")) + aes(pp5) + geom_bar() +xlim(c(-40,10))
  + labs(title = "ssu 5' end", x = "5' end position") )
(p_p5_rs <- p%+%(filter(data, id == "ssu")) + aes(pp3) + geom_bar() + xlim(c(0,40))
  + labs(title = "ssu 5' end", x = "5' end position") )
(p_p5_rs <- p%+%(filter(data, id == "ssu")) + aes(pp3) + geom_bar() + xlim(c(-400,40))
  + labs(title = "ssu pp3' end", x = "pp3' end position") )

# mrna
(p_lenfra_mrna <- p%+%(filter(data, id == "rmrna")) + 
    aes(length, fill = factor(frame)) + geom_bar(position = "dodge")
 + labs(title = "rmrna length frame distribution"))
(p_lenfra_mrna <- p%+%(filter(data, id == "smrna")) + 
    aes(length, fill = factor(frame)) + geom_bar(position = "dodge") 
     + labs(title = "smrna length frame distribution"))

(p_lenfra_mrna <- p%+%(filter(data, id == "rmrna")) + 
    aes(length) + geom_bar()
  + labs(title = "rmrna length frame distribution"))
(p_lenfra_mrna <- p%+%(filter(data, id == "smrna")) + 
    aes(length) + geom_bar() 
  + labs(title = "smrna length frame distribution"))

########## syo?????? #################
head(data)
filter(data,id == c("syossu","ssu")) %>% group_by(id,gene_name) %>% summarise(num = n()) -> data_test
head(data_test)
arrange(data_test,desc(num)) -> data_test
head(data_test)

data_y <- filter(data, gene_name == "YJR045C", id == c("syossu","ssu") )
data_y <- filter(data, gene_name == "YFL039C", id == c("syossu","ssu") )
data_y <- filter(data, gene_name == "YER165W", id == c("syossu","ssu") )
data_y <- filter(data, gene_name == "YLR044C", id == c("syossu","ssu") )
head(data_y)
length(filter(data_y, id == "ssu")$id)
length(filter(data_y, id == "syossu")$id)
length(filter(data_y, id == "ssu", p5 <=2 & p3 >= 0)$id)
length(filter(data_y, id == "syossu", p5 <=2 & p3 >= 0)$id)

data %>% filter(id == c("syossu","ssu"),p5 <=2 & p3 >= 0) -> data_aug;length(data_aug$id)
length(filter(data_aug, id == "syossu")$id)
length(filter(data_aug, id == "ssu")$id)

(plot_lenp5_ssu <- p%+%(filter(data_y, id == "ssu")) + aes(p5) + 
    geom_bar() + xlim(c(-20,8)) + labs(title = "ssu 5' end distribution"))
(plot_lenp5_ssu <- p%+%(filter(data_y, id == "syossu")) + aes(p5) + 
    geom_bar() + xlim(c(-20,8)) + labs(title = "syossu 5' end distribution"))

(p_p3_ssu <- p%+%(filter(data_y, id == "ssu")) + aes(x = p3) + 
    geom_bar() + xlim(c(-5,40)) + labs(title = "ssu 3' end" , x = "3' end position"))
(p_p3_ssu <- p%+%(filter(data_y, id == "syossu")) + aes(x = p3) + 
    geom_bar() + xlim(c(-5,40)) + labs(title = "syossu 3' end" , x = "3' end position"))
(plot_lenp5_ssu <- p%+%(filter(data_y, id == "ssu", p5 == (-12))) + aes(length) + 
    geom_bar() + xlim(c(16,50)) + labs(title = "ssu 5' end == -12 length-distribution"))
(plot_lenp5_ssu <- p%+%(filter(data_y, id == "syossu", p5 == (-12))) + aes(length) + 
    geom_bar() + xlim(c(16,50)) + labs(title = "syossu 5' end == -12 length-distribution"))

# #### p5 p3 
# ( plot_19_p5 <- p%+%(filter(data, id == "ssu", length == 19)) + aes(p5) + geom_bar() + xlim(c(-50,10))+  labs(title = "ssu 19nt 5' end", x = "5' end position"))
# ( plot_19_p3 <- p%+%(filter(data, id == "ssu", length == 19)) + aes(p3) + geom_bar() + xlim(c(-5,30))+  labs(title = "ssu 19nt 3' end", x = "3' end position"))
# 
# ( plot_29_p5 <- p%+%(filter(data, id == "ssu", length == 29)) + aes(p5) + geom_bar() + xlim(c(-50,10))+  labs(title = "ssu 29nt 5' end", x = "5' end position"))
# ( plot_29_p3 <- p%+%(filter(data, id == "ssu", length == 29)) + aes(p3) + geom_bar() + xlim(c(-5,40))+  labs(title = "ssu 29nt 3' end", x = "3' end position"))
# 
# ( plot_37_p5 <- p%+%(filter(data, id == "ssu", length == 37)) + aes(p5) + geom_bar() + xlim(c(-50,10))+  labs(title = "ssu 37nt 5' end", x = "5' end position"))
# ( plot_37_p3 <- p%+%(filter(data, id == "ssu", length == 37)) + aes(p3) + geom_bar() + xlim(c(-5,40))+  labs(title = "ssu 37nt 3' end", x = "3' end position"))

######## ##########################热图
head(data)
#### 分组摘要 
heat_map <- data %>% group_by(id,p5,length) %>% summarise(count = n())
head(heat_map)
# + scale_fill_gradient(low = "white", high = "red")
h <- ggplot(heat_map) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), 
    legend.position = "right",legend.text =  element_text(face="bold", size=20),
    axis.text.x = element_text(size = 26,face = "bold"), 
    axis.text.y = element_text(size = 26,face = "bold"), axis.title = element_text(size = 26, face = "bold"))
# ,panel.grid =element_blank()
( h_point <- h%+%(filter(heat_map, id == "ssu")) + aes(x = p5, y = length , fill = count ) 
  + geom_tile() + scale_fill_gradientn(colours = terrain.colors(10)) 
  + scale_x_continuous(limits = c(-40,30), breaks = seq(-40,30,10)) + ylim(c(17,100))
  + labs(title = "ssu 5' end position"))

( h_point <- h%+%(filter(heat_map, id == "rs")) + aes(x = p5, y = length , fill = count ) 
  + geom_tile() + scale_fill_gradientn(colours = terrain.colors(10)) 
  + scale_x_continuous(limits = c(-40,30), breaks = seq(-40,30,10)) + ylim(c(17,100))
  + labs(title = "ssu 5' end position"))

( h_point <- h%+%(filter(heat_map, id == "syors")) + aes(x = p5, y = length , fill = count ) 
  + geom_tile() + scale_fill_gradientn(colours = terrain.colors(10)) 
  + scale_x_continuous(limits = c(-40,30), breaks = seq(-40,30,10)) + ylim(c(17,100))
  + labs(title = "ssu 5' end position"))

( h_point <- h%+%(filter(heat_map, id == "mrna")) + aes(x = p5, y = length , fill = count ) 
  + geom_tile() + scale_fill_gradientn(colours = terrain.colors(10), limits =c(0,150)) 
  + scale_x_continuous(limits = c(-40,30), breaks = seq(-40,30,10)) + ylim(c(17,100))
  + labs(title = "ssu 5' end position"))

( h_point <- h%+%(filter(heat_map, id == "syossu")) + aes(x = p5, y = length , fill = count ) 
  + geom_tile() + scale_fill_gradientn(colours = terrain.colors(10)) 
  + scale_x_continuous(limits = c(-40,0), breaks = seq(-40,0,10)) + ylim(c(17,100))
  + labs(title = "ssu 5' end position"))


heat_map <- data %>% group_by(id,p3,length) %>% summarise(count = n())
head(heat_map,100)
?scale_fill_gradientn
( h_point <- h%+%(filter(heat_map, id == "ssu")) + aes(x = p3, y = length , fill = count ) 
  + geom_tile() + scale_fill_gradientn(colours = terrain.colors(10)) 
  + scale_x_continuous(limits = c(-10,40), breaks = seq(-10,40,10))
  + labs(title = "ssu 3' position")
  + scale_y_continuous(limits = c(17,100)))

( h_point <- h%+%(filter(heat_map, id == "syossu")) + aes(x = p3, y = length , fill = count ) 
  + geom_tile() + scale_fill_gradientn(colours = terrain.colors(10)) 
  + scale_x_continuous(limits = c(-10,40), breaks = seq(-10,40,10))
  + labs(title = "ssu 3' position")
  + scale_y_continuous(limits = c(17,100)))

# ggplotly()
( plot_point2_p3 <- p%+%(filter(data, id == "rs")) + aes(x = p5, y = length) + geom_hex() + xlim(c(-80,40)))

############   算 score chunk， 已核对是准确的  ##################
##注意取的是 12个碱基的 要稍修改下
position_score<-read_csv('G:/kozak/iniation_score/position_score_21nt.csv',col_names = T)
position_score
# 因为 19 29 37 最短是19 
position_score <- position_score[,3:21]
position_score
colnames(position_score) <- c(-12:-1,1:7)
position_score <- as.data.frame(position_score)
position_score
#计算每条序列的score, 值一一替换进去，每行每行的来
sum2 <- function(x){
  as.numeric(x) ->x
  result <- sum(x)-6
  return(result)
}
get_score_value <- function(use2_data,name){
  data_matrix <- str_split(use2_data$gene_read,"", simplify = T)
  for (i in 1:nrow(data_matrix)){
    for (ii in 1:ncol(data_matrix)){
      if (data_matrix[i,ii] == 'A'){
        data_matrix[i,ii] <- position_score[1,ii]
      }else if (data_matrix[i,ii] == 'G'){
        data_matrix[i,ii] <- position_score[3,ii]
      }else if (data_matrix[i,ii] == 'C'){
        data_matrix[i,ii] <- position_score[2,ii]
      }else if (data_matrix[i,ii] == 'TRUE'){
        data_matrix[i,ii] <- position_score[4,ii]
      }else if (data_matrix[i,ii] == 'T'){
        data_matrix[i,ii] <- position_score[4,ii]
      }
    }
  }
  ####计算每行的score之和
  apply(data_matrix, 1, sum2) -> score_value
  as.data.frame(score_value) -> score_value
  ####把算好的 score_value 合并进 data里
  id <- rep(name,length(score_value$score_value))
  cbind(id,score_value,reads = use2_data$gene_read, gene_id = use2_data$gene_id) -> md
  return(md)
}
# test<- data.frame(x = c(1:10, 1:3), y=1:13)
# test
# duplicated(test$x)
# test[!duplicated(test$x), ]
#######  get value unique
getuniqScorevalue <- function(use_data,name,long){
  ####  这一步是数据框去重
  uni_data <- use_data[!duplicated(use_data$gene_id),]
  uni_data$gene_read <- str_sub(uni_data$gene_read,start = -long, end = (-long)+18)
  # print(head(use_data));print(head(uni_data))
  print(length(use_data$gene_read))
  print(length(uni_data$gene_read))  # 
  get_score_value(uni_data,name) -> md
  return(md)
}
# 不 unique 的
getScorevalue <- function(use_data,name,long){
  use_data$gene_read <- str_sub(use_data$gene_read,start = -long, end = (-long)+18)
  # print(head(use_data))
  # uni_data <- unique(sub_data)
  print(length(use_data$gene_read));
  print(length(unique(use_data$gene_read)))  # 
  get_score_value(use_data,name) -> md
  return(md)
}

# #### rpkm
# head(data_rpkm);head(data_gene)
# data_gene2 <- select(data_gene, gene_id, seq); head(data_gene2)
# length(data_gene$gene_id);length(data_gene2$gene_id)  # 5884
# ### change seq to 12 nt
# str_sub(data_gene2$seq, start = 2, end = 20) -> data_gene2$seq ;head(data_gene2);length(data_gene2$gene_id) # 5884

# density plot ------------------------------------------------------------
head(data)
use_s19 <- filter(data, id == "ssu", p3 == 6, length == 19)
# use_s19 <- filter(data, id == "syossu", p3 == 6, length == 19)
head(use_s19)
getuniqScorevalue(use_s19,"ssu19",19) -> uni19 
head(uni19);length(uni19$id)
getScorevalue(use_s19,"ssu19",19) -> ssu19 ; head(ssu19);length(ssu19$id);length(unique(ssu19$reads))
ssu19 <- cbind(ssu19,freq = rep(1/length(ssu19$id),length(ssu19$id)) )
head(ssu19);length(ssu19$id)
new_ssu19 <- ssu19 %>% group_by(id,score_value,reads,gene_id) %>% summarise(count = n())
head(new_ssu19);length(new_ssu19$id)
#这句 control 了 mrna 
head(data_rpkm)
rpkm_ssu19 <- merge(new_ssu19,data_rpkm,by = "gene_id")
# head(rpkm_ssu19)
rpkm_ssu19 <- cbind(rpkm_ssu19, freq = rpkm_ssu19$count/log2(rpkm_ssu19$rpkm1))
rpkm_ssu19 <- cbind(rpkm_ssu19, freq2 = rpkm_ssu19$count/rpkm_ssu19$rpkm1)
rpkm_ssu19$freq <- rpkm_ssu19$freq/sum(rpkm_ssu19$freq)
rpkm_ssu19$freq2 <- rpkm_ssu19$freq2/sum(rpkm_ssu19$freq2)
head(rpkm_ssu19);length(rpkm_ssu19$gene_id)
# sum(rpkm_ssu19$freq)
# sum(rpkm_ssu19$freq2)
### 这个 把 new 进行加权, 求和为1
new_ssu19$freq <- (new_ssu19$count)/(length(ssu19$id))
head(new_ssu19);length(new_ssu19$id)
sum(new_ssu19$count);sum(new_ssu19$freq)

####
use_s29 <- filter(data, id == "ssu", p3 == 16, length == 29 )
# use_s29 <- filter(data, id == "syossu", p3 == 16, length == 29 )
getuniqScorevalue(use_s29,"ssu29",29) -> uni29
getScorevalue(use_s29,"ssu29",29) -> ssu29
ssu29 <- cbind(ssu29,freq = rep(1/length(ssu29$id),length(ssu29$id)) )
head(ssu29);length(ssu29$id);length(unique(ssu29$reads))
new_ssu29 <- ssu29 %>% group_by(id,score_value,reads,gene_id) %>% summarise(count = n())
head(new_ssu29);length(new_ssu29$id)
#这句 control 了 mrna 
head(data_rpkm)
rpkm_ssu29 <- merge(new_ssu29,data_rpkm,by = "gene_id")
# head(rpkm_ssu29)
rpkm_ssu29 <- cbind(rpkm_ssu29, freq = rpkm_ssu29$count/log2(rpkm_ssu29$rpkm1))
rpkm_ssu29 <- cbind(rpkm_ssu29, freq2 = rpkm_ssu29$count/rpkm_ssu29$rpkm1)
rpkm_ssu29$freq <- rpkm_ssu29$freq/sum(rpkm_ssu29$freq)
rpkm_ssu29$freq2 <- rpkm_ssu29$freq2/sum(rpkm_ssu29$freq2)
head(rpkm_ssu29);length(rpkm_ssu29$gene_id)
sum(rpkm_ssu29$freq)
### 这个 把 new 进行加权, 求和为1
new_ssu29$freq <- (new_ssu29$count)/(length(ssu29$id))
head(new_ssu29);length(new_ssu29$id)
sum(new_ssu29$count);sum(new_ssu29$freq)

####
use_s47 <- filter(data, id == "ssu", p3 == 16, length == 47 )
# use_s47 <- filter(data, id == "syossu", p3 == 16, length == 47 )
getuniqScorevalue(use_s47,"ssu47",29) -> uni47
getScorevalue(use_s47,"ssu47",29) -> ssu47
ssu47 <- cbind(ssu47,freq = rep(1/length(ssu47$id),length(ssu47$id)) )
head(ssu47);length(unique(ssu47$reads))
new_ssu47 <- ssu47 %>% group_by(id,score_value,reads,gene_id) %>% summarise(count = n())
head(new_ssu47)
#这句 control 了 mrna 
head(data_rpkm)
rpkm_ssu47 <- merge(new_ssu47,data_rpkm,by = "gene_id")
# head(rpkm_ssu47)
rpkm_ssu47 <- cbind(rpkm_ssu47, freq = rpkm_ssu47$count/log2(rpkm_ssu47$rpkm1))
rpkm_ssu47 <- cbind(rpkm_ssu47, freq2 = rpkm_ssu47$count/rpkm_ssu47$rpkm1)
rpkm_ssu47$freq <- rpkm_ssu47$freq/sum(rpkm_ssu47$freq)
rpkm_ssu47$freq2 <- rpkm_ssu47$freq2/sum(rpkm_ssu47$freq2)
head(rpkm_ssu47);length(rpkm_ssu47$gene_id)
sum(rpkm_ssu47$freq)
### 这个 把 new 进行加权, 求和为1
new_ssu47$freq <- (new_ssu47$count)/(length(ssu47$id))
head(new_ssu47);length(new_ssu47$id)
sum(new_ssu47$freq)

# suse_s19 <- filter(data, id == "syossu", p3 == 6, length == 19)
# suse_s47 <- filter(data, id == "syossu", p3 == 16, length == 47 )

####
use_s37 <- filter(data, id == "ssu", p3 == 24, length == 37 )
# use_s37 <- filter(data, id == "syossu", p3 == 24, length == 37 )
getuniqScorevalue(use_s37,"ssu37",37) -> uni37
getScorevalue(use_s37,"ssu37",37) -> ssu37
ssu37 <- cbind(ssu37,freq = rep(1/length(ssu37$id),length(ssu37$id)) )
head(ssu37);length(unique(ssu37$reads))
new_ssu37 <- ssu37 %>% group_by(id,score_value,reads,gene_id) %>% summarise(count = n())
#这句 control 了 mrna 
head(data_rpkm)
rpkm_ssu37 <- merge(new_ssu37,data_rpkm,by = "gene_id")
# head(rpkm_ssu37)
rpkm_ssu37 <- cbind(rpkm_ssu37, freq = rpkm_ssu37$count/log2(rpkm_ssu37$rpkm1))
rpkm_ssu37 <- cbind(rpkm_ssu37, freq2 = rpkm_ssu37$count/rpkm_ssu37$rpkm1)
rpkm_ssu37$freq <- rpkm_ssu37$freq/sum(rpkm_ssu37$freq)
rpkm_ssu37$freq2 <- rpkm_ssu37$freq2/sum(rpkm_ssu37$freq2)
head(rpkm_ssu37);length(rpkm_ssu37$gene_id)
sum(rpkm_ssu37$freq)
### 这个 把 new 进行加权, 求和为1
new_ssu37$freq <- (new_ssu37$count)/(length(ssu37$id))
sum(new_ssu37$count);sum(new_ssu37$freq)
head(new_ssu37)
###

use_s55 <- filter(data, id == "ssu", p3 == 24, length == 55 )
use_s55 <- filter(data, id == "syossu", p3 == 24, length == 55 )
getuniqScorevalue(use_s55,"ssu55",37) -> uni55
getScorevalue(use_s55,"ssu55",37) -> ssu55
ssu55 <- cbind(ssu55,freq = rep(1/length(ssu55$id),length(ssu55$id)) )
head(ssu55);length(unique(ssu55$reads))
new_ssu55 <- ssu55 %>% group_by(id,score_value,reads,gene_id) %>% summarise(count = n())
head(new_ssu55)
#这句 control 了 mrna 
head(data_rpkm)
rpkm_ssu55 <- merge(new_ssu55,data_rpkm,by = "gene_id")
# head(rpkm_ssu55)
rpkm_ssu55 <- cbind(rpkm_ssu55, freq = rpkm_ssu55$count/log2(rpkm_ssu55$rpkm1))
rpkm_ssu55 <- cbind(rpkm_ssu55, freq2 = rpkm_ssu55$count/rpkm_ssu55$rpkm1)
rpkm_ssu55$freq <- rpkm_ssu55$freq/sum(rpkm_ssu55$freq)
rpkm_ssu55$freq2 <- rpkm_ssu55$freq2/sum(rpkm_ssu55$freq2)
head(rpkm_ssu55);length(rpkm_ssu55$gene_id)
sum(rpkm_ssu55$freq)
### 这个 把 new 进行加权, 求和为1
new_ssu55$freq <- (new_ssu55$count)/(length(ssu55$id))
sum(new_ssu55$count);sum(new_ssu55$freq)

### 无意义
# use_ram <- head(data3$ram_aug,2000 )
# head(use_ram)
# getuniqScorevalue(use_ram,"ram",20) -> uniram
# getScorevalue(use_ram,"ram",20) -> ram
# head(ram);length(unique(ram$reads))
# new_ram <- ram %>% group_by(id,score_value,reads) %>% summarise(count = n())
# new_ram$count <- (new_ram$count)/(length(ram$id))
# # sum(new_ram$count)
# head(new_ram);length(new_ram$id)


## 所有的 gene
# use_gene <- str_sub(data_gene$seq,start = 2, end = 20)
# head(use_gene)
# getuniqScorevalue(use_gene,"gene",19) -> gene
# head(gene)

##### gene 的 control 了 mrna
###########因为 取了 sample 确保每次都一样
# # 选取覆盖 aug 的 mrna-seq reads
# head(data)
# data_m_aug <- filter(data, id == "mrna", (p5 <= 2 & p3 >= 0) )
# head(data_m_aug);length(data_m_aug$read_id) #32349
# head(data_gene)
# mrna_gene <- data_m_aug$gene_name;head(mrna_gene);length(mrna_gene) # 32349
# ### 其实 选取覆盖 aug的 基因，其实就已经只有
# # uniq_mrna_gene <- unique(mrna_gene);length(uniq_mrna_gene) # 3183
# # sample_mrna_gene <- sample(mrna_gene,5000)
# # uniq_sample_mrna_gene <- unique(sample_mrna_gene);length(uniq_sample_mrna_gene) # 886
# # 更正 选择 1500 条， 目的 variant数目与 实验的差不多
# sample_mrna_gene <- sample(mrna_gene,3000)
# sample_mrna_gene <- sample(mrna_gene,10000)
# sample_mrna_gene <- sample(mrna_gene,500)
# uniq_sample_mrna_gene <- unique(sample_mrna_gene);length(uniq_sample_mrna_gene) # 579
# sort(uniq_sample_mrna_gene) -> jj
# write.table(sample_mrna_gene,"G:/40s/sample_mrna_gene_500.csv", quote = F ,sep = ",",eol = "\n",col.names = F, row.names = F)


## control 了的 gene , 如果直接用 reads 用下面那个
#取了 随机的 3000
# sample_mrna_gene <- read_csv("G:/40s/sample_mrna_gene.csv", col_names = F )
# head(sample_mrna_gene);length(sample_mrna_gene$X1) # 3000
# # length(sample_mrna_gene$X1);length(unique(sample_mrna_gene$X1)) #886
# sample_data_gene <- filter(data_gene, gene_id %in% sample_mrna_gene$X1 )
# head(sample_data_gene);length(sample_data_gene$gene_id) # 875
# # length(unique(sample_data_gene$gene_id))
# reads_sam_gene <- str_sub(sample_data_gene$seq,start = 2, end = 20)
# head(reads_sam_gene)
# unirea_sam_gene <- unique(reads_sam_gene)
# length(unirea_sam_gene);length(unirea_sam_gene)
# getuniqScorevalue(unirea_sam_gene,"gene",19) -> gene
# head(gene)

##基因如果直接用reads 的话 用这个
# sample_mrna_gene <- read_csv("G:/40s/sample_mrna_gene_500.csv", col_names = F )
# sample_mrna_gene <- read_csv("G:/40s/sample_mrna_gene_1500.csv", col_names = F )
sample_mrna_gene <- read_csv("G:/40s/sample_mrna_gene_3000.csv", col_names = F )
# sample_mrna_gene <- read_csv("G:/40s/sample_mrna_gene_10000.csv", col_names = F )
head(sample_mrna_gene);length(sample_mrna_gene$X1)
head(d21nt);nrow(d21nt)## data_gene 即是d21nt
sample_gene <- merge(d21nt,sample_mrna_gene, by.x = "gene_id", by.y = "X1") %>% select(-"score_21nt") 
head(sample_gene);nrow(sample_gene) # 2985
##为了寻找差异的原因
# which(is.na(sample_gene))
# write.csv(sample_gene,"G:\\difference.csv")
# which(is.na(sample_gene))
# length(unique(sample_gene$seq)) # 875
sample_gene$gene_read <- str_sub(sample_gene$reads_21nt,start = 2, end = 20)
head(sample_gene)
use_gene <- select(sample_gene,-"reads_21nt")
head(use_gene)
# head(reads_sam_gene)
# unirea_sam_gene <- unique(reads_sam_gene)
# length(unirea_sam_gene);length(unirea_sam_gene)
getuniqScorevalue(use_gene,"gene",19) -> unigene
getScorevalue(use_gene,"gene",19) -> gene
gene <- cbind(gene,freq = rep(1/length(gene$id),length(gene$id)) )
head(gene);length(unique(gene$reads))
new_gene <- gene %>% group_by(id,score_value,reads,gene_id) %>% summarise(count = n())
head(new_gene)
#这句 control 了 mrna 
head(data_rpkm)
rpkm_gene <- merge(new_gene,data_rpkm,by = "gene_id")
# head(rpkm_gene)
rpkm_gene <- cbind(rpkm_gene, freq = rpkm_gene$count/log2(rpkm_gene$rpkm1))
rpkm_gene <- cbind(rpkm_gene, freq2 = rpkm_gene$count/rpkm_gene$rpkm1)
rpkm_gene$freq <- rpkm_gene$freq/sum(rpkm_gene$freq)
rpkm_gene$freq2 <- rpkm_gene$freq2/sum(rpkm_gene$freq2)
head(rpkm_gene);length(rpkm_gene$gene_id)
sum(rpkm_gene$freq)
### 这个 把 new 进行加权, 求和为1
new_gene$freq <- (new_gene$count)/(length(gene$id))
sum(new_gene$count);sum(new_gene$freq)
head(new_gene)

####
use_r28 <- filter(data, id == "rrs", p5 == (-12) , length == 28)
head(use_r28);length(use_r28$read_id);
getuniqScorevalue(use_r28,"rs28",28) -> uni28
getScorevalue(use_r28,"rs28",28) -> rs28
rs28 <- cbind(rs28,freq = rep(1/length(rs28$id),length(rs28$id)) )
head(rs28);length(unique(rs28$reads))
new_rs28 <- rs28 %>% group_by(id,score_value,reads,gene_id) %>% summarise(count = n())
#这句 control 了 mrna 
head(data_rpkm)
rpkm_rs28 <- merge(new_rs28,data_rpkm,by = "gene_id")
# head(rpkm_rs28)
rpkm_rs28 <- cbind(rpkm_rs28, freq = rpkm_rs28$count/log2(rpkm_rs28$rpkm1))
rpkm_rs28 <- cbind(rpkm_rs28, freq2 = rpkm_rs28$count/rpkm_rs28$rpkm1)
rpkm_rs28$freq <- rpkm_rs28$freq/sum(rpkm_rs28$freq)
rpkm_rs28$freq2 <- rpkm_rs28$freq2/sum(rpkm_rs28$freq2)
head(rpkm_rs28);length(rpkm_rs28$gene_id)
sum(rpkm_rs28$freq)
### 这个 把 new 进行加权, 求和为1
new_rs28$freq <- (new_rs28$count)/(length(rs28$id))
sum(new_rs28$count);sum(new_rs28$freq)
# sum(new_ssu19$count)
# sum(is.na(ssu19))
# sum(is.na(ssu37))
# sum(is.na(ram))
# sum(is.na(gene))
# sum(is.na(rs28))
# which(is.na(ssu19))
# head(ssu19)

#### 注意！！！ 不想写代码这里 我就只是 换个变量名跑而已
# syossu_score <- rbind(ssu19,ssu29, ssu37, gene,rs28,ssu47,ssu55)

### new score 是 score 的压缩  uni score 等于 newscore（除了 count 那一列）
# ssu_score <-rbind(ssu19,ssu29, ssu37, gene,rs28,ssu47,ssu55)
ssu_score <-rbind(ssu19,ssu29, ssu37, gene,rs28)
# ssu_score <- cbind(ssu_score,count = rep(1,length(ssu_score$id)))
head(ssu_score);length(ssu_score$id)
# new_score <-rbind(as.data.frame(new_ssu19),as.data.frame(new_ssu29), as.data.frame(new_ssu37), as.data.frame(new_gene)
            # ,as.data.frame(new_rs28),as.data.frame(new_ssu47),as.data.frame(new_ssu55)) 
new_score <-rbind(as.data.frame(new_ssu19),as.data.frame(new_ssu29), as.data.frame(new_ssu37), as.data.frame(new_gene)
                  ,as.data.frame(new_rs28)) 
length(new_score$id)
# uni_score <- rbind(uni19,uni29, uni37, unigene,uni28,uni47,uni55) 
uni_score <- rbind(uni19,uni29, uni37, unigene,uni28) 
length(uni_score$id)
# summary(new_ssu19$count);summary(new_ssu19$freq);
rpkm_score <-rbind(rpkm_ssu19,rpkm_ssu29, rpkm_ssu37, rpkm_gene,rpkm_rs28, rpkm_ssu47, rpkm_ssu55) 
head(rpkm_score)
head(ssu_score)

# 保存下来，下次不用再跑
write_csv(ssu_score,"G:/kozak/ssu/ssu_score.csv")
write_csv(new_score,"G:/kozak/ssu/new_score.csv")
write_csv(uni_score,"G:/kozak/ssu/uni_score.csv")
# write.csv(filter(ssu_score, id %in% c("ssu19","ssu37","gene")),"G:\\reads.csv")
# write.csv(filter(new_score, id %in% c("ssu19","ssu37","gene")),"G:\\uniq_reads.csv")
#################
ppp <- ggplot(rpkm_score)
(box <- ppp + aes(x = id , y= log2(rpkm1) , fill = id) + geom_boxplot() )
write.csv(rpkm_ssu55$gene_id,"G:\\ssu55.csv")
write.csv(rpkm_ssu37$gene_id,"G:\\ssu37.csv")
write.csv(rpkm_ssu47$gene_id,"G:\\ssu47.csv")
write.csv(rpkm_ssu29$gene_id,"G:\\ssu29.csv")
# rpkm_score <-rbind(rpkm_ssu19,rpkm_ssu29, rpkm_ssu37, ram, rpkm_gene,rpkm_rs28,ssu47,ssu55) 
?sort
# sort(new_score$reads) -> score_1
# sort(uni_score$reads) -> score_2
# head(score_1,20);head(score_2,20)

# new_score <-rbind(new_ssu19,new_ssu29, new_ssu37, new_ram, new_gene,new_rs28) 
# warnings()
# head(score);length(score$id)
# head(new_score);length(new_score$id)
# length(uni_score$id)
# new_score <- score %>% group_by(id,score_value,reads) %>% summarise(count = n())
# fi <- "G:/40s/score_reads.csv"
# write.table(score,fi, quote = F ,sep = ",",eol = "\n",col.names = T, row.names = F)
####### 画图
pp <- ggplot() + theme_bw() +theme(plot.title = element_text(hjust = 0.5, size = 30), 
  legend.position = c(0.9,0.85),legend.title = element_text(face="bold", size=10),legend.text = element_text(face = "bold", size = 23),
  axis.text.x = element_text(size = 26,face = "bold"), 
  axis.text.y = element_text(size = 26,face = "bold"), axis.title = element_text(size = 26, face = "bold"))

# geom_rug
######### 自动 clor的 不能多于 5个 否则报错
###测试
# head(nnew_score);head(uni)
# nnew_score <- new_score[,-4]

# head(ssu_score)
# ( p_score <- pp%+%(filter(ssu_score, id == c("ssu29","ssu19","ssu37","gene","rs29"))) +
#     aes(score_value, color = id,weights = count) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)  
#   + labs(title = "initiation score density"))
# lines(density(ssu19$score_value,kernel = c("gaussian")),col = 1, lwd = 3, lty = 1)
# lines(density(ssu29$score_value,kernel = c("gaussian")), col = 2, lwd = 3, lty = 2)
# plot(density(ssu37$score_value,kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
# lines(density(gene$score_value,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
# lines(density(rs28$score_value,kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
# legend("topright",c("ssu19","ssu29","ssu37","gene","rs19"), col = 1:5, lty = 1:5,lwd = 3,bty = "n")
# 
# lines(density(ssu19$score_value,kernel = c("gaussian"),weights = ssu19$count),col = 1, lwd = 3, lty = 1)
# lines(density(ssu29$score_value,kernel = c("gaussian"),weights = ssu29$count), col = 2, lwd = 3, lty = 2)
# plot(density(ssu37$score_value,kernel = c("gaussian"),weights = ssu37$count),col = 3, lwd = 3, lty = 3)
# lines(density(gene$score_value,kernel = c("gaussian"),weights = gene$count),col = 4, lwd = 3, lty = 4)
# lines(density(rs28$score_value,kernel = c("gaussian"),weights = rs28$count),col = 5, lwd = 3, lty = 5)
# # lines(density(new_ram$score_value,kernel = c("gaussian"),weights = new_ram$count),col = 6, lwd = 3, lty = 6)
# legend("topright",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n")
# 
# lines(density(ssu19$score_value,kernel = c("gaussian"),weights = ssu19$freq),col = 1, lwd = 3, lty = 1)
# lines(density(ssu29$score_value,kernel = c("gaussian"),weights = ssu29$freq), col = 2, lwd = 3, lty = 2)
# plot(density(ssu37$score_value,kernel = c("gaussian"),weights = ssu37$freq),col = 3, lwd = 3, lty = 3)
# lines(density(gene$score_value,kernel = c("gaussian"),weights = gene$freq),col = 4, lwd = 3, lty = 4)
# lines(density(rs28$score_value,kernel = c("gaussian"),weights = rs28$freq),col = 5, lwd = 3, lty = 5)
# # lines(density(new_ram$score_value,kernel = c("gaussian"),weights = new_ram$count),col = 6, lwd = 3, lty = 6)
# legend("topright",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n")
# 
# ( p_score <- pp%+%(filter(uni_score, id == c("ssu29","ssu19","ssu37","gene","ram","rs29"))) +
#     aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)  
#   + labs(title = "initiation score density") + xlim(c(0.25,0.9)))
# 
# plot(density(uni19$score_value,kernel = c("gaussian")),col = 1, lwd = 3, lty = 1)
# lines(density(uni29$score_value,kernel = c("gaussian")), col = 2, lwd = 3, lty = 2)
# lines(density(uni37$score_value,kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
# lines(density(unigene$score_value,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
# lines(density(uni28$score_value,kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
# legend("topright",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n")


# < > 19 29 37 nt的分析########## file 加载 -------------------------------------------------
ssu_score <- read_csv("G:/kozak/ssu/ssu_score.csv",col_names = T)  # 全部没有 count
new_score <- read_csv("G:/kozak/ssu/new_score.csv",col_names = T)  #已经count
uni_score <- read_csv("G:/kozak/ssu/uni_score.csv",col_names = T)  

new_ssu19 <- filter(new_score,id == "ssu19")
new_ssu29 <- filter(new_score,id == "ssu29")
new_ssu37 <- filter(new_score,id == "ssu37")
new_gene <- filter(new_score,id == "gene")

ssu19 <- filter(ssu_score,id == "ssu19")
ssu29 <- filter(ssu_score,id == "ssu29")
ssu37 <- filter(ssu_score,id == "ssu37")
gene <- filter(ssu_score,id == "gene")
head(ssu19)

# < > 都能测到 19 29 37 nt的分析  ------------------------------------------------
#### 看 19 29 37 共同基因的递进 以及gene
####试一下 merge 找 共同的 基因 是否值越来越高，加权重
head(new_ssu19);head(new_ssu29);head(new_ssu37);head(new_gene)
# head(uni19);head(uni29);head(uni37)
# new19 <- ungroup(new_ssu19) %>% select(-c("id","reads","freq")) %>% rename(score_19 = score_value,count_19 = count)
new19 <- ungroup(new_ssu19) %>% select(-c("id","reads","freq")) %>% rename(score_19 = score_value,count_19 = count)
new29 <- ungroup(new_ssu29) %>%select(-c("id","reads","freq")) %>% rename(score_29 = score_value,count_29 = count)
new37 <- ungroup(new_ssu37) %>%select(-c("id","reads","freq")) %>% rename(score_37 = score_value,count_37 = count)
newgene <- ungroup(new_gene) %>%select(-c("id","reads","freq")) %>% rename(count_gene = count)

new192937 <- merge(new19,new29,by = "gene_id") %>% merge(new37,by = "gene_id") 
head(new192937);nrow(new192937) # 57个基因
(num1 <- sum(new192937$count_19)); (num2 <- sum(new192937$count_29)); (num3 <- sum(new192937$count_37))
# 458 1265 1013 

newall <- mutate(new192937,ssu_19 = count_19/num1,ssu_29 = count_29/num2,ssu_37 = count_37/num3) %>%
  select(-c("score_29","score_37","count_19","count_29","count_37"))%>%
  rename(initiation_score = score_19)
head(newall);nrow(newall) # 57
sum(newall$ssu_19);sum(newall$ssu_29);

lines(density(newall$score_19,kernel = c("gaussian")),col = 2, lwd = 3, lty = 2)
lines(density(newall$score_29, kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
plot(density(newall$score_37,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
lines(density(new_gene$score_value, kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex =1.7)

lines(density(newall$initiation_score,kernel = c("gaussian"),weights = newall$ssu_19),col = 2, lwd = 3, lty = 2)
lines(density(newall$initiation_score, kernel = c("gaussian"),weights = newall$ssu_29),col = 3, lwd = 3, lty = 3)
plot(density(newall$initiation_score,kernel = c("gaussian"),weights = newall$ssu_37),col = 4, lwd = 3, lty = 4)
lines(density(new_gene$score_value, kernel = c("gaussian"),weights = new_gene$freq),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex =1.7)
### 试一下没有 count 的，exart 的是否也一样 调整bw   #有点不太一样， 
lines(density(ssu19$score_value,kernel = c("gaussian"),bw =0.02),col = 2, lwd = 3, lty = 2)
lines(density(ssu29$score_value, kernel = c("gaussian"),bw =0.02),col = 3, lwd = 3, lty = 3)
plot(density(ssu37$score_value,kernel = c("gaussian"),bw =0.02),col = 4, lwd = 3, lty = 4)
lines(density(gene$score_value, kernel = c("gaussian"),bw =0.02),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex =1.7)

head(newall)
newall2 <- gather(newall,key = "kind",value = "frequency",-c("initiation_score","gene_id"))
head(newall2)

(p00 <- p%+%newall2 + aes(x = initiation_score, color = factor(kind), weights = frequency) +
    geom_line(stat = "density",size = 1.2) +
    xlim(c(0.4,0.85)) 
  # + scale_fill_discrete(name="factor (kind)",breaks=c("freq_19", "freq_29", "freq_37"),labels=c("ssu_19", "ssu_29", "ssu_37"))
)


# 没有加权重，错误
plot(ecdf(newall$score_19),col="red")
plot(ecdf(ssu19$score_value),col="red")
lines(ecdf(ssu29$score_value),col="red")
lines(ecdf(ssu37$score_value),col="red")

#加上权重
#同时有 19 29 37 reads 的基因
head(newall);nrow(newall) # 57 个基因
# head(ssu19)
head(ssu_score)
cuma <- select(ssu_score,c("id","score_value","gene_id","freq")) %>%
  filter(gene_id %in% newall$gene_id)
head(cuma) # 3847
(p00 <- p%+%cuma + aes(y = score_value, fill = factor(id), x= factor(id)) 
  + geom_boxplot(notch = F, varwidth = T)
  # + 
    # xlim(c(0.4,0.85)) 
  # + scale_fill_discrete(name="factor (kind)",breaks=c("freq_19", "freq_29", "freq_37"),labels=c("ssu_19", "ssu_29", "ssu_37"))
)
wilcox.test(cuma[cuma$id == "ssu19",]$score_value,cuma[cuma$id == "ssu37",]$score_value)
wilcox.test(cuma[cuma$id == "ssu29",]$score_value,cuma[cuma$id == "ssu37",]$score_value)
wilcox.test(cuma[cuma$id == "ssu19",]$score_value,cuma[cuma$id == "ssu29",]$score_value)
wilcox.test(cuma[cuma$id == "ssu19",]$score_value,cuma[cuma$id == "gene",]$score_value)
wilcox.test(cuma[cuma$id == "ssu29",]$score_value,cuma[cuma$id == "gene",]$score_value)
wilcox.test(cuma[cuma$id == "ssu37",]$score_value,cuma[cuma$id == "gene",]$score_value)

# rm(cuma)
plot(ecdf(filter(cuma, id == "ssu19")$score_value),col=1, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(filter(cuma, id == "ssu29")$score_value),col=2, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(filter(cuma, id == "ssu37")$score_value),col=3, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(gene$score_value),col=4, lwd = 1,do.points=FALSE, verticals=TRUE)
legend("topleft",c("ssu 19","ssu 29","ssu 37","gene"), col = 1:4,lwd = 1,bty = "n",cex =1.7)
# 用不了 freq 因为freq 都一样呀，都是 1/sum 得来的

# plot(ecdf(1/filter(cuma, id == "ssu19")$score_value),col=1, lwd = 1,do.points=FALSE, verticals=TRUE)
# lines(ecdf(1/filter(cuma, id == "ssu29")$score_value),col=2, lwd = 1,do.points=FALSE, verticals=TRUE)
# lines(ecdf(1/filter(cuma, id == "ssu37")$score_value),col=3, lwd = 1,do.points=FALSE, verticals=TRUE)
# lines(ecdf(1/gene$score_value),col=4, lwd = 1,do.points=FALSE, verticals=TRUE)
# legend("topleft",c("ssu 19","ssu 29","ssu 37","gene"), col = 1:4,lwd = 1,bty = "n",cex =1.7)

# install.packages("sROC")
library(sROC)
?kCDF
head(cuma)
plot(kCDF(filter(cuma, id == "ssu19")$score_value,bw = 0.01),col=1, lwd = 1)
plot(kCDF(filter(cuma, id == "ssu29")$score_value,bw = 0.01),col=2, lwd = 1)
plot(kCDF(filter(cuma, id == "ssu37")$score_value,bw = 0.01),col=3, lwd = 1)
plot(kCDF(gene$score_value,bw = 0.01),col=4, lwd = 1)
warnings()
# 用不了 freq 因为样本不够

plot(kCDF(filter(cuma, id == "ssu19")$freq,bw = 0.005),col=1, lwd = 2)
plot(kCDF(filter(cuma, id == "ssu29")$freq,bw = 0.01),col=2, lwd = 1)

?kCDF
aaa <- kCDF(filter(cuma, id == "ssu29")$freq,bw = 0.005)
lines(aaa$x,aaa$Fhat, col=2, lwd = 2)
bbb <- kCDF(filter(cuma, id == "ssu37")$freq,bw = 0.005)
lines(bbb$x,bbb$Fhat, col=4, lwd = 2)


plot(kCDF(filter(cuma, id == "ssu37")$freq,bw = 0.01),col=3, lwd = 1)
plot(kCDF(gene$score_value,bw = 0.01),col=4, lwd = 1)

# plot(density(new_ssu192937$score_value,kernel = c("gaussian")),col = 1, lwd = 3, lty = 1)
# lines(density(new_gene$score_value,kernel = c("gaussian")),col = 2, lwd = 3, lty = 2)
# legend("topleft",c("ssu 19/29/37","gene"), col = 1:2, lty = 1:2,lwd = 3,bty = "n",cex =1.7)

# boxplot  不是这个
head(cuma)
cuma2 <- filter(cuma, id != "rs28")
(p99 <- pp%+%cuma2 + aes(x = factor(id), y = score_value) +
    geom_boxplot(notch = T))

head(newall2)
newall2 <- mutate(newall2,fc = initiation_score/frequency)

(p44 <- p%+%newall2 + aes(x = factor(kind), y= fc)
  + geom_boxplot(notch = T,varwidth = T)
  + labs(x = NULL, y= "initiation_score/frequency"))
(wilcox.test(newall2[newall2$kind == "ssu_19",]$fc,newall2[newall2$kind == "ssu_29",]$fc))
(wilcox.test(newall2[newall2$kind == "ssu_29",]$fc,newall2[newall2$kind == "ssu_37",]$fc))
(wilcox.test(newall2[newall2$kind == "ssu_19",]$fc,newall2[newall2$kind == "ssu_37",]$fc))

# (p3 <- p%+%newall2 + aes(x= log(fc),fill = factor(kind),color = factor(kind))
#   + geom_line(stat = "density",alpha = 0.4)
#   )

head(newall2)
plot(ecdf(filter(newall2, kind == "ssu_19")$frequency),col=1, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(filter(newall2, kind == "ssu_29")$frequency),col=2, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(filter(newall2, kind == "ssu_37")$frequency),col=3, lwd = 1,do.points=FALSE, verticals=TRUE)
# lines(ecdf(gene$score_value),col=4, lwd = 1,do.points=FALSE, verticals=TRUE)
legend("right",c("ssu 19","ssu 29","ssu 37"), col = 1:3,lwd = 1,bty = "n",cex =1.7)

library(sROC)
plot(kCDF(filter(newall2, kind == "ssu19")$frequency),col=1, lwd = 1)
plot(kCDF(filter(newall2, kind == "ssu29")$frequency,bw = 0.01),col=2, lwd = 1)
plot(kCDF(filter(newall2, kind == "ssu37")$frequency,bw = 0.01),col=3, lwd = 1)
# plot(kCDF(gene$score_value,bw = 0.01),col=4, lwd = 1)
warnings()



# 画钱老师那个图 -----------------------------------------------------------------
pdf("G:/p1h.pdf",width = 7,height = 7)
lines(density(newall$initiation_score,kernel = c("gaussian"),weights = newall$ssu_19,bw=0.025),col = 1, lwd = 3, lty = 1)
lines(density(newall$initiation_score, kernel = c("gaussian"),weights = newall$ssu_29,bw=0.025),col = 2, lwd = 3, lty = 1)
plot(density(newall$initiation_score,kernel = c("gaussian"),weights = newall$ssu_37,bw=0.025),col = 3, lwd = 3, lty = 1)
lines(density(new_gene$score_value, kernel = c("gaussian"),weights = new_gene$freq,bw=0.025),col = 4, lwd = 3, lty = 1)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 1:4, lty = 1,lwd = 3,bty = "n",cex =1.7)
dev.off()

d37 <- density(newall$initiation_score,kernel = c("gaussian"),weights = newall$ssu_37,bw=0.025)
d37xy <- cbind(as.data.frame(d37$x),as.data.frame(d37$y))  %>% rename(x ="d37$x",y="d37$y" )
(max37 <- max(d37xy$y)) #10.684
(mix37 <- max(d37xy[d37xy$x < 0.6,]$y)) #1.7379

d29 <- density(newall$initiation_score, kernel = c("gaussian"),weights = newall$ssu_29,bw=0.025)
d29xy <- cbind(as.data.frame(d29$x),as.data.frame(d29$y))  %>% rename(x ="d29$x",y="d29$y" )
(max29 <- max(d29xy$y)) #7.661
(mix29 <- max(d29xy[d29xy$x < 0.6,]$y)) # 1.1409

d19 <- density(newall$initiation_score,kernel = c("gaussian"),weights = newall$ssu_19,bw=0.025)
d19xy <- cbind(as.data.frame(d19$x),as.data.frame(d19$y))  %>% rename(x ="d19$x",y="d19$y" )
(max19 <- max(d19xy$y)) # 7.1177
(mix19 <- max(d19xy[d19xy$x < 0.6,]$y)) # 1.9319

# 包括29
max_value <- as.numeric(c(max19,max29,max37))
mix_value <- as.numeric(c(mix19,mix29,mix37))
不包括 29
max_value <- as.numeric(c(max19,max37))
mix_value <- as.numeric(c(mix19,mix37))
type <- c("ssu19","ssu29","ssu37")
type <- c("ssu19","ssu37")
value <- data.frame(type,max_value,mix_value) %>% as.tibble()#%>% as.data.frame() ; #colnames(value) <- name;# rownames(value) <- c("max_val","mix_val")
head(value)
value_ra <- mutate(value, ratio = max_value/mix_value)
value_ga <- gather(value,key = "subtype",value = "value",-"type")


#
head(value_ga)
factor(value_ga$subtype)
(p1a_1 <- p%+%value_ga + aes(fill = rev(factor(subtype)),y = value,x=factor(type))
  + geom_bar(stat = "identity",position = "dodge2",width = .5)
  # + ylim(c(0,10))
  + labs(x= NULL,y="peak value")
  + scale_fill_discrete(limits =c("max_value","mix_value"),labels=c("low initiation score","high initiation score"))
)

(p1a_2 <- p%+%value_ga + aes(x = rev(factor(subtype)),y = value,fill=factor(type))
  + geom_bar(stat = "identity",position = "dodge2",width = .5)
  # + ylim(c(0,10))
  + labs(x= NULL,y="peak value")
  + scale_x_discrete(limits =c("max_value","mix_value"),labels=c("low initiation score","high initiation score"))
)


head(value_ra)
(p1a_3 <- p%+%value_ra + aes(x = factor(type),y = ratio)
  + geom_bar(stat = "identity",width = .3,fill ="#F8766D")
  + ylim(c(0,10))
  + labs(x= NULL,y="peak ratio (high/low)")
)


# 保存 ----------------------------------------------------------------------
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(8,10))) ####将页面分成2*2矩阵
print(p1a_1, vp = vplayout(2:4,2:7))   ###将(2,1)的位置画图b
print(p1a_2, vp = vplayout(5:7,2:6))
print(p1a_3, vp = vplayout(5:7,7:9))

dev.off()
# 23 16



# < > 所有 19 29 37 nt的分析，并不一定是共同基因  ---------------------------------------------------
head(ssu19);head(ssu29);head(ssu37);head(gene)
1/887
# head(uni19);head(uni29);head(uni37)
new19 <- ungroup(ssu19) %>% select(-c("id","reads")) %>% mutate(id = "ssu19")
new29 <- ungroup(ssu29) %>%select(-c("id","reads")) %>% mutate(id = "ssu29")
new37 <- ungroup(ssu37) %>%select(-c("id","reads")) %>% mutate(id = "ssu37")
newgene <- ungroup(gene) %>%select(-c("id","reads")) %>% mutate(id = "gene")
head(new19);length(unique(new19$gene_id)) # 233
newall <- rbind(new19,new29,new37,newgene)
head(newall);sum(newall$freq)
head(newall)
?geom_boxplot

# 画个 boxplot 看 score ------------------------------------------------------
(p9_1a <- p%+%newall + aes(x = factor(id), y= score_value)
  + geom_boxplot(notch = T,varwidth = F, fill = "grey")
 + labs(x ="ribosome small subunit length", y ="initiation score")
 + ylim(c(min(newall$score_value),0.9))
 + scale_x_discrete(limits=c("ssu19","ssu29","ssu37","gene"))
 
  # + ylim(c(0.65,0.8))
  # + geom_jitter()
)
(p9_2a <- p%+%newall + aes(x = factor(id), y= score_value)
  + geom_boxplot(notch = T,varwidth = F,outlier.shape = NA, fill = "grey")
  + labs(x ="ribosome small subunit length", y ="initiation score")
  + ylim(c(min(newall$score_value),0.9))
  + scale_x_discrete(limits=c("ssu19","ssu29","ssu37","gene"))
  # + geom_jitter()
)
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(4,11))) ####将页面分成2*2矩阵
print(p9_1a, vp = vplayout(2:3,2:5))
print(p9_2a, vp = vplayout(2:3,7:10))
dev.off()
# 22 15 英尺

wilcox.test(new19$score_value,new37$score_value)# 0.01
wilcox.test(new29$score_value,new37$score_value) # 8.6e-05
wilcox.test(new19$score_value,new29$score_value) # 0.55
wilcox.test(new19$score_value,newgene$score_value) # 7.2e-05
wilcox.test(new29$score_value,newgene$score_value) # 9.2e-09
wilcox.test(newgene$score_value,new37$score_value) # 0.13
2^(-6)

# density plot ------------------------------------------------------------
# # 看默认的 bw 是多少
# t <- density(new37$score_37,kernel = c("gaussian"))
# str(t) #0.00805
head(new19)
pdf("G:/p1a.pdf",width = 7,height = 7)
lines(density(new19$score_value,kernel = c("gaussian"), bw = 0.025),col = 1, lwd = 2, lty = 1)
lines(density(new29$score_value, kernel = c("gaussian"), bw = 0.025),col = 2, lwd = 2, lty = 1)
plot(density(new37$score_value,kernel = c("gaussian"), bw = 0.025),col = 3, lwd = 2, lty = 1)
lines(density(newgene$score_value, kernel = c("gaussian"), bw = 0.025),col = 4, lwd = 2, lty = 1)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 1:4, lty = 1,lwd = 2,bty = "n",cex =1.7)
dev.off()

# 画钱老师那个图 -----------------------------------------------------------------
d37 <- density(new37$score_value,kernel = c("gaussian"), bw = 0.025)
d37xy <- cbind(as.data.frame(d37$x),as.data.frame(d37$y))  %>% rename(x ="d37$x",y="d37$y" )
(max37 <- max(d37xy$y)) #8.9658
(mix37 <- max(d37xy[d37xy$x < 0.6,]$y)) #1.6543

d29 <- density(new29$score_value,kernel = c("gaussian"), bw = 0.025)
d29xy <- cbind(as.data.frame(d29$x),as.data.frame(d29$y))  %>% rename(x ="d29$x",y="d29$y" )
(max29 <- max(d29xy$y)) #6.6935
(mix29 <- max(d29xy[d29xy$x < 0.6,]$y)) # 1.6129

d19 <- density(new19$score_value,kernel = c("gaussian"), bw = 0.025)
d19xy <- cbind(as.data.frame(d19$x),as.data.frame(d19$y))  %>% rename(x ="d19$x",y="d19$y" )
(max19 <- max(d19xy$y)) # 7.4871
(mix19 <- max(d19xy[d19xy$x < 0.6,]$y)) # 1.2084

# 包括19
max_value <- as.numeric(c(max19,max29,max37))
mix_value <- as.numeric(c(mix19,mix29,mix37))
不包括 19
max_value <- as.numeric(c(max29,max37))
mix_value <- as.numeric(c(mix29,mix37))
type <- c("ssu19","ssu29","ssu37")
type <- c("ssu29","ssu37")
value <- data.frame(type,max_value,mix_value) %>% as.tibble()#%>% as.data.frame() ; #colnames(value) <- name;# rownames(value) <- c("max_val","mix_val")
head(value)
value_ra <- mutate(value, ratio = max_value/mix_value)
value_ga <- gather(value,key = "subtype",value = "value",-"type")


#
head(value_ga)
factor(value_ga$subtype)
(p1a_1 <- p%+%value_ga + aes(fill = rev(factor(subtype)),y = value,x=factor(type))
  + geom_bar(stat = "identity",position = "dodge2",width = .5)
  + ylim(c(0,10))
  + labs(x= NULL,y="peak value")
  + scale_fill_discrete(limits =c("max_value","mix_value"),labels=c("low initiation score","high initiation score"))
)

(p1a_2 <- p%+%value_ga + aes(x = rev(factor(subtype)),y = value,fill=factor(type))
  + geom_bar(stat = "identity",position = "dodge2",width = .5)
  + ylim(c(0,10))
  + labs(x= NULL,y="peak value")
  + scale_x_discrete(limits =c("max_value","mix_value"),labels=c("low initiation score","high initiation score"))
)


head(value_ra)
(p1a_3 <- p%+%value_ra + aes(x = factor(type),y = ratio)
  + geom_bar(stat = "identity",width = .3,fill ="#F8766D")
  + ylim(c(0,10))
  + labs(x= NULL,y="peak ratio (high/low)")
)


# 保存 ----------------------------------------------------------------------
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(8,10))) ####将页面分成2*2矩阵
print(p1a_1, vp = vplayout(2:4,2:7))   ###将(2,1)的位置画图b
print(p1a_2, vp = vplayout(5:7,2:6))
print(p1a_3, vp = vplayout(5:7,7:9))

dev.off()
# 23 16

# 其他 ----------------------------------------------------------------------
#与上面的图一样
lines(density(new19$score_value,kernel = c("gaussian"), bw = 0.025,weights = new19$freq),col = 2, lwd = 3, lty = 2)
lines(density(new29$score_value, kernel = c("gaussian"), bw = 0.025,weights = new29$freq),col = 3, lwd = 3, lty = 3)
plot(density(new37$score_value,kernel = c("gaussian"), bw = 0.025,weights = new37$freq),col = 4, lwd = 3, lty = 4)
lines(density(newgene$score_value, kernel = c("gaussian"), bw = 0.025,weights = newgene$freq),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex =1.7)

# 插曲 试试 new_ssu 是否跟 ssu一样
head(new_ssu19);sum(new_ssu19$count)
new19 <- ungroup(new_ssu19) %>% select(-c("id","reads","count")) %>% mutate(id = "ssu19")
new29 <- ungroup(new_ssu29) %>%select(-c("id","reads","count")) %>% mutate(id = "ssu29")
new37 <- ungroup(new_ssu37) %>%select(-c("id","reads","count")) %>% mutate(id = "ssu37")
newgene <- ungroup(new_gene) %>%select(-c("id","reads","count")) %>% mutate(id = "gene")
head(new19);length(unique(new19$gene_id)) # 233
newall <- rbind(new19,new29,new37,newgene)
head(newall);sum(newall$freq)

# #与上面图不一样，因为没加权重，所以尽管加权重方式不一样，但只要数据一样，设置的 bw 一样就OK。
lines(density(new19$score_value,kernel = c("gaussian"), bw = 0.025),col = 2, lwd = 3, lty = 2)
lines(density(new29$score_value, kernel = c("gaussian"), bw = 0.025),col = 3, lwd = 3, lty = 3)
plot(density(new37$score_value,kernel = c("gaussian"), bw = 0.025),col = 4, lwd = 3, lty = 4)
lines(density(newgene$score_value, kernel = c("gaussian"), bw = 0.025),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex =1.7)
# #与上面的图一样
lines(density(new19$score_value,kernel = c("gaussian"), bw = 0.025,weights = new19$freq),col = 2, lwd = 3, lty = 2)
lines(density(new29$score_value, kernel = c("gaussian"), bw = 0.025,weights = new29$freq),col = 3, lwd = 3, lty = 3)
plot(density(new37$score_value,kernel = c("gaussian"), bw = 0.025,weights = new37$freq),col = 4, lwd = 3, lty = 4)
lines(density(newgene$score_value, kernel = c("gaussian"), bw = 0.025,weights = newgene$freq),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex =1.7)
# ##已经test完了，所以把上面的注释了，以免 newxx 搞混了

# 用不了 score 因为没有加权重
# 用不了 freq 因为都一样

## 0115 感觉用 freq 是错的

# 要用 new_ssu
head(new19)
pdf("G:/p3a.pdf",width = 7,height = 7)
plot(ecdf(new19$freq),col=1, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(new29$freq),col=2, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(new37$freq),col=3, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(newgene$freq),col=4, lwd = 1,do.points=FALSE, verticals=TRUE)
legend("topright",c("ssu 19","ssu 29","ssu 37","gene"), col = 1:4,lwd = 1,bty = "n",cex =1.7)
dev.off()

# 要用 ssu
pdf("G:/p2a.pdf",width = 7,height = 7)
plot(ecdf(new19$score_value),col=1, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(new29$score_value),col=2, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(new37$score_value),col=3, lwd = 1,do.points=FALSE, verticals=TRUE)
lines(ecdf(newgene$score_value),col=4, lwd = 1,do.points=FALSE, verticals=TRUE)
legend("topleft",c("ssu 19","ssu 29","ssu 37","gene"), col = 1:4,lwd = 1,bty = "n",cex =1.7)
dev.off()

# plot(ecdf(1/filter(cuma, id == "ssu19")$score_value),col=1, lwd = 1,do.points=FALSE, verticals=TRUE)
# lines(ecdf(1/filter(cuma, id == "ssu29")$score_value),col=2, lwd = 1,do.points=FALSE, verticals=TRUE)
# lines(ecdf(1/filter(cuma, id == "ssu37")$score_value),col=3, lwd = 1,do.points=FALSE, verticals=TRUE)
# lines(ecdf(1/gene$score_value),col=4, lwd = 1,do.points=FALSE, verticals=TRUE)
# legend("topleft",c("ssu 19","ssu 29","ssu 37","gene"), col = 1:4,lwd = 1,bty = "n",cex =1.7)

# install.packages("sROC")
library(sROC)
# 要用 new_ssu
# pdf("G:/p2a.pdf",width = 7,height = 7)
# b <- 0.015
# k <- "epanechnikov"
# # b <- 0.02
# # k <- "normal"
# a1 <- kCDF(filter(newall, id == "ssu19")$freq,bw = b,kernel = k)
# a2 <- kCDF(filter(newall, id == "ssu29")$freq,bw = b,kernel = k)
# a3 <- kCDF(filter(newall, id == "ssu37")$freq,bw = b,kernel = k)
# a4 <- kCDF(filter(newall, id == "gene")$freq,bw = b,kernel = k)
# plot(a1$x,a1$Fhat,col=1, lwd = 2,type = "l")
# lines(a2$x,a2$Fhat,col=2, lwd = 2)
# lines(a3$x,a3$Fhat,col=3, lwd = 2)
# lines(a4$x,a4$Fhat,col=4, lwd = 2)
# legend("topleft",c("ssu 19","ssu 29","ssu 37","gene"), col = 1:4,lwd = 2,bty = "n",cex =1.7)
# dev.off()

?kCDF
# 要用 ssu 
pdf("G:/p4a.pdf",width = 7,height = 7)
b <- 0.015
k <- "epanechnikov"
# b <- 0.02
# k <- "normal"
a1 <- kCDF(filter(newall, id == "ssu19")$score_value,bw = b,kernel = k)
a2 <- kCDF(filter(newall, id == "ssu29")$score_value,bw = b,kernel = k)
a3 <- kCDF(filter(newall, id == "ssu37")$score_value,bw = b,kernel = k)
a4 <- kCDF(filter(newall, id == "gene")$score_value,bw = b,kernel = k)
plot(a1$x,a1$Fhat,col=1, lwd = 2,type = "l")
lines(a2$x,a2$Fhat,col=2, lwd = 2)
lines(a3$x,a3$Fhat,col=3, lwd = 2)
lines(a4$x,a4$Fhat,col=4, lwd = 2)
legend("topleft",c("ssu 19","ssu 29","ssu 37","gene"), col = 1:4,lwd = 2,bty = "n",cex =1.7)
dev.off()


# boxplot  不是这个
# head(newall)
# (p99 <- p%+%newall + aes(x = factor(id), y = score_value) +
#     geom_boxplot(notch = T,varwidth = T))
# 
# head(newall2)
# newall2 <- mutate(newall2,fc = initiation_score/frequency)
# 
# (p44 <- p%+%newall2 + aes(x = factor(kind), y= fc)
#   + geom_boxplot(notch = T,varwidth = T)
#   + labs(x = NULL, y= "initiation_score/frequency"))
# (wilcox.test(newall2[newall2$kind == "ssu_19",]$fc,newall2[newall2$kind == "ssu_29",]$fc))
# (wilcox.test(newall2[newall2$kind == "ssu_29",]$fc,newall2[newall2$kind == "ssu_37",]$fc))
# (wilcox.test(newall2[newall2$kind == "ssu_19",]$fc,newall2[newall2$kind == "ssu_37",]$fc))


# 按照 score 分bin 用new_ssu ---想要说明，score 越高，reads的占比越大--------------------------------------------------------
head(new_ssu19);sum(new_ssu19$count)
new19 <- ungroup(new_ssu19) %>% select(-c("id","reads","count")) %>% mutate(id = "ssu19")
new29 <- ungroup(new_ssu29) %>%select(-c("id","reads","count")) %>% mutate(id = "ssu29")
new37 <- ungroup(new_ssu37) %>%select(-c("id","reads","count")) %>% mutate(id = "ssu37")
newgene <- ungroup(new_gene) %>%select(-c("id","reads","count")) %>% mutate(id = "gene")
head(new19);length(unique(new19$gene_id)) # 233
newall <- rbind(new19,new29,new37,newgene)
head(newall);sum(newall$freq)

head(newall) # 
head(new19)
new192 <- new19
new292 <- new29
new372 <- new37

bin <- 8
as.numeric(cut2(new19$score_value,g = bin)) -> new19$bin_score
as.numeric(cut2(new37$score_value,g = bin)) -> new37$bin_score
as.numeric(cut2(new29$score_value,g = bin)) -> new29$bin_score
as.numeric(cut2(newgene$score_value,g = bin)) -> newgene$bin_score
# tapply(log2(newall$ssu_19), INDEX = newall$bin_score,FUN = mean) -> mean_19
# tapply(log2(newall$ssu_29), INDEX = newall$bin_score,FUN = mean) -> mean_29
# tapply(log2(newall$ssu_37), INDEX = newall$bin_score,FUN = mean) -> mean_37
head(new19)
tapply(new19$freq, INDEX = new19$bin_score,FUN = sum) -> mean_19
tapply(new29$freq, INDEX = new29$bin_score,FUN = sum) -> mean_29
tapply(new37$freq, INDEX = new37$bin_score,FUN = sum) -> mean_37
tapply(newgene$freq, INDEX = newgene$bin_score,FUN = sum) -> mean_gene

bar1 <- tapply(new19$freq, INDEX = new19$bin_score,FUN = sd)/sqrt(
  tapply(new19$freq, INDEX = new19$bin_score,FUN = length))
bar2 <- tapply(new29$freq, INDEX = new29$bin_score,FUN = sd)/sqrt(
  tapply(new29$freq, INDEX = new29$bin_score,FUN = length))
bar3 <- tapply(new37$freq, INDEX = new37$bin_score,FUN = sd)/sqrt(
  tapply(new37$freq, INDEX = new37$bin_score,FUN = length))
bar4 <- tapply(newgene$freq, INDEX = newgene$bin_score,FUN = sd)/sqrt(
  tapply(newgene$freq, INDEX = newgene$bin_score,FUN = length))
bar1 # max =

bins <- seq(1,bin)
dbase <- cbind(as.data.frame(bins),as.data.frame(mean_19),
               as.data.frame(mean_29),as.data.frame(mean_37),as.data.frame(mean_gene)) %>%
  gather(key = "type",value = "freq",-"bins")
head(dbase)
mean_19
# dplot <- rbind(new19,new29,new37)
# head(dplot)
# ?geom_bar
(p5a <- p%+%dbase + aes(x = factor(bins), y= freq,fill = factor(type)) 
  # + geom_point(size = 5,alpha = .7)
  + geom_bar(position = "dodge",stat = "identity")
  # + theme(legend.title=element_blank()) 
  + labs(x = "initiation socor bins",y = "mean of frequency of ssu_reads ")
  + geom_errorbar(aes(ymax= c(mean_19+bar1,mean_29+bar2,mean_37+bar3,mean_gene+bar4),
                      ymin = c(mean_19-bar1,mean_29-bar2,mean_37-bar3,mean_gene-bar4)),
                  width =.3,size =1.5,position=position_dodge(.9))
)
pdf("G:/p5a.pdf",width = 9,height = 9)
print(p5a)
dev.off()

# # 尝试用 mean --------------------------------------------------------------
tapply(new19$freq, INDEX = new19$bin_score,FUN = mean) -> mean_19
tapply(new29$freq, INDEX = new29$bin_score,FUN = mean) -> mean_29
tapply(new37$freq, INDEX = new37$bin_score,FUN = mean) -> mean_37
tapply(newgene$freq, INDEX = newgene$bin_score,FUN = mean) -> mean_gene

bar1 <- tapply(new19$freq, INDEX = new19$bin_score,FUN = sd)/sqrt(
  tapply(new19$freq, INDEX = new19$bin_score,FUN = length))
bar2 <- tapply(new29$freq, INDEX = new29$bin_score,FUN = sd)/sqrt(
  tapply(new29$freq, INDEX = new29$bin_score,FUN = length))
bar3 <- tapply(new37$freq, INDEX = new37$bin_score,FUN = sd)/sqrt(
  tapply(new37$freq, INDEX = new37$bin_score,FUN = length))
bar4 <- tapply(newgene$freq, INDEX = newgene$bin_score,FUN = sd)/sqrt(
  tapply(newgene$freq, INDEX = newgene$bin_score,FUN = length))
bar1 # max =

bins <- seq(1,bin)
dbase2 <- cbind(as.data.frame(bins),as.data.frame(mean_19),
               as.data.frame(mean_29),as.data.frame(mean_37),as.data.frame(mean_gene)) %>%
  gather(key = "type",value = "freq",-"bins")
head(dbase)
mean_19
# dplot <- rbind(new19,new29,new37)
# head(dplot)
# ?geom_bar
(p6a <- p%+%dbase2 + aes(x = factor(bins), y= freq,fill = factor(type)) 
  # + geom_point(size = 5,alpha = .7)
  + geom_bar(position = "dodge",stat = "identity")
  # + theme(legend.title=element_blank()) 
  + labs(x = "initiation socor bins",y = "mean of frequency of ssu_reads ")
  + geom_errorbar(aes(ymax= c(mean_19+bar1,mean_29+bar2,mean_37+bar3,mean_gene+bar4),
                      ymin = c(mean_19-bar1,mean_29-bar2,mean_37-bar3,mean_gene-bar4)),
                  width =.2,size =1,position=position_dodge(.9))
)
pdf("G:/p6a.pdf",width = 9,height = 9)
print(p6a)
dev.off()

# 同上 尝试用 freq 分bin 用 new ssu --期望看到 freq 越高， socore 越高-不能用 sum,score 很大-----------------------------------------------
head(new_ssu19);sum(new_ssu19$count)
new19 <- ungroup(new_ssu19) %>% select(-c("id","reads","count")) %>% mutate(id = "ssu19")
new29 <- ungroup(new_ssu29) %>%select(-c("id","reads","count")) %>% mutate(id = "ssu29")
new37 <- ungroup(new_ssu37) %>%select(-c("id","reads","count")) %>% mutate(id = "ssu37")
newgene <- ungroup(new_gene) %>%select(-c("id","reads","count")) %>% mutate(id = "gene")
head(new19);length(unique(new19$gene_id)) # 233
newall <- rbind(new19,new29,new37,newgene)
head(newall);sum(newall$freq)

head(newall) # 
head(new19)
new192 <- new19
new292 <- new29
new372 <- new37

bin <- 5
as.numeric(cut2(new19$freq,g = bin)) -> new19$bin_freq
as.numeric(cut2(new37$freq,g = bin)) -> new37$bin_freq
as.numeric(cut2(new29$freq,g = bin)) -> new29$bin_freq
as.numeric(cut2(newgene$freq,g = 6)) -> newgene$bin_freq
# tapply(log2(newall$ssu_19), INDEX = newall$bin_score,FUN = mean) -> mean_19
# tapply(log2(newall$ssu_29), INDEX = newall$bin_score,FUN = mean) -> mean_29
# tapply(log2(newall$ssu_37), INDEX = newall$bin_score,FUN = mean) -> mean_37
head(new19)
tapply(new19$score_value, INDEX = new19$bin_freq,FUN = mean) -> mean_19
tapply(new29$score_value, INDEX = new29$bin_freq,FUN = mean) -> mean_29
tapply(new37$score_value, INDEX = new37$bin_freq,FUN = mean) -> mean_37
tapply(newgene$score_value, INDEX = newgene$bin_freq,FUN = mean) -> mean_gene

bar1 <- tapply(new19$score_value, INDEX = new19$bin_freq,FUN = sd)/sqrt(
  tapply(new19$score_value, INDEX = new19$bin_freq,FUN = length))
bar2 <- tapply(new29$score_value, INDEX = new29$bin_freq,FUN = sd)/sqrt(
  tapply(new29$score_value, INDEX = new29$bin_freq,FUN = length))
bar3 <- tapply(new37$score_value, INDEX = new37$bin_freq,FUN = sd)/sqrt(
  tapply(new37$score_value, INDEX = new37$bin_freq,FUN = length))
bar4 <- tapply(newgene$score_value, INDEX = newgene$bin_freq,FUN = sd)/sqrt(
  tapply(newgene$score_value, INDEX = newgene$bin_freq,FUN = length))
bar1 # max =

bins <- seq(1,bin-1)
dbase <- cbind(as.data.frame(bins),as.data.frame(mean_19),
               as.data.frame(mean_29),as.data.frame(mean_37),as.data.frame(mean_gene)) %>%
  gather(key = "type",value = "freq",-"bins")
head(dbase)
mean_19
# dplot <- rbind(new19,new29,new37)
# head(dplot)
# ?geom_bar
(p7a <- p%+%dbase + aes(x = factor(bins), y= freq,fill = factor(type)) 
  # + geom_point(size = 5,alpha = .7)
  + geom_bar(position = "dodge",stat = "identity")
  + ylim(c(0,1))
  + labs(x = "frequency of ssu_reads bins",y = "mean of initiation socor ")
  + geom_errorbar(aes(ymax= c(mean_19+bar1,mean_29+bar2,mean_37+bar3,mean_gene+bar4),
                      ymin = c(mean_19-bar1,mean_29-bar2,mean_37-bar3,mean_gene-bar4)),
                  width =.3,size =1.5,position=position_dodge(.9))
)
pdf("G:/p7a.pdf",width = 9,height = 9)
print(p7a)
dev.off()

# <> 尝试用ssu的，同样做，---但用不了 mean 和 freq 因为都一样-------------------------------------------------------------------
head(ssu_score);head(ssu19)
new19 <- ungroup(ssu19) %>% select(-c("id","reads")) %>% mutate(id = "ssu19")
new29 <- ungroup(ssu29) %>%select(-c("id","reads")) %>% mutate(id = "ssu29")
new37 <- ungroup(ssu37) %>%select(-c("id","reads")) %>% mutate(id = "ssu37")
newgene <- ungroup(gene) %>%select(-c("id","reads")) %>% mutate(id = "gene")
head(new19);length(unique(new19$gene_id)) # 233
newall <- rbind(new19,new29,new37,newgene)
head(newall);sum(newall$freq)

head(newall) # 
head(new19)
new192 <- new19
new292 <- new29
new372 <- new37

bin <- 8
as.numeric(cut2(new19$score_value,g = bin)) -> new19$bin_score
as.numeric(cut2(new37$score_value,g = bin)) -> new37$bin_score
as.numeric(cut2(new29$score_value,g = bin)) -> new29$bin_score
as.numeric(cut2(newgene$score_value,g = bin)) -> newgene$bin_score
# tapply(log2(newall$ssu_19), INDEX = newall$bin_score,FUN = mean) -> mean_19
# tapply(log2(newall$ssu_29), INDEX = newall$bin_score,FUN = mean) -> mean_29
# tapply(log2(newall$ssu_37), INDEX = newall$bin_score,FUN = mean) -> mean_37
head(new19)
tapply(new19$freq, INDEX = new19$bin_score,FUN = sum) -> mean_19
tapply(new29$freq, INDEX = new29$bin_score,FUN = sum) -> mean_29
tapply(new37$freq, INDEX = new37$bin_score,FUN = sum) -> mean_37
tapply(newgene$freq, INDEX = newgene$bin_score,FUN = sum) -> mean_gene

bar1 <- tapply(new19$freq, INDEX = new19$bin_score,FUN = sd)/sqrt(
  tapply(new19$freq, INDEX = new19$bin_score,FUN = length))
bar2 <- tapply(new29$freq, INDEX = new29$bin_score,FUN = sd)/sqrt(
  tapply(new29$freq, INDEX = new29$bin_score,FUN = length))
bar3 <- tapply(new37$freq, INDEX = new37$bin_score,FUN = sd)/sqrt(
  tapply(new37$freq, INDEX = new37$bin_score,FUN = length))
bar4 <- tapply(newgene$freq, INDEX = newgene$bin_score,FUN = sd)/sqrt(
  tapply(newgene$freq, INDEX = newgene$bin_score,FUN = length))
bar1 # max =

bins <- seq(1,bin)
dbase <- cbind(as.data.frame(bins),as.data.frame(mean_19),
               as.data.frame(mean_29),as.data.frame(mean_37),as.data.frame(mean_gene)) %>%
  gather(key = "type",value = "freq",-"bins")
head(dbase)
mean_19
# dplot <- rbind(new19,new29,new37)
# head(dplot)
# ?geom_bar
(p5a <- p%+%dbase + aes(x = factor(bins), y= freq,fill = factor(type)) 
  # + geom_point(size = 5,alpha = .7)
  + geom_bar(position = "dodge",stat = "identity")
  # + theme(legend.title=element_blank()) 
  + labs(x = "initiation socor bins",y = "sum of frequency of ssu_reads ")
  + geom_errorbar(aes(ymax= c(mean_19+bar1,mean_29+bar2,mean_37+bar3,mean_gene+bar4),
                      ymin = c(mean_19-bar1,mean_29-bar2,mean_37-bar3,mean_gene-bar4)),
                  width =.3,size =1.5,position=position_dodge(.9))
)
pdf("G:/p5a.pdf",width = 9,height = 9)
print(p5a)
dev.off()

# ## 为了说明分更高的，更多进入 37，即37的freq更高，reads 占比更大 ## 想了想 用sum freq 是合理的 ---------

# ## 尝试 把 19 29 37 所有score 分bin，用 ssu 再分别算 ? --------------------------------------
head(ssu19);head(ssu29);head(ssu37);head(gene)
1/887
# head(uni19);head(uni29);head(uni37)
new19 <- ungroup(ssu19) %>% select(-c("id","reads")) %>% mutate(id = "ssu19")
new29 <- ungroup(ssu29) %>%select(-c("id","reads")) %>% mutate(id = "ssu29")
new37 <- ungroup(ssu37) %>%select(-c("id","reads")) %>% mutate(id = "ssu37")
newgene <- ungroup(gene) %>%select(-c("id","reads")) %>% mutate(id = "gene")
head(new19);length(unique(new19$gene_id)) # 233
newall <- rbind(new19,new29,new37,newgene)
head(newall);sum(newall$freq)
head(newall) # 
head(new19)
new192 <- new19
new292 <- new29
new372 <- new37

bin <- 8
head(newall)
as.numeric(cut2(newall$score_value,g = bin)) -> newall$bin_score
# tapply(log2(newall$ssu_19), INDEX = newall$bin_score,FUN = mean) -> mean_19
# tapply(log2(newall$ssu_29), INDEX = newall$bin_score,FUN = mean) -> mean_29
# tapply(log2(newall$ssu_37), INDEX = newall$bin_score,FUN = mean) -> mean_37
head(newall)
tapply(newall[newall$id == "ssu19",]$freq, INDEX = newall[newall$id == "ssu19",]$bin_score,FUN = sum) -> mean_19
tapply(newall[newall$id == "ssu37",]$freq, INDEX = newall[newall$id == "ssu37",]$bin_score,FUN = sum) -> mean_37
tapply(newall[newall$id == "ssu29",]$freq, INDEX = newall[newall$id == "ssu29",]$bin_score,FUN = sum) -> mean_29
tapply(newall[newall$id == "gene",]$freq, INDEX = newall[newall$id == "gene",]$bin_score,FUN = sum) -> mean_gene
## ssu 数据只能用 sum 不能用mean

bar1 <- tapply(newall[newall$id == "ssu19",]$freq, INDEX = newall[newall$id == "ssu19",]$bin_score,FUN = sd)/sqrt(
  tapply(newall[newall$id == "ssu19",]$freq, INDEX = newall[newall$id == "ssu19",]$bin_score,FUN = length))
bar2 <- tapply(newall[newall$id == "ssu29",]$freq, INDEX = newall[newall$id == "ssu29",]$bin_score,FUN = sd)/sqrt(
  tapply(newall[newall$id == "ssu29",]$freq, INDEX = newall[newall$id == "ssu29",]$bin_score,FUN = length))
bar3 <- tapply(newall[newall$id == "ssu37",]$freq, INDEX = newall[newall$id == "ssu37",]$bin_score,FUN = sd)/sqrt(
  tapply(newall[newall$id == "ssu37",]$freq, INDEX = newall[newall$id == "ssu37",]$bin_score,FUN = length))
bar4 <- tapply(newall[newall$id == "gene",]$freq, INDEX = newall[newall$id == "gene",]$bin_score,FUN = sd)/sqrt(
  tapply(newall[newall$id == "gene",]$freq, INDEX = newall[newall$id == "gene",]$bin_score,FUN = length))


bar1 # max =
bar2
bins <- seq(1,bin)
dbase <- cbind(as.data.frame(bins),as.data.frame(mean_19),
               as.data.frame(mean_29),as.data.frame(mean_37),as.data.frame(mean_gene)) %>%
  gather(key = "type",value = "freq",-"bins")
head(dbase)
mean_19
# dplot <- rbind(new19,new29,new37)
# head(dplot)

(p8a <- p%+%dbase + aes(x = factor(bins), y= freq,fill = factor(type)) 
  # + geom_point(size = 5,alpha = .7)
  + geom_bar(position = "dodge",stat = "identity")
  # + theme(legend.title=element_blank()) 
  + labs(x = "initiation socor bins(all)",y = "sum of frequency of ssu_reads ")
  + geom_errorbar(aes(ymax= c(mean_19+bar1,mean_29+bar2,mean_37+bar3,mean_gene+bar4),
                      ymin = c(mean_19-bar1,mean_29-bar2,mean_37-bar3,mean_gene-bar4)),
                  width =.3,size =1.5,position=position_dodge(.9))
)
pdf("G:/p8a.pdf",width = 9,height = 9)
print(p8a)
dev.off()

wilcox.test(mean_37[5],mean_29[5])
# 以前的 ---------------------------------------------------------------------
rm(dbin)
initiation_score_bins <- seq(1,10)
dbin <- cbind(as.data.frame(initiation_score_bins),mean_19,mean_29,mean_37) %>%
  gather(key = "kind", value = "score_frequency",-"initiation_score_bins")
head(dbin)
(p3 <- p%+%dbin + aes(x =initiation_score_bins ,y= score_frequency,color= factor(kind) )
  + geom_point(size = 2.5) 
  + geom_errorbar(aes(ymax= c(mean_19+bar1,mean_29+bar2,mean_37+bar3),ymin = c(mean_19-bar1,mean_29-bar2,mean_37-bar3)),
                              width =.1,size =.5)
)

(p3 <- p%+%dbin + aes(x =initiation_score_bins ,y= log(score_weight),color= factor(kind) )
  + geom_line(size = 2.5) 
    )
?kCDF

# 体现 热图
head(newall)
newall2 <- select(newall,-c("score_29","score_37","count_19","count_29","count_37")) %>%
  rename(initiation_score = score_19) %>%
  gather(key = "kind",value = "frequency",-c("initiation_score","gene_id"))
head(newall2)
(p999 <- pp%+%newall2 + aes(x = initiation_score, y = frequency) +
    geom_point() + facet_wrap(~kind))

(p999 <- pp%+%newall2 + aes(x = initiation_score, y = frequency,fill = frequency) +
    geom_tile() + facet_wrap(~kind,nrow = 3))
# geom_point 用于离散

??cdf
cuma2 <- filter(cuma,id != "rs28")
(p99 <- pp%+%cuma2 + aes(score_value, y = ..y..,color = factor(id)) +
    geom_step(stat = "ecdf"))

# lines(ecdf(newall$))
# a <- ecdf(newall$score_19)
# str(a)
# ()
# ?ecdf

new_ssu192937 <-  merge(new_ssu19[,1:5],new_ssu29[,1:5], by = c("gene_id","reads","score_value"))
head(new_ssu192937);nrow(new_ssu192937) # 114
new_ssu192937 <- merge(new_ssu192937, new_ssu37[,1:5],by = c("gene_id","reads","score_value"))
head(new_ssu192937);nrow(new_ssu192937) # 57
(s19 <- sum(new_ssu192937$count.x));(s29 <- sum(new_ssu192937$count.y));(s37 <- sum(new_ssu192937$count))
# 458 1265 1013 
new_ssu192937 <- new_ssu192937 %>% mutate(freq19 = count.x/s19) %>% mutate(freq29 = count.y/s29) %>% mutate(freq37 = count/s37)
head(new_ssu192937);nrow(new_ssu192937) # 57

nrow(new_gene);sum(new_gene$count)
nrow(new_ssu192937);
plot(density(new_ssu192937$score_value,kernel = c("gaussian")),col = 1, lwd = 3, lty = 1)
lines(density(new_gene$score_value,kernel = c("gaussian")),col = 2, lwd = 3, lty = 2)
legend("topleft",c("ssu 19/29/37","gene"), col = 1:2, lty = 1:2,lwd = 3,bty = "n",cex =1.7)
head(new_ssu192937)

lines(density(new_ssu192937$score_value,kernel = c("gaussian")),col = 2, lwd = 3, lty = 2)
lines(density(new_ssu192937$score_value, kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
plot(density(new_ssu192937$score_value,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
lines(density(new_gene$score_value, kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex =1.7)

lines(density(new_ssu192937$score_value,weights = new_ssu192937$freq19,kernel = c("gaussian")),col = 2, lwd = 3, lty = 2)
lines(density(new_ssu192937$score_value,weights = new_ssu192937$freq29, kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
plot(density(new_ssu192937$score_value,weights = new_ssu192937$freq37,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
lines(density(new_gene$score_value,weights = new_gene$freq, kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex =1.7)

lines(density(new_ssu192937$score_value,weights = new_ssu192937$freq19,kernel = c("gaussian"),bw = 0.015),col = 2, lwd = 3, lty = 2)
lines(density(new_ssu192937$score_value,weights = new_ssu192937$freq29, kernel = c("gaussian"),bw = 0.015),col = 3, lwd = 3, lty = 3)
plot(density(new_ssu192937$score_value,weights = new_ssu192937$freq37,kernel = c("gaussian"),bw = 0.015),col = 4, lwd = 3, lty = 4)
lines(density(new_gene$score_value,weights = new_gene$freq, kernel = c("gaussian"),bw = 0.015),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19","ssu_29","ssu_37","gene"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex =1.7)


new_ssu19_37 <- select(new_ssu192937,c("freq19","freq29","freq37"))
d <- c("freq19","freq29","freq37")
new_ssu19_37 <- select(new_ssu192937,c("count.x","count.y","count"))
d <- c("count.x","count.y","count")
for (i in d){
  new_ssu19_37[,i] <- log(new_ssu19_37[,i])
}

value <- corr.test(new_ssu19_37, method = "s")
corrplot(value$r, method = "number",type="lower" )
corrplot(value$p, method = "number",type="lower" )
library(PerformanceAnalytics)#加载包
chart.Correlation(new_ssu19_37, histogram =T, pch=19, method = "spearman")

##除以 rpkm_ssu
head(new_ssu192937);nrow(new_ssu192937) # 57
file_rpkm <- "G:/40s/mrna_gene_rpkm.csv"
data_rpkm <- read_csv(file_rpkm, col_names = T)
head(data_rpkm)
new_new <- merge(new_ssu192937,data_rpkm, by = "gene_id") %>% select(gene_id,reads,score_value,count.x,
                                                                     count.y,count,rpkm_ssu)
head(new_new)
new_new <- new_new %>% mutate(freq19 = count.x/rpkm_ssu) %>% mutate(freq29 = count.y/rpkm_ssu
) %>% mutate(freq37 = count/rpkm_ssu)
head(new_new)
(s19 <- sum(new_new$freq19));(s29 <- sum(new_new$freq29));(s37 <- sum(new_new$freq37))
# 458 1265 1013 
new_new <- new_new %>% mutate(freq19 = freq19/s19) %>% mutate(freq29 = freq29/s29) %>% mutate(freq37 = freq37/s37)
head(new_new)

lines(density(new_new$score_value,weights = new_new$freq19,kernel = c("gaussian")),col = 2, lwd = 3, lty = 2)
lines(density(new_new$score_value,weights = new_new$freq29, kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
plot(density(new_new$score_value,weights = new_new$freq37,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
legend("topleft",c("ssu_19","ssu_29","ssu_37"), col = 2:4, lty = 2:4,lwd = 3,bty = "n")



#### 看看 把 19 29 37 合在一起的reads
head(new_score)
new_3_score <- filter(new_score, id %in% c("ssu19","ssu29","ssu37")) 
head(new_3_score);nrow(new_3_score)  # 984
(sum11 <- sum(new_3_score$count)) # 5430
new_3_score <- mutate(new_3_score, freq3 = count/sum11)
(sum(new_3_score$freq3))
# plot(density(new_3_score$score_value,kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
# lines(density(new_gene$score_value,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
# lines(density(new_rs28$score_value,kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
# legend("topleft",c("ssu_19+29+37","gene","rs28"), col = 3:5, lty = 3:5,lwd = 3,bty = "n")

#### 看看 把 19 29 37 47 55合在一起的reads
head(new_score)
new_5_score <- filter(new_score, id %in% c("ssu19","ssu29","ssu37","ssu47","ssu55")) 
head(new_5_score);nrow(new_5_score)  # 1144
(sum111 <- sum(new_5_score$count)) # 5970
new_5_score <- mutate(new_5_score, freq5 = count/sum111)
(sum(new_5_score$freq5))

plot(density(new_5_score$score_value,kernel = c("gaussian")),col = 2, lwd = 3, lty = 2)
lines(density(new_3_score$score_value,kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
lines(density(new_gene$score_value,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
lines(density(new_rs28$score_value,kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19+29+37+47+55","ssu_19+29+37","gene","rs28"), col = 2:5, lty = 2:5,lwd = 5,bty = "n",cex = 1.7)
## 加上权重
head(new_3_score);nrow(new_3_score);sum(new_3_score$count) # 984 5430
head(new_5_score);nrow(new_5_score);sum(new_5_score$count)  # 1144 5970
plot(density(new_5_score$score_value,weights = new_5_score$freq5,kernel = c("gaussian")),col = 2, lwd = 3, lty = 2)
lines(density(new_3_score$score_value,weights = new_3_score$freq3, kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
lines(density(new_gene$score_value,weights = new_gene$freq,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
lines(density(new_rs28$score_value,weights = new_rs28$freq,kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19+29+37+47+55","ssu_19+29+37","gene","rs28"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex=1.7)


plot(density(new_5_score$score_value,kernel = c("gaussian"),bw= 0.017),col = 2, lwd = 3, lty = 2)
lines(density(new_3_score$score_value,kernel = c("gaussian"),bw= 0.017),col = 3, lwd = 3, lty = 3)
lines(density(new_gene$score_value,kernel = c("gaussian"),bw= 0.017),col = 4, lwd = 3, lty = 4)
lines(density(new_rs28$score_value,kernel = c("gaussian"),bw= 0.017),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19+29+37+47+55","ssu_19+29+37","gene","rs28"), col = 2:5, lty = 2:5,lwd = 5,bty = "n",cex = 1.7)
## 加上权重
head(new_3_score);nrow(new_3_score);sum(new_3_score$count) # 984 5430
head(new_5_score);nrow(new_5_score);sum(new_5_score$count)  # 1144 5970
plot(density(new_5_score$score_value,weights = new_5_score$freq5,kernel = c("gaussian"),bw= 0.015),col = 2, lwd = 3, lty = 2)
lines(density(new_3_score$score_value,weights = new_3_score$freq3, kernel = c("gaussian"),bw= 0.015),col = 3, lwd = 3, lty = 3)
lines(density(new_gene$score_value,weights = new_gene$freq,kernel = c("gaussian"),bw= 0.015),col = 4, lwd = 3, lty = 4)
lines(density(new_rs28$score_value,weights = new_rs28$freq,kernel = c("gaussian"),bw= 0.015),col = 5, lwd = 3, lty = 5)
legend("topleft",c("ssu_19+29+37+47+55","ssu_19+29+37","gene","rs28"), col = 2:5, lty = 2:5,lwd = 3,bty = "n",cex=1.7)

#
head(new_score)
( p_score <- pp%+%(filter(new_score, id == c("ssu29","ssu19","ssu37","gene","ram","rs29"))) +
    aes(score_value, color = id, weights = count) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)  
  + labs(title = "initiation score density")+ xlim(c(0.25,0.9)))
( p_score <- pp%+%(filter(new_score, id == c("ssu29","ssu19","ssu37","gene","ram","rs29"))) +
    aes(score_value, color = id, weights = freq) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)  
  + labs(title = "initiation score density")+ xlim(c(0.25,0.9)))
plot(density(new_ssu19$score_value,kernel = c("gaussian")),col = 1, lwd = 3, lty = 1)
lines(density(new_ssu29$score_value,kernel = c("gaussian")), col = 2, lwd = 3, lty = 2)
lines(density(new_ssu37$score_value,kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
lines(density(new_gene$score_value,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
lines(density(new_rs28$score_value,kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
# lines(density(new_ram$score_value,kernel = c("gaussian")),col = 6, lwd = 3, lty = 6)
# legend("topright",c("ssu19","ssu29","ssu37","gene","rs19","ram"), col = 1:6, lty = 1:6,lwd = 3,bty = "n")
legend("topleft",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n",cex =1.7)

# lines(density(new_ssu19$score_value,kernel = c("gaussian"),weights = new_ssu19$count),col = 1, lwd = 3, lty = 1)
# plot(density(new_ssu29$score_value,kernel = c("gaussian"),weights = new_ssu29$count), col = 2, lwd = 3, lty = 2)
# lines(density(new_ssu37$score_value,kernel = c("gaussian"),weights = new_ssu37$count),col = 3, lwd = 3, lty = 3)
# lines(density(new_gene$score_value,kernel = c("gaussian"),weights = new_gene$count),col = 4, lwd = 3, lty = 4)
# lines(density(new_rs28$score_value,kernel = c("gaussian"),weights = new_rs28$count),col = 5, lwd = 3, lty = 5)
# # lines(density(new_ram$score_value,kernel = c("gaussian"),weights = new_ram$count),col = 6, lwd = 3, lty = 6)
# legend("topright",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n")

lines(density(new_ssu19$score_value,kernel = c("gaussian"),weights = new_ssu19$freq),col = 1, lwd = 3, lty = 1)
lines(density(new_ssu29$score_value,kernel = c("gaussian"),weights = new_ssu29$freq), col = 2, lwd = 3, lty = 2)
plot(density(new_ssu37$score_value,kernel = c("gaussian"),weights = new_ssu37$freq),col = 3, lwd = 3, lty = 3)
lines(density(new_gene$score_value,kernel = c("gaussian"),weights = new_gene$freq),col = 4, lwd = 3, lty = 4)
lines(density(new_rs28$score_value,kernel = c("gaussian"),weights = new_rs28$freq),col = 5, lwd = 3, lty = 5)
# lines(density(new_ram$score_value,kernel = c("gaussian"),weights = new_ram$freq),col = 6, lwd = 3, lty = 6)
legend("topleft",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n",cex =1.7)

### 加了 bw
lines(density(new_ssu19$score_value,kernel = c("gaussian"),weights = new_ssu19$freq, bw = 0.017),col = 1, lwd = 3, lty = 1)
lines(density(new_ssu29$score_value,kernel = c("gaussian"),weights = new_ssu29$freq, bw = 0.017), col = 2, lwd = 3, lty = 2)
plot(density(new_ssu37$score_value,kernel = c("gaussian"),weights = new_ssu37$freq, bw = 0.017),col = 3, lwd = 3, lty = 3)
lines(density(new_gene$score_value,kernel = c("gaussian"),weights = new_gene$freq, bw = 0.017),col = 4, lwd = 3, lty = 4)
lines(density(new_rs28$score_value,kernel = c("gaussian"),weights = new_rs28$freq, bw = 0.017),col = 5, lwd = 3, lty = 5)
# lines(density(new_ram$score_value,kernel = c("gaussian"),weights = new_ram$freq),col = 6, lwd = 3, lty = 6)
legend("topleft",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n",cex =1.7)

head(rpkm_score)
# write.csv(filter(rpkm_score, log(rpkm1) <= 0),"G:\\rpkm 小于1.csv")
lines(density((filter(rpkm_ssu19,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_ssu19,freq >0))$freq ),col = 1, lwd = 3, lty = 1)
lines(density((filter(rpkm_ssu29,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_ssu29,freq >0))$freq ), col = 2, lwd = 3, lty = 2)
lines(density((filter(rpkm_ssu37,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_ssu37,freq >0))$freq ),col = 3, lwd = 3, lty = 3)
lines(density((filter(rpkm_gene,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_gene,freq >0))$freq ),col = 4, lwd = 3, lty = 4)
plot(density((filter(rpkm_rs28,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_rs28,freq >0))$freq ),col = 5, lwd = 3, lty = 5)
# lines(density(rpkm_ram$score_value,kernel = c("gaussian"),weights = rpkm_ram$freq),col = 6, lwd = 3, lty = 6)
legend("topright",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n")

lines(density(rpkm_ssu19$score_value,kernel = c("gaussian"),weights = rpkm_ssu19$freq),col = 1, lwd = 3, lty = 1)
lines(density(rpkm_ssu29$score_value,kernel = c("gaussian"),weights = rpkm_ssu29$freq), col = 2, lwd = 3, lty = 2)
plot(density(rpkm_ssu37$score_value,kernel = c("gaussian"),weights = rpkm_ssu37$freq),col = 3, lwd = 3, lty = 3)
lines(density(rpkm_gene$score_value,kernel = c("gaussian"),weights = rpkm_gene$freq),col = 4, lwd = 3, lty = 4)
lines(density(rpkm_rs28$score_value,kernel = c("gaussian"),weights = rpkm_rs28$freq),col = 5, lwd = 3, lty = 5)
# lines(density(rpkm_ram$score_value,kernel = c("gaussian"),weights = rpkm_ram$freq),col = 6, lwd = 3, lty = 6)
legend("topright",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n")

lines(density(rpkm_ssu19$score_value,kernel = c("gaussian"),weights = rpkm_ssu19$freq2),col = 1, lwd = 3, lty = 1)
lines(density(rpkm_ssu29$score_value,kernel = c("gaussian"),weights = rpkm_ssu29$freq2), col = 2, lwd = 3, lty = 2)
lines(density(rpkm_ssu37$score_value,kernel = c("gaussian"),weights = rpkm_ssu37$freq2),col = 3, lwd = 3, lty = 3)
lines(density(rpkm_gene$score_value,kernel = c("gaussian"),weights = rpkm_gene$freq2),col = 4, lwd = 3, lty = 4)
plot(density(rpkm_rs28$score_value,kernel = c("gaussian"),weights = rpkm_rs28$freq2),col = 5, lwd = 3, lty = 5)
# lines(density(rpkm_ram$score_value,kernel = c("gaussian"),weights = rpkm_ram$freq),col = 6, lwd = 3, lty = 6)
legend("topright",c("ssu19","ssu29","ssu37","gene","rs28"), col = 1:5, lty = 1:5,lwd = 3,bty = "n")

# 看 27 47 37 55 的区别 ###
( p_score <- pp%+%(filter(ssu_score, id == c("ssu29","ssu47","ssu37","gene","ssu55","ssu19"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)  + labs(title = "initiation score density"))

lines(density(ssu19$score_value,kernel = c("gaussian")),col = 1, lwd = 3, lty = 1)
lines(density(ssu29$score_value,kernel = c("gaussian")), col = 2, lwd = 3, lty = 2)
lines(density(ssu37$score_value,kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
lines(density(gene$score_value,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
lines(density(ssu47$score_value,kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
plot(density(ssu55$score_value,kernel = c("gaussian")),col = 6, lwd = 3, lty = 6)
legend("topleft",c("ssu19","ssu29","ssu37","gene","ssu47","ssu55"), col = 1:6, lty = 1:6,lwd = 3,bty = "n",cex=1.5)

( p_score <- pp%+%(filter(uni_score, id == c("ssu29","ssu47","ssu37","gene","ssu55","ssu19"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)  
  + labs(title = "initiation score density") + theme(legend.position = c(0.4,0.8) ) + xlim (c(0.3,0.9)))

lines(density(uni19$score_value,kernel = c("gaussian")),col = 1, lwd = 3, lty = 1)
lines(density(uni29$score_value,kernel = c("gaussian")), col = 2, lwd = 3, lty = 2)
lines(density(uni37$score_value,kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
lines(density(unigene$score_value,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
plot(density(uni47$score_value,kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
lines(density(uni55$score_value,kernel = c("gaussian")),col = 6, lwd = 3, lty = 6)
legend("topleft",c("ssu19","ssu29","ssu37","gene","uni47","uni55"), col = 1:6, lty = 1:6,lwd = 3,bty = "n",cex =1.5)

( p_score <- pp%+%(filter(new_score,id == c("ssu29","ssu47","ssu37","gene","ssu55","ssu19"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)  
  + labs(title = "initiation score density") + theme(legend.position = c(0.4,0.8))+ xlim (c(0.3,0.9)))

lines(density(new_ssu19$score_value,kernel = c("gaussian")),col = 1, lwd = 3, lty = 1)
lines(density(new_ssu29$score_value,kernel = c("gaussian")), col = 2, lwd = 3, lty = 2)
lines(density(new_ssu37$score_value,kernel = c("gaussian")),col = 3, lwd = 3, lty = 3)
lines(density(new_gene$score_value,kernel = c("gaussian")),col = 4, lwd = 3, lty = 4)
plot(density(new_ssu47$score_value,kernel = c("gaussian")),col = 5, lwd = 3, lty = 5)
lines(density(new_ssu55$score_value,kernel = c("gaussian")),col = 6, lwd = 3, lty = 6)
legend("topleft",c("ssu19","ssu29","ssu37","gene","ssu47","ssu55"), col = 1:6, lty = 1:6,lwd = 3,bty = "n", cex = 1.7)

# lines(density(new_ssu19$score_value,kernel = c("gaussian"),weights = new_ssu19$count),col = 1, lwd = 3, lty = 1)
# plot(density(new_ssu29$score_value,kernel = c("gaussian"),weights = new_ssu29$count), col = 2, lwd = 3, lty = 2)
# lines(density(new_ssu37$score_value,kernel = c("gaussian"),weights = new_ssu37$count),col = 3, lwd = 3, lty = 3)
# lines(density(new_gene$score_value,kernel = c("gaussian"),weights = new_gene$count),col = 4, lwd = 3, lty = 4)
# lines(density(new_ssu47$score_value,kernel = c("gaussian"),weights = new_ssu47$count),col = 5, lwd = 3, lty = 5)
# plot(density(new_ssu55$score_value,kernel = c("gaussian"),weights = new_ssu55$count),col = 6, lwd = 3, lty = 6)
# legend("topright",c("ssu19","ssu29","ssu37","gene","uni47","uni55"), col = 1:6, lty = 1:6,lwd = 3,bty = "n")

lines(density(new_ssu19$score_value,kernel = c("gaussian"),weights = new_ssu19$freq),col = 1, lwd = 3, lty = 1)
lines(density(new_ssu29$score_value,kernel = c("gaussian"),weights = new_ssu29$freq), col = 2, lwd = 3, lty = 2)
lines(density(new_ssu37$score_value,kernel = c("gaussian"),weights = new_ssu37$freq),col = 3, lwd = 3, lty = 3)
lines(density(new_gene$score_value,kernel = c("gaussian"),weights = new_gene$freq),col = 4, lwd = 3, lty = 4)
lines(density(new_ssu47$score_value,kernel = c("gaussian"),weights = new_ssu47$freq),col = 5, lwd = 3, lty = 5)
plot(density(new_ssu55$score_value,kernel = c("gaussian"),weights = new_ssu55$freq),col = 6, lwd = 3, lty = 6)
legend("topleft",c("ssu19","ssu29","ssu37","gene","ssu47","ssu55"), col = 1:6, lty = 1:6,lwd = 3,bty = "n", cex = 1.7)

head(rpkm_score)

lines(density(rpkm_ssu19$score_value,kernel = c("gaussian"),weights = rpkm_ssu19$freq),col = 1, lwd = 3, lty = 1)
lines(density(rpkm_ssu29$score_value,kernel = c("gaussian"),weights = rpkm_ssu29$freq), col = 2, lwd = 3, lty = 2)
lines(density(rpkm_ssu37$score_value,kernel = c("gaussian"),weights = rpkm_ssu37$freq),col = 3, lwd = 3, lty = 3)
lines(density(rpkm_gene$score_value,kernel = c("gaussian"),weights = rpkm_gene$freq),col = 4, lwd = 3, lty = 4)
lines(density(rpkm_ssu47$score_value,kernel = c("gaussian"),weights = rpkm_ssu47$freq),col = 5, lwd = 3, lty = 5)
plot(density(rpkm_ssu55$score_value,kernel = c("gaussian"),weights = rpkm_ssu55$freq),col = 6, lwd = 3, lty = 6)
legend("topright",c("ssu19","ssu29","ssu37","gene","uni47","uni55"), col = 1:6, lty = 1:6,lwd = 3,bty = "n")

lines(density((filter(rpkm_ssu19,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_ssu19,freq >0))$freq),col = 1, lwd = 3, lty = 1)
lines(density((filter(rpkm_ssu29,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_ssu29,freq >0))$freq), col = 2, lwd = 3, lty = 2)
lines(density((filter(rpkm_ssu37,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_ssu37,freq >0))$freq),col = 3, lwd = 3, lty = 3)
lines(density((filter(rpkm_gene,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_gene,freq >0))$freq),col = 4, lwd = 3, lty = 4)
plot(density((filter(rpkm_ssu47,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_ssu47,freq >0))$freq),col = 5, lwd = 3, lty = 5)
lines(density((filter(rpkm_ssu55,freq >0))$score_value,kernel = c("gaussian"),weights = (filter(rpkm_ssu55,freq >0))$freq),col = 6, lwd = 3, lty = 6)
legend("topright",c("ssu19","ssu29","ssu37","gene","uni47","uni55"), col = 1:6, lty = 1:6,lwd = 3,bty = "n")

( p_score <- pp%+%(filter(new_score, id == c("ssu29","ssu19","ssu37","gene","ram","rs29"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)  + labs(title = "initiation score density"))
( p_score <- pp%+%(filter(new_score, id == c("ssu29","ssu19","ssu37","gene","ram","rs29"))) +
    aes(score_value, color = id, group = id, weights = count) + geom_density(kernel = c("gaussian"),size = 0.9)  + labs(title = "initiation score density"))

# filter(new_score,id == "rs28") -> hh
# sum(hh$count)
# ?write.csv
# write.csv(new_score,"G:\\new.csv")
# write.csv(uni_score,"G:\\uni.csv")
##
( p_score <- pp%+%(filter(score, id == c("ssu37","ssu19"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9) + xlim(c(0.2,0.9)) + labs(title = "initiation score density"))

( p_score <- pp%+%(filter(score, id == c("ssu29","ssu19","ssu37","gene","rs28"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9) + xlim(c(0.2,0.9)) + labs(title = "initiation score density"))

( p_score <- pp%+%(filter(new_score, id == c("ssu29","gene","ssu19","ssu37","rs28"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9) 
  + xlim(c(0.2,0.9)) + labs(title = "initiation score density"))

( p_score <- pp%+%(filter(new_score, id == c("gene","ssu19"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9) 
  + xlim(c(0.2,0.9)) + labs(title = "initiation score density"))

( p_score <- pp%+%(filter(score, id == c("gene","ssu37"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)+ labs(title = "initiation score density"))

( p_score <- pp%+%(filter(score, id == c("rs19","gene","ssu29"))) +
    aes(score_value, color = id) + geom_line(stat="density",kernel = c("gaussian"),size = 0.9)+ labs(title = "initiation score density"))

# head(ssu19)
# , xlab = "ssu19:800 ssu29:501 ssu37:324 /n ram:2000 gene:5715 rs19:234"
?density

#### 画 counts
filescorecount <- "G:/40s/scorecounts.csv"
counts <- read_delim(filescorecount,col_names = T, delim = ",")
head(counts)
ppp <- ggplot(counts) + theme_bw()  +theme(plot.title = element_text(hjust = 0.5, size = 30), 
                                           legend.position = c(0.9,0.85),legend.title = element_text(face="bold", size=10),legend.text = element_text(face = "bold", size = 20),
                                           axis.text.x = element_text(size = 23,face = "bold"), 
                                           axis.text.y = element_text(size = 23,face = "bold"), axis.title = element_text(size = 26, face = "bold"))
(p_score_counts <- ppp + aes(x= id, y =num, fill = kind) + geom_bar(stat= "identity", position = "dodge"))
?geom_density




########### position 的 覆盖度  peak coverage  ################
getFreqposition <- function(distance,data){
  freq <-vector(mode="numeric",length=0)
  len <- seq(0,distance)
  for (i in len ){
    # dataf <- filter(data1, p5 <= i & p3 >= i )
    dataf <- filter(data, p5 <= i, p3 >= i )
    freq[i+1] <- length(dataf$read_id)
  }
  freq2 <- data.frame(len,freq)
  return(freq2)
}
freq1 <- getFreqposition(500,data1)
freq2 <- getFreqposition(500,data2)
freq1 <- cbind(freq1, id = rep("ribo", length(freq1$len)))
freq2 <- cbind(freq2, id = rep("mRNA", length(freq2$len)))
freq <- rbind(freq1, freq2)
head(freq)
pp <- ggplot(freq) +theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 29), 
    legend.position = c(0.8,0.8),legend.text =  element_text(face="bold", size=20),
    axis.text.x = element_text(size = 23,face = "bold"), axis.text.y = element_text(size = 23,face = "bold"), 
    axis.title = element_text(size = 26, face = "bold"), panel.grid =element_blank())
( p_freq <- pp + aes(x = len, y = freq/(max(freq)), color = id) + geom_line(size = 1) 
  + theme(legend.text = element_text(face = "bold",size = 30)) + xlab("position[condon]") + ylab("read count"))

########韦恩图######
# id == c("ssu29","ssu37") 的用法是有误的 不等于 id == "ssu19" | id == "ssu29" 神奇
s19_29 <- filter(score,id == "ssu19" | id == "ssu29")
# tail(s19_29)
s29_37 <- filter(score,id == "ssu29" | id == "ssu37")
s19_37 <- filter(score,id == "ssu19" | id == "ssu37")
s19_29_37 <- filter(score,id == "ssu19" | id == "ssu29" | id == "ssu37")
sall <- filter(score,id == "ssu19" | id == "ssu29" | id == "ssu37" | id == "gene")
length(s19_29$reads);length(s29_37$reads);length(s19_37$reads);length(s19_29_37$reads);length(sall$reads)
length(unique(s19_29$reads))
length(unique(s29_37$reads))
length(unique(s19_37$reads))
length(unique(s19_29_37$reads))
length(filter(score, id =="gene")$reads)
length(unique(sall$reads))
length(unique(score$reads))
1301-1060
1625-1173
#### 计算之后画 韦恩图

draw.triple.venn(area1=800, area2=501, area3=324
                 ,n12=241, n23=143, n13=164, n123=96
                 ,category = c('19nt','29nt','37nt')
                 ,col=c('red','green','blue'),fill=c('red','green','blue')
                 ,cat.col=c('red','green','blue')
                 ,reverse = FALSE)


################# weblogo ###################
?ggseqlogo

reads_s19 <- filter(data, id == "ssu", p3 == 6, length >= 19 & length < 75 )$gene_read
reads_s19 <- filter(data, id == "ssu", p3 == 6, length == 19)$gene_read
reads_s19 <- str_sub(reads_s19,start = -19, end = -1)
unirea_s19 <- unique(reads_s19)
length(reads_s19);length(unirea_s19)
head(unirea_s19)

reads_s29 <- filter(data, id == "ssu", p3 == 16, length == 29 | length == 47 )$gene_read
reads_s29 <- filter(data, id == "ssu", p3 == 16, length == 29 )$gene_read
reads_s29 <- str_sub(reads_s29,start = -29, end = -1)
unirea_s29 <- unique(reads_s29)
length(reads_s29);length(unirea_s29)
head(unirea_s29)

reads_s37 <- filter(data, id == "ssu", p3 == 24, length == 37 | length == 55 )$gene_read
reads_s37 <- filter(data, id == "ssu", p3 == 24, length == 37 )$gene_read
reads_s37 <- str_sub(reads_s37,start = -37, end = -1)
unirea_s37 <- unique(reads_s37)
length(reads_s37);length(unirea_s37)
head(unirea_s37)

reads_ram <- head(data3$ram_aug,2000 )
head(reads_ram)
unirea_ram <- unique(reads_ram)
length(reads_ram);length(unirea_ram)

#### gene 的没有control mrna
reads_gene <- str_sub(data_gene$seq,start = 2, end = 20)
head(reads_gene)
unirea_gene <- unique(reads_gene)
length(reads_gene);length(unirea_gene)

# control 了的 mrna
sample_mrna_gene <- read_csv("G:/40s/sample_mrna_gene.csv", col_names = F )
head(sample_mrna_gene);length(sample_mrna_gene$X1) # 3000
# length(sample_mrna_gene$X1);length(unique(sample_mrna_gene$X1)) #886
sample_data_gene <- filter(data_gene, gene_id %in% sample_mrna_gene$X1 )
head(sample_data_gene);length(sample_data_gene$gene_id) # 875
# length(unique(sample_data_gene$gene_id))
reads_sam_gene <- str_sub(sample_data_gene$seq,start = 2, end = 20)
head(reads_sam_gene)
unirea_sam_gene <- unique(reads_sam_gene)
length(unirea_sam_gene);length(unirea_sam_gene)

####用的 reads
##基因如果直接用reads 的话 用这个
sample_mrna_gene <- read_csv("G:/40s/sample_mrna_gene.csv", col_names = F )
head(sample_mrna_gene);length(sample_mrna_gene$X1)
head(data_gene)
sample_gene <- merge(data_gene,sample_mrna_gene, by.x = "gene_id", by.y = "X1") 
head(sample_gene);length(sample_gene$seq)
# which(is.na(sample_gene))
length(unique(sample_gene$seq)) # 875

reads_sam_gene <- str_sub(sample_gene$seq,start = 2, end = 20)
head(reads_sam_gene)
# unirea_sam_gene <- unique(reads_sam_gene)
# length(unirea_sam_gene);length(unirea_sam_gene)

#######??????????? 解决不了这个问题 
# ?data.frame
# data.frame(x1 = jj, id = seq(1,886)) -> x1
# data.frame(x2 = sample_data_gene$gene_id, id = seq(1,875)) -> x2
# head(x1);head(x2)
# hhhh <- merge(x1,x2 , by = "id", all =T)
# which(is.na(hhhh))

##### 用的 我们的 ribo
# reads_r19 <- filter(data, id == "rs", p5 == (-12) , length >= 19)$gene_read

reads_r19 <- filter(data, id == "rs", p5 == (-12) , length == 28)$gene_read
reads_r19 <- str_sub(reads_r19,start = 1, end = 19)
print(head(reads_r19))
unirea_r19 <- unique(reads_r19)
print(length(reads_r19));
print(length(unirea_r19))
head(unirea_r19)


#################### 前面都是废的  只有这里开始
## 直接借用这里的  ssu_score <-rbind(ssu19,ssu29, ssu37, gene,rs28,ssu47,ssu55)
head(ssu19)
( p_s19_logo <- ggseqlogo(as.character(ssu19$reads),method='bits', seq_type = "dna") +theme_logo(base_size = 20)
  + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 19nt logo"))
( p_s29_logo <- ggseqlogo(as.character(ssu29$reads),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
  + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 29nt logo"))
( p_s37_logo <- ggseqlogo(as.character(ssu37$reads),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
  + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 37nt logo"))
# ( p_gene_logo <- ggseqlogo(as.character(unirea_gene),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
# + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="gene logo"))
# ( p_sam_gene_logo <- ggseqlogo(as.character(unirea_sam_gene),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="gene logo"))
# ( p_ram_logo <- ggseqlogo(as.character(unirea_ram),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
# + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ramdon logo"))
( p_gene_logo <- ggseqlogo(as.character(gene$reads),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
  + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="gene 19nt logo"))
( p_rs28_logo <- ggseqlogo(as.character(rs28$reads),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
  + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="rs 28nt logo"))
( p_ssu47_logo <- ggseqlogo(as.character(ssu47$reads),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
  + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 47nt logo"))
( p_ssu55_logo <- ggseqlogo(as.character(ssu55$reads),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
  + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 55nt logo"))



#### reads 的用这个
# ( p_s19_logo <- ggseqlogo(as.character(reads_s19),method='bits', seq_type = "dna") +theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 19nt logo"))
# ( p_s29_logo <- ggseqlogo(as.character(reads_s29),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 29nt logo"))
# ( p_s37_logo <- ggseqlogo(as.character(reads_s37),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 37nt logo"))
# # ( p_gene_logo <- ggseqlogo(as.character(unirea_gene),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
# # + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="gene logo"))
# ( p_sam_gene_logo <- ggseqlogo(as.character(reads_sam_gene),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="gene logo"))
# # ( p_ram_logo <- ggseqlogo(as.character(unirea_ram),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
# # + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ramdon logo"))
# ( p_rs19_logo <- ggseqlogo(as.character(reads_r19),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="rs 19nt logo"))
# ( p_s19_logo <- ggseqlogo(as.character(unirea_s19),method='bits', seq_type = "dna") +theme_logo(base_size = 20) 
#      + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 19nt logo"))
# ( p_s29_logo <- ggseqlogo(as.character(unirea_s29),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 29nt logo"))
# ( p_s37_logo <- ggseqlogo(as.character(unirea_s37),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 37nt logo"))
# # ( p_gene_logo <- ggseqlogo(as.character(unirea_gene),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   # + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="gene logo"))
# ( p_sam_gene_logo <- ggseqlogo(as.character(unirea_sam_gene),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="gene logo"))
# # ( p_ram_logo <- ggseqlogo(as.character(unirea_ram),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   # + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ramdon logo"))
# ( p_rs19_logo <- ggseqlogo(as.character(unirea_r19),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="rs 19nt logo"))
# 
# #### reads 的用这个
# ( p_s19_logo <- ggseqlogo(as.character(reads_s19),method='bits', seq_type = "dna") +theme_logo(base_size = 20) 
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 19nt logo"))
# ( p_s29_logo <- ggseqlogo(as.character(reads_s29),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 29nt logo"))
# ( p_s37_logo <- ggseqlogo(as.character(reads_s37),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ssu 37nt logo"))
# # ( p_gene_logo <- ggseqlogo(as.character(unirea_gene),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
# # + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="gene logo"))
# ( p_sam_gene_logo <- ggseqlogo(as.character(reads_sam_gene),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="gene logo"))
# # ( p_ram_logo <- ggseqlogo(as.character(unirea_ram),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
# # + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="ramdon logo"))
# ( p_rs19_logo <- ggseqlogo(as.character(reads_r19),method='bits', seq_type = "dna")+theme_logo(base_size = 20)
#   + theme(plot.title = element_text(hjust = 0.5, size = 30)) + labs(title ="rs 19nt logo"))

filescorecount <- "G:/40s/weblogocounts.csv"
counts <- read_delim(filescorecount,col_names = T, delim = ",")
head(counts)
ppp <- ggplot(counts) + theme_bw()  +theme(plot.title = element_text(hjust = 0.5, size = 30), 
                                           legend.position = c(0.9,0.85),legend.title = element_text(face="bold", size=10),legend.text = element_text(face = "bold", size = 20),
                                           axis.text.x = element_text(size = 23,face = "bold"), 
                                           axis.text.y = element_text(size = 23,face = "bold"), axis.title = element_text(size = 26, face = "bold"))
(p_score_counts <- ppp + aes(x= id, y =num, fill = kind) + geom_bar(stat= "identity", position = "dodge"))

########### save weblogo  ##########
###新建画图页面
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

# print(p_s19_logo, vp = vplayout(1,1:3))   ###将（1,1)和(1,2)的位置画图c
# print(p_s29_logo, vp = vplayout(2,1:5))   ###将(2,1)的位置画图b
# print(p_s37_logo, vp = vplayout(3,1:5))  ###将（2,2)的位置画图a
# # print(p_gene_logo, vp = vplayout(4,1:3))  ###将（2,2)的位置画图a
# print(p_gene_logo, vp = vplayout(4,1:3))
# # print(p_ram_logo, vp = vplayout(5,1:3))  ###将（2,2)的位置画图a
# print(p_rs28_logo, vp = vplayout(5,1:3))  ###将（2,2)的位置画图a
# print(p_ssu47_logo, vp = vplayout(5,1:3))
# print(p_ssu55_logo, vp = vplayout(5,1:3))

dev.off() ##画下一幅图，记得关闭窗口


# 备份的 之前的 -----------------------------------------------------------------
#### 统计reads数目
r_c <- length((filter(data, id == "ribo"))$chr)
r_c_28 <- length((filter(data, id == "ribo", length == 28))$chr)
m_c <- length((filter(data, id == "mRNA"))$chr)
m_c_28 <- length((filter(data, id == "mRNA", length == 28))$chr)
count <- data.frame(num = c(r_c, r_c_28, m_c, m_c_28),  id = c("ribo","ribo_28","mRNA","mRNA_28"), x = c(0,1,2,3))
pp <- ggplot(count) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 29), 
                                         legend.position = c(0.8,0.8),legend.text =  element_text(face="bold", size=20),
                                         axis.text.x = element_text(size = 23,face = "bold"), axis.text.y = element_text(size = 23,face = "bold"), 
                                         axis.title = element_text(size = 26, face = "bold"), panel.grid =element_blank())
count 
(p_count <- pp + aes(x = id , y= num)+ geom_bar(stat = "identity",position = "identity") + labs(title = "read counts"))
# position 的 覆盖度  ##
getFreqposition <- function(distance,data){
  freq <-vector(mode="numeric",length=0)
  len <- seq(0,distance)
  for (i in len ){
    # dataf <- filter(data1, p5 <= i & p3 >= i )
    dataf <- filter(data, p5 <= i, p3 >= i )
    freq[i+1] <- length(dataf$read_id)
  }
  freq2 <- data.frame(len,freq)
  return(freq2)
}
freq1 <- getFreqposition(500,data1)
freq2 <- getFreqposition(500,data2)
freq1 <- cbind(freq1, id = rep("ribo", length(freq1$len)))
freq2 <- cbind(freq2, id = rep("mRNA", length(freq2$len)))
freq <- rbind(freq1, freq2)
head(freq)
pp <- ggplot(freq) +theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 29), 
                                       legend.position = c(0.8,0.8),legend.text =  element_text(face="bold", size=20),
                                       axis.text.x = element_text(size = 23,face = "bold"), axis.text.y = element_text(size = 23,face = "bold"), 
                                       axis.title = element_text(size = 26, face = "bold"), panel.grid =element_blank())
( p_freq <- pp + aes(x = len, y = freq/(max(freq)), color = id) + geom_line(size = 1) 
  + theme(legend.text = element_text(face = "bold",size = 30)) + xlab("position[condon]") + ylab("read count"))
#  weblogo ？ ##
### ？？？？？  有什么用 以及为什么要做？？ 可以推导出什么之类的
unique((filter(data, id == "ssu", p5 == (-12), length == 29))$gene_read) -> ssu29_reads
unique((filter(data, id == "ssu", p5 == (-12), length == 37))$gene_read) -> ssu37_reads
unique((filter(data, id == "rs", p5 == (-12), length == 29))$gene_read) -> rs29_reads
# length((filter(data, id == "rs", p5 == (-12), length == 29)$gene_read))
unique(data3$ram_aug) -> ram_reads; sample(ram_reads,421) -> ramreads
head(ssu_reads)
head(ramreads)
length(ssu29_reads) # kind   count
length(ssu37_reads) # kind    count
length(ramreads) #  kind    count
length(rs29_reads) #  kind    count
( plot_ram_logo <- ggseqlogo(ramreads,method='bits', seq_type = "dna"))
( plot_ssu29_logo <- ggseqlogo(ssu29_reads,method='bits', seq_type = "dna"))
( plot_ssu37_logo <- ggseqlogo(ssu37_reads,method='probability', seq_type = "dna"))
( plot_rs29_logo <- ggseqlogo(rs29_reads,method='bits', seq_type = "dna"))


# RPKM   样本的 相关性 ##
?read_delim
file8 <- "G:/my-riboseq/gene_28_rpkm.csv"
# file9 <- "G:/my-riboseq/gene_2_28_rpkm.csv"
# file10 <- "G:/my-riboseq/gene_length.csv"
# file1 <- "G:/igolia/last_fp_1_final.csv"
# file2 <- "G:/igolia/last_Ho-1-T-final.csv"
data8 <- read_delim(file8, col_names = T, delim = ",")
head(data8); length(data8$gene_id) #439283  # 与原始的进行比对
p10 <- ggplot(data8) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 29), 
                                          legend.position = c(0.8,0.8),legend.text =  element_text(face="bold", size=20),
                                          axis.text.x = element_text(size = 23,face = "bold"), axis.text.y = element_text(size = 23,face = "bold"), 
                                          axis.title = element_text(size = 26, face = "bold"), panel.grid =element_blank())
(p_relate <- p10 + aes(log2(rpkm1),log2(rpkm2) ) + geom_point() 
  + xlab ("replicate #1 [rpkM]") +ylab ("replicate #2 [rpkM]")
  + geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x))
# ?cor.test
cor.test(x = log2(data8$rpkm1), y = log2(data8$rpkm2))
# Pearson's product-moment correlation
# 
# data:  log2(data8$rpkm1) and log2(data8$rpkm2)
# t = 155.25, df = 5280, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9007348 0.9104335
# sample estimates:
#       cor
# 0.9057026
###新建画图页面
grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(3,8))) ####将页面分成2*2矩阵
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p_ri_len_fra, vp = vplayout(1,1:3))   ###将（1,1)和(1,2)的位置画图c
print(p_mR_len_fra, vp = vplayout(1,4:6))   ###将(2,1)的位置画图b
print(p_count, vp = vplayout(1,7))

print(p_p5_rs, vp = vplayout(2,1:2))   ###将（1,1)和(1,2)的位置画图c
print(p_p5_mR, vp = vplayout(2,3:4))   ###将(2,1)的位置画图b
print(p_28ntp5_rs, vp = vplayout(2,5:6))
print(p_28ntp5_mR, vp = vplayout(2,7:8))

print(p_p3_rs, vp = vplayout(3,1:2))   ###将（1,1)和(1,2)的位置画图c
print(p_p3_mR, vp = vplayout(3,3:4))   ###将(2,1)的位置画图b
print(p_28ntp3_rs, vp = vplayout(3,5:6))
print(p_28ntp3_mR, vp = vplayout(3,7:8))

dev.off() ##画下一幅图，记得关闭窗口




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
options(digits=4)
options(tibble.width = Inf) # 表示 tibble 总是打印所有列 ，比如 使用 head 等函数的时候
#last_plot()
p <- ggplot() + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30), 
        legend.position = "right",legend.text =  element_text(face="bold", size=23),
        legend.title = element_blank(),
        # axis.text.x = element_text(size = 26,face = "bold"),
        # axis.ticks = element_line(),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(.20,"cm"),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(size = 25,face = "bold"),  
        axis.title = element_text(size = 28, face = "bold"))
pp <- ggplot() + theme_classic()
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
a






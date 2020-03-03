#  IR analysis 
#  input file 
#             -----IRTools: ****.quant.IRI.genes.txt
#           
# IRI : CIR_RPKM/CER_RPKM

## --------20190814-------PAY ATTENTION:  ( trouble in ggplot2  not in the base plot function)
# pdf works in the console but not in a function or when you source this from a file.
#  it will produce a corrupt pdf file and the way to fix it us use
#  print(p) within a function.-----but print can only save one figure  use gridExtra????
#   print(get(p))
# In a console. "p" is automatically printed but not in a function or when you source the file.


## 20191024 try to find obvious IRI changes use new criteria

path <- "/Users/wangkai/Dropbox/QingChen_Peng_Lab/00_Projects/Seto/"
setwd(path)

# library(gridExtra)-------
library(plyr)
library(ggplot2)
library("biomaRt")
library(grid)
library(futile.logger)
library(VennDiagram)


# pdf(paste0(timepoint,".pdf"))
#-----------------------------------------------------------------
#    1.  data preparation 
#    1.1  read the data of IRTools  
#---------------------------------------------------------
#   list.files(path)[grep("WT",(list.files(path)))]
IR_WT <- read.table(list.files(path)[grep("WT.quant.IRI.genes.",(list.files(path)))],header = T)
Group <- rep("WT",dim(IR_WT)[1])
IR_WT <- cbind(Group,IR_WT)
# colnames(IR_WT)[2]<-"gene"
IR_KD <- read.table(list.files(path)[grep("KD.quant.IRI.genes.",(list.files(path)))],header = T)
Group <- rep("KD",dim(IR_KD)[1])
IR_KD <- cbind(Group,IR_KD)
# colnames(IR_KD)[2]<-"gene"
head(IR_KD)
##  data pretreated : set threshold
#   (1) ______ IRI remove IRI > 1 (including Inf ones)
#   (2)________ gene_CER_RPKM > 1 

## 20191024  only remove both IRI >1 later in merged files

# IR_WT <- IR_WT[which(0<IR_WT$gene_IRI &IR_WT$gene_IRI<=1&IR_WT$gene_CER_RPKM>1 ),]
# dim(IR_WT)  
# IR_KD <- IR_KD[which(0<IR_KD$gene_IRI&IR_KD$gene_IRI<=1&IR_KD$gene_CER_RPKM>1 ),]
# dim(IR_KD)  

IR_WT <- IR_WT[which(0<IR_WT$gene_IRI&IR_WT$gene_CER_RPKM>1 ),]
dim(IR_WT)  
IR_KD <- IR_KD[which(0<IR_KD$gene_IRI&IR_KD$gene_CER_RPKM>1 ),]
dim(IR_KD) 
## rbind ( group difference)
IR_WT_IR_KD_rbind<-rbind(IR_WT,IR_KD) 
summary(IR_WT_IR_KD_rbind$gene_CER_RPKM)
plot(density(log10(IR_WT_IR_KD_rbind$gene_CER_RPKM)))
## merge / join  ( each IR change in different group)  common ones
IR_WT_IR_KD_merge <- merge(IR_WT,IR_KD,by = "gene_id")
dim(IR_WT_IR_KD_merge)  

head(IR_WT_IR_KD_merge)
IR_changes <- IR_WT_IR_KD_merge$gene_IRI.y/IR_WT_IR_KD_merge$gene_IRI.x
# summary(IR_changes)
plot(density(IR_changes))
IR_WT_IR_KD_merge <- cbind(IR_WT_IR_KD_merge,IR_changes)

###  20191024  only remove both IRI >1 later in merged files
both_high1 <- which(IR_WT_IR_KD_merge$gene_IRI.x>1&IR_WT_IR_KD_merge$gene_IRI.y>1)
IR_WT_IR_KD_merge <- IR_WT_IR_KD_merge[-both_high1,]

# ----extract the top changes genes : add some thresholds:--------
plot(density(IR_WT_IR_KD_merge$IR_changes))

IRI_threshold <- quantile(c(IR_WT_IR_KD_merge$gene_IRI.x,IR_WT_IR_KD_merge$gene_IRI.y),0.25)

# up_index<-which(IR_WT_IR_KD_merge$IR_changes >1 & IR_WT_IR_KD_merge$gene_IRI.y>IRI_threshold)
# down_index <-which(IR_WT_IR_KD_merge$IR_changes <1 & IR_WT_IR_KD_merge$gene_IRI.x>IRI_threshold)
# IR_WT_IR_KD_merge <-IR_WT_IR_KD_merge[c(up_index,down_index),]


# 20191024  1.   add CER  > 10
up_index<-which(IR_WT_IR_KD_merge$IR_changes >1 & IR_WT_IR_KD_merge$gene_IRI.y>IRI_threshold &IR_WT_IR_KD_merge$gene_CER_RPKM.x >10)
down_index <-which(IR_WT_IR_KD_merge$IR_changes <1 & IR_WT_IR_KD_merge$gene_IRI.x>IRI_threshold & IR_WT_IR_KD_merge$gene_CER_RPKM.y >10)
write.csv(IR_WT_IR_KD_merge[up_index,],row.names =FALSE ,quote = FALSE,file = "20191024_IR_changes_increase_CER_10.csv")
write.csv(IR_WT_IR_KD_merge[down_index,],row.names =FALSE ,quote = FALSE,file = "20191024_IR_changes_decrease_CER_10.csv")
length(up_index)
length(down_index)
# 20191024  2.   add CER  > 50
up_index<-which(IR_WT_IR_KD_merge$IR_changes >1 & IR_WT_IR_KD_merge$gene_IRI.y>IRI_threshold &IR_WT_IR_KD_merge$gene_CER_RPKM.x >50)
down_index <-which(IR_WT_IR_KD_merge$IR_changes <1 & IR_WT_IR_KD_merge$gene_IRI.x>IRI_threshold & IR_WT_IR_KD_merge$gene_CER_RPKM.y >50)
write.csv(IR_WT_IR_KD_merge[up_index,],row.names =FALSE ,quote = FALSE,file = "20191024_IR_changes_increase_CER_50.csv")

write.csv(IR_WT_IR_KD_merge[down_index,],row.names =FALSE ,quote = FALSE,file = "20191024_IR_changes_decrease_CER_50.csv")

length(up_index)
length(down_index)

# 20191024  3.  add CER  > 10 and IRI changge > 2 or  < 0.5
up_index<-which(IR_WT_IR_KD_merge$IR_changes >2 & IR_WT_IR_KD_merge$gene_IRI.y>IRI_threshold &IR_WT_IR_KD_merge$gene_CER_RPKM.x >10)
down_index <-which(IR_WT_IR_KD_merge$IR_changes < 0.5  & IR_WT_IR_KD_merge$gene_IRI.x>IRI_threshold & IR_WT_IR_KD_merge$gene_CER_RPKM.y >10)
write.csv(IR_WT_IR_KD_merge[up_index,],row.names =FALSE ,quote = FALSE,file = "20191024_IR_changes_increase_2_CER_10_.csv")
write.csv(IR_WT_IR_KD_merge[down_index,],row.names =FALSE ,quote = FALSE,file = "20191024_IR_changes_decrease_0dot5_CER_10.csv")
length(up_index)
length(down_index)

#---------------------------------------------------------
# 2 --------draw figure ----based on rbind date
#   -----------see group difference
#---------------------------------------------------------
# 
# ### gene_IRI
# p_gene_IRI_1<-ggplot(IR_WT_IR_KD_rbind ,aes(x = gene_IRI,fill = Group))+geom_density(alpha = .3)+ 
#   labs(title = "gene_IRI (all the filtered gene)")+ theme(plot.title = element_text(hjust = 0.5))
# 
# p_gene_IRI_2<-ggplot(IR_WT_IR_KD_rbind,aes(x = Group,y=gene_IRI))+geom_violin()+
#   geom_boxplot(width = 0.1,fill = "black",outlier.colour =NA)+stat_summary(fun.y = "median",shape =21,geom = "point",fill = "white",size =1)+
#   labs(title = "gene_IRI (all the filtered gene)") +  theme(plot.title = element_text(hjust = 0.5))
# 
# 
# ### gene_IRI
# ggplot(IR_WT_IR_KD_rbind ,aes(x = log2(gene_IRI),fill = Group))+geom_density(alpha = .3)+ 
#   labs(title = "gene_IRI (all the filtered gene)")+ theme(plot.title = element_text(hjust = 0.5))
# 
# 2.1  gene_CIR_length
p1<- ggplot(IR_WT_IR_KD_rbind ,aes(x = gene_CIR_length,fill = Group))+geom_density(alpha = .3) + labs(title = "gene_CIR_length(all the filtered gene) \n xlim(0,250K)")+xlim(0,250000)
p2<-ggplot(IR_WT_IR_KD_rbind_top1000 ,aes(x = gene_CIR_length,fill = Group))+geom_density(alpha = .3) + labs(title = "gene_CIR_length(top1000)  \n xlim(0,250K)")+xlim(0,250000)
p3<-ggplot(IR_WT_IR_KD_rbind_bottom1000 ,aes(x = gene_CIR_length,fill = Group))+geom_density(alpha = .3) + labs(title = "gene_CIR_length(bottom1000) \n xlim(0,250K)")+xlim(0,250000)
p1
p2
p3
## 2.2  gene_CIR_RPKM
p4<-ggplot(IR_WT_IR_KD_rbind ,aes(x = gene_CIR_RPKM,fill = Group))+geom_density(alpha = .3)+xlim(0,10) + labs(title = "gene_CIR_RPKM (all the filtered gene) \n xlim(0,10)")
p5<-ggplot(IR_WT_IR_KD_rbind_top1000 ,aes(x = gene_CIR_RPKM,fill = Group))+geom_density(alpha = .3)+xlim(0,10)  + labs(title = "gene_CIR_RPKM (top1000) \n xlim(0,10)")
p6<-ggplot(IR_WT_IR_KD_rbind_bottom1000 ,aes(x = gene_CIR_RPKM,fill = Group))+geom_density(alpha = .3)+xlim(0,10)  + labs(title = "gene_CIR_RPKM (bottom1000) \n xlim(0,10)")
p4
p5
p6
### 2.3  gene_IRI
p7<-ggplot(IR_WT_IR_KD_rbind ,aes(x = gene_IRI,fill = Group))+geom_density(alpha = .3)+ 
  labs(title = "gene_IRI (all the filtered gene)")+ theme(plot.title = element_text(hjust = 0.5))
p8<-ggplot(IR_WT_IR_KD_rbind_top1000 ,aes(x = gene_IRI,fill = Group))+geom_density(alpha = .3)  + 
  labs(title = "gene_IRI (top1000)")+ theme(plot.title = element_text(hjust = 0.5))
p9<-ggplot(IR_WT_IR_KD_rbind_bottom1000 ,aes(x = gene_IRI,fill = Group))+geom_density(alpha = .3)  + 
  labs(title = "gene_IRI (bottom1000)")+ theme(plot.title = element_text(hjust = 0.5))
p7
p8
p9
## log scale
p7<-ggplot(IR_WT_IR_KD_rbind ,aes(x =log2(gene_IRI),fill = Group))+geom_density(alpha = .3)+ 
  labs(title = "gene_IRI:log scale (all the filtered gene)")+ theme(plot.title = element_text(hjust = 0.5))
p8<-ggplot(IR_WT_IR_KD_rbind_top1000 ,aes(x = log2(gene_IRI),fill = Group))+geom_density(alpha = .3)  + 
  labs(title = "gene_IRI:log scale (top1000)")+ theme(plot.title = element_text(hjust = 0.5))
p9<-ggplot(IR_WT_IR_KD_rbind_bottom1000 ,aes(x =log2(gene_IRI),fill = Group))+geom_density(alpha = .3)  + 
  labs(title = "gene_IRI:log scale (bottom1000)")+ theme(plot.title = element_text(hjust = 0.5))

# ggplot(IR_WT_IR_KD_rbind,aes(x = Group,y=gene_IRI))+geom_violin()+
#   geom_boxplot(width = 0.1,fill = "black",outlier.colour =NA)+stat_summary(fun.y = "median",shape =21,geom = "point",fill = "white",size =1)+
#   labs(title = "gene_IRI (all the filtered gene)") +  theme(plot.title = element_text(hjust = 0.5))
# 
# ggplot(IR_WT_IR_KD_rbind_top1000,aes(x = Group,y=gene_IRI))+geom_violin()+
#   geom_boxplot(width = 0.1,fill = "black",outlier.colour =NA)+stat_summary(fun.y = "median",shape =21,geom = "point",fill = "white",size =2.5)+
#   labs(title = "gene_IRI (top1000)") +  theme(plot.title = element_text(hjust = 0.5))
# ggplot(IR_WT_IR_KD_rbind_bottom1000,aes(x = Group,y=gene_IRI))+geom_violin()+
#   geom_boxplot(width = 0.1,fill = "black",outlier.colour =NA)+stat_summary(fun.y = "median",shape =21,geom = "point",fill = "white",size =2.5)+
#   labs(title = "gene_IRI (bottom1000)") +  theme(plot.title = element_text(hjust = 0.5))
# 
## 2.4  host gene FPKM
p10 <- ggplot(IR_WT_IR_KD_rbind ,aes(x = log2(gene_CER_RPKM+1),fill = Group))+geom_density(alpha = .3)+ 
  labs(title = "log2(gene_CER_RPKM +1) (all the filtered gene)")+  theme(plot.title = element_text(hjust = 0.5))
p11<-ggplot(IR_WT_IR_KD_rbind_top1000 ,aes(x = log2(gene_CER_RPKM+1),fill = Group))+geom_density(alpha = .3)+
  labs(title = "log2(gene_CER_RPKM +1) (top1000)")+  theme(plot.title = element_text(hjust = 0.5))
p12<-ggplot(IR_WT_IR_KD_rbind_bottom1000 ,aes(x = log2(gene_CER_RPKM+1),fill = Group))+geom_density(alpha = .3)+ 
  labs(title = "log2(gene_CER_RPKM +1) (bottom1000)")+  theme(plot.title = element_text(hjust = 0.5))
p10
p11
p12
# ggplot(IR_WT_IR_KD_rbind,aes(x = Group,y=log2(gene_CER_RPKM+1)))+geom_violin()+
#   geom_boxplot(width = 0.1,fill = "black",outlier.colour =NA)+stat_summary(fun.y = "median",shape =21,geom = "point",fill = "white",size =2.5)+
#   labs(title = "log2(gene_CER_RPKM +1) (all the filtered gene)") +  theme(plot.title = element_text(hjust = 0.5))
# 
# ggplot(IR_WT_IR_KD_rbind_top1000,aes(x = Group,y=log2(gene_CER_RPKM+1)))+geom_violin()+
#   geom_boxplot(width = 0.1,fill = "black",outlier.colour =NA)+stat_summary(fun.y = "median",shape =21,geom = "point",fill = "white",size =2.5)+
#   labs(title = "log2(gene_CER_RPKM +1) (top1000)") +  theme(plot.title = element_text(hjust = 0.5))
# 
# ggplot(IR_WT_IR_KD_rbind_bottom1000,aes(x = Group,y=log2(gene_CER_RPKM+1)))+geom_violin()+
#   geom_boxplot(width = 0.1,fill = "black",outlier.colour =NA)+stat_summary(fun.y = "median",shape =21,geom = "point",fill = "white",size =2.5)+
#   labs(title = "log2(gene_CER_RPKM +1) (bottom1000)") +  theme(plot.title = element_text(hjust = 0.5))


# 3 ----------draw figure---based on merged date----------------------------------------------------
my_theme_for_density_2d <- theme(panel.grid.major = element_line(colour="black"), panel.grid.minor = element_line(colour="black",linetype="dashed", size=0.2) )+theme_bw()

head(IR_WT_IR_KD_merge)

p13<-ggplot(IR_WT_IR_KD_merge, aes(gene_IRI.x,gene_IRI.y )) +
  # stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE,alpha = 0.9) +
  labs(y = "KD", x = "WT",title = "gene_IRI") +
  # scale_fill_continuous(low = "white",high = "black") +
  # coord_cartesian(expand = FALSE) +
  # my_theme_for_density_2d +
  geom_point(shape = '.', col = 'black')+ theme(plot.title = element_text(hjust = 0.5))
  # +geom_smooth(method = "auto")
## log scale
p13<-ggplot(IR_WT_IR_KD_merge, aes(log2(gene_IRI.x),log2(gene_IRI.y) )) +
  # stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE,alpha = 0.9) +
  labs(y = "KD", x = "WT",title = "log2(gene_IRI)") +
  # scale_fill_continuous(low = "white",high = "black") +
  # coord_cartesian(expand = FALSE) +
  # my_theme_for_density_2d +
  geom_point(shape = '.', col = 'black')+ theme(plot.title = element_text(hjust = 0.5))
# +geom_smooth(method = "auto")

ggplot(IR_WT_IR_KD_merge, aes(gene_IRI.x,gene_IRI.y))+ labs(y = "KD", x = "WT",title = "gene_IRI")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point(shape = '.',col = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,'red',"black"),
             size = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1),
             alpha = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1))

 
ggplot(IR_WT_IR_KD_merge, aes(log2(gene_IRI.x),log2(gene_IRI.y)))+ labs(y = "KD", x = "WT",title = "log2(gene_IRI)")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point(col = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,'red',"black"),
             size = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1),
             alpha = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1))+
  geom_text(aes(label=ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,as.character(gene_id),'')),position=position_jitter())

# Install ggrepel package if needed
# install.packages("devtools")
# devtools::install_github("slowkow/ggrepel")
library(ggrepel)

# down
p_down<-ggplot(IR_WT_IR_KD_merge, aes(log2(gene_IRI.x),log2(gene_IRI.y)))+ labs(y = "KD", x = "WT",title = "log2(gene_IRI)")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point(col = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,'red',"black"),
             size = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1),
             alpha = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1))+
  geom_text_repel(data =IR_WT_IR_KD_merge[down_index,],aes(label = gene_id) )
# up

up <- IR_WT_IR_KD_merge[up_index,]
up_order <- up[order(up$IR_changes,decreasing = T),]
up_order$IR_changes

p_up_1<-ggplot(IR_WT_IR_KD_merge, aes(log2(gene_IRI.x),log2(gene_IRI.y)))+ labs(y = "KD", x = "WT",title = "log2(gene_IRI)")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point(col = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,'red',"black"),
             size = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1),
             alpha = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1))+
  geom_text_repel(data =up_order[1:100,],aes(label = gene_id) ,segment.size = 0.3,force =10 )  


p_up_2<-ggplot(IR_WT_IR_KD_merge, aes(log2(gene_IRI.x),log2(gene_IRI.y)))+ labs(y = "KD", x = "WT",title = "log2(gene_IRI)")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point(col = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,'red',"black"),
             size = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1),
             alpha = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1))+
  geom_text_repel(data =up_order[101:200,],aes(label = gene_id),segment.size = 0.3,force =10 )          

p_up_3<-ggplot(IR_WT_IR_KD_merge, aes(log2(gene_IRI.x),log2(gene_IRI.y)))+ labs(y = "KD", x = "WT",title = "log2(gene_IRI)")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point(col = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,'red',"black"),
             size = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1),
             alpha = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1))+
  geom_text_repel(data =up_order[201:300,],aes(label = gene_id),segment.size = 0.3,force =10 )          

p_up_4<-ggplot(IR_WT_IR_KD_merge, aes(log2(gene_IRI.x),log2(gene_IRI.y)))+ labs(y = "KD", x = "WT",title = "log2(gene_IRI)")+ theme(plot.title = element_text(hjust = 0.5))+
  geom_point(col = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,'red',"black"),
             size = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1),
             alpha = ifelse(IR_WT_IR_KD_merge$gene_id%in% IR_WT_IR_KD_merge[c(up_index,down_index),]$gene_id,1,0.1))+
  geom_text_repel(data =up_order[300:dim(up_order)[1],],aes(label = gene_id) ,segment.size = 0.3,force =10)          





##############
p14<-ggplot(IR_WT_IR_KD_merge, aes(log2(gene_CIR_RPKM.x+1),log2(gene_CIR_RPKM.y+1) )) +
  stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE,alpha = 0.9) +
  labs(y = "KD", x = "WT",title = "log2(gene_CIR_RPKM+1)") +
  scale_fill_continuous(low = "white",high = "black") +
  coord_cartesian(expand = FALSE) +
  my_theme_for_density_2d +
  geom_point(shape = '.', col = 'red')+ theme(plot.title = element_text(hjust = 0.5))

p15<-ggplot(IR_WT_IR_KD_merge, aes(log2(gene_CER_RPKM.x+1),log2(gene_CER_RPKM.y+1) )) +
  stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE,alpha = 0.9) +
  labs(y = "KD", x = "WT",title = "log2(gene_CER_RPKM+1)") +
  scale_fill_continuous(low = "white",high = "black") +
  coord_cartesian(expand = FALSE) +
  xlim(0,10)+ylim(0,10)+
  my_theme_for_density_2d +
  geom_point(shape = '.', col = 'red')+ theme(plot.title = element_text(hjust = 0.5))
p13
p14
p15

plot(density(log10(IR_WT_IR_KD_rbind$gene_CER_RPKM)+1),main ="log10(gene_CER_RPKM +1)" )


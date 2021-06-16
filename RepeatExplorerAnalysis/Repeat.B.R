if (!requireNamespace('devtools', quietly = TRUE))
       install.packages('devtools')
devtools::install_github('kevinblighe/PCAtools')

install.packages(c("BiocManager","devtools","ggplot2","tidyr","scales","cowplot","ggfortify","corrplot","tidyverse","ggrepel","ggdendro","caret","MASS","reshape2","ROCR"))

#library(PCAtools)
library(scales)
library(cowplot)
library(ggfortify)
#library(corrplot)
library(tidyverse)
library(ggrepel)
#library(ggdendro)
#library(caret)
#library(MASS)
library(reshape2)
#library(ROCR)
#library(gridExtra)
#library(testit)
#library(data.table)
#library(factoextra)

sessionInfo()
setwd('Y:/corrinne/BEgenomes/Bmanuscript/repeats')

# import file with all counts
# evaluate if 0.01% cutoff is reasonable (here, cluster 274)

# Columns added in order to rename for code development 
data <- read.table("B.clusters.counts.txt", header = T, sep="\t")
data$annotation <- data$annotation %>% replace_na("none")

# Add size/percent columns
data$size <- as.numeric(rowSums(data[,-(1:2)]))
data$percent <- cumsum(data$size)/sum(data$size)
which(data$size>sum(data$size)*0.0001)
# clusters 1-345
which(rowMeans(data[,3:18])>10)
# clusters 1-363

data$clusternum <- as.numeric(row.names(data))

# 345 clusters each contain >= 0.01% of reads
png("cotton.cutoff.png", 5000, 5000, pointsize=12, res=600)
ggplot(data, aes(x=clusternum, y=percent)) + geom_line(size=1) + 
    geom_vline(xintercept=345, color='yellow3', size=1) + 
    geom_vline(xintercept=363, color='yellow3', size=1) + 
    scale_x_log10(labels=comma) + scale_y_log10() + 
    geom_vline(xintercept=0, color="grey")+ geom_hline(yintercept=0, color="grey")
dev.off()



##########################################################
################### PCA of scaled data ###################
##########################################################

clusters <- data[c(1:345),c(2:(ncol(data)-3))]

# make table for PCA
counts <- as.data.frame(t(clusters [,(2:(ncol(clusters)))]))
counts$species <- as.factor(row.names(counts))

# just a few stats, in case we want to see them
summary(counts)
ncol(counts)

##################### PCA with prcomp
pr.counts <- prcomp(counts[,-346], scale=TRUE)
summary(pr.counts)

### how many PC capture the information

pveC <- 100*pr.counts$sdev^2/sum(pr.counts$sdev^2)
par(mfrow=c(1,2))
plot(pveC, type ="o", ylab="PVE", xlab="Principal Component", col ="blue")
plot(cumsum(pveC), type="o", ylab ="Cumulative PVE", xlab="Principal Component", col ="brown3")
# most information in the first 5


countpca <- data.frame('species' = counts$species, pr.counts$x[,1:5]) # first five components

countpca$continent <- c("Americas","Americas","Americas","Americas","Americas",
    "Africa","Africa","Africa","Africa","Africa","Australia","Australia","Australia","Australia","Australia","Australia")

countpcaFig <- ggplot(data = countpca, aes(x = PC1, y = PC2, col = continent)) + 
    geom_point() + geom_text_repel(aes(label=species)) +
    scale_color_manual(values=c("Americas"="darkgreen","Africa"="darkred","Australia"="darkblue")) +
    labs(title = "", x="PC1 (11.6%)", y="PC2 (8.1%)")


ggsave(filename="Figure_nicePCA.jpg",plot=countpcaFig,device="jpeg",width=8,height=8,units="in",dpi=400)




########### are cluster abundances correlated
correlations <- cor(counts[,1:100])
res1 <- cor.mtest(counts[,1:100], conf.level = .99)
png("correlation.plot.png", 10000, 10000, pointsize=12, res=600)
corrplot(correlations, p.mat = res1$p, insig = "blank",type="upper",tl.pos = "td", tl.cex = 0.5,diag=F)
dev.off()




########### characterize composition ###########

# 9 multiplier represents # reads (x) * 90nt/read * 1 kb/1000nt * 100% = # reads * 9 = # Kb in entire genome for that cluster 
Kb <- data.frame(clusters[1], apply(clusters[2:ncol(clusters)], 2, function (x) x*9))
Kbsum <- aggregate(. ~annotation, data=Kb, FUN=sum)
Kbsum <- Kbsum[-c(2,9:11),]
Kbsum$annotation <- gsub("DNA.*","DNA",Kbsum$annotation)
Kbsum[1,1] <- "unknown"
Kbsum[7,1] <- "unknown"
Kbsum <- aggregate(. ~annotation, data=Kbsum, FUN=sum)


Kbsum$D <- rowMeans(Kbsum[,2:6])
Kbsum$A <- rowMeans(Kbsum[,7:8])
Kbsum$F <- Kbsum[,9]
Kbsum$B <- Kbsum[,10]
Kbsum$E <- Kbsum[,11]
Kbsum$C <- rowMeans(Kbsum[,12:13])
Kbsum$G <- rowMeans(Kbsum[,14:16])
Kbsum$K <- Kbsum[,17]


Kbsum$Dmin <- apply(Kbsum[,2:6], 1, min)
Kbsum$Amin <- apply(Kbsum[,7:8], 1, min)
Kbsum$Fmin <- Kbsum[,9]
Kbsum$Bmin <- Kbsum[,10]
Kbsum$Emin <- Kbsum[,11]
Kbsum$Cmin <- apply(Kbsum[,12:13], 1, min)
Kbsum$Gmin <- apply(Kbsum[,14:16], 1, min)
Kbsum$Kmin <- Kbsum[,17]


Kbsum$Dmax <- apply(Kbsum[,2:6], 1, max)
Kbsum$Amax <- apply(Kbsum[,7:8], 1, max)
Kbsum$Fmax <- Kbsum[,9]
Kbsum$Bmax <- Kbsum[,10]
Kbsum$Emax <- Kbsum[,11]
Kbsum$Cmax <- apply(Kbsum[,12:13], 1, max)
Kbsum$Gmax <- apply(Kbsum[,14:16], 1, max)
Kbsum$Kmax <- Kbsum[,17]


Kbsum <- Kbsum[,-(2:17)]
Kbm <- melt(Kbsum[,-(10:25)])
min <- c(Kbsum$Dmin,
    Kbsum$Amin,
    Kbsum$Fmin,
    Kbsum$Bmin,
    Kbsum$Emin,
    Kbsum$Cmin,
    Kbsum$Gmin,
    Kbsum$Kmin)
max <- c(Kbsum$Dmax,
    Kbsum$Amax,
    Kbsum$Fmax,
    Kbsum$Bmax,
    Kbsum$Emax,
    Kbsum$Cmax,
    Kbsum$Gmax,
    Kbsum$Kmax)
Kbm$min <- min
Kbm$max <- max

limits <- aes(ymax=Kbm$max, ymin=Kbm$min)
dodge <- position_dodge(width=0.9)

TEamount <- Kbm %>% 
    mutate(variable=fct_relevel(variable,"D","F","B","E","A","G","C","K")) %>% 
    ggplot(aes(x=annotation, y=value, fill = variable)) + 
    geom_bar(stat = "identity",position = dodge) + 
    scale_y_continuous(labels=comma) + 
    scale_fill_manual(values=c("darkgreen",  "#67000d","#fed976","#a63603","#7f2704", 
        "#6baed6","#3182bd","#08519c"), labels=c("D-genome 885 Mbp",
            "F-genome 1311 Mbp","B-genome 1350 Mbp","E-genome 1560 Mbp","A-genome 1697 Mbp",
            "G-genome 1785 Mbp","C-genome 1980 Mbp","K-genome 2572 Mbp")) +
    geom_errorbar(limits, position = dodge) + 
    labs(title = "", x="", y="Amount in genome (kb)") + 
    theme_set(theme_grey(base_size=12)) + 
    theme(axis.text = element_text(size = rel(1.5)), 
        plot.title=element_text(face="bold", hjust=0.5),
        axis.title.x = element_text(face="bold", hjust=0.5, vjust=9), 
        axis.title.y = element_text(face="bold", vjust=2),
        axis.text.x = element_text(size = 12, angle = 45, hjust=0.7),
        axis.text.y = element_text(size = 12, angle = 45, vjust=0.1)) + 
    theme(legend.title=element_blank(), legend.position=c(0.15,0.8)) + 
    theme(legend.margin = margin(0, 5, 5, 5))


ggsave(filename="TEamounts.jpg",plot=TEamount,device="jpeg",width=10,height=10,units="in",dpi=400)


sum(Kbsum$D)/1000 # convert to Mb
# [1] 286.6
sum(Kbsum$A)/1000 # convert to Mb
# [1] 992.9
sum(Kbsum$F)/1000 # convert to Mb
# [1] 607.4
sum(Kbsum$B)/1000 # convert to Mb
# [1] 627.1
sum(Kbsum$E)/1000 # convert to Mb
# [1] 567.5
sum(Kbsum$C)/1000 # convert to Mb
# [1] 1253.0
sum(Kbsum$G)/1000 # convert to Mb
# [1] 1022.0
sum(Kbsum$K)/1000 # convert to Mb
# [1] 1617.6


sum(Kbsum$D)/10/885 # convert to % genome
# [1] 32.4
sum(Kbsum$A)/10/1697 # convert to % genome
# [1] 58.5
sum(Kbsum$F)/10/1311 # convert to % genome
# [1] 46.3
sum(Kbsum$B)/10/1350 # convert to % genome
# [1] 46.5
sum(Kbsum$E)/10/1560 # convert to % genome
# [1] 36.4
sum(Kbsum$C)/10/1980 # convert to % genome
# [1] 63.3
sum(Kbsum$G)/10/1785 # convert to % genome
# [1] 57.3
sum(Kbsum$K)/10/2572 # convert to % genome
# [1] 62.9

########### compare A1 vs A2 ###########

### plot a 1:1 line and display cluster comparisons relative to 1:1 expection
# remember, these have nearly equivalent genome sizes

TEkeep <- c("LTR/Gypsy","LTR","LTR/Copia","DNA/MULE-MuDR","DNA/TcMar-Mariner")

# mean
bef <- clusters[,c(1,9:11)]
bef <- bef[bef$annotation %in% TEkeep,]

bef$signbe <- {ifelse((bef$B01 > bef$E02), "positive", "negative")}
bef$percbe <- bef$B01/bef$E02
bef <- bef %>% mutate(highlightbe=case_when(percbe > 1.25 | percbe < 0.75 ~ annotation, FALSE ~ ""))
bef <- bef %>% mutate(clusterbe=case_when(percbe > 1.25 | percbe < 0.75 ~ row.names(bef), FALSE ~ ""))
bef[is.na(bef)] <- ""
bef$labelbe <- paste0(bef$highlightbe,bef$clusterbe)
bef$labelbe <- gsub("LTR/","",bef$labelbe)

BvE <- ggplot(bef, aes(x=B01, y=E02, shape=signbe, color=signbe)) + 
    geom_point(size=2) + geom_abline(intercept=0, slope=1) + 
    scale_color_manual(values=c("positive"="darkgoldenrod2", "negative"="darkred")) +  
    scale_x_continuous(expand = c(0, 0), limits=c(0,1500)) + 
    scale_y_continuous(expand = c(0, 0), limits=c(0,1500)) + 
    theme(legend.position="none", axis.title.x = element_text(face="italic", hjust=0.5), 
        axis.title.y = element_text(face="italic", vjust=0.5)) + 
    geom_text_repel(box.padding = 0.5,aes(label=highlightbe)) + 
    labs(title = "", x="G.anomalum", y="G.somalense")


bef$signbf <- {ifelse((bef$B01 > bef$F01), "positive", "negative")}
bef$percbf <- bef$B01/bef$F01
bef <- bef %>% mutate(highlightbf=case_when(percbf > 1.25 | percbf < 0.75 ~ annotation, FALSE ~ ""))
bef <- bef %>% mutate(clusterbf=case_when(percbf > 1.25 | percbf < 0.75 ~ row.names(bef), FALSE ~ ""))
bef[is.na(bef)] <- ""
bef$labelbf <- paste0(bef$highlightbf,bef$clusterbf)
bef$labelbf <- gsub("LTR/","",bef$labelbf)

BvF <- ggplot(bef, aes(x=B01, y=F01, shape=signbf, color=signbf)) + 
    geom_point(size=2) + geom_abline(intercept=0, slope=1) + 
    scale_color_manual(values=c("positive"="darkgoldenrod2", "negative"="darkred")) +  
    scale_x_continuous(expand = c(0, 0), limits=c(0,1500)) + 
    scale_y_continuous(expand = c(0, 0), limits=c(0,1500)) + 
    theme(legend.position="none", axis.title.x = element_text(face="italic", hjust=0.5), 
        axis.title.y = element_text(face="italic", vjust=0.5)) + 
    geom_text_repel(box.padding = 0.5,aes(label=highlightbf)) + 
    labs(title = "", x="G.anomalum", y="G.longicalyx")


FigPlot <- plot_grid(BvF,BvE)

ggsave(filename="Figure_TEtwoPanel.jpg",plot=FigPlot,device="jpeg",width=15,height=9,units="in",dpi=400)



# table(A1A2$annotation)
# DNA/MULE-MuDR           LTR     LTR/Copia     LTR/Gypsy 
#             2            15            17           141 

# table(sigClus$annotation)
# LTR LTR/Copia LTR/Gypsy 
#  3         7        54 

# are the repetitive values scaling with genome size
x <- c(885,1311,1350,1560,1697,1785,1980,2572)
y <- c(286.6,607.4,627.1,567.5,992.9,1022,1253,1617.6)

z <- cbind(x,y)

xyz <- ggplot(z, aes(x,y)) + geom_line() + geom_point()+ 
    geom_smooth(method='lm', formula= y~x)+geom_text_repel(aes(label=c("D","F","B","E","A","G","C","K"))) +
    labs(title = "", x="genome size", y="Mb repeats in genome")

ggsave(filename="SuppFigureTEexpAmt.jpg",plot=xyz,device="jpeg",width=15,height=9,units="in",dpi=400)

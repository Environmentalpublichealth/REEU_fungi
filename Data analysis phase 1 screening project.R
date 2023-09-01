setwd("~/Downloads/")
Data <- read.csv("Screening project phase 1 sheet new measurements.csv")

library(dplyr)
library(reshape2)

# remove rows with na

Data[is.na(Data)] <- 0

Data <- Data[-c(58,61,134,303),]

# remove duplicate rows, first duplicate removed and second kept

table(duplicated(Data[,1:2]))

uniqueData <- Data[rownames(unique(Data[,1:2], fromLast = TRUE)),]
  
names(Data)

W1R1 <- dcast(uniqueData, Strain ~ Treatment, value.var = c("W1...R1"))
       
W1R2 <- dcast(uniqueData, Strain ~ Treatment, value.var = c("W1...R2"))
              
W2R1 <- dcast(uniqueData, Strain ~ Treatment, value.var = c("W2...R1"))
              
W2R2 <- dcast(uniqueData, Strain ~ Treatment, value.var = c("W2...R2"))

Density <- dcast(uniqueData, Strain ~ Treatment, value.var = c("Density.ranking..1...most.dense..3...least.dense."))

# PFOA calculation
PFOA <- data.frame(Strain=W1R1$Strain, Diff1 = W1R1$PFOA - W1R1$Control,
                   Diff2=W1R2$PFOA - W1R2$Control,
                   Diff3=W2R1$PFOA - W2R1$Control,
                   Diff4=W2R2$PFOA - W2R2$Control,
                   Dense=Density$PFOA/Density$Control)
PFOA$R_mean <- rowMeans(PFOA[,2:5])


# PFOS calculation
PFOS <- data.frame(Strain=W1R1$Strain, Diff1 = W1R1$PFOS - W1R1$Control,
                   Diff2=W1R2$PFOS - W1R2$Control,
                   Diff3=W2R1$PFOS - W2R1$Control,
                   Diff4=W2R2$PFOS - W2R2$Control,
                   Dense=Density$PFOS/Density$Control)
PFOS$R_mean <- rowMeans(PFOS[,2:5])

View(PFOA)
#final isolate selections:

PFOA[PFOA$R_mean > 0 & PFOA$Dense < 1,]

PFOS[PFOS$R_mean > 0 & PFOS$Dense < 1,]

PFOA[PFOA$R_mean > 0 | PFOA$Dense < 1,]

good_PFOS <- PFOS[PFOS$R_mean > 0 | PFOS$Dense < 1,]
good_PFOS[order(good_PFOS$R_mean, decreasing = T),]

good_PFOS <- PFOS[PFOS$R_mean > 0,]
good_PFOA <- PFOA[PFOA$R_mean > 0,]

intersect(good_PFOA$Strain, good_PFOS$Strain)

#get p-values from data
pA <- apply(scale(PFOA[,-1]), 1, function(x) t.test(x[2:5], mu = 0)$p.value)

PFOA$Pvalue <- pA

library(ggplot2)

PFOA$Pvaluelog <- -log10(PFOA$Pvalue)

plotA <- ggplot(data=PFOA, aes(x=R_mean, y=Pvaluelog)) + geom_point() + 
  geom_text(aes(label=ifelse(R_mean>.5,as.character(Strain),'')),hjust=0,vjust=0, size = 6, position = position_jitter(h=.1,w=.5)) + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + theme_minimal(base_size = 16) + ggtitle("Volcano Plot of Average Difference in Radial Growth Between PFOA and Control") +
  labs(y= "-log10(p-value)", x = "Average difference in radial growth in mm (PFOA - Control)")

plotA

pS <- apply(scale(PFOS[,-1]), 1, function(x) t.test(x[2:5], mu = 0)$p.value)

PFOS$Pvalue <- pS

library(ggplot2)

PFOS$Pvaluelog <- -log10(PFOS$Pvalue)

plotS <- ggplot(data=PFOS, aes(x=R_mean, y=Pvaluelog)) + geom_point() + 
  geom_text(aes(label=ifelse(R_mean>.5,as.character(Strain),'')), hjust=0,vjust=0, size = 6, position = position_jitter(h=.1,w=.5)) + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + theme_minimal(base_size = 16) + ggtitle("Volcano Plot of Average Difference in Radial Growth Between PFOS and Control") +
  labs(y= "-log10(p-value)", x = "Average difference in radial growth in mm (PFOS - Control)")

plotS





# TAG lenght&saturation analysis
library("stringr")
library("ggplot2")
library("plyr")
library("reshape2")
library("pheatmap")

# indicate path to .csv file with average values for each replicate
todo_treatments <- read.delim(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Beatson_MS.txt"), sep = " ")
sample_info <- read.delim(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Beatson_sample_info.txt"), sep = " ")

# MSTUS NORMALIZATION (normalized to sum of peak areas)
MSTUS_norm <- sweep(todo_treatments[,3:ncol(todo_treatments)] ,2,colSums(todo_treatments[,3:ncol(todo_treatments)])/100000000,`/`)

# LOG2 scale (plotting, Limma)
log_MSTUS_norm <- log2(MSTUS_norm)

### PCA
pca <- prcomp(t(na.omit(log_MSTUS_norm)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info$group)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

#PCA by repeat
ggplot(PCAi , aes(PC1, PC2, col= as.character(sample_info$replicate))) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

# probando

inh_1 <- log_MSTUS_norm$Average.DGAT1 - log_MSTUS_norm$Average.DMSO1
inh_2 <- log_MSTUS_norm$Average.DGAT2 - log_MSTUS_norm$Average.DMSO2
inh_3 <- log_MSTUS_norm$Average.DGAT3 - log_MSTUS_norm$Average.DMSO3

MKC8866_1 <- log_MSTUS_norm$Average.MKC1 - log_MSTUS_norm$Average.DMSO1
MKC8866_2 <- log_MSTUS_norm$Average.MKC2 - log_MSTUS_norm$Average.DMSO1
MKC8866_3 <- log_MSTUS_norm$Average.MKC3 - log_MSTUS_norm$Average.DMSO1

cot_1 <- log_MSTUS_norm$Average.co.treat.1 - log_MSTUS_norm$Average.DMSO1
cot_2 <- log_MSTUS_norm$Average.co.treat.2 - log_MSTUS_norm$Average.DMSO1
cot_3 <- log_MSTUS_norm$Average.co.treat.3 - log_MSTUS_norm$Average.DMSO1

final_FCs <- cbind(todo_treatments[,1:2],MKC8866_1, MKC8866_2, MKC8866_3, inh_1, inh_2, inh_3, cot_1, cot_2, cot_3)

# TAGs
TG <- subset(final_FCs, final_FCs$C..Lipid.Class == "TG")

# duplicated entries
TG <- ddply(TG,2,numcolwise(mean))

# anotation. Saturation & lenght
TG$lenght <- sapply(str_extract_all(TG$C..Identification, "[0-9]{2}(?:)"), function(x) sum(as.numeric(x)))
TG$saturation <- sapply(str_extract_all(TG$C..Identification, "(?<=:)[0-9]+"), function(x) sum(as.numeric(x)))
numeros_magicos <- cumsum(rev(table(TG$saturation))) + 0.5

TG_avg <- as.data.frame(cbind( ID = TG$C..Identification, 
                               MKC8866 = rowMeans(TG[,2:4]), 
                               inh = rowMeans(TG[,5:7]), 
                               cot = rowMeans(TG[,8:10]),
                               MKC8866.sd = apply(TG[,2:4], 1, function(x) sd(as.numeric(x))) / sqrt(3), 
                               inh.sd = apply(TG[,5:7], 1, function(x) sd(as.numeric(x))) / sqrt(3), 
                               cot.sd = apply(TG[,8:10], 1, function(x) sd(as.numeric(x))) / sqrt(3)))

# format for plotting. Means, SDs, indiv values
melt1 <- melt(TG, id.vars = 1, measure.vars = 2:10)
melt_avg <- melt(TG_avg, id.vars = 1, measure.vars = 2:4)
melt_sd <- melt(TG_avg, id.vars = 1, measure.vars = 5:7)
melt_sd$variable <- melt_avg$variable

## formatting, lenght & saturation
melt1$lenght <- sapply(str_extract_all(melt1$C..Identification, "[0-9]{2}(?:)"), function(x) sum(as.numeric(x)))
melt1$saturation <- sapply(str_extract_all(melt1$C..Identification, "(?<=:)[0-9]+"), function(x) sum(as.numeric(x)))
melt1$group <- gsub("\\_.*","",melt1$variable)

melt1 <- melt1[order(melt1$lenght, decreasing = TRUE),]
melt1 <- melt1[order(melt1$saturation, decreasing = TRUE),]
melt1$C..Identification <- factor(melt1$C..Identification, levels = unique(melt1$C..Identification))
melt1$variable <- factor(melt1$variable, levels = c("cot", "inh", "MKC8866"))

melt_avg$ID <- factor(melt_avg$ID, levels = levels(melt1$C..Identification))
melt_avg$variable <- factor(melt_avg$variable, levels = levels(melt1$variable))

melt_sd$ID <- factor(melt_sd$ID, levels = levels(melt1$C..Identification))
melt_sd$variable <- factor(melt_sd$variable, levels = levels(melt1$variable))

ggplot() +
  geom_bar(data=melt_avg, aes(fill= `variable`, y= as.numeric(melt_avg$value), x= `ID`),position="dodge", stat="identity", color="black") +
  geom_errorbar(data=melt_sd, aes(fill = `variable`, y=as.numeric(melt_avg$value) ,x=`ID`, ymin=as.numeric(melt_avg$value)-as.numeric(melt_sd$value),ymax=as.numeric(melt_avg$value)+as.numeric(melt_sd$value)), position=position_jitterdodge(jitter.width=0.1, dodge.width=0.9),stat="identity") +
  geom_point(data=melt1 ,pch=21, colour="black", aes(fill = `variable`,  group= `group`, y=melt1$value,x=melt1$C..Identification), position=position_jitterdodge(jitter.width=0.1, dodge.width=0.9)) +
  labs(x = "number of double bonds", y = "log2 FC (MKC8866 vs control)", element_text(face = "bold", angle = 0)) +
  geom_hline(yintercept = 0, size = 1) +
  coord_flip() +
  geom_vline(xintercept = numeros_magicos[-length(numeros_magicos)] ,linetype = 2, size = 1, color = "black") 







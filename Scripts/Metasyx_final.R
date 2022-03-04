### METASYX FINAL
library("ggplot2")
library("sva")
library("colorRamps")
library("pheatmap")
library("stringr")
library("reshape2")
library("plyr")

# indicate path to .csv file with average values for each replicate
# Load data and metadata 
# data: samples in columns, features in rows
todo_treatments <- read.csv(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Metasyx_MD.csv"))
Sample_Info <- read.csv(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Metasyx_sample_info.csv"))
anotacion <- todo_treatments[,c(1,2,10)]

# MSTUS NORMALIZATION (normalized to sum of peak areas)
MSTUS_norm <- sweep(todo_treatments[,12:ncol(todo_treatments)] ,2,colSums(todo_treatments[,12:ncol(todo_treatments)])/100000000,`/`)
row.names(MSTUS_norm) <- todo_treatments$Peak_ID

# LOG2 scale 
log_MSTUS_norm <- log2(MSTUS_norm)

### PCA
pca <- prcomp(t(na.omit(log_MSTUS_norm)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= Sample_Info$Group)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

#PCA by repeat
ggplot(PCAi , aes(PC1, PC2, col= as.factor(Sample_Info$replicate))) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

# repeat 1 looks particularly ugly...
# BATCH EFFECT CORRECTION WITH SVA
mod1 <- model.matrix(~0 + Group, Sample_Info)
mod0 <- model.matrix(~1, Sample_Info)
svobj <- svaseq(as.matrix(MSTUS_norm), mod1, mod0) 

# "Clean" dataset. regress identified sva from data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_count <- cleanY(log_MSTUS_norm, mod1, svobj$sv[,1:3]) # how many sv to use (the more, more transformation)

### PCA after sva
pca <- prcomp(t(na.omit(cleaned_count)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2

variance <- eigs*100/sum(eigs)
cumvar <- cumsum(variance)
eig.veamos <- data.frame(eig = eigs, variance = variance, cumvariance = cumvar)
variances <- barplot(eig.veamos[, 2], names.arg=1:nrow(eig.veamos),main = "Variances", xlab = "Principal Components",ylab = "Percentage of variances", col ="steelblue")
lines(x = variances, eig.veamos[, 2], type="b", pch=19, col = "red")

PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= Sample_Info$Group)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

# fold changes for each replicate individually
ind <- which(Sample_Info$Treatment == "DMSO")

for (i in list(ind)){
  ooo <- cleaned_count[,i+5] - cleaned_count[,i]
}

FCs <- merge(anotacion, ooo, by.x = 1, by.y = 0)

FCs_avg <- as.data.frame(cbind("24h" = rowMeans(FCs[,4:8]), 
                               "48h" = rowMeans(FCs[,9:13]), 
                               "72h" = rowMeans(FCs[,14:18]),
                               "24h.sd" = apply(FCs[,4:8], 1, function(x) sd(as.numeric(x))) / sqrt(5), 
                               "48h.sd" = apply(FCs[,9:13], 1, function(x) sd(as.numeric(x))) / sqrt(5), 
                               "72h.sd" = apply(FCs[,14:18], 1, function(x) sd(as.numeric(x))) / sqrt(5)))

FCs_avg <- cbind(FCs[,1:3], FCs_avg)

##lenght & saturation
FCs$lenght <- sapply(str_extract_all(FCs$Compound_Name, "[0-9]{2}(?=:)"), function(x) sum(as.numeric(x)))
FCs$saturation <- sapply(str_extract_all(FCs$Compound_Name, "(?<=:)[0-9]+"), function(x) sum(as.numeric(x)))

FCs_avg$lenght <- sapply(str_extract_all(FCs$Compound_Name, "[0-9]{2}(?=:)"), function(x) sum(as.numeric(x)))
FCs_avg$saturation <- sapply(str_extract_all(FCs$Compound_Name, "(?<=:)[0-9]+"), function(x) sum(as.numeric(x)))

# HEATMAP OF FOLD CHANGES FOR ALL LIPIDS THROUGH 3 TIME POINTS
annotation_row <- as.data.frame(FCs_avg$Lipid_Class)
row.names(annotation_row) <- FCs_avg$Peak_ID
row.names(FCs_avg) <- FCs_avg$Peak_ID

newCols <- rev(primary.colors(length(unique(FCs_avg$Lipid_Class)))) #reverse just cause I didnt like colors
names(newCols) <- unique(FCs_avg$Lipid_Class)
newCols <- list("FCs_avg$Lipid_Class" = newCols)

my_palette <- colorRampPalette(c("navy","white","red"))(n = 20) # or whatever colors and number
myBreaks <- c(seq(-1, 0, length.out=ceiling(length(my_palette)/2) + 1), 
              seq(1/length(my_palette), 1, length.out=floor(length(my_palette)/2))) # from -1 to 1 log2FC 

pheatmap(na.omit(FCs_avg[,4:6]), 
         cluster_cols=FALSE, 
         color = my_palette, 
         border_color = NA, 
         fontsize_col = 14, 
         show_rownames = FALSE,
         annotation_row = annotation_row,
         annotation_colors = newCols,
         breaks = myBreaks)


# explore each lipid group changes through time
timeplot = function(lipido, color) {
  ooo <- subset(FCs_avg, FCs_avg$Lipid_Class == lipido)
  melt_ooo <- melt(ooo, id.vars = 1, measure.vars = 4:6)
  p <- ggplot(melt_ooo, aes(melt_ooo$variable, melt_ooo$value))
  p + geom_boxplot(outlier.shape = NA, size =2, alpha = 0.5) + 
    geom_point(data=melt_ooo,pch=21, size=3, colour="black", aes(fill= lipido, group= `variable`, y= `value`,x= `variable`), position=position_jitterdodge(jitter.width=0.4, dodge.width=0.9)) +
    theme_classic() +
    theme(text = element_text(size=20), axis.title.y = element_text(size=16), legend.position = "none") +
    ggtitle(lipido) +
    scale_fill_manual(values = color) +
    geom_hline(yintercept = 0 ,linetype = 2, size = 1, color = "black") +
    ylim(-1, 1.2) +
    labs(x="", y = "Log2 FC (MKC8866/DMSO)", element_text(face = "bold", angle = 0,)) + 
    annotate('text', x = c("24h", "48h", "72h"), y = 1.1 ,label = c(paste("P-value==", t.test(ooo$`24h`, mu = 0, alternative = "two.sided")$p.value),
                                                                    paste("P-value==", t.test(ooo$`48h`, mu = 0, alternative = "two.sided")$p.value),
                                                                    paste("P-value==", t.test(ooo$`72h`, mu = 0, alternative = "two.sided")$p.value)) ,parse = TRUE, size=3)
}

# check what lipid groups there are
unique(FCs_avg$Lipid_Class)
# example, lets plot "Triacylglyceride"
timeplot("Triacylglyceride", "red")


# plot all lipid groups. using same colors as in heatmap
for (i in unique(FCs_avg$Lipid_Class)) {
  g <- timeplot(i, newCols$`FCs_avg$Lipid_Class`[i])
  print(g) 
  #ggsave( paste( "time", i, ".pdf"), plot=last_plot(), width=8, height=5)
}

#or plot them all together to evaluate in 1 look
library(gridExtra)
p <- list()
for (i in unique(FCs_avg$Lipid_Class)) {
  p[[i]] <- timeplot(i, newCols$`FCs_avg$Lipid_Class`[i])
}
do.call(grid.arrange,p)

### saturation & elongation indexes 
# plot by saturation
# divide in groups: 0,1,2,3+ saturated bonds
satuplot = function(lipido, color) {
  ooo <- subset(FCs_avg, FCs_avg$Lipid_Class == lipido)
  ooo$saturange <- NA
  ooo$saturange[ooo$saturation == 0] <-  "0"
  ooo$saturange[ooo$saturation == 1] <-  "1"
  ooo$saturange[ooo$saturation == 2] <-  "2"
  ooo$saturange[ooo$saturation >= 3] <-  "3+"
  average_tri <- ddply(ooo,12,numcolwise(mean))
  average_tri$se24 <- c(sd(ooo$`24h`[ooo$saturange == "0"]), sd(ooo$`24h`[ooo$saturange == "1"]), sd(ooo$`24h`[ooo$saturange == "2"]), sd(ooo$`24h`[ooo$saturange == "3+"]))
  average_tri$se48 <- c(sd(ooo$`48h`[ooo$saturange == "0"]), sd(ooo$`48h`[ooo$saturange == "1"]), sd(ooo$`48h`[ooo$saturange == "2"]), sd(ooo$`48h`[ooo$saturange == "3+"]))
  average_tri$se72 <- c(sd(ooo$`72h`[ooo$saturange == "0"]), sd(ooo$`72h`[ooo$saturange == "1"]), sd(ooo$`72h`[ooo$saturange == "2"]), sd(ooo$`72h`[ooo$saturange == "3+"]))
  melt2 <- melt(average_tri, id.vars = 1, measure.vars = 2:4)
  melt3 <- melt(average_tri, id.vars = 1, measure.vars = 10:12)
  melt_todo <- cbind(melt2, melt3[,3])
  melt_dots <- melt(ooo, id.vars = 12, measure.vars = 4:6)
  ggplot() +
    geom_bar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$saturange,ymin=melt_todo$value-melt_todo$`melt3[, 3]`,ymax=melt_todo$value+melt_todo$`melt3[, 3]`), position="dodge", stat="identity", color="black") +
    geom_errorbar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$saturange,ymin=melt_todo$value-melt_todo$`melt3[, 3]`,ymax=melt_todo$value+melt_todo$`melt3[, 3]`), position="dodge",stat="identity") +
    geom_point(data=melt_dots,pch=21, colour="black", aes(fill=melt_dots$variable ,group=melt_dots$variable, y=melt_dots$value,x=melt_dots$saturange), position=position_jitterdodge(jitter.width=0.1, dodge.width=0.9)) +
    scale_fill_manual(values= c("white", "grey50", color)) +
    labs(x = "number of double bonds", y = "log2 FC (MKC8866 vs control)", element_text(face = "bold", angle = 0)) +
    ylim(-1, 1) +
    geom_hline(yintercept = 0, linetype="dashed", size = 2)
}

satuplot("Triacylglyceride"  , "green")

# elongation
# plot by total chain lenght
# divided into ranges matching the 4 quartiles
lenghtplot = function(lipido, color) {
  ooo <- subset(FCs_avg, FCs_avg$Lipid_Class == lipido)
  uuu <- summary(ooo$lenght)
  ooo$saturange <- NA
  ooo$saturange[ooo$lenght <= round(uuu[[2]])] <- paste( "<" , round(uuu[[2]]))
  ooo$saturange[ooo$lenght > round(uuu[[2]])] <- paste(round(uuu[[2]]),"-", round(uuu[[3]]))
  ooo$saturange[ooo$lenght > round(uuu[[3]])] <- paste(round(uuu[[3]]),"-", round(uuu[[5]]))
  ooo$saturange[ooo$lenght >= round(uuu[[5]])] <- paste( ">" , round(uuu[[5]]))
  average_tri <- ddply(ooo,12,numcolwise(mean))
  average_tri$se24 <- c(sd(ooo$`24h`[ooo$saturange == paste( "<" , round(uuu[[2]]))]), sd(ooo$`24h`[ooo$saturange == paste(round(uuu[[2]]),"-", round(uuu[[3]]))]), sd(ooo$`24h`[ooo$saturange == paste(round(uuu[[3]]),"-", round(uuu[[5]]))]), sd(ooo$`24h`[ooo$saturange == paste( ">" , round(uuu[[5]]))]))
  average_tri$se48 <- c(sd(ooo$`48h`[ooo$saturange == paste( "<" , round(uuu[[2]]))]), sd(ooo$`48h`[ooo$saturange == paste(round(uuu[[2]]),"-", round(uuu[[3]]))]), sd(ooo$`48h`[ooo$saturange == paste(round(uuu[[3]]),"-", round(uuu[[5]]))]), sd(ooo$`48h`[ooo$saturange == paste( ">" , round(uuu[[5]]))]))
  average_tri$se72 <- c(sd(ooo$`72h`[ooo$saturange == paste( "<" , round(uuu[[2]]))]), sd(ooo$`72h`[ooo$saturange == paste(round(uuu[[2]]),"-", round(uuu[[3]]))]), sd(ooo$`72h`[ooo$saturange == paste(round(uuu[[3]]),"-", round(uuu[[5]]))]), sd(ooo$`72h`[ooo$saturange == paste( ">" , round(uuu[[5]]))]))
  melt2 <- melt(average_tri, id.vars = 1, measure.vars = 2:4)
  melt3 <- melt(average_tri, id.vars = 1, measure.vars = 10:12)
  melt_todo <- cbind(melt2, melt3[,3])
  melt_dots <- melt(ooo, id.vars = 12, measure.vars = 4:6)
  melt_todo$saturange <- factor(melt_todo$saturange, levels = c(paste( "<" , round(uuu[[2]])), paste(round(uuu[[2]]),"-", round(uuu[[3]])), paste(round(uuu[[3]]),"-", round(uuu[[5]])), paste( ">" , round(uuu[[5]]))))
  ggplot() +
    geom_bar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$saturange,ymin=melt_todo$value-melt_todo$`melt3[, 3]`,ymax=melt_todo$value+melt_todo$`melt3[, 3]`), position="dodge", stat="identity", color="black") +
    geom_errorbar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$saturange,ymin=melt_todo$value-melt_todo$`melt3[, 3]`,ymax=melt_todo$value+melt_todo$`melt3[, 3]`), position="dodge",stat="identity") +
    geom_point(data=melt_dots,pch=21, colour="black", aes(fill=melt_dots$variable ,group=melt_dots$variable, y=melt_dots$value,x=melt_dots$saturange), position=position_jitterdodge(jitter.width=0.1, dodge.width=0.9)) +
    scale_fill_manual(values= c("white", "grey50", color)) +
    labs(x = "number of carbons", y = "log2 FC (MKC8866 vs control)", element_text(face = "bold", angle = 0)) +
    ylim(-1, 1) +
    geom_hline(yintercept = 0, linetype="dashed", size = 2)
}

# example
unique(FCs_avg$Lipid_Class)

satuplot("Sphyngomyelin"  , "blue")
lenghtplot("Triacylglyceride"  , "magenta")

# extra bomboclap
# consolidate duplicated entries (mean)
FCs_o <- ddply(FCs ,2,numcolwise(mean))
FCs_avg_o <- ddply(FCs_avg ,2,numcolwise(mean))

#anotar 
anot <- unique(anotacion[,2:3])

FCs_o <- merge(anot, FCs_o, by.x = 1, by.y = 1)
FCs_avg_o <- merge(anot, FCs_avg_o, by.x = 1, by.y = 1)

# plot
probandou = function(lipido) {
  ooo <- subset(FCs_o, FCs_o$Lipid_Class == lipido )
  uuu <- subset(FCs_avg_o, FCs_avg_o$Lipid_Class == lipido)
  ooo <- ooo[order(ooo$lenght, decreasing = TRUE),]
  ooo <- ooo[order(ooo$saturation, decreasing = TRUE),]
  niveles <- unique(ooo$Compound_Name)
  melt_ooo <- melt(ooo, id.vars = 1, measure.vars = 13:17)
  melt_ooo$Compound_Name <- factor(melt_ooo$Compound_Name, levels = niveles)
  uuu$Compound_Name <- factor(uuu$Compound_Name, levels = niveles)
  numeros_magicos <- cumsum(rev(table(ooo$saturation))) + 0.5
  numero1 <- round(max(melt_ooo$value + 0.1), digits = 1)
  numero2 <- round(min(melt_ooo$value - 0.1), digits = 1)
  ggplot() +
    geom_bar(data=uuu, aes( y= as.numeric(uuu$`72h`), x= uuu$Compound_Name),position="dodge", stat="identity", color="black") +
    geom_errorbar(data=uuu, aes(y=as.numeric(uuu$`72h`) ,x= uuu$Compound_Name, ymin=as.numeric(uuu$`72h`)-as.numeric(uuu$`72h.sd`),ymax=as.numeric(uuu$`72h`)+as.numeric(uuu$`72h.sd`), position= "dodge",stat="identity")) +
    geom_point(data=melt_ooo ,pch=21, fill ="black", colour="black", aes( y=as.numeric(melt_ooo$value), x=melt_ooo$Compound_Name)) +
    labs(x = "number of double bonds", y = "log2 FC (MKC8866 vs control)", element_text(face = "bold", angle = 0)) +
    geom_hline(yintercept = 0, size = 1) +
    ylim(numero2, numero1) +
    coord_flip() +
    geom_vline(xintercept = numeros_magicos[-length(numeros_magicos)] ,linetype = 2, size = 1, color = "black") 
}

probandou("Lyso Phosphatidylcholine")  

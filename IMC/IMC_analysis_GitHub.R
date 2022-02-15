options(stringsAsFactors = F)
library(umap)
library(Hmisc)
library(corrplot)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(fractal)
library(fields)
library(scales)
library(plotrix)
library(raster)
library(gridExtra)
library(ggpubr)
library(plot3D)
# library(rgl)
library(psych)
library(plotly)
library(akima)
library(lattice)
library(ellipse)
library(corrplot)
library(ggrepel)
library(Rtsne)

c25 <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

all_plot <- read.table("all_plot_data.txt", sep="\t", header = T, row.names = F)
all_plot[1:10,1:4]
# check some stuff ####
ind_test <- which(all_plot$ID == "003" & all_plot$DAY == "D00" & all_plot$ROI == "ROI 02")
test <- all_plot[ind_test,]
plot(test$IDO, test$IDO_pred)
test$new_IDO_pred <- ifelse(test$IDO >= 11, 1, 0)
plot(test$IDO, test$new_IDO_pred)
table(test$IDO_pred, test$new_IDO_pred)
glm(test$new_IDO_pred ~ test$IDO)
plot(test$X, 1000-(test$Y), pch=16, col=factor(test$new_IDO_pred))
test$IDO_pred <- test$new_IDO_pred
test$new_IDO_pred <- NULL

all_plot[ind_test,] <- test



# -- >> Figure 5C <<< -----
# count fraction of CD8 cells per total number of cells ####
pp <- unique(all_plot$ID) 
# pp <- pp[-which(pp %in% c("001", "002"))]
all_plot$ID_DAY <- paste0(all_plot$ID, "_", all_plot$DAY)
dd <- unique(all_plot$ID_DAY)
for(ii in 1:length(pp)){
  if(ii == 1){
    cd8_fraction_0 <- c()
    foxp3_fraction_0 <- c()
    IDO_fraction_0 <- c()
    cd8_fraction_15 <- c()
    foxp3_fraction_15 <- c()
    IDO_fraction_15 <- c()
  }
  
  temp <- all_plot[which(all_plot$ID == pp[ii]),]
  temp_0 <- temp[which(temp$DAY == "D00"),]
  temp_15 <- temp[which(temp$DAY == "D15"),]
  cd8_fraction_0 <- c(cd8_fraction_0, sum(temp_0$CD8) / ( sum(temp_0$Tumor) + sum(temp_0$Stroma) ))
  foxp3_fraction_0 <- c(foxp3_fraction_0, sum(temp_0$FOXP3_pred) / ( sum(temp_0$Tumor) + sum(temp_0$Stroma) ))
  IDO_fraction_0 <- c(IDO_fraction_0, sum(temp_0$IDO_pred) / ( sum(temp_0$Tumor) + sum(temp_0$Stroma) ))
  cd8_fraction_15 <- c(cd8_fraction_15, sum(temp_15$CD8) / ( sum(temp_15$Tumor) + sum(temp_15$Stroma) ))
  foxp3_fraction_15 <- c(foxp3_fraction_15, sum(temp_15$FOXP3_pred) / ( sum(temp_15$Tumor) + sum(temp_15$Stroma) ))
  IDO_fraction_15 <- c(IDO_fraction_15, sum(temp_15$IDO_pred) / ( sum(temp_15$Tumor) + sum(temp_15$Stroma) ))
  
  if(ii == length(pp)){
    names(cd8_fraction_0) <- pp
    cd8_fraction_0 <- cd8_fraction_0[order(cd8_fraction_0, decreasing = T)]
    names(cd8_fraction_15) <- pp
    cd8_fraction_15 <- cd8_fraction_15[order(cd8_fraction_15, decreasing = T)]
    names(foxp3_fraction_0) <- pp
    foxp3_fraction_0 <- foxp3_fraction_0[order(foxp3_fraction_0, decreasing = T)]
    names(foxp3_fraction_15) <- pp
    foxp3_fraction_15 <- foxp3_fraction_15[order(foxp3_fraction_15, decreasing = T)]
    names(IDO_fraction_0) <- pp
    IDO_fraction_0 <- IDO_fraction_0[order(IDO_fraction_0, decreasing = T)]
    names(IDO_fraction_15) <- pp
    IDO_fraction_15 <- IDO_fraction_15[order(IDO_fraction_15, decreasing = T)]
  }
}
par(mfrow=c(1,3))
barplot(c(cd8_fraction_0, cd8_fraction_15),
        col=rep(c("blue", "red"), each=7), 
        las=2, main="Ratio CD8+/Tumor cells")
barplot(c(foxp3_fraction_0, foxp3_fraction_15),
        col=rep(c("blue", "red"), each=7), 
        las=2, main="Ratio FOXP3+/Tumor cells")
barplot(c(IDO_fraction_0, IDO_fraction_15),
        col=rep(c("blue", "red"), each=7), 
        las=2, main="Ratio IDO+/Tumor cells")

fraction_all <- data.frame(CD8 = c(cd8_fraction_15[pp],cd8_fraction_0[pp]),
                           FOXP3 = c(foxp3_fraction_15[pp],foxp3_fraction_0[pp]),
                           IDO = c(IDO_fraction_15[pp],IDO_fraction_0[pp]),
                           ID = rep(pp, 2),
                           DAY = rep(c("D15", "D00"), each = 7))
# fraction_all <- fraction_all[order(fraction_all$CD8, decreasing = F),]

# --- >> Figure S8B <<< ----
# make paired boxplots for the fraction ####
paired_fraction <- data.frame(cbind(fraction_all[1:7,], 
                                    fraction_all[8:14,]))
colnames(paired_fraction) <- gsub("\\.1","_D00",colnames(paired_fraction))
colnames(paired_fraction)[1:5] <- paste0(colnames(paired_fraction)[1:5],"_D15")
gg1 <- ggpaired(paired_fraction, cond1 = "CD8_D00", cond2 = "CD8_D15", fill=c("blue", "red")) + 
  ggtitle("CD8") + ylim(0,0.04) + theme(axis.text.x = element_text(angle=90)) +
  stat_compare_means(method="t.test", paired = T)
gg2 <- ggpaired(paired_fraction, cond1 = "IDO_D00", cond2 = "IDO_D15", fill=c("blue", "red")) + 
  ggtitle("IDO")  + ylim(0,0.5) + theme(axis.text.x = element_text(angle=90)) +
  stat_compare_means(method="t.test", paired = T)
gg3 <- ggpaired(paired_fraction, cond1 = "FOXP3_D00", cond2 = "FOXP3_D15", fill=c("blue", "red")) + 
  ggtitle("FOXP3") + ylim(0,0.03) + theme(axis.text.x = element_text(angle=90)) +
  stat_compare_means(method="t.test", paired = T)

grid.arrange(gg1, gg2, gg3, ncol=3)
t.test(paired_fraction$CD8_D15, paired_fraction$CD8_D00, paired = T)
t.test(paired_fraction$CD8_D15, paired_fraction$CD8_D00, paired = T)$stderr
t.test(paired_fraction$IDO_D15, paired_fraction$IDO_D00, paired = T)
t.test(paired_fraction$IDO_D15, paired_fraction$IDO_D00, paired = T)$stderr
t.test(paired_fraction$FOXP3_D15, paired_fraction$FOXP3_D00, paired = T)
t.test(paired_fraction$FOXP3_D15, paired_fraction$FOXP3_D00, paired = T)$stderr


# calculate distances new ####
all_plot$ID_ROI_DAY <- paste0(all_plot$ID, "_", all_plot$ROI, "_", all_plot$DAY)
table(all_plot$ID_ROI_DAY)
uid <- unique(all_plot$ID_ROI)

for(ii in 1:length(uid)){
  if(ii == 1){
    dist_ido_tumor <- c()
    dist_cd8_tumor <- c()
    dist_ido_cd8 <- c()
    
  }
  if(substr(uid[ii], 1, 3) %in% c("001", "002")){
    message(paste(" !!! >> 1708 skipping", ii))
    next
  }
  # create temp data
  temp <- all_plot[which(all_plot$ID_ROI_DAY == uid[ii]),]
  temp_ido <- temp[which(temp$IDO_pred == 1), c(17:18, 27, 29)]
  temp_cd8 <- temp[which(temp$CD8 == 1), c(17:18, 27, 29)]
  temp_tumor <- temp[which(temp$Tumor == 1), c(17:18, 27, 29)]
  
  # calculate distances
  if(length(temp_tumor$X) == 0){
    message(paste(" !!! >> 1720 skipping", ii))
    next
  } else {
    if(length(temp_ido$X) == 0){
      message(paste(" !!! >> 1724 skipping", ii))
      next
    } else {
      rr_ido_tumor <- data.frame(DIST = as.numeric(unlist(rdist(x1 = temp_ido[,1:2],
                                                                x2 = temp_tumor[,1:2],
                                                                compact = T))),
                                 ID = temp_ido[1,3],
                                 DAY = temp_ido[1,4])
      rr_ido_tumor$SCORE <- tryCatch(mean(rr_ido_tumor$DIST[which(rr_ido_tumor$DIST < 100)]) / nrow(rr_ido_tumor), error=function(ee) return(NA))
    }
    if(length(temp_cd8$X) == 0){
      message(paste(" !!! >> 1735 skipping", ii))
      next
    } else {
      if(ii %in% 81:83){
        message(paste("> working with",uid[ii]))
      }
      rr_cd8_tumor <- data.frame(DIST = as.numeric(unlist(rdist(x1 = temp_cd8[,1:2],
                                                                x2 = temp_tumor[,1:2],
                                                                compact = T))),
                                 ID = temp_cd8[1,3],
                                 DAY = temp_cd8[1,4])
      rr_cd8_tumor$SCORE <- tryCatch(mean(rr_cd8_tumor$DIST[which(rr_cd8_tumor$DIST < 100)]) / nrow(rr_cd8_tumor), error=function(ee) return(NA))
      if(ii %in% 81:83){
        message(paste("    > nrow",nrow(rr_cd8_tumor)))
      }
    }
    
  }
  if(length(temp_cd8$X) == 0 & length(temp_ido$X) == 0){
    message(paste(" !!! >> 1764 skipping", ii))
    next
  } else {
    rr_ido_cd8 <- data.frame(DIST = as.numeric(unlist(rdist(x1 = temp_ido[,1:2],
                                                            x2 = temp_cd8[,1:2],
                                                            compact = T))),
                             ID = temp_ido[1,3],
                             DAY = temp_ido[1,4])
    rr_ido_cd8$SCORE <- tryCatch(mean(rr_ido_cd8$DIST[which(rr_ido_cd8$DIST < 100)]) / nrow(rr_ido_cd8), error=function(ee) return(NA))
  }
  
  if(length(temp_cd8$X) == 0 & length(temp_foxp3$X) == 0){
    message(paste(" !!! >> 1786 skipping", ii))
    next
  } else {
    rr_cd8_foxp3 <- data.frame(DIST = as.numeric(unlist(rdist(x1 = temp_cd8[,1:2],
                                                              x2 = temp_foxp3[,1:2],
                                                              compact = T))),
                               ID = temp_cd8[1,3],
                               DAY = temp_cd8[1,4])
    rr_cd8_foxp3$SCORE <- tryCatch(mean(rr_cd8_foxp3$DIST[which(rr_cd8_foxp3$DIST < 100)]) / length(which(rr_cd8_foxp3$DIST < 100)), error=function(ee) return(NA))
  }
  
  # save distances
  ind <- which(rr_ido_tumor$DIST < 100)
  if(length(ind) > 0){
    dist_ido_tumor <- data.frame(rbind(dist_ido_tumor, 
                                       rr_ido_tumor[ind,]))
  }
  ind <- which(rr_cd8_tumor$DIST < 100)
  message(ii)
  if(ii %in% 81:83){
    message(paste("   > less than 100",length(ind)))
  }
  if(length(ind) > 0)
    dist_cd8_tumor <- data.frame(rbind(dist_cd8_tumor, 
                                       rr_cd8_tumor[ind,]))
  
  ind <- which(rr_ido_cd8$DIST < 100)
  if(length(ind > 0))
    dist_ido_cd8 <- data.frame(rbind(dist_ido_cd8, 
                                     rr_ido_cd8[ind,]))

}


head(dist_cd8_tumor)
# --- >>> Figure 5F <<< ----
mm <- dist_ido_tumor[order(dist_ido_tumor$SCORE, decreasing = T),]
mm$ID_DAY <- paste0(mm$ID, "_", mm$DAY)
agg <- aggregate(mm$SCORE, by=list(mm$ID_DAY), FUN=median)
agg <- agg[order(agg$x, decreasing = F),]
agg
table(mm$ID_DAY)
ggplot(mm, aes(x=factor(ID_DAY, levels = rev(agg$Group.1)), y=log2(SCORE), fill=DAY)) + 
  ylab("log2(SCORE)") + xlab("") + theme_bw() +
  theme(axis.text.y = element_text(size=12)) + 
  geom_boxplot(col="black") + 
  coord_flip() + 
  scale_fill_manual(values=c("D00" = "blue", "D15" = "red")) 

# --- >>> Figure 5G <<< ----
mm <- dist_cd8_tumor[order(dist_cd8_tumor$SCORE, decreasing = T),]
mm$ID_DAY <- paste0(mm$ID, "_", mm$DAY)
agg <- aggregate(mm$SCORE, by=list(mm$ID_DAY), FUN=median)
agg <- agg[order(agg$x, decreasing = F),]
agg
table(mm$ID_DAY)
ggplot(mm, aes(x=factor(ID_DAY, levels = rev(agg$Group.1)), y=log2(SCORE), fill=DAY)) + 
  ylab("log2(SCORE)") + xlab("") + theme_bw() +
  theme(axis.text.y = element_text(size=12)) + 
  geom_boxplot(col="black") + 
  coord_flip() + 
  scale_fill_manual(values=c("D00" = "blue", "D15" = "red")) 

# --- >>> Figure S8C <<< ----
mm <- dist_foxp3_tumor[order(dist_foxp3_tumor$SCORE, decreasing = T),]
mm$ID_DAY <- paste0(mm$ID, "_", mm$DAY)
agg <- aggregate(mm$SCORE, by=list(mm$ID_DAY), FUN=median)
agg <- agg[order(agg$x, decreasing = F),]
agg
table(mm$ID_DAY)
ggplot(mm, aes(x=factor(ID_DAY, levels = rev(agg$Group.1)), y=log2(SCORE), fill=DAY)) + 
  ylab("log2(SCORE)") + xlab("") + theme_bw() +
  theme(axis.text.y = element_text(size=12)) + 
  geom_boxplot(col="black") + 
  coord_flip() + 
  scale_fill_manual(values=c("D00" = "blue", "D15" = "red")) 



# --- >>> Figure S9B <<< ----
# check markers abundance and cell types ####
colnames(all_plot)
colSums(table(all_plot$ID, all_plot$DAY))
hist(all_plot$CD8)
plot(sort(all_plot$CD4))
quantile(all_plot$CD3, seq(0,1,length.out=20)) #34
all_plot$CD3_plot <- ifelse(all_plot$CD3 > 34, all_plot$CD3, 0)
quantile(all_plot$CD68, seq(0,1,length.out=20)) #35
all_plot$CD68_plot <- ifelse(all_plot$CD68 > 35, all_plot$CD68, 0)
quantile(all_plot$CD8a, seq(0,1,length.out=20)) #15
all_plot$CD8a_plot <- ifelse(all_plot$CD8a > 15, all_plot$CD8a, 0)
quantile(all_plot$CD4, seq(0,1,length.out=20)) #43
all_plot$CD4_plot <- ifelse(all_plot$CD4 > 43, all_plot$CD4, 0)
quantile(all_plot$Pan.Keratin, seq(0,1,length.out=20)) #5
all_plot$PK_plot <- ifelse(all_plot$Pan.Keratin > 5, all_plot$Pan.Keratin, 0)
quantile(all_plot$Collagen.Type.1, seq(0,1,length.out=20)) #5
all_plot$CL_plot <- ifelse(all_plot$Collagen.Type.1 > 5, all_plot$Collagen.Type.1, 0)
plot(all_plot$CD4,
     all_plot$FOXP3)
all_plot$CellType <- ifelse(all_plot$CD8a_plot > 0 & all_plot$CD3_plot > 0, "CD8a", 
                            ifelse(all_plot$CD4_plot > 0 & all_plot$CD3_plot > 0, "CD4",
                                   ifelse(all_plot$CD68_plot > 0 & all_plot$CD8a_plot == 0 & all_plot$CD4_plot == 0, "Macrophage",
                                          ifelse(all_plot$PK_plot > 0 & all_plot$CL_plot == 0, "PK",
                                                 ifelse(all_plot$CL_plot > 0 & all_plot$PK_plot == 0, "CL","Other")))))
table(all_plot$CellType)
cbind(round(table(ifelse(all_plot$CD8a_plot == 0, 0, 1)[all_plot$DAY == "D00" & all_plot$ID == "026"]) / 79368, 3),
      round(table(ifelse(all_plot$CD8a_plot == 0, 0, 1)[all_plot$DAY == "D15" & all_plot$ID == "026"]) / 314621, 3))


ggplot(all_plot[which(all_plot$ID == "003"),] %>% filter(CellType == "CD8a"), 
       aes(x=log10(CD3 * CD8a * GranzymeB),
           ..scaled..,
           col=DAY,
           fill=DAY)) +
  theme_bw() +
  geom_density(alpha=0.25) +
  scale_color_manual(values = c("D00" = "blue", "D15" = "red")) +
  scale_fill_manual(values = c("D00" = "blue", "D15" = "red"))

ggplot(all_plot[which(all_plot$ID == "003"),] %>% filter(CellType == "CD8a"), 
       aes(x=log10(CD3 * CD8a * CD27),
           ..scaled..,
           col=DAY,
           fill=DAY)) +
  theme_bw() +
  geom_density(alpha=0.25) +
  scale_color_manual(values = c("D00" = "blue", "D15" = "red")) +
  scale_fill_manual(values = c("D00" = "blue", "D15" = "red"))


ggplot(all_plot[which(all_plot$ID == "008"),] %>% filter(CellType == "CD8a"), 
       aes(x=log10(CD3 * CD8a * GranzymeB),
           ..scaled..,
           col=DAY,
           fill=DAY)) +
  theme_bw() +
  geom_density(alpha=0.25) +
  scale_color_manual(values = c("D00" = "blue", "D15" = "red")) +
  scale_fill_manual(values = c("D00" = "blue", "D15" = "red"))

ggplot(all_plot[which(all_plot$ID == "008"),] %>% filter(CellType == "CD8a"), 
       aes(x=log10(CD3 * CD8a * CD27),
           ..scaled..,
           col=DAY,
           fill=DAY)) +
  theme_bw() +
  geom_density(alpha=0.25) +
  scale_color_manual(values = c("D00" = "blue", "D15" = "red")) +
  scale_fill_manual(values = c("D00" = "blue", "D15" = "red"))

head(all_plot[which(all_plot$ID == "008"),] %>% filter(CellType == "CD8a"))
to_compare <- data.frame(DENS = log10(all_plot$CD3 * all_plot$CD8a * all_plot$GranzymeB),
                         DAY = all_plot$DAY,
                         PT = all_plot$ID)
to_compare <- to_compare[-which(is.infinite(to_compare$DENS)),]
head(to_compare)
ks.test(to_compare$DENS[to_compare$DAY == "D00" & to_compare$PT == "008"],
        to_compare$DENS[to_compare$DAY == "D15"& to_compare$PT == "008"])


head(all_plot[which(all_plot$ID == "008"),] %>% filter(CellType == "CD8a"))
to_compare <- data.frame(DENS = log10(all_plot$CD3 * all_plot$CD8a * all_plot$CD27),
                         DAY = all_plot$DAY,
                         PT = all_plot$ID)
to_compare <- to_compare[-which(is.infinite(to_compare$DENS)),]
head(to_compare)
ks.test(to_compare$DENS[to_compare$DAY == "D00" & to_compare$PT == "008"],
        to_compare$DENS[to_compare$DAY == "D15"& to_compare$PT == "008"])

head(all_plot[which(all_plot$ID == "003"),] %>% filter(CellType == "CD8a"))
to_compare <- data.frame(DENS = log10(all_plot$CD3 * all_plot$CD8a * all_plot$GranzymeB),
                         DAY = all_plot$DAY,
                         PT = all_plot$ID)
to_compare <- to_compare[-which(is.infinite(to_compare$DENS)),]
head(to_compare)
ks.test(to_compare$DENS[to_compare$DAY == "D00" & to_compare$PT == "003"],
        to_compare$DENS[to_compare$DAY == "D15"& to_compare$PT == "003"])


head(all_plot[which(all_plot$ID == "003"),] %>% filter(CellType == "CD8a"))
to_compare <- data.frame(DENS = log10(all_plot$CD3 * all_plot$CD8a * all_plot$CD27),
                         DAY = all_plot$DAY,
                         PT = all_plot$ID)
to_compare <- to_compare[-which(is.infinite(to_compare$DENS)),]
head(to_compare)
ks.test(to_compare$DENS[to_compare$DAY == "D00" & to_compare$PT == "003"],
        to_compare$DENS[to_compare$DAY == "D15"& to_compare$PT == "003"])





options(stringsAsFactors = F)
library(amap)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(fields)
library(scales)
library(plotrix)
library(gridExtra)
library(Rtsne)
library(flowCore)
library(matrixStats)
library(citrus)
library(ggpubr)
library(ConsensusClusterPlus)
library(ggridges)
library(FlowSOM)
library(ggrepel)
library(ggpp)
library(ggpmisc)
library(dplyr)

c25 <- c("dodgerblue2",
         "#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2",
         "#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")


setwd("CyTOF")

files_input <- list.files(path = "/Volumes/CYTOF-IMC/PIOdunsi/CITN/TIL_122920/FCS files/", full.names = T, pattern = "*.fcs")[-1]
flow_data <- read.flowSet(files = files_input, 
                          transformation = FALSE,
                          truncate_max_range = F,
                          verbose=T)

panel_fcs <- flowCore::pData(flowCore::parameters(flow_data[[1]]))
panel_fcs
panel_fcs$desc <- as.character(sapply(panel_fcs$desc, function(xx) strsplit(xx,"_")[[1]][2]))
panel_fcs

# transform data archsin 5 ####
fcs_transformed <- fsApply(flow_data, function(x){
  ex <- exprs(x)
  message(paste("number of cells:", nrow(ex), "for", identifier(x)))
  # archsin transform with cofactor 5
  cofactor = 5
  ex <- asinh(ex / cofactor)
  exprs(x) <- ex
  x
})


expr <- fsApply(fcs_transformed, function(xx){
  if(nrow(exprs(xx)) > 0){
    ex <- data.frame(exprs(xx))
    ex$File <- identifier(xx)
    ex$PatientID <- gsub("TIL","",strsplit(strsplit(identifier(xx), "_")[[1]][3], " ")[[1]][2])
    ex$DAY <- sapply(ex$PatientID, function(xx) strsplit(xx, "-")[[1]][4])
    message(paste("Patient", unique(ex$PatientID), "-", nrow(ex), "cells"))
    return(ex)
  }
})
expr <- do.call("rbind", expr)
colnames(expr) <- c(panel_fcs$desc, "File", "PatientID", "DAY")
expr <- expr[,-which(is.na(colnames(expr)))]
table(expr$DAY)
colnames(expr)
# alternatively read in results
# expr_sub <- read.table("CITN_CyTOF_GitHub_processed.txt", header = T, row.names = 1, sep='\t')

expr_sub <- expr[,-which(colnames(expr) %in% c("IL-17a", "IDO1", "IL-4", "Arg1", "Foxp3", "Eomes"))] # they have aspecific signal
expr_sub <- expr_sub[c(sample(which(expr_sub$DAY == "D000"), 1625, replace = F),
                   sample(which(expr_sub$DAY == "D015"), 1625, replace = F)),]
table(expr_sub$DAY)

#### read in processed files ####
rng <- colQuantiles(as.matrix(expr_sub[,1:30]), probs = c(0.01, 0.99))
expr01 <- t((t(expr_sub[,1:30]) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

colnames(expr_sub)
tSNE_new <- Rtsne(expr_sub[,1:30], 
                  verbose=T,
                  check_duplicates = FALSE)

# note - no setseed() was used, embedding will most likely differ as well as clustering with Kmeans #
tsne_dr <- data.frame(tSNE1=tSNE_new$Y[,1],
                      tSNE2=tSNE_new$Y[,2],
                      expr01,
                      PatientID = expr_sub$PatientID,
                      DAY=expr_sub$DAY)
tsne_dr$CLUSTERS <- kmeans(expr_sub[,1:30], centers = 12, iter.max = 1000)$cluster

# --- >> Figure 6A <<< -----
ggplot(tsne_dr, aes(x=tSNE1, y=tSNE2, col=factor(CLUSTERS))) + 
  geom_point(size=1.2) +
  theme_bw() +
  scale_color_manual(values = c25) +
  guides(colour = guide_legend(title = "CLUSTER", override.aes = list(size=6), ncol = 2), legend = element_text(size = 20)) 


# --- >> Figure 6B <<< ----
agg <- aggregate(expr_sub[,1:30], by=list(tsne_dr$CLUSTER), FUN=median)
rws <- colnames(agg)[-1]
cls <- agg$Group.1
agg <- agg[,-1]
agg <- t(apply(agg,2,scale))
rownames(agg) <- rws
colnames(agg) <- cls
Heatmap(agg, 
        clustering_distance_rows = "pearson", 
        cluster_columns = F,
        clustering_distance_columns = "pearson",
        col = colorRamp2(c(-2, -0.5,0.5, 2), c("#437EA4","white","white","#E41A1C")), 
        name = "Intensity")

# --- >> Figure 6C <<< ----
ggplot(tsne_dr, aes(x=tSNE1, y=tSNE2, col=DAY)) + 
  geom_point(size=1) +
  theme_bw() +
  scale_color_manual(values = c("D000"="blue", "D015"="red")) +
  guides(colour = guide_legend(override.aes = list(size=6)), legend = element_text(size = 20)) +
facet_wrap(~ DAY)

# --- >> Figure 6D <<< ----
# get number of cells per each cluster and compare d0 and d15
n_cells_cluster <- as.data.frame.matrix(table(tsne_dr$DAY, tsne_dr$CLUSTER))
n_cells_cluster <- melt(t(n_cells_cluster))
colnames(n_cells_cluster) <- c("CLUSTER", "DAY", "nCells")
head(n_cells_cluster)
ggplot(n_cells_cluster, aes(x=factor(CLUSTER), y=nCells, fill=DAY)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("D000" = "blue", "D015" = "red")) +
  xlab("CLUSTER") +
  theme(axis.text.x = element_text(size=15), 
        axis.title.x = element_text(size = 15))


# --- >> Figure S6E <<< ----
# plot all tSNEs
colnames(tsne_dr)
for(ii in 3:33){
  if(ii == 33){
    gg <- ggplot(tsne_dr, aes(x=tSNE1, y=tSNE2, col=DAY)) + 
      geom_point(size=0.75) +
      theme_bw() +
      scale_color_manual(values = c("D000"="blue", "D015"="red")) +
      facet_wrap(~ DAY)
    print(gg)
    
    gg <- ggplot(tsne_dr, aes(x=tSNE1, y=tSNE2, col=factor(CLUSTER))) + 
      geom_point(size=0.5) +
      theme_bw() +
      scale_color_manual(values = c("1"=c25[1], "2"=c25[2], "3"=c25[3], "4"=c25[4], 
                                    "5"=c25[5], "6"=c25[6], "7"=c25[7],"8"=c25[8],
                                    "9"=c25[9], "10"=c25[10], "11"=c25[11], "12"=c25[12])) +
      guides(colour = guide_legend(override.aes = list(size=6)), legend = element_text(size = 20)) +
      facet_wrap(~ DAY)
    print(gg)
  } else {
    gg <- ggplot(tsne_dr, aes(x=tSNE1, y=tSNE2, col=tsne_dr[,ii])) + 
      geom_point(size=0.45) +
      scale_color_gradientn(name = colnames(tsne_dr)[ii],
                            colors = c("purple", alpha("purple",0.5), alpha("green", 0.5), "green")) +
      theme_bw() + 
      theme(plot.tag = element_text(size=12)) +
      facet_wrap(~ DAY)
    print(gg)
  }
  message(paste("Done for", colnames(tsne_dr)[ii]))
}

# --- >> Figure S6A<<< ----
# calculat expression of IFNg in CD8, CD4 and NK ####
ex <- expr_sub
ex[1:10,1:2]
colnames(ex)
ex[,1:30] <- apply(ex[,1:30], 2, function(xx) log1p(xx))
ex <- data.frame(merge(tsne_dr[,c(1,2,35)], ex, by="row.names"))

ex$CD8_T <- ifelse(ex$CD8a > 1 & ex$CD3 > 1, "CD8_T", "")
ex$CD4_T <- ifelse(ex$CD4 > 1 & ex$CD3 > 1 & ex$CD8a < 1, "CD4_T", "")
ex$NK <- ifelse(ex$CD56 > 1 & ex$CD3 < 1 & ex$CD326 < 1, "NK", "")
ex$CellType <- apply(ex[,38:40], 1, function(xx) paste0(xx, collapse = ""))
table(ex$CellType)

# calculate number of cells per time point
tt <- round((table(ex$CellType, ex$DAY) / nrow(ex)) * 100, 3)
d0_d15_cd8cd4cd56_cytotoxic_df <- data.frame(DAY = c(rep(c("D000", "D015"), 3)),
                                             nCells = c(13, 90, 7, 33, 5, 33),
                                             CellType = c(rep(c("CD8 T", "CD4 T", "NK"), each = 2)))
ggplot(d0_d15_cd8cd4cd56_cytotoxic_df, aes(x=DAY, y=nCells, fill=DAY)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_bw() +
  facet_wrap(~ factor(CellType, levels = c("CD8 T", "CD4 T", "NK"))) +
  scale_fill_manual(values = c("D000" = "blue", "D015"= "red")) +
  guides(name="DAY")

d0_d15_cd8cd4cd56_cytotoxic_df

prop_cd8 <- prop.test(c(13,90), c(108, 225),conf.level = .95, alternative = "less")
prop_cd4 <- prop.test(c(7,33), c(111, 188),conf.level = .95, alternative = "less")
prop_nk <- prop.test(c(5,33), c(247, 810),conf.level = .95, alternative = "less")


ggplot(ex %>% 
         filter(CellType != ""), 
       aes(x=DAY, y=IFNg, fill=DAY)) +
  geom_boxplot(col="black") +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values = c("D000" = "blue", "D015" = "red")) +
  facet_wrap(~ factor(CellType, levels=c("CD8_T", "CD4_T", "NK")), ncol=3) +
  stat_compare_means(method = "t.test")



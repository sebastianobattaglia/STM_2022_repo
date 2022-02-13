options(stringsAsFactors = FALSE)
setwd("RNASeq")
library(AMR)
library(gplots)
library(tximport)
library(readr)
library(limma)
library(edgeR)
library(biomaRt)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(reshape)
library(qusage)
library(circlize)
library(ggpubr)
library(clusterProfiler)
library(gage)
library(GSVA)
library(DOSE)
library(enrichplot)
library(RColorBrewer)
library(amap)
library(dendextend)
library(gridExtra)
library(survival)
library(plot3D)
library(matrixStats)
library(survminer)
library(matrixStats)


mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = 'www.ensembl.org')
t2g <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name","transcript_biotype","entrezgene"), mart = mart)
t2g <- t2g[which(t2g$transcript_biotype == "protein_coding"),]
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name, entrezid = entrezgene)
head(t2g)


# genes in the IDO pathway ####
genes <- c("IDO1","TDO2","AFMID","KMO","KYAT1","KYAT3","KYNU","HAAO","QPRT")
d2 <- read.table("Input_files/RNASeq_normalized.txt",header = T, row.names = 1)
rna <- d2;dim(rna);head(rna)
rna <- rna[,c(9:16,1:8)]
dim(rna)
hist(data.matrix(rna), main="", xlab="")

# update gene names to current 
rownames(rna)[which(rownames(rna) == "SDCCAG3")] <- "ENTR1"
rownames(rna)[which(rownames(rna) == "ATP5E")] <- "ATP5F1E"
rownames(rna)[which(rownames(rna) == "USMG5")] <- "ATP5MD"
rownames(rna)[which(rownames(rna) == "CYR61")] <- "CCN1"
rownames(rna)[which(rownames(rna) == "C5orf42")] <- "CPLANE1"
rownames(rna)[which(rownames(rna) == "TMEM110")] <- "STIMATE"
rownames(rna)[which(rownames(rna) == "SLC35E2")] <- "SLC35E2A"
rownames(rna)[which(rownames(rna) == "C11orf70")] <- "CFAP300"
rownames(rna)[which(rownames(rna) == "FAM46C")] <- "TENT5C"
rownames(rna)[which(rownames(rna) == "C17orf105")] <- "CFAP97D1"
rownames(rna)[which(rownames(rna) == "FAM109A")] <- "PHETA1"
rownames(rna)[which(rownames(rna) == "NOTCH2NL")] <- "NOTCH2NLA"


# hist(rna)
# ---- >> Figure 2A << ----
pp$variance
p <- prcomp(t(rna))
summary(p)
screeplot(p, npcs = min(10, length(p$sdev)), type = "line")
ss <- unlist(sapply(colnames(rna),function(xx) strsplit(xx,"_")[[1]][1]))
pat <- unlist(sapply(colnames(rna),function(xx) strsplit(xx,"_")[[1]][2]))
p1 <- data.frame(PC1=p$x[,1],
                 PC2=p$x[,2],
                 SAMPLE=ss,
                 PATIENT=pat,
                 GROUP=ifelse(pat == "22.022", "22",
                              ifelse(pat == "22.026", "26",
                                     ifelse(pat == "22.029", "29",
                                            ifelse(pat == "58.001", "01",
                                                   ifelse(pat == "58.006", "06",
                                                          ifelse(pat == "58.007", "07",
                                                                 ifelse(pat == "58.008", "08",
                                                                        ifelse(pat == "58.003", "03","")))))))))

ggplot(p1, aes(x=as.numeric(PC1), y=as.numeric(PC2), col=SAMPLE, group=GROUP)) + 
  geom_line(col="gray", size=0.6) +
  geom_point(size=6) + 
  theme_bw() + xlab("PC1 - 20% Variance") + ylab("PC2 = 12% Variance") +
  geom_label_repel(show.legend = F, aes(label = PATIENT), col="black", force = 4) +
  theme(axis.title = element_text(size=18),legend.text = element_text(size=15), legend.title = element_text(size=18)) +
  scale_color_manual(values = c(Pre = "blue", Post = "red"))  

# prepare design matrix for DEG analysis analysis ####
patient <- sapply(colnames(rna), function(xx) strsplit(xx, "_")[[1]][[2]])
group <- sapply(colnames(rna), function(xx) strsplit(xx, "_")[[1]][[1]])
dm <- model.matrix(~patient+group)
rownames(dm) <- gsub("X","",colnames(rna))
dm[,"groupPre"] <- c(rep(1,8), rep(0,8))

# do DEG analysis with limma ####
fit2 <- lmFit(rna, dm, method = "robust")
fit2 <- eBayes(fit2, trend = T)
res <- topTable(fit2, coef = 9, n=Inf); sum(res2$adj.P.Val <= 0.05)
head(res)

# --- >> Figure 2B << ----
res_sig <- res[which(res$adj.P.Val <= 0.01 & abs(res$logFC) >= 0.5),];res_sig <- res_sig[order(res_sig$adj.P.Val,decreasing = F),];head(res_sig)
heat <- rna[which(rownames(rna) %in% rownames(res_sig)),]
rownames(heat) <- res_sig$ext_gene
dim(heat)
heat2 <- t(apply(heat,1,scale));rownames(heat2) <- rownames(heat); colnames(heat2) <- gsub("X","",colnames(heat))
head(heat2)
Heatmap(heat2,
        show_row_names = F, 
        clustering_distance_rows = "euclidean", 
        clustering_distance_columns = "euclidean", 
        col=colorRamp2(colors = c(alpha("cyan", 0.5), "white", alpha("red", 0.8)), breaks = c(-1,0,1)), 
        name = "Expression")

# -- >> FIGURE 2C << ----
# check Tryptophan pathways ####
path <- read.gmt("Input_files/serotonin_tryptophan_nicotinamide_pathways.txt") 
ath_tryptophan <- data.frame(ont=path$ont[grep("tryp", path$ont, ignore.case = T)],
                              gene=path$gene[grep("tryp", path$ont, ignore.case = T)])
path_tryptophan <- path_tryptophan[-which(duplicated(path_tryptophan$gene)),]
path_tryptophan <- data.frame(ont="TRIPTOPHAN_PATHWAY",
                              gene=path_tryptophan$gene)

set.seed(13333)
gsea_path <- GSEA(ranked_genes, 
                  TERM2GENE = path_tryptophan, 
                  pvalueCutoff = 1);gsea_path@result
gseaplot2(gsea_path,geneSetID = 1, pvalue_table = F, base_size = 15, rel_heights = c(1, 0.25, 0.75))


# --- >> Figure 2D << ----
# evaluate IDO genes 
genes <- c("IDO1", "IDO2","TDO2","AFMID","KMO","KYAT1","KYAT3","KYNU","HAAO","QPRT", "NADK", "AOC1", "ALDH1B1", "ALDH7A1", "OGDHL", "IL4I1", "CD38", "CAT", "CYP1B1", "NUDT12", "NT5C1B", "MAOA", "NAMPT", "AOX1")
res_ido <- merge(data.frame(genes), res2, by.x="genes", by.y="row.names")
res_ido <- res_ido[order(res_ido$logFC, decreasing = T),]; head(res_ido)

ind_FC <- t(apply(rna[which(rownames(rna) %in% genes),],1,function(xx){
  x001 <- as.numeric(xx[which(colnames(rna) == "Post_58.001")])-as.numeric(xx[which(colnames(rna) == "Pre_58.001")])
  x003 <- as.numeric(xx[which(colnames(rna) == "Post_58.003")])-as.numeric(xx[which(colnames(rna) == "Pre_58.003")])
  x006 <- as.numeric(xx[which(colnames(rna) == "Post_58.006")])-as.numeric(xx[which(colnames(rna) == "Pre_58.006")])
  x007 <- as.numeric(xx[which(colnames(rna) == "Post_58.007")])-as.numeric(xx[which(colnames(rna) == "Pre_58.007")])
  x008 <- as.numeric(xx[which(colnames(rna) == "Post_58.008")])-as.numeric(xx[which(colnames(rna) == "Pre_58.008")])
  x022 <- as.numeric(xx[which(colnames(rna) == "Post_22.022")])-as.numeric(xx[which(colnames(rna) == "Pre_22.022")])
  x026 <- as.numeric(xx[which(colnames(rna) == "Post_22.026")])-as.numeric(xx[which(colnames(rna) == "Pre_22.026")])
  x029 <- as.numeric(xx[which(colnames(rna) == "Post_22.029")])-as.numeric(xx[which(colnames(rna) == "Pre_22.029")])
  all <- cbind(x001,x003,x006,x007,x008,x022,x026,x029)
  all
})); ind_FC <- data.frame(ind_FC); colnames(ind_FC) <- c("x001","x003","x006","x007","x008","x022","x026","x029")
ind_FC$Mean <- rowMeans(ind_FC)
ind_FC <- ind_FC[order(ind_FC$Mean, decreasing = T),]
order_genes  <- rownames(ind_FC)
res2[which(rownames(res2) %in% order_genes),]
significant <- ifelse(rownames(ind_FC) %in% rownames(res2)[which(res2$P.Value <= 0.1 & res2$logFC > 0)], "UP", 
                      ifelse(rownames(ind_FC) %in% rownames(res2)[which(res2$P.Value <= 0.1 & res2$logFC < 0)], "DOWN", "NO")); names(significant) <- rownames(ind_FC)
ind_FC  <- melt(t(ind_FC[,1:7]), id.vars = rownames(ind_FC))
ind_FC$SIGNIFICANT <- significant[match(ind_FC$X2, names(significant))]

ggplot(ind_FC, aes(y=value, x=factor(X2, levels = rev(order_genes)), fill=SIGNIFICANT)) + geom_hline(yintercept = 0, linetype="dashed", col="blue") + geom_boxplot(col="darkblue") + 
  xlab("IDO Genes") + ylab("Log2FC") + ylim(-7,+7) + coord_flip() + theme_bw() +
  scale_fill_manual("Groups\np.value < 0.1", values = c("UP" = scales::alpha("red",0.5), "DOWN" = scales::alpha("blue",0.3), "NO" = scales::alpha("lightgray",0.5))) +
  theme(axis.text = element_text(size=15),axis.title=element_text(size=14))



# --- >> Figure S5A << ----
# do REACTOME PA analysis
rna_gsva <- rna[which(rownames(rna) %in% rownames(res_sig)),]
rownames(rna_gsva) <- gsub("\\.","-",rownames(rna_gsva))
rna_gsva2 <- merge(rna_gsva, t2g, by.x="row.names", by.y="Gene.name"); rna_gsva2 <- rna_gsva2[-which(duplicated(rna_gsva2$Row.names)),];head(rna_gsva2)
rownames(res_sig)[-which(gsub("\\.","-",rownames(res_sig)) %in% rna_gsva2$Row.names)]
rownames(rna_gsva2) <- rna_gsva2$NCBI.gene.ID
rna_gsva2 <- rna_gsva2[,2:17]

genes_for_ppp <- rownames(rna_gsva2)
ppp <- ReactomePA::enrichPathway(genes_for_ppp, organism = "human", 
                                 pvalueCutoff = 0.25, qvalueCutoff = 0.25, 
                                 universe = as.character(unique(t2g$NCBI.gene.ID)))
top <- data.frame(ppp[1:30,2:6])
top$GeneRatio2 <- unlist(sapply(top$GeneRatio, function(xx) as.numeric(strsplit(xx, "/")[[1]][[1]]) / as.numeric(strsplit(xx, "/")[[1]][[2]]) ))
top <- top[order(top$GeneRatio2, decreasing = T ),]
top <- top[1:24,]
ggplot(top, aes(x=factor(Description, levels=rev(as.character(top$Description))), y=GeneRatio2, color=-log10(p.adjust))) + ylim(0, 0.07) +
  geom_point(size=5.5) + coord_flip() + theme_bw() + scale_color_gradient(low="blue", high="red") + xlab("Ranked Pathways") + ylab("GeneRatio") +
  theme(axis.text = element_text(size=12), 
        axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), 
        axis.text.y=element_blank())


# -- >> Figure S5B << ----
hlas <- res[grep("^HLA", rownames(res)),]
hlas$Significant <- ifelse(hlas$logFC < 0 & hlas$P.Value < 0.1, "DOWN", 
                           ifelse(hlas$logFC > 0 & hlas$P.Val < 0.1, "UP", "NO"))
hlas  <- hlas[order(hlas$logFC, decreasing = T), ]
ggplot(hlas, aes(x=factor(rownames(hlas), levels = rownames(hlas)), y=logFC, fill=Significant, col=Significant)) + geom_bar(stat = "identity") +
  scale_fill_manual("Groups\np < 0.1", values = c("UP" = scales::alpha("red",0.5), "DOWN" = scales::alpha("blue",0.3), "NO" = scales::alpha("lightgray",0.5))) + 
  scale_color_manual("Groups\np < 0.1", values = c("UP" = "red", "DOWN" = "blue", "NO" = "gray")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 15), 
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size=18)) + 
  xlab("Genes")

hist_genes <- res[c(grep("^hdac",rownames(res), ignore.case = T),
                    grep("^kmt",rownames(res), ignore.case = T),
                    grep("^kdm",rownames(res), ignore.case = T),
                    grep("^prmt",rownames(res), ignore.case = T),
                    grep("^sirt",rownames(res), ignore.case = T),
                    grep("^suv39h",rownames(res), ignore.case = T),
                    grep("ezh2",rownames(res), ignore.case = T)),]
hist_genes$Significant <- ifelse(hist_genes$logFC < 0 & hist_genes$P.Value < 0.1, "DOWN", 
                           ifelse(hist_genes$logFC > 0 & hist_genes$P.Val < 0.1, "UP", "NO"))
hist_genes  <- hist_genes[order(hist_genes$logFC, decreasing = T), ]
ggplot(hist_genes, aes(x=factor(rownames(hist_genes), levels = rownames(hist_genes)), y=logFC, fill=Significant, col=Significant)) + geom_bar(stat = "identity") +
  scale_fill_manual("Groups\np < 0.1", values = c("UP" = scales::alpha("red",0.5), "DOWN" = scales::alpha("blue",0.3), "NO" = scales::alpha("lightgray",0.5))) + 
  scale_color_manual("Groups\np < 0.1", values = c("UP" = "red", "DOWN" = "blue", "NO" = "gray")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 15), 
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size=18)) + 
  xlab("Genes")

# --- >> Figure S5C << ----
rna_ifn <- merge(ifn_pathways, rna, by.x="gene", by.y="row.names")
rna_ifn <- rna_ifn[-which(duplicated(rna_ifn$gene)),]
tt <- as.data.frame.matrix(table(as.character(rna_ifn$gene), as.character(rna_ifn$term)))
rownames(tt)
rna_ifn <- merge(rna_ifn, tt, by.x="gene", by.y="row.names")
rna_ifn$PATHWAY <- ifelse(rna_ifn$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING == 1 & rna_ifn$REACTOME_INTERFERON_GAMMA_SIGNALING == 1, "BOTH", 
                          ifelse(rna_ifn$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING == 1 & rna_ifn$REACTOME_INTERFERON_GAMMA_SIGNALING == 0, "ALPHA_BETA",
                                 ifelse(rna_ifn$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING == 0 & rna_ifn$REACTOME_INTERFERON_GAMMA_SIGNALING == 1, "GAMMA", "IFN")))

res_ind_FC <- t(apply(rna_ifn[,3:18],1,function(xx){
  x022 <- xx[1]-xx[9]
  x026 <- xx[2]-xx[10]
  x029 <- xx[3]-xx[11]
  x001 <- xx[4]-xx[12]
  x006 <- xx[5]-xx[13]
  x007 <- xx[6]-xx[14]
  x008 <- xx[7]-xx[15]
  x003 <- xx[8]-xx[16]
  all <- cbind(x001,x003,x006,x007,x008,x022,x026,x029)
  # all$mean <- rowMeans(all)
  all
})); res_ind_FC <- data.frame(res_ind_FC); res_ind_FC$Mean <- rowMeans(res_ind_FC); rownames(res_ind_FC) <- rna_ifn$gene
colnames(res_ind_FC) <- c(sort(unique(as.numeric(sapply(colnames(rna_ifn)[3:10], function(xx) strsplit(xx, "\\.")[[1]][2])))), "Mean")
ind_sig <- which(rownames(res_ind_FC) %in% rownames(genes_ranked)[which(genes_ranked$adj.P.Val <= 0.05 & abs(genes_ranked$logFC) >= 0.25)])

mm <- melt(t(res_ind_FC[ind_sig,1:8])); head(mm)
mm$COLOR <- rna_ifn$PATHWAY[match(mm$X2, rna_ifn$gene)]
ggplot(mm, aes(x=factor(X2, levels = rownames(res_ind_FC)[order(res_ind_FC$Mean, decreasing = F)]), y=value, fill=COLOR)) + 
  geom_hline(yintercept = 0, linetype="dashed", col="blue") + ylab("Log2FC") + xlab("GENE") +
  geom_boxplot(col="black") + coord_flip() + 
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15)) + theme_bw()



# --- >> Figure 2E << ---- 
genes_ranked_2 <- ranked_genes
names(genes_ranked_2) <- rownames(genes_ranked)[order(genes_ranked$logFC, decreasing = T)]
df_ifn <- ifn_pathways
gsea_path <- GSEA(genes_ranked_2, 
                  TERM2GENE = ifn_pathways,
                  nPerm = 1000,
                  minGSSize = 5, 
                  pvalueCutoff = 1);gsea_path@result[,2:8]
gseaplot2(gsea_path,geneSetID = 1, pvalue_table = F, base_size = 15, rel_heights = c(1, 0.25, 0.75))

# --- >> Figure 2F << ---- 
ranked_genes <- res$logFC
names(ranked_genes) <- rownames(res)
ranked_genes <- ranked_genes[order(ranked_genes, decreasing = T)]

gmt_custom <- qusage::read.gmt("Input_files/c2.cp.reactome.v7.2.symbols.gmt")
gmt_custom <- qusage::read.gmt("Input_files/h.all.v7.2.symbols.gmt")
gmt_custom <- data.frame(ont = rep(names(gmt_custom), as.numeric(lapply(gmt_custom, length))),
                         gene = unlist(gmt_custom))
gmt_custom %>%
  filter(grepl("ifn", ont, ignore.case = T)) %>%
  unique(ont)
gsea_custom <- clusterProfiler::GSEA(ranked_genes, 
                                     TERM2GENE = gmt_custom, 
                                     pvalueCutoff = 0.05); nrow(gsea_custom@result)
View(gsea_custom@result)
sets <- c("HALLMARK_WNT_BETA_CATENIN_SIGNALING",
          "HALLMARK_TGF_BETA_SIGNALING",
          "HALLMARK_MYC_TARGETS_V2")

set_IDs <- which(gsea_custom@result$Description %in% sets)
gseaplot2(gsea_custom, 
          color = c("brown", "purple", "black"),base_size = 15,
          geneSetID = set_IDs)

# -- >> Figure S5E << ----
# confirm with Nanostring.
nano <- data.frame(t(read.table("Nanostring/Nanostring_normalized_data.txt", sep='\t', header = F, row.names = 1)))[,-1]
nano[1:10,1:10]
nano$Patient <- sapply(nano$ID, function(xx) strsplit(strsplit(xx, ":")[[1]][1], "_")[[1]][1])
nano$Group <- sapply(nano$ID, function(xx) strsplit(strsplit(xx, ":")[[1]][1], "_")[[1]][2])
hist(log2(data.matrix(nano[2:nrow(nano), 3:786])))
nano[2:nrow(nano), 3:786] <- log2(data.matrix(nano[2:nrow(nano), 3:786]))
# keep same genes as above
nano_ifn <- nano[-1,c(which(colnames(nano) %in% rownames(res_ind_FC)[ind_sig]), 787, 788)]
nano_ind_FC <- t(apply(data.matrix(nano_ifn[,-c(19:20)]),2,function(xx){
  x024 <- xx[1]-xx[13]
  x025 <- xx[2]-xx[22]
  x001 <- xx[3]-xx[18]
  x003 <- xx[4]-xx[19]
  x004 <- xx[5]-xx[12]
  x005 <- xx[6]-xx[20]
  x023 <- xx[7]-xx[21]
  x007 <- xx[8]-xx[23]
  x026 <- xx[9]-xx[24]
  x002 <- xx[10]-xx[11]
  x029 <- xx[15]-xx[14]
  x008 <- xx[17]-xx[16]
  all <- cbind(x024,x025,x001,x003,x004,x005,x023,x007,x026,x002,x029,x008)
  # all$mean <- rowMeans(all)
  all
})); nano_ind_FC <- data.frame(nano_ind_FC); nano_ind_FC$Mean <- rowMeans(nano_ind_FC)
nano_ind_FC$SIGN <- round(apply(nano_ifn[,-c(19,20)], 2, function(xx) t.test(as.numeric(xx) ~ nano_ifn$Group)$p.value), 3)
nano_ind_FC$SIGN <- ifelse(nano_ind_FC$SIGN <= 0.1, "YES", "NO")

mm <- melt(t(nano_ind_FC[,c(1:12)])); head(mm)
mm$SIGN <- nano_ind_FC$SIGN[match(mm$X2, rownames(nano_ind_FC))]
ggplot(mm, aes(x=factor(X2, levels = rownames(nano_ind_FC)[order(nano_ind_FC$Mean, decreasing = F)]), y=value, fill=SIGN)) + 
  geom_hline(yintercept = 0, linetype="dashed", col="maroon") + ylab("Log2FC") + xlab("GENE") +
  geom_boxplot(col="black") + coord_flip() + theme_bw() +
  scale_fill_manual(values = c(NO = "lightgray", YES = alpha("green", 0.4)), name="p < 0.1") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15))



# --- >>> Figure S9 <<< ----
# run GSVA and take info at the patient level ####
res_gsva <- gsva(data.matrix(rna), 
                 qusage::read.gmt("~/OneDrive - Roswell Park Comprehensive Cancer Center/GSEA_signatures/h.all.v7.2.symbols.gmt"))
rownames(res_gsva) <- gsub("HALLMARK_", "",rownames(res_gsva))


df_gsva_paired <- data.frame(Pathway = rep(rownames(res_gsva),each=8),
                             Patient = rep(sapply(colnames(res_gsva)[1:8], function(xx) paste0("Pt.",strsplit(xx, "\\.")[[1]][2])),nrow(res_gsva)),
                             after = melt(t(res_gsva[,1:8]))$value,
                             before = melt(t(res_gsva[,9:16]))$value)

ggpaired(df_gsva_paired[which(df_gsva_paired$Pathway == "MYC_TARGETS_V2"),], 
         label = "Patient",
         cond1 = "before",
         cond2 = "after", 
         color = "condition", 
         fill = "condition",
         repel = T, label.rectangle = T) +
  scale_color_manual(values = c("before" = "blue", "after" = "red")) +
  scale_fill_manual(values = c("before" = scales::alpha("blue", 0.5), "after" = scales::alpha("red", 0.5)))

ggpaired(df_gsva_paired[which(df_gsva_paired$Pathway == "WNT_BETA_CATENIN_SIGNALING"),], 
         label = "Patient",
         cond1 = "before",
         cond2 = "after", 
         color = "condition", 
         fill = "condition",
         repel = T, label.rectangle = T) +
  scale_color_manual(values = c("before" = "blue", "after" = "red")) +
  scale_fill_manual(values = c("before" = scales::alpha("blue", 0.5), "after" = scales::alpha("red", 0.5)))

ggpaired(df_gsva_paired[which(df_gsva_paired$Pathway == "TGF_BETA_SIGNALING"),], 
         label = "Patient",
         cond1 = "before",
         cond2 = "after", 
         color = "condition", 
         fill = "condition",
         repel = T, label.rectangle = T) +
  scale_color_manual(values = c("before" = "blue", "after" = "red")) +
  scale_fill_manual(values = c("before" = scales::alpha("blue", 0.5), "after" = scales::alpha("red", 0.5)))




# --- >> Table S5 <<< -----
# use caret for RF classification ####
library(caret)
library(randomForest)
library(e1071)
library(pROC)
library(MLeval)
jci <- read.table("Verhaak_subClass_JCI_genes.txt", header = T, sep='\t')
tcga <- read.table("~/Downloads/ov_tcga/data_expression_all_sample_Zscores.txt", sep='\t', header = T)
dim(tcga)
tcga[1:10,1:6]
jci$Gene_Symbol[-which(jci$Gene_Symbol %in% toupper(tcga$Hugo_Symbol))]
tcga <- tcga[which(toupper(tcga$Hugo_Symbol) %in% jci$Gene_Symbol),]
tcga <- tcga[-1,]
rownames(tcga) <- tcga$Hugo_Symbol
tcga <- tcga[,-c(1:2)]
colnames(tcga) <- gsub("\\.","-",substr(colnames(tcga),1,12))
tcga <- t(tcga)
tcga[1:10,1:6]
classes <- read.table("tcga_subclass_4ML.txt")
colnames(classes) <- c("TCGA_ID", "GROUP", "CLASS")
head(classes)
mm  <- data.frame(cbind(classes, tcga[match(classes$TCGA_ID, rownames(tcga)),]))
mm <- mm[-which(mm$TCGA_ID == "TCGA-13-0760"),]
set.seed(1234)
ind_train <- sample(1:nrow(mm), nrow(mm)*.9,replace = F)
mm_train <- mm[ind_train,]
mm_test <- mm[-ind_train,]
set.seed(1111)
rf_model<-train(CLASS~., data=mm_train[,-c(1,2)],
                method="rf",
                trControl=trainControl(method="cv",number=10),
                prox=TRUE,
                allowParallel=TRUE, 
                ntree=500)
varImp(rf_model)
print(rf_model)
set.seed(2222)
new_pred <- predict(rf_model, newdata = mm_test)
cc <- confusionMatrix(new_pred, factor(mm_test$CLASS))

jci$Gene_Symbol[-which(jci$Gene_Symbol %in% rownames(rna))]
rownames(rna)[grep("NREP",rownames(rna),ignore.case = T)]
genes <- c(jci$Gene_Symbol, "AGFG2", "C1orf116", "C2CD2L", "TYMP","ACKR3","PMEPA1","SUSD6","HSPA13","NREP")
heat <- rna[which(rownames(rna) %in% genes),]
rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat, 1, scale))
colnames(heat) <- cls

to_test <- t(heat)
colnames(to_test)[which(colnames(to_test) == "AGFG2")] <- "HRBL"
colnames(to_test)[which(colnames(to_test) == "C2CD2L")] <- "TMEM24"
colnames(to_test)[which(colnames(to_test) == "TYMP")] <- "ECGF1"
colnames(to_test)[which(colnames(to_test) == "ACKR3")] <- "CXCR7"
colnames(to_test)[which(colnames(to_test) == "PMEPA1")] <- "TMEPAI"
colnames(to_test)[which(colnames(to_test) == "SUSD6")] <- "KIAA0247"
colnames(to_test)[which(colnames(to_test) == "HSPA13")] <- "STCH"
colnames(to_test)[which(colnames(to_test) == "NREP")] <- "C5orf13"

head(to_test)
set.seed(3333)
new_pred_validation <- predict(rf_model, newdata = to_test)
res_pred <- data.frame(cbind(as.character(new_pred_validation), rownames(to_test)))
colnames(res_pred) <- c("Class", "ID")
res_pred$Time <- sapply(res_pred$ID, function(xx) strsplit(as.character(xx), "_")[[1]][1])
res_pred$Sample <- sapply(res_pred$ID, function(xx) strsplit(as.character(xx), "_")[[1]][2])
table(res_pred$Class, res_pred$Time)


# --- >>> Figure S11A <<< ----
# use results from  IDO FK866 #
# setwd("~/OneDrive - Roswell Park Comprehensive Cancer Center/CITN_RNASeq")
rna_new <- read.table("Input_files/CountsSubread.csv", sep=",", header = T)
head(rna_new)
meta_rna_new <- read.table("Input_files/myPhenoFQ.csv", sep=",", header = T, skip = 1)
colnames(meta_rna_new)[c(1,4)] <- c("ID","Treatment")
meta_rna_new$ID <- gsub("-","\\.",meta_rna_new$ID)
colnames(rna_new)[3:19] <- meta_rna_new$Treatment

# read orthologous genes ####
orth_genes <- read.table("Input_files/Ortholog_Human_Mouse_genes.txt", sep = "\t", header = T)

# calculate TPM for new RNAseq ####
head(rna_new)
rna_cpm <- merge(orth_genes, rna_new, by.x="mouse", by.y="Gene")
rna_cpm <- rna_cpm[-which(rowSums(rna_cpm[,-c(1:4)]) < 10),]
rna_cpm <- aggregate(rna_cpm[,5:21], by=list(rna_cpm$human), FUN=mean)
rownames(rna_cpm) <- rna_cpm$Group.1
rna_cpm$Group.1 <- NULL
colnames(rna_cpm) <- gsub("\\.[0-9]","",colnames(rna_cpm))

# run GSVA ###
gmt_custom <- list(GO_NAD_METABOLIC_PROCESS = qusage::read.gmt("Input_files/c5.go.bp.v7.5.1.symbols.gmt")[["GO_NAD_METABOLIC_PROCESS"]])

gsva_res <- gsva(log2(data.matrix(rna_cpm)+1), 
                 kcdf="Gaussian",
                 gset.idx.list = gmt_custom,
                 method="gsva")
colnames(gsva_res) <- gsub("\\.[0-9]","",colnames(gsva_res))

TukeyHSD(aov(as.numeric(gsva_res) ~ factor(colnames(gsva_res), levels = c("vehicle", "Incyte", "FK866", "FK866_Incyte"))))
####
#                             diff        lwr       upr     p adj
# Incyte-vehicle        0.173829572 -0.4430779 0.7907371 0.8407378
# FK866-vehicle        -0.044398373 -0.6342748 0.5454780 0.9959997
# FK866_Incyte-vehicle -0.008699934 -0.5985763 0.5811765 0.9999692
# FK866-Incyte         -0.218227945 -0.7600640 0.3236081 0.6479778
# FK866_Incyte-Incyte  -0.182529506 -0.7243656 0.3593066 0.7582376
# FK866_Incyte-FK866    0.035698438 -0.4751495 0.5465464 0.9967887
####
boxplot(as.numeric(gsva_res) ~ factor(colnames(gsva_res), levels = c("vehicle", "Incyte", "FK866", "FK866_Incyte")), 
        main="GO_NAD_METABOLIC_PROCESS",
        las=2, 
        col=c("gray", "maroon", "blue", "purple"),
        ylab="", 
        xlab="")


# --- >>> Figure 7C <<< -----
# compare DEGs across groups ####
library(UpSetR)
library(DOSE)
library(enrichplot)
# library(pathview)

# note major updates in GSEA GO BP signatures
# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/MSigDB_v7.0_Release_Notes#C5_.28Gene_Ontology_collection.29_-_Major_overhaul

gmt_custom <- qusage::read.gmt("Input_files/c5.go.bp.v7.5.1.symbols.gmt")
gmt_custom <- data.frame(ont = rep(names(gmt_custom), as.numeric(lapply(gmt_custom, length))),
                         gene = unlist(gmt_custom))

ab_vehicle <- read.table("Input_files/IDOA2a_v_IDOVehicle.csv", sep=',', header = T)
aa_vehicle <- read.table("Input_files/IDOA2b_v_IDOVehicle.csv", sep=',', header = T)
incyte_vehicle <- read.table("Input_files/IDOIncyte_v_IDOVehicle.csv", sep=',', header = T)
triple_vehicle <- read.table("Input_files/IDOIncyteA2aA2b_v_IDOVehicle.csv", header = T, sep=",")

# run GSEA for triple
ranked_genes <- triple_vehicle$log2FoldChange
names(ranked_genes) <- triple_vehicle$Gene
ranked_genes <- sort(ranked_genes, decreasing = T)
head(ranked_genes);tail(ranked_genes)

gsea_custom <- clusterProfiler::GSEA(ranked_genes, 
                                     TERM2GENE = gmt_custom, 
                                     # nPerm = 1500,
                                     minGSSize = 5,
                                     pvalueCutoff = 0.05)
gsea_triple_gobp <- gsea_custom

# run GSEA for incyte
ranked_genes <- incyte_vehicle$log2FoldChange
names(ranked_genes) <- incyte_vehicle$Gene
ranked_genes <- sort(ranked_genes, decreasing = T)
head(ranked_genes);tail(ranked_genes)

gsea_custom <- clusterProfiler::GSEA(ranked_genes, 
                                     TERM2GENE = gmt_custom, 
                                     # nPerm = 1500,
                                     minGSSize = 5,
                                     pvalueCutoff = 0.05)
gsea_incyte_gobp <- gsea_custom

# run GSEA for A2a
ranked_genes <- aa_vehicle$log2FoldChange
names(ranked_genes) <- aa_vehicle$Gene
ranked_genes <- sort(ranked_genes, decreasing = T)
head(ranked_genes);tail(ranked_genes)

gsea_custom <- clusterProfiler::GSEA(ranked_genes, 
                                     TERM2GENE = gmt_custom, 
                                     # nPerm = 1500,
                                     minGSSize = 5,
                                     pvalueCutoff = 0.05)
gsea_A2a_gobp <- gsea_custom

# run GSEA for A2b
ranked_genes <- ab_vehicle$log2FoldChange
names(ranked_genes) <- ab_vehicle$Gene
ranked_genes <- sort(ranked_genes, decreasing = T)
head(ranked_genes);tail(ranked_genes)

gsea_custom <- clusterProfiler::GSEA(ranked_genes, 
                                     TERM2GENE = gmt_custom, 
                                     # nPerm = 1500,
                                     minGSSize = 5,
                                     pvalueCutoff = 0.05)
gsea_A2b_gobp <- gsea_custom

# gather info on immune pathways from GSEA results ####
grep_immune <- function(xx){
  ind_keep <- c(grep("immun", xx@result$Description, ignore.case = T),
                grep("cytok", xx@result$Description, ignore.case = T),
                grep("chemok", xx@result$Description, ignore.case = T),
                grep("tcr", xx@result$Description, ignore.case = T),
                grep("_il", xx@result$Description, ignore.case = T),
                grep("_t_cell", xx@result$Description, ignore.case = T),
                grep("_B_cell", xx@result$Description, ignore.case = T),
                grep("lympho", xx@result$Description, ignore.case = T),
                grep("leuko", xx@result$Description, ignore.case = T),
                grep("macropha", xx@result$Description, ignore.case = T),
                grep("nad", xx@result$Description, ignore.case = T),
                grep("chemotax", xx@result$Description, ignore.case = T),
                grep("interferon", xx@result$Description, ignore.case = T),
                grep("interleuk", xx@result$Description, ignore.case = T),
                grep("tgf", xx@result$Description, ignore.case = T),
                grep("antigen", xx@result$Description, ignore.case = T),
                grep("mhc", xx@result$Description, ignore.case = T),
                grep("cd4_", xx@result$Description, ignore.case = T),
                grep("cd8_", xx@result$Description, ignore.case = T),
                grep("tumor_necrosis_factor", xx@result$Description, ignore.case = T),
                grep("transforming_growth_factor", xx@result$Description, ignore.case = T))
  ind_keep <- unique(ind_keep)
  data_keep <- data.frame(xx@result[ind_keep,c(2,3,5,6,7,11)])
  data_keep <- data_keep[which(data_keep$p.adjust < 0.05),]
  return(data_keep)
}

imm_A2a <- grep_immune(gsea_A2a_gobp)
imm_A2b <- grep_immune(gsea_A2b_gobp)
imm_incyte <- grep_immune(gsea_incyte_gobp)
imm_triple <- grep_immune(gsea_triple_gobp)

gplots::venn(list(A2a=imm_A2a$Description,
                  A2b=imm_A2b$Description,
                  Triple=imm_triple$Description))

imm_common <- data.frame(rbind(imm_triple[intersect(imm_A2b$Description, imm_triple$Description),c(3,5)],
                               imm_A2a[intersect(imm_A2a$Description, imm_triple$Description),c(3,5)]))

# do heatmap with leading edge genes from immune pathways ####
genes_common <- unique(gmt_custom[which(gmt_custom$ont %in% rownames(imm_common)),2])
deg_files <- list.files(path = "Input_files", pattern = "*.csv", full.names = T)
for(ii in 1:4){
  if(ii == 1){
    temp_deg_files <- deg_files[grep("_v_IDOVehicle.csv", deg_files)]
  }
  temp_deg <- read.table(temp_deg_files[ii], sep=",", header = T)
  temp_fc <- ifelse(temp_deg$padj < 0.05 & abs(temp_deg$log2FoldChange) > 1, temp_deg$log2FoldChange, 0)
  names(temp_fc) <- temp_deg$Gene
  
  
  if(ii == 1){
    deg_heat <- data.frame(GeneName=names(temp_fc),
                           temp_fc)
  } else {
    deg_heat <- data.frame(cbind(deg_heat,
                                 temp_fc[match(deg_heat$GeneName, names(temp_fc))]))
    if(sum(is.na(deg_heat[,ii+1])) > 0) { 
      deg_heat <- deg_heat[-which(is.na(deg_heat[,ii+1])),]
    } else {
      deg_heat <- deg_heat
    }
  }
  colnames(deg_heat)[ii+1] <- strsplit(temp_deg_files[ii], "_")[[1]][1]
}

heat <- deg_heat[which(deg_heat$GeneName %in% genes_common),]
rownames(heat) <- heat$GeneName
heat$GeneName <- NULL
cls <- colnames(heat)
rws <- rownames(heat)
rownames(heat) <- rws
colnames(heat) <- c("A2a", "A2b", "EPA", "EPA_A2a_A2b")
heat <- heat[which(apply(heat, 1, function(x) sum(x == 0)) < 4),]
dim(heat)

breaks <- seq(-2,2,length.out=50)
clrs <- colorRampPalette(colors = c("orange", "white", "purple"))(length(breaks))
pheatmap(heat,
         # scale = "cols",
         show_rownames = T,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows = T,
         border_color = "lightgray",
         na_col = "#F5F5F5",
         fontsize_col = 15,
         fontsize_row = 15,
         color = clrs, 
         name = "Log2FC",
         breaks = breaks)




############################################################################################
## The code below is for the analyses of the metabolism data for the manuscript:
## "IDO1 inhibition and subsequent metabolic adaptations constrains anti-tumor immune 
## responses in the tumor microenvironment of patients with ovarian cancer"
############################################################################################

# rm(list = ls())
library("doBy")
library("ggplot2")
library("ggsci")
library("RColorBrewer")
library("DESeq2")
library("igraph")
library("graphite")
library("psych")
library("pheatmap")
library("ggpubr")
library("scales")
library("reshape2")
library("Partiallyoverlapping")

##############################################
## Comparison of plasma Kyn and Trp between time points
##############################################

plasma <- read.csv(file = "./data_metabolism/Metabolite_plasma.csv", stringsAsFactors = FALSE)
plasma$Day <- as.character(plasma$Day)
plasma$ID <- gsub("\\..*", "", plasma$ID)

df <- plasma
df_0 <- subset(df, Day == "0")
df_14 <- subset(df, Day == "14")
df_15 <- subset(df, Day == "15")
rownames(df_0) <- df_0$ID
rownames(df_14) <- df_14$ID
rownames(df_15) <- df_15$ID

patients <- unique(plasma$ID)
plasma_all <- plasma

vvec <- c("Kyn", "Trp", "KT_ratio")
for (i in 1:length(vvec)) {
  v <- vvec[i]
  d0 <- df_0[patients,v]
  d14 <- df_14[patients,v]
  d15 <- df_15[patients,v]
  t12 <- t.test(d14, d0, paired = TRUE)
  t23 <- t.test(d15, d14, paired = TRUE)
  t13 <- t.test(d15, d0, paired = TRUE)
  output_0 <- rbind(unlist(t12), unlist(t23), unlist(t13))
  colnames(output_0)
  output <- output_0[, c(6, 8, 4, 5, 1, 2, 3, 10)]
  colnames(output) <- c("Estimate", "SE", "CI_lower", "CI_upper", "t", "df", "p value", "comparison")
  fname <- paste0("Plasma_", v, "_2.csv")
  # write.csv(output, file = fname, quote = FALSE, row.names = FALSE)
}

plasma_0 <- subset(plasma_all, Day == 0)
plasma_14 <- subset(plasma_all, Day == 14)
plasma_15 <- subset(plasma_all, Day == 15)
rownames(plasma_0) <- plasma_0$ID
rownames(plasma_14) <- plasma_14$ID
rownames(plasma_15) <- plasma_15$ID

kyn_0 <- plasma_0[patients,]$Kyn
kyn_14 <- plasma_14[patients,]$Kyn
kyn_15 <- plasma_15[patients,]$Kyn
trp_0 <- plasma_0[patients,]$Trp
trp_14 <- plasma_14[patients,]$Trp
trp_15 <- plasma_15[patients,]$Trp
kt_0 <- plasma_0[patients,]$KT_ratio
kt_14 <- plasma_14[patients,]$KT_ratio
kt_15 <- plasma_15[patients,]$KT_ratio

my_comparisons <- list(c("0", "14"), c("14", "15"), c("0", "15"))
ggpaired(plasma, x="Day", y="Kyn", color = "Day", id = "ID", xlab = "Day", ylab = "Plasma Kyn(uM)",
         line.color = "gray", line.size = 0.4, palette = "aaas", point.size = 1.5) + 
  theme_bw() +
  stat_compare_means(comparisons=my_comparisons, paired = TRUE, method = "t.test", na.rm = TRUE, label = "p.format") +
  guides(fill=FALSE, color=FALSE)

ggpaired(plasma, x="Day", y="Trp", color = "Day", id = "ID", xlab = "Day", ylab = "Plasma Trp(uM)",
         line.color = "gray", line.size = 0.4, palette = "aaas", point.size = 1.5) + 
  theme_bw() +
  stat_compare_means(comparisons=my_comparisons, paired = TRUE, method = "t.test") +
  guides(fill=FALSE, color=FALSE)

ggpaired(plasma, x="Day", y="KT_ratio", color = "Day", id = "ID", xlab = "Day", ylab = "Plasma Kyn/Trp",
         line.color = "gray", line.size = 0.4, palette = "aaas", point.size = 1.5) + 
  theme_bw() +
  stat_compare_means(comparisons=my_comparisons, paired = TRUE, method = "t.test") +
  guides(fill=FALSE, color=FALSE)

##############################################
## Comparison of ascites Kyn and Trp between time points
##############################################

df <- read.csv(file = "./data_metabolism/Metabolite_ascites.csv", stringsAsFactors = FALSE)
df$Day <- as.character(df$Day)
df$ID <- gsub("\\..*", "", df$ID)

df_0 <- subset(df, Day == "0")
df_14 <- subset(df, Day == "14")
df_15 <- subset(df, Day == "15")
rownames(df_0) <- df_0$ID
rownames(df_14) <- df_14$ID
rownames(df_15) <- df_15$ID
patients <- unique(df$ID)

my_comparisons <- list(c("0", "15"))
df.sub <- subset(df, Day!="14")

vvec <- c("Kyn", "Trp", "KT_ratio")

## Kyn
v <- vvec[1]
d0 <- df_0[patients,v]
d15 <- df_15[patients,v]
t13 <- Partover.test(d15, d0, var.equal = TRUE, stacked = TRUE, conf.level = 0.95)

stat.test <- compare_means(Kyn ~ Day, data = df.sub, method = "t.test")
stat.test$p[1] <- round(t13$p.value, 4) # label the p-value from the partially paired test
stat.test$y.position <- 6.7
ggpaired(df.sub, x="Day", y="Kyn", color = "Day", id = "ID", xlab = "Day", ylab = "Ascites Kyn(uM)",
         line.color = "gray", line.size = 0.4, palette = pal_aaas("default")(3)[c(1,3)], point.size = 1.5) + 
  stat_pvalue_manual(stat.test, label = "p") +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)

## Trp
v <- vvec[2]
d0 <- df_0[patients,v]
d15 <- df_15[patients,v]
t13 <- Partover.test(d15, d0, var.equal = TRUE, stacked = TRUE, conf.level = 0.95)

stat.test <- compare_means(Trp ~ Day, data = df.sub, method = "t.test")
stat.test$p[1] <- round(t13$p.value, 4) # label the p-value from the partially paired test
max(df.sub$Trp, na.rm = TRUE)
stat.test$y.position <- 55

ggpaired(df.sub, x="Day", y="Trp", color = "Day", id = "ID", xlab = "Day", ylab = "Ascites Trp(uM)",
         line.color = "gray", line.size = 0.4, palette = pal_aaas("default")(3)[c(1,3)], point.size = 1.5) + 
  stat_pvalue_manual(stat.test, label = "p") +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)

## KT ratio
v <- vvec[3]
d0 <- df_0[patients,v]
d15 <- df_15[patients,v]
t13 <- Partover.test(d15, d0, var.equal = TRUE, stacked = TRUE, conf.level = 0.95)

stat.test <- compare_means(KT_ratio ~ Day, data = df.sub, method = "t.test")
stat.test$p[1] <- round(t13$p.value, 4) # label the p-value from the partially paired test
max(df.sub$KT_ratio, na.rm = TRUE)
stat.test$y.position <- 0.28

ggpaired(df.sub, x="Day", y="KT_ratio", color = "Day", id = "ID", xlab = "Day", ylab = "Ascites Kyn/Trp",
         line.color = "gray", line.size = 0.4, palette = pal_aaas("default")(3)[c(1,3)], point.size = 1.5) + 
  stat_pvalue_manual(stat.test, label = "p") +
  theme_bw() +
  guides(fill=FALSE, color=FALSE)

##############################################
# correlation between plasma and ascites kyn and trp levels
##############################################

### Plasma data
plasma <- read.csv(file = "./data_metabolism/Metabolite_plasma.csv", stringsAsFactors = FALSE)
plasma$Day <- as.character(plasma$Day)
plasma$ID <- gsub("\\..*", "", plasma$ID)
df <- plasma

df_0 <- subset(df, Day == "0")
df_14 <- subset(df, Day == "14")
df_15 <- subset(df, Day == "15")
rownames(df_0) <- df_0$ID
rownames(df_14) <- df_14$ID
rownames(df_15) <- df_15$ID

### Ascites sata
dfa <- read.csv(file = "./data_metabolism/Metabolite_ascites.csv", stringsAsFactors = FALSE)
dfa$Day <- as.character(dfa$Day)
dfa$ID <- gsub("\\..*", "", dfa$ID)

dfa_0 <- subset(dfa, Day == "0")
dfa_15 <- subset(dfa, Day == "15")
rownames(dfa_0) <- dfa_0$ID
rownames(dfa_15) <- dfa_15$ID

patients <- unique(plasma$ID)
plasma_all <- plasma

## Kyn
v <- "Kyn"
p.d0 <- df_0[patients,v]
p.d15 <- df_15[patients,v]
a.d0 <- dfa_0[patients,v]
a.d15 <- dfa_15[patients,v]

df.plot <- data.frame(plasma = p.d0, ascites = a.d0)
lm.fit <- lm(plasma~ascites, df.plot)
pval <- coef(summary(lm.fit))[2,4]
ggplot(df.plot, aes(x = ascites, y = plasma)) + 
  geom_point() +
  stat_smooth(method = "lm", col = pal_aaas("default")(1)[1], se = FALSE, size = 0.5) +
  xlab(paste0("Ascites ", v)) + ylab(paste0("Plasma ", v)) +
  theme_bw() +
  ggtitle(paste0("Day 0, p=", round(pval, 3)))

df.plot <- data.frame(plasma = p.d15, ascites = a.d15)
lm.fit <- lm(plasma~ascites, df.plot)
pval <- coef(summary(lm.fit))[2,4]
ggplot(df.plot, aes(x = ascites, y = plasma)) + 
  geom_point() +
  stat_smooth(method = "lm", col = pal_aaas("default")(1)[1], se = FALSE, size = 0.5) +
  xlab(paste0("Ascites ", v)) + ylab(paste0("Plasma ", v)) +
  theme_bw() +
  ggtitle(paste0("Day 15, p=", round(pval, 3)))
  
## Kyn/Trp ratio
v <- "KT_ratio"
p.d0 <- df_0[patients,v]
p.d15 <- df_15[patients,v]
a.d0 <- dfa_0[patients,v]
a.d15 <- dfa_15[patients,v]

df.plot <- data.frame(plasma = p.d0, ascites = a.d0)
lm.fit <- lm(plasma~ascites, df.plot)
pval <- coef(summary(lm.fit))[2,4]

ggplot(df.plot, aes(x = ascites, y = plasma)) + 
  geom_point() +
  stat_smooth(method = "lm", col = pal_aaas("default")(1)[1], se = FALSE, size = 0.5) +
  xlab(paste0("Ascites ", v)) + ylab(paste0("Plasma ", v)) +
  theme_bw() +
  ggtitle(paste0("Day 0, p=", round(pval, 3)))

df.plot <- data.frame(plasma = p.d15, ascites = a.d15)
lm.fit <- lm(plasma~ascites, df.plot)
pval <- coef(summary(lm.fit))[2,4]

ggplot(df.plot, aes(x = ascites, y = plasma)) + 
  geom_point() +
  stat_smooth(method = "lm", col = pal_aaas("default")(1)[1], se = FALSE, size = 0.5) +
  xlab(paste0("Ascites ", v)) + ylab(paste0("Plasma ", v)) +
  theme_bw() +
  ggtitle(paste0("Day 15, p=", round(pval, 3)))

##############################################
# Compare tumor CD3 and CD8:CD3 at day 0 and 15
##############################################

df.comp <- read.csv(file = "./data_metabolism/TIL.csv", header = TRUE)
df.comp$id <- paste0(df.comp$SubjectID, "-", df.comp$Timepoint)
df.comp.d0 <- subset(df.comp, Timepoint == 0)
df.comp.d14 <- subset(df.comp, Timepoint == 14)
df.comp.d15 <- subset(df.comp, Timepoint == 15)
rownames(df.comp.d0) <- df.comp.d0$SubjectID
rownames(df.comp.d14) <- df.comp.d14$SubjectID
rownames(df.comp.d15) <- df.comp.d15$SubjectID

## CD8
x1 <- df.comp.d0$IT.IHC.CD8; x15 <- df.comp.d15$IT.IHC.CD8
df.compare <- data.frame(ID = c(as.character(df.comp.d0$SubjectID), as.character(df.comp.d15$SubjectID)),
                         Day = rep(c(0, 15), each = 17),
                         level = c(x1, x15))

ggpaired(df.compare, x="Day", y="level", color = "Day", id = "ID", xlab = "Day", ylab = "CD3",
         line.color = "gray", line.size = 0.4, palette = "aaas", point.size = 1.5) + 
  theme_bw() +
  guides(fill=FALSE, color=FALSE)

## CD8/CD3
x1 <- df.comp.d0$IT.IHC.CD3.8; x15 <- df.comp.d15$IT.IHC.CD3.8
df.compare <- data.frame(ID = c(as.character(df.comp.d0$SubjectID), as.character(df.comp.d15$SubjectID)),
                         Day = rep(c(0, 15), each = 17),
                         level = c(x1, x15))

ggpaired(df.compare, x="Day", y="level", color = "Day", id = "ID", xlab = "Day", ylab = "CD3",
         line.color = "gray", line.size = 0.4, palette = "aaas", point.size = 1.5) + 
  theme_bw() +
  guides(fill=FALSE, color=FALSE)

##############################################
# process the first metabolomics data set
##############################################

meta_tumor <- read.csv("./data_metabolism/metabolomics_tumor_1.csv")
meta_info_tumor <- read.csv("./data_metabolism/sample_info_meta_tumor.csv")
id_tumor <- apply(meta_info_tumor, 1, paste0, collapse = "-")
id_tumor <- gsub(" ", "", id_tumor)

rownames(meta_tumor) <- meta_tumor[,1]
meta_tumor <- meta_tumor[,-1]
colnames(meta_tumor) <- id_tumor

df_meta_t_duplicate <- t(meta_tumor)
meta_info_tumor_unique <- unique(meta_info_tumor)
id_tumor_unique <- apply(meta_info_tumor_unique, 1, paste0, collapse = "-")
id_tumor_unique <- as.vector(gsub(" ", "", id_tumor_unique))

df_meta_t <- matrix(NA, nrow = length(id_tumor_unique), ncol = ncol(df_meta_t_duplicate))
for (i in 1:length(id_tumor_unique)) {
  # i <- 1
  this.ind <- which(id_tumor == id_tumor_unique[i])
  if (length(this.ind) == 1) {
    df_meta_t[i, ] <- df_meta_t_duplicate[this.ind, ]
    next
  }
  df.sub <- df_meta_t_duplicate[this.ind,]
  vec.temp <- colMeans(df.sub)
  df_meta_t[i, ] <- vec.temp
}

rownames(df_meta_t) <- id_tumor_unique
colnames(df_meta_t) <- rownames(meta_tumor)

df_meta_t_0 <- df_meta_t[2:15, ]
df_meta_t_15 <- df_meta_t[16:29, ]
idmet <- substr(rownames(df_meta_t_0), 1, 6)
df_meta_t_norm_0 <- df_meta_t_0[, 2:16]
df_meta_t_norm_15 <- df_meta_t_15[, 2:16]
for(i in 1:15) {
  df_meta_t_norm_0[,i] <- df_meta_t_norm_0[,i] - df_meta_t_0[,1]
  df_meta_t_norm_15[,i] <- df_meta_t_norm_15[,i] - df_meta_t_15[,1]
}

df_meta_t_norm_diff <- df_meta_t_norm_15 - df_meta_t_norm_0
rownames(df_meta_t_norm_diff) <- idmet

# process the second metabolite data set
meta_tumor_2_d0 <- read.csv("./data_metabolism/Metabolite_tumor_2_D0.csv", stringsAsFactors = FALSE)
meta_tumor_2_d15 <- read.csv("./data_metabolism/Metabolite_tumor_2_D15.csv", stringsAsFactors = FALSE)
meta_tumor_2_d0_mean <- summaryBy(. ~ Sample, data = meta_tumor_2_d0, FUN = "mean", keep.names = TRUE)
meta_tumor_2_d15_mean <- summaryBy(. ~ Sample, data = meta_tumor_2_d15, FUN = "mean", keep.names = TRUE)
meta_2_d0 <- meta_tumor_2_d0_mean[,-1]
meta_2_d15 <- meta_tumor_2_d15_mean[,-1]
rownames(meta_2_d0) <- meta_tumor_2_d0_mean$Sample
rownames(meta_2_d15) <- meta_tumor_2_d15_mean$Sample
df_meta_2_d0 <- meta_2_d0[rownames(df_meta_t_norm_diff), ]
df_meta_2_d15 <- meta_2_d15[rownames(df_meta_t_norm_diff), ]
trp_d0 <- df_meta_t_0[, "Tryptophan"]
trp_d15 <- df_meta_t_15[, "Tryptophan"]
df_meta_2_norm_d0 <- df_meta_2_d0 - trp_d0
df_meta_2_norm_d15 <- df_meta_2_d15 - trp_d15
df_meta_2_norm_diff <- df_meta_2_norm_d15 - df_meta_2_norm_d0
# metabolites in trp and nicotinamide pathways
df_meta_2_norm_diff_sub <- df_meta_2_norm_diff[, c(2, 10, 20, 29, 30, 31, 32)]

##############################################
# build correlation network
##############################################

rna_info <- read.csv(file = "./data_metabolism/sample_info_rna.csv")
id_rna <- apply(rna_info, 1, paste0, collapse = "-")
id_rna <- gsub(" ", "", id_rna)
rna.readin <- read.delim(file = "./data_metabolism/CITN_prepost_rna_counts.txt")
rna.0 <- as.matrix(rna.readin[, -1])
rna.1 <- rlogTransformation(round(rna.0))
rna.readin[,1] <- as.character(rna.readin[,1])
rm.ind <- which(grepl("^[0-9]{1,2}-", rna.readin[,1]),1)
rna.2 <- rna.1[-rm.ind,]
rownames(rna.2) <- gene.names <- as.character(rna.readin[-rm.ind, 1])
df_rna <- t(rna.2)
rownames(df_rna) <- id_rna
inter.2 <- intersect(id_rna, id_tumor_unique)
df_rna_in <- df_rna[inter.2, ]
df_rna_all_0 <- df_rna_in[9:16, ]
df_rna_all_15 <- df_rna_in[1:8, ]
rownames(df_rna_all_0) <- substr(rownames(df_rna_all_0), 1, 6)
rownames(df_rna_all_15) <- substr(rownames(df_rna_all_15), 1, 6)
df_rna_all_0 <- df_rna_all_0[rownames(df_rna_all_15), ]
df_rna_diff <- df_rna_all_15 - df_rna_all_0

## The code below can be used to load KEGG tryptophan and nicotinamide metabolism pathways through "graphite" package
# humanKegg <- pathways("hsapiens", "kegg")
# kegg.entry <- humanKegg@entries
# which(grepl("tryp", names(kegg.entry), ignore.case = TRUE)) 
# pSymbol_1 <- convertIdentifiers(kegg.entry[[30]], "SYMBOL")
# pSymbol_2 <- convertIdentifiers(kegg.entry[[66]], "SYMBOL")

# load KEGG tryptophan and nicotinamide metabolism pathways
load("./data_metabolism/pathways_trp_nic.rda")
p.graph_1 <- pathwayGraph(pSymbol_1, which = "mixed")
p.graph_2 <- pathwayGraph(pSymbol_2, which = "mixed")
p.igraph_1 <- graph_from_graphnel(p.graph_1)
p.igraph_2 <- graph_from_graphnel(p.graph_2)
comp.nodes_1 <- which(grepl("KEGGCOMP:", V(p.igraph_1)$name))
comp.nodes_2 <- which(grepl("KEGGCOMP:", V(p.igraph_2)$name))
V(p.igraph_1)$name <- gsub("SYMBOL:", "", V(p.igraph_1)$name)
V(p.igraph_1)$name <- gsub("KEGGCOMP:", "", V(p.igraph_1)$name)
V(p.igraph_2)$name <- gsub("SYMBOL:", "", V(p.igraph_2)$name)
V(p.igraph_2)$name <- gsub("KEGGCOMP:", "", V(p.igraph_2)$name)

# map compound-id to compound names
keggcomp.df <- read.csv(file = "./data_metabolism/kegg_compound_id.csv", header = FALSE)
keggcomp.vec <- as.character(keggcomp.df[,2])
names(keggcomp.vec) <- keggcomp.df[,1]

V(p.igraph_1)$name[comp.nodes_1] <- keggcomp.vec[V(p.igraph_1)$name[comp.nodes_1]]
V(p.igraph_2)$name[comp.nodes_2] <- keggcomp.vec[V(p.igraph_2)$name[comp.nodes_2]]
p.igraph <- igraph::union(p.igraph_1, p.igraph_2)

node_names <- V(p.igraph)$name
type <- rep("gene", length(node_names))
type[node_names %in% keggcomp.vec] <- "metabolite"
df.node <- data.frame(node_names, type)
node_names_in <- intersect(colnames(df_rna_diff), V(p.igraph_1)$name)
df_rna_trp <- df_rna_diff[, node_names_in]

gene.s <- read.csv(file = "./data_metabolism/genes_correlation.csv", header = FALSE)[,1]
sub.net.nodes <- intersect(gene.s, colnames(df_rna_trp))

df_trp <- data.frame(df_rna_trp[, sub.net.nodes], 
                     df_meta_t_norm_diff[rownames(df_rna_trp),], 
                     df_meta_2_norm_diff_sub[rownames(df_rna_trp),])
cor.mat <- cor(df_trp, method = "pearson") 

df_trp_met <- data.frame(df_meta_t_norm_diff, df_meta_2_norm_diff_sub)
# the correlation network between metabolites used info from 14 patients having metabolomics data
cor.met.2 <- cor(df_trp_met, method = "pearson") 
delta <- 0.6
Amet.2 <- cor.met.2
Amet.2[abs(cor.met.2)>delta] <- 1
Amet.2[abs(cor.met.2)<=delta] <- 0

A <- matrix(0, ncol = ncol(cor.mat), nrow = nrow(cor.mat))
A[abs(cor.mat)>delta] <- 1 
colnames(A) <- rownames(A) <- colnames(df_trp)
# the adjacency matrix among metabolites was replaced by the results from a larger group obtained above
A[13:34, 13:34] <- Amet.2 
diag(A) <- 0
diag(Amet.2) <- 0
this.memb.2 <- colnames(A) <- rownames(A) <- rownames(cor.mat)
netgraph <- graph_from_adjacency_matrix(A, mode = "undirected")

comm <- fastgreedy.community(netgraph)
V(netgraph)$size <- 3
V(netgraph)$color <- comm$membership
plot(netgraph) # visualize the correlation network

## obtain the cluster scores
clust <- c("Kynurenine", "Serotonin", "MAOA", "Nicotinamide", "None")
df_memb <- data.frame(
  name = V(netgraph)$name,
  clustr_membership = clust[comm$membership]
)
metabolites <- V(netgraph)$name[13:34]
df_meta_t_diff <- cbind(df_meta_t_norm_diff, df_meta_2_norm_diff)
colnames(df_meta_t_diff) <- gsub("-", "\\.", colnames(df_meta_t_diff))

# get all members in the same cluster as kynurenine
m.kyn <- V(netgraph)$name[comm$membership == comm$membership[13]]
m.kyn <- intersect(m.kyn, metabolites)
m.kyn <- gsub("X", "", m.kyn)
df_meta_kyn <- df_meta_t_diff[, m.kyn]
pc.fit.kyn <- prcomp(df_meta_kyn, scale. = TRUE)
pc.kyn <- pc.fit.kyn$x[,1]

# get all members in the same cluster as serotonin
m.ser <- V(netgraph)$name[comm$membership == comm$membership[16]]
m.ser <- intersect(m.ser, metabolites)
m.ser <- gsub("X", "", m.ser)
df_meta_ser <- df_meta_t_diff[, m.ser]
pc.fit.ser <- prcomp(df_meta_ser, scale. = TRUE)
pc.ser <- pc.fit.ser$x[,1]

# get all members in the same cluster as nicotinamide
m.nic <- V(netgraph)$name[comm$membership == comm$membership[15]]
m.nic <- intersect(m.nic, metabolites)
df_meta_nic <- df_meta_t_diff[, m.nic]
pc.fit.nic <- prcomp(df_meta_nic, scale. = TRUE)
pc.nic <- pc.fit.nic$x[,1]

SAM_SAH <- df_meta_2_norm_diff$SAM - df_meta_2_norm_diff$SAH

cor.test(SAM_SAH, pc.kyn)
df.plot <- data.frame(x = pc.kyn, y = SAM_SAH)
ggplot(df.plot, aes(x = x, y = y)) + 
  geom_point() +
  stat_smooth(method = "lm", col = pal_aaas("default")(1)[1], se = FALSE, size = 0.5) +
  ylab("SAM/SAH ratio D15 vs. D00") +
  xlab("Kyn reduction score") +
  theme_bw() +
  # the correlation coefficient and 95% CI is from the cor.test function above
  ggtitle(paste0("r = -0.79 (-0.93, -0.46)")) +
  theme(plot.title = element_text(hjust = 1, size = 10, colour = "darkgreen"),
        axis.title.x = element_text(colour = "darkgreen"))

cor.test(SAM_SAH, -pc.ser)
df.plot <- data.frame(x = -pc.ser, y = SAM_SAH)
ggplot(df.plot, aes(x = x, y = y)) + 
  geom_point() +
  stat_smooth(method = "lm", col = pal_aaas("default")(1)[1], se = FALSE, size = 0.5) +
  ylab("SAM/SAH ratio D15 vs. D00") +
  xlab("Ser elevation score") +
  theme_bw() +
  ggtitle(paste0("r = 0.28 (-0.29, 0.71)")) +
  theme(plot.title = element_text(hjust = 1, size = 10, colour = "orangered"),
        axis.title.x = element_text(colour = "orangered"))

cor.test(SAM_SAH, pc.nic)
df.plot <- data.frame(x = pc.nic, y = SAM_SAH)
ggplot(df.plot, aes(x = x, y = y)) + 
  geom_point() +
  stat_smooth(method = "lm", col = pal_aaas("default")(1)[1], se = FALSE, size = 0.5) +
  ylab("SAM/SAH ratio D15 vs. D00") +
  xlab("Nic elevation score") +
  theme_bw() +
  ggtitle(paste0("r = -0.66 (-0.88, -0.20)")) +
  theme(plot.title = element_text(hjust = 1, size = 10, colour = "darkcyan"),
        axis.title.x = element_text(colour = "darkcyan"))

##############################################
# PCA of the metabolomics data
##############################################

df_meta_trp <- df_meta_t[, 2:16]
for (i in 1:ncol(df_meta_trp)) {
  df_meta_trp[,i] <- df_meta_trp[,i] - df_meta_t[,1]
}

df_meta_2_norm <- rbind(df_meta_2_norm_d0, df_meta_2_norm_d15)
df_meta_trp_all <- cbind(df_meta_trp[-1,], df_meta_2_norm)

meta_map <- read.csv("./data_metabolism/metabolite_info.csv", stringsAsFactors = FALSE)
meta_map$SMILES <- NULL
pc.fit.met <- prcomp(df_meta_trp_all, scale. = TRUE)
df.projection <- data.frame(met = names(pc.fit.met$rotation[,1]), 
                            score = -pc.fit.met$rotation[,1], 
                            Pathway = meta_map$Pathway)
df.projection$colour <- "orangered"
df.projection <- orderBy(~score, df.projection)
df.projection$met <- factor(df.projection$met, levels = df.projection$met)
df.projection$hjust <- ifelse(df.projection$score > 0, 1.2, -0.2)
df.projection$colour <- ifelse(df.projection$score < 0, "firebrick1", "steelblue")

pal <- pal_nejm()(7)
pal[1] <- "chocolate4"
  
ggplot(df.projection, aes(met, score, hjust = hjust, label = met)) +
  # geom_text(aes(y = 0), size = 3, colour = "black") +
  geom_bar(stat = "identity", aes(fill = Pathway)) +
  theme_bw() +
  coord_flip() + labs(x = "", y = "- PC1 projection") +
  scale_fill_manual(values = pal) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.key.size = unit(0.1, "in"),
        legend.position="bottom")

plot(pc.fit.met$x[,1], pc.fit.met$x[,2], pch=19, col = rep(c("dodgerblue", "orangered"), each=14),
     xlab = "PC1", ylab="PC2")
legend("bottomright", legend = c("Day0", "Day15"), col = c("dodgerblue", "orangered"), pch=19)



##############################################
# Compare all tumor metaolite levels at day 0 and 15
##############################################

df_meta_t_norm_0_temp <-  data.frame(Day = 0, ID = substr(rownames(df_meta_t_norm_0), 1, 6), df_meta_t_norm_0)
df_meta_t_norm_15_temp <-  data.frame(Day = 15, ID = substr(rownames(df_meta_t_norm_15), 1, 6), df_meta_t_norm_15)
df_meta_t_norm_temp <- rbind(df_meta_t_norm_0_temp, df_meta_t_norm_15_temp)
df_meta_melt_t <- melt(df_meta_t_norm_temp, id = c("Day", "ID"))

df_meta_2_norm_d0_temp <- data.frame(Day = 0, ID = rownames(df_meta_2_norm_d0), df_meta_2_norm_d0)
df_meta_2_norm_d15_temp <- data.frame(Day = 15, ID = rownames(df_meta_2_norm_d15), df_meta_2_norm_d15)
df_meta_2_norm_temp <- rbind(df_meta_2_norm_d0_temp, df_meta_2_norm_d15_temp)

df_meta_melt <- melt(df_meta_2_norm_temp, id = c("Day", "ID"))

df_meta_melt_all <- rbind(df_meta_melt_t, df_meta_melt)
df_meta_melt_all$Day <- factor(df_meta_melt_all$Day)
df_meta_melt_all$variable <- as.character(df_meta_melt_all$variable)

df_meta_melt_all$variable[df_meta_melt_all$variable == "Nicotinic.acid.ribonucleotide"] <- "NARN"
df_meta_melt_all$variable[df_meta_melt_all$variable == "X5.Hydroxytryptophan"] <- "5-HydroxyTrp"
df_meta_melt_all$variable[df_meta_melt_all$variable == "X5.Methoxyindole.3.acetic.acid"] <- "5-MI-3-AA"
df_meta_melt_all$variable[df_meta_melt_all$variable == "X5.Hydroxyindoleacetic.acid"] <- "5-HIAA"
df_meta_melt_all$variable[df_meta_melt_all$variable == "X3hydroxyanthranilic.acid"] <- "3-HAnthranilate"
df_meta_melt_all$variable[df_meta_melt_all$variable == "Methyl.2.pyridone.5.carboxamide"] <- "2-PY"
df_meta_melt_all$variable[df_meta_melt_all$variable == "N.methylnicotinamide"] <- "MNA"
df_meta_melt_all$variable[df_meta_melt_all$variable == "Hydroxy.kynurenine"] <- "HKYN"
df_meta_melt_all$variable[df_meta_melt_all$variable == "Uridine.monophosphate"] <- "UMP"
df_meta_melt_all$variable[df_meta_melt_all$variable == "Inosine.monophosphate"] <- "IMP"
df_meta_melt_all$variable[df_meta_melt_all$variable == "Adenosine.monophosphate"] <- "AMP"
df_meta_melt_all$variable[df_meta_melt_all$variable == "Guanosine.monophosphate"] <- "GMP"
df_meta_melt_all$variable[df_meta_melt_all$variable == "Uridine.monophosphate"] <- "UMP"
df_meta_melt_all$variable[df_meta_melt_all$variable == "Hydroxy.glutarate"] <- "HGlutarate"
df_meta_melt_all$variable <- gsub("\\.", " ", df_meta_melt_all$variable)

ggplot(df_meta_melt_all) +
  geom_boxplot(aes(x = Day, y = value, group = Day))+
  geom_point(aes(x = Day, y = value), color = "gray70", size = 0.5) +
  geom_line(aes(x  = Day, y = value, group = ID), color = "gray70", lwd = 0.5) +
  # scale_x_continuous(breaks = c(1,2), labels = c("No Treatment", "Treatment"))+
  xlab("Day") + ylab("log(level)") + 
  facet_wrap(~variable, scales = "free", ncol = 6) +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"))

##############################################
# Correlation of the cluster scores with the CyTOF
##############################################

df_cytof_0 <- read.csv(file = "./data_metabolism/CyTOF.csv", stringsAsFactors = FALSE)
df_cytof_0$X <- NULL
df_cytof_0$File <- NULL
id_cytof <- df_cytof_0$PatientID

id_sample_cytof <- substr(id_cytof, 5, 10)
id_day_cytof <- df_cytof_0$DAY

df_cytof <- df_cytof_0
df_cytof$PatientID <- id_sample_cytof
df_cytof_s <- summaryBy(. ~ PatientID + DAY, df_cytof, FUN = "mean", na.rm = TRUE, keep.names = TRUE)

df_cytof_s$PatientID[df_cytof_s$PatientID == "22-21-"] <- "22-021"
df_cytof_s$PatientID[df_cytof_s$PatientID == "22-22-"] <- "22-022"
df_cytof_s$PatientID[df_cytof_s$PatientID == "22-23-"] <- "22-023"

df_cytof_d0 <- subset(df_cytof_s, DAY == "D000")
df_cytof_d0_2 <- df_cytof_d0[, -(1:2)]
rownames(df_cytof_d0_2) <- df_cytof_d0$PatientID

df_cytof_d15 <- subset(df_cytof_s, DAY == "D015")
df_cytof_d15_2 <- df_cytof_d15[, -(1:2)]
rownames(df_cytof_d15_2) <- df_cytof_d15$PatientID

id.in.cytof <- intersect(df_cytof_d0$PatientID, df_cytof_d15$PatientID)
id.in.cytof <- intersect(id.in.cytof, names(pc.kyn))

df_cytof_diff <- (df_cytof_d15_2[id.in.cytof,] - df_cytof_d0_2[id.in.cytof,])

pc.kyn.cytof <- pc.kyn[id.in.cytof]
pc.ser.cytof <- pc.ser[id.in.cytof]
pc.nic.cytof <- pc.nic[id.in.cytof]

df.cluster.cytof <- data.frame(Kyn = pc.kyn.cytof,
                               Ser = pc.ser.cytof,
                               Nic = pc.nic.cytof)
colnames(df_cytof_diff)[colnames(df_cytof_diff) == "IL.2"] <- "IL-2"

cor.cytof <- corr.test(df.cluster.cytof, df_cytof_diff, method = "spearman", adjust = "none")
p.cor.cytof <- cor.cytof$p
r.cor.cytof <- cor.cytof$r
  cytof_ci <- cor.cytof$ci

fdr <- p.adjust(cytof_ci$p, method = "BH")
cytof_ci$fdr <- fdr
p_adj_corr <- p.adjust(cytof_ci$p, method = "BH")

paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "orangered"))(paletteLength)
myBreaks <- c(seq(min(r.cor.cytof), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(r.cor.cytof)/paletteLength, max(r.cor.cytof), length.out=floor(paletteLength/2)))

colnames(r.cor.cytof) <- gsub("\\.", "-", colnames(r.cor.cytof))

pheatmap(r.cor.cytof, scale = 'none', cluster_rows = FALSE,
         treeheight_col = 0, main = "CyTOF",
         color = myColor, breaks = myBreaks)


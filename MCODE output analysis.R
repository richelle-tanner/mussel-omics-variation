#process MCODE cluster analysis
#uses node table that is created after using MCODE, not the export from MCODE!
#MCODE paramters: node cluster score = 0.5, k core = 3
#import edge table as well to look at connections 

#first priority is figuring out the gene ontology of the genes within each cluster
#can also look at edges only within the cluster, and then where they reach outside of the cluster
#eigenvalue analysis for the list of clustered genes within each cluster


library(ggplot2)
library(evolqg)
library(vegan)
library(dplyr)
library(reshape2)
library(RVAideMemoire)
library(stringr)
library(effsize)
library(GOfuncR)

####import and cleaning data####
#import data
ddr <- "~/Documents/Cytoscape_projects/M_Cali/"
protein_node <- read.csv(paste(ddr, "allTreat_proteinIG default node_6.16.20.csv", sep=''))
gene_node <- read.csv(paste(ddr, "allTreat_geneIG default node_6.16.20.csv", sep=''))
gene_matrix <- read.delim("~/Documents/R_Data/musselRNASeq/final_files/gene_matrixllsImpute_6.8.20.txt")
protein_matrix <- read.delim("~/Documents/R_Data/musselRNASeq/final_files/protein_matrixllsImpute_6.8.20.txt") #use lls impute method 
annotation <- read.csv("~/Documents/R_Data/musselRNASeq/final_files/Blast2GO_table_sepGOIDs.csv", header=TRUE)
KEGG_annotation <- read.table("~/Documents/R_Data/musselRNASeq/reduced_genes_proteins_pathways.txt", fill=TRUE, header=TRUE)


# #get rid of redundant columns
protein_node <- protein_node[,c(4,6,7,10)]
gene_node <- gene_node[,c(3,5,6,9)]

gene_node <- gene_node[,c(1,4)]
protein_node <- protein_node[,c(1,4)]

#names of genes
names <- data.frame(gene=row.names(gene_matrix))

#reduce gene_node and protein_node to the 1519 genes
colnames(gene_node) <- c("MCODE_Cluster", "gene")
colnames(protein_node) <- c("MCODE_Cluster", "gene")

#scale the data here
gene_matrix_scaled <- scale(t(gene_matrix))
protein_matrix_scaled <- scale(t(protein_matrix))
#flip it back for now so that everything downstream works
gene_matrix <- t(gene_matrix_scaled)
protein_matrix <- t(protein_matrix_scaled)

#data frames for each cluster
protein_cluster_names  <- row.names(table(protein_node$MCODE_Cluster))
protein_cluster_names <- protein_cluster_names[-1]

protein_cluster_list <- setNames(replicate(length(protein_cluster_names),data.frame()),protein_cluster_names)
for (i in 1:length(protein_cluster_names)){
  protein_cluster_list[[i]] <- protein_node[which(protein_node$MCODE_Cluster == protein_cluster_names[i]),]
}

gene_cluster_names  <- row.names(table(gene_node$MCODE_Cluster))
gene_cluster_names <- gene_cluster_names[-1]

gene_cluster_list <- setNames(replicate(length(gene_cluster_names),data.frame()),gene_cluster_names)
for (i in 1:length(gene_cluster_names)){
  gene_cluster_list[[i]] <- gene_node[which(gene_node$MCODE_Cluster == gene_cluster_names[i]),]
}

#create new DFs for each cluster for expression values
gene_cluster_names_d <- c(paste(gene_cluster_names, "gene", sep='_'), paste(gene_cluster_names, "protein", sep='_'))
gene_expMat_cluster_list <- setNames(replicate(length(gene_cluster_names_d),data.frame()),gene_cluster_names_d)
for (i in 1:length(gene_cluster_names)){
  transcripts <- gene_cluster_list[[i]]$gene
  gene_expMat_cluster_list[[i]] <- gene_matrix[which(row.names(gene_matrix) %in% transcripts),]
  gene_expMat_cluster_list[[i+6]] <- protein_matrix[which(row.names(protein_matrix) %in% transcripts),]
}

protein_cluster_names_d <- c(paste(protein_cluster_names, "gene", sep='_'), paste(protein_cluster_names, "protein", sep='_'))
protein_expMat_cluster_list <- setNames(replicate(length(protein_cluster_names_d),data.frame()),protein_cluster_names_d)
for (i in 1:length(protein_cluster_names)){
  transcripts <- protein_cluster_list[[i]]$gene
  protein_expMat_cluster_list[[i]] <- gene_matrix[which(row.names(gene_matrix) %in% transcripts),]
  protein_expMat_cluster_list[[i+7]] <- protein_matrix[which(row.names(protein_matrix) %in% transcripts),]
}


#create new DFs for each treatment for each cluster
ind_treatments <- c("FAE", "POL", "FAE", "FAP", "POH", "CGP", "POH", "FAE", "FAP", "FAE", "POL", "FAP", "FAP", "FAE",
                    "FAP", "POH", "CGP", "CGP", "POL", "FAE", "FAP", "POH", "POH", "POL", "CGP", "POL", "POL", "POH",
                    "POH", "POH", "POL", "POL", "FAE", "FAP", "POH", "CGP", "FAE", "FAP", "POH", "FAE", "FAP")
CGP_ind <- which(ind_treatments %in% "CGP")
FAE_ind <- which(ind_treatments %in% "FAE")
FAP_ind <- which(ind_treatments %in% "FAP")
POH_ind <- which(ind_treatments %in% "POH")
POL_ind <- which(ind_treatments %in% "POL")
treatments <- c("CGP", "FAE", "FAP", "POH", "POL")

gene_treat_cluster <- expand.grid(treatments, gene_cluster_names_d)
gene_treat_cluster$names <- paste(gene_treat_cluster$Var1, "_", gene_treat_cluster$Var2, "_gene", sep='')
gene_treat_cluster <- gene_treat_cluster$names

counter <- 0
gene_expMat_cluster_treat_list <- setNames(replicate(length(gene_treat_cluster), data.frame()), gene_treat_cluster)
for (i in 1:length(gene_cluster_names_d)){
  gene_expMat_cluster_treat_list[[i+counter]] <- gene_expMat_cluster_list[[i]][,CGP_ind]
  gene_expMat_cluster_treat_list[[(i+counter+1)]] <- gene_expMat_cluster_list[[i]][,FAE_ind]
  gene_expMat_cluster_treat_list[[(i+counter+2)]] <- gene_expMat_cluster_list[[i]][,FAP_ind]
  gene_expMat_cluster_treat_list[[(i+counter+3)]] <- gene_expMat_cluster_list[[i]][,POH_ind]
  gene_expMat_cluster_treat_list[[(i+counter+4)]] <- gene_expMat_cluster_list[[i]][,POL_ind]
  counter <- counter+4
}

protein_treat_cluster <- expand.grid(treatments, protein_cluster_names_d)
protein_treat_cluster$names <- paste(protein_treat_cluster$Var1, "_", protein_treat_cluster$Var2, "_protein", sep='')
protein_treat_cluster <- protein_treat_cluster$names

counter <- 0
protein_expMat_cluster_treat_list <- setNames(replicate(length(protein_treat_cluster), data.frame()), protein_treat_cluster)
for (i in 1:length(protein_cluster_names_d)){
  protein_expMat_cluster_treat_list[[i+counter]] <- protein_expMat_cluster_list[[i]][,CGP_ind]
  protein_expMat_cluster_treat_list[[(i+counter+1)]] <- protein_expMat_cluster_list[[i]][,FAE_ind]
  protein_expMat_cluster_treat_list[[(i+counter+2)]] <- protein_expMat_cluster_list[[i]][,FAP_ind]
  protein_expMat_cluster_treat_list[[(i+counter+3)]] <- protein_expMat_cluster_list[[i]][,POH_ind]
  protein_expMat_cluster_treat_list[[(i+counter+4)]] <- protein_expMat_cluster_list[[i]][,POL_ind]
  counter <- counter+4
}

####calculate MAD and distribution of of MAD within genes in each cluster####
MAD_function <- function(data_input) {
  MAD <- vector(length=length(data_input[1,]))
  #calculate MAD per gene/protein
  for (i in 1:length(data_input[1,])) {
    sub <- as.numeric(as.character(data_input[,i]))
    #find median for each gene
    median <- median(sub)
    #calculate deviations
    deviations <- sub - median
    #calculate and store MAD
    MAD[i] <- median(abs(deviations))
  }
  return(MAD)
}

#this is genes (rows) x treatment (columns)
counter <- 0
gene_cluster_MAD_list <- setNames(replicate(length(gene_cluster_names_d),data.frame()),paste("MAD_", gene_cluster_names_d, sep=''))
for (i in 1:length(gene_cluster_MAD_list)){
  CGP_gene_MAD <- MAD_function(t(gene_expMat_cluster_treat_list[[i+counter]])) #got rid of scale here, did it above on the whole dataset 2.21.20
  FAE_gene_MAD <- MAD_function(t(gene_expMat_cluster_treat_list[[i+counter+1]]))
  FAP_gene_MAD <- MAD_function(t(gene_expMat_cluster_treat_list[[i+counter+2]]))
  POH_gene_MAD <- MAD_function(t(gene_expMat_cluster_treat_list[[i+counter+3]]))
  POL_gene_MAD <- MAD_function(t(gene_expMat_cluster_treat_list[[i+counter+4]]))
  
  MAD_gene_df <- data.frame(cbind(CGP_gene_MAD, FAE_gene_MAD, FAP_gene_MAD, POH_gene_MAD, POL_gene_MAD))
  gene_cluster_MAD_list[[i]] <- MAD_gene_df
  row.names(gene_cluster_MAD_list[[i]]) <- row.names(gene_expMat_cluster_list[[i]])
  
  counter <- counter+4
}

counter <- 0
protein_cluster_MAD_list <- setNames(replicate(length(protein_cluster_names_d),data.frame()),paste("MAD_", protein_cluster_names_d, sep=''))
for (i in 1:length(protein_cluster_MAD_list)){
  CGP_protein_MAD <- MAD_function(t(protein_expMat_cluster_treat_list[[i+counter]]))
  FAE_protein_MAD <- MAD_function(t(protein_expMat_cluster_treat_list[[i+counter+1]]))
  FAP_protein_MAD <- MAD_function(t(protein_expMat_cluster_treat_list[[i+counter+2]]))
  POH_protein_MAD <- MAD_function(t(protein_expMat_cluster_treat_list[[i+counter+3]]))
  POL_protein_MAD <- MAD_function(t(protein_expMat_cluster_treat_list[[i+counter+4]]))
  
  MAD_protein_df <- data.frame(cbind(CGP_protein_MAD, FAE_protein_MAD, FAP_protein_MAD, POH_protein_MAD, POL_protein_MAD))
  protein_cluster_MAD_list[[i]] <- MAD_protein_df
  row.names(protein_cluster_MAD_list[[i]]) <- row.names(protein_expMat_cluster_list[[i]])
  
  counter <- counter+4
}

#this prints the protein ones as clusters 8-14, need to change it manually
#also it goes to the Cytoscape folder
for (i in 1:length(gene_cluster_MAD_list)) {
  data <- gene_cluster_MAD_list[[i]]
  plot <- ggplot(data=data, aes(x=data[,1])) +
    geom_density(fill='blue', aes(x=data[,1]), alpha=.5) +
    geom_density(fill='red', aes(x=data[,2]), alpha=.5) +
    geom_density(fill='green', aes(x=data[,3]), alpha=.5) +
    geom_density(fill='purple', aes(x=data[,4]), alpha=.5) +
    geom_density(fill='yellow', aes(x=data[,5]), alpha=.5) +
    geom_rug() +
    labs(x="MAD score") +
    theme_bw()
  filename <- paste(ddr, "gene_density_cluster", i, "plot_6.17.20.pdf", sep='')
  ggsave(filename, plot=plot)
}


for (i in 1:length(protein_cluster_MAD_list)) {
  data <- protein_cluster_MAD_list[[i]]
  plot <- ggplot(data=data, aes(x=data[,1])) +
    geom_density(fill='blue', aes(x=data[,1]), alpha=.5) +
    geom_density(fill='red', aes(x=data[,2]), alpha=.5) +
    geom_density(fill='green', aes(x=data[,3]), alpha=.5) +
    geom_density(fill='purple', aes(x=data[,4]), alpha=.5) +
    geom_density(fill='yellow', aes(x=data[,5]), alpha=.5) +
    geom_rug() +
    labs(x="MAD score") +
    theme_bw()
  filename <- paste(ddr, "protein_density_cluster", i, "plot_6.17.20.pdf", sep='')
  ggsave(filename, plot=plot)
}

####eigenvalue calculation####
#make covariance matrices
#need gene/protein matrices that are subsetted to the specific clusters
#use protein_expMat_cluster_treat_list and gene_expMat_cluster_treat_list
library(evolqg)

gene_eigen <- setNames(replicate(length(gene_treat_cluster),data.frame()),paste("eigenvals_", gene_treat_cluster, sep=''))
gene_eigenStat_cluster <- data.frame(treatment=rep(c("CGP", "FAE", "FAP", "POH", "POL"), 12),
                                ICV=rep(0, length(gene_eigen)), SD=rep(0, length(gene_eigen)))

for (i in 1:length(gene_eigen)){
  df <- t(gene_expMat_cluster_treat_list[[i]])           
  for (j in 1:length(df[1,])) {
    if (sum(df[,j]) == 0 | sum(df[,j]) == "NaN") {
      df[,j] <- NA
    }
  }
  dropset <- colnames(df)[which(colSums(is.na(df)) != 0)]
  if (length(dropset) > 0) {
    dropset_colnum <- which(colnames(df) %in% dropset)
    dropset_vals <- data.frame(num=dropset_colnum, val=rep("NA", length(dropset_colnum)))
    keepset_colnum <- colnames(df)[-dropset_colnum]
    df <- df[ , apply(df, 2, function(x) !any(is.na(x)))]
    #eigenvalu calcs
    covar <- cov(df)
    
    eigenvals <- eigen(covar)
    #insert NAs where we had to take columns out
    keepset_vals <- data.frame(num=keepset_colnum, val=eigenvals$values)
    eigenvals_full <- rbind(keepset_vals, dropset_vals)
    gene_eigen[[i]] <- eigenvals_full$val
  }
  if (length(dropset) == 0) {
    covar <- cov(df)
    eigenvals <- eigen(covar)
    gene_eigen[[i]] <- eigenvals$values
  }
  gene_eigenStat_cluster[i,2] <- CalcICV(covar)
  gene_eigenStat_cluster[i,3] <- CalcEigenVar(covar, sd=TRUE, sample=length(eigenvals$values))
}

protein_eigen <- setNames(replicate(length(protein_treat_cluster),data.frame()),paste("eigenvals_", protein_treat_cluster, sep=''))
protein_eigenStat_cluster <- data.frame(treatment=rep(c("CGP", "FAE", "FAP", "POH", "POL"), 14),
                                     ICV=rep(0, length(protein_eigen)), SD=rep(0, length(protein_eigen)))

for (i in 1:length(protein_eigen)){
  df <- t(protein_expMat_cluster_treat_list[[i]])
  for (j in 1:length(df[1,])) {
    if (sum(df[,j]) == 0 | sum(df[,j] == "NaN")) {
      df[,j] <- NA
    }
  }
  #keep track of which columns those are, and take them out so the rest will run
  dropset <- colnames(df)[which(colSums(is.na(df)) != 0)]
  if (length(dropset) > 0) {
    dropset_colnum <- which(colnames(df) %in% dropset)
    dropset_vals <- data.frame(num=dropset_colnum, val=rep("NA", length(dropset_colnum)))
    keepset_colnum <- colnames(df)[-dropset_colnum]
    df <- df[ , apply(df, 2, function(x) !any(is.na(x)))]
    #eigenvalu calcs
    covar <- cov(df)
    
    eigenvals <- eigen(covar)
    #insert NAs where we had to take columns out
    keepset_vals <- data.frame(num=keepset_colnum, val=eigenvals$values)
    eigenvals_full <- rbind(keepset_vals, dropset_vals)
    protein_eigen[[i]] <- eigenvals_full$val
  }
  if (length(dropset) == 0) {
    covar <- cov(df)
    eigenvals <- eigen(covar)
    protein_eigen[[i]] <- eigenvals$values
  }
  protein_eigenStat_cluster[i,2] <- CalcICV(covar)
  protein_eigenStat_cluster[i,3] <- CalcEigenVar(covar, sd=TRUE, sample=length(eigenvals$values))
}

protein_eigenStat_cluster$cluster <- rep(seq(1,7), each=10)
protein_eigenStat_cluster$type <- rep(rep(c("gene", "protein"), each=5), 7)

gene_eigenStat_cluster$cluster <- rep(seq(1,6), each=10)
gene_eigenStat_cluster$type <- rep(rep(c("gene", "protein"), each=5), 6)

#export SD and ICV stats
write.table(gene_eigenStat_cluster, "~/geneClusters_eigenStat.txt", sep="\t", row.names = F, quote = F)
write.table(protein_eigenStat_cluster, "~/Eigenvalue plots/proteinClusters_eigenStat.txt", sep="\t", row.names = F, quote = F)

#calculate eigenvalue % variance explained
perc_var_gene <- vector()
for (i in 1:length(gene_eigen)) {
  gene_eigen[[i]] <- as.numeric(gene_eigen[[i]])
  perc_var_gene[i] <- gene_eigen[[i]][1]/sum(na.omit(gene_eigen[[i]]))
}
perc_var_gene_table <- data.frame(perc_var=perc_var_gene, treatment=rep(c("CGP", "FAE", "FAP", "POH", "POL"),12), 
                                  type=rep(c("gene", "protein"),each=30), cluster=rep(rep(c(1:6), each=5), 2))

perc_var_protein <- vector()
for (i in 1:length(protein_eigen)) {
  protein_eigen[[i]] <- as.numeric(protein_eigen[[i]])
  perc_var_protein[i] <- protein_eigen[[i]][1]/sum(na.omit(protein_eigen[[i]]))
}
perc_var_protein_table <- data.frame(perc_var=perc_var_protein, treatment=rep(c("CGP", "FAE", "FAP", "POH", "POL"),14), 
                                  type=rep(c("gene", "protein"),each=35), cluster=rep(rep(c(1:7), each=5), 2))

#put the list elements into one data frame for plotting
plot_gene_eigen <- melt(gene_eigen, level=1)
plot_protein_eigen <- melt(protein_eigen, level=1)

#extra columns for labeling
colnames(plot_gene_eigen) <- c("eigenvalue", "list_label")
colnames(plot_protein_eigen) <- c("eigenvalue", "list_label")
plot_gene_eigen$list_label1 <- str_sub(plot_gene_eigen$list_label, start=11)
plot_gene_eigen$treatment <- str_sub(plot_gene_eigen$list_label1, 1, 3)
plot_gene_eigen$cluster <- str_sub(plot_gene_eigen$list_label1, 13, 13)
plot_protein_eigen$list_label1 <- str_sub(plot_protein_eigen$list_label, start=11)
plot_protein_eigen$treatment <- str_sub(plot_protein_eigen$list_label1, 1, 3)
plot_protein_eigen$cluster <- str_sub(plot_protein_eigen$list_label1, 13, 13)
plot_gene_eigen$data_type <- str_sub(plot_gene_eigen$list_label1, 15, 18)
plot_protein_eigen$data_type <- str_sub(plot_protein_eigen$list_label1, 15, 18)
plot_gene_eigen$eigenvalue <- as.numeric(plot_gene_eigen$eigenvalue)
plot_protein_eigen$eigenvalue <- as.numeric(plot_protein_eigen$eigenvalue)

#separate out into clusters again, and assign the x value

gene_eigenPlot_cluster1 <- subset(plot_gene_eigen, cluster=="1")
gene_eigenPlot_cluster1$x <- rep(c(1:410), 10)

gene_eigenPlot_cluster2 <- subset(plot_gene_eigen, cluster=="2")
gene_eigenPlot_cluster2$x <- rep(c(1:321), 10)

gene_eigenPlot_cluster3 <- subset(plot_gene_eigen, cluster=="3")
gene_eigenPlot_cluster3$x <- rep(c(1:126), 10)

gene_eigenPlot_cluster4 <- subset(plot_gene_eigen, cluster=="4")
gene_eigenPlot_cluster4$x <- rep(c(1:91), 10)

gene_eigenPlot_cluster5 <- subset(plot_gene_eigen, cluster=="5")
gene_eigenPlot_cluster5$x <- rep(c(1:6), 10)

gene_eigenPlot_cluster6 <- subset(plot_gene_eigen, cluster=="6")
gene_eigenPlot_cluster6$x <- rep(c(1:5), 10)

#plot
#colors
eigenPlot_colorMan <- c("blue", "blue", "green", "green", "orange", "orange", "red", "red", "purple", "purple")

#only show the first 15 values
gene_eigenPlot_cluster1$plot_label <- paste(gene_eigenPlot_cluster1$treatment, gene_eigenPlot_cluster1$data_type, sep="_")
ggplot(gene_eigenPlot_cluster1, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_cluster1.pdf")

gene_eigenPlot_cluster2$plot_label <- paste(gene_eigenPlot_cluster2$treatment, gene_eigenPlot_cluster2$data_type, sep="_")
ggplot(gene_eigenPlot_cluster2, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_cluster2.pdf")

gene_eigenPlot_cluster3$plot_label <- paste(gene_eigenPlot_cluster3$treatment, gene_eigenPlot_cluster3$data_type, sep="_")
ggplot(gene_eigenPlot_cluster3, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_cluster3.pdf")

gene_eigenPlot_cluster4$plot_label <- paste(gene_eigenPlot_cluster4$treatment, gene_eigenPlot_cluster4$data_type, sep="_")
ggplot(gene_eigenPlot_cluster4, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_cluster4.pdf")

gene_eigenPlot_cluster5$plot_label <- paste(gene_eigenPlot_cluster5$treatment, gene_eigenPlot_cluster5$data_type, sep="_")
ggplot(gene_eigenPlot_cluster5, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,6) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_cluster5.pdf")

gene_eigenPlot_cluster6$plot_label <- paste(gene_eigenPlot_cluster6$treatment, gene_eigenPlot_cluster6$data_type, sep="_")
ggplot(gene_eigenPlot_cluster6, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,5) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_cluster6.pdf")



#separate out into clusters again, and assign the x value
#42, 11, 78, 8, 6, 6, 61, 11

protein_eigenPlot_cluster1 <- subset(plot_protein_eigen, cluster=="1")
protein_eigenPlot_cluster1$x <- rep(c(1:52), 10)

protein_eigenPlot_cluster2 <- subset(plot_protein_eigen, cluster=="2")
protein_eigenPlot_cluster2$x <- rep(c(1:11), 10)

protein_eigenPlot_cluster3 <- subset(plot_protein_eigen, cluster=="3")
protein_eigenPlot_cluster3$x <- rep(c(1:16), 10)

protein_eigenPlot_cluster4 <- subset(plot_protein_eigen, cluster=="4")
protein_eigenPlot_cluster4$x <- rep(c(1:8), 10)

protein_eigenPlot_cluster5 <- subset(plot_protein_eigen, cluster=="5")
protein_eigenPlot_cluster5$x <- rep(c(1:55), 10)

protein_eigenPlot_cluster6 <- subset(plot_protein_eigen, cluster=="6")
protein_eigenPlot_cluster6$x <- rep(c(1:4), 10)

protein_eigenPlot_cluster7 <- subset(plot_protein_eigen, cluster=="7")
protein_eigenPlot_cluster7$x <- rep(c(1:30), 10)

#plot
#only show the first 15 values
protein_eigenPlot_cluster1$plot_label <- paste(protein_eigenPlot_cluster1$treatment, protein_eigenPlot_cluster1$data_type, sep="_")
ggplot(protein_eigenPlot_cluster1, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_protein_cluster1.pdf")

protein_eigenPlot_cluster2$plot_label <- paste(protein_eigenPlot_cluster2$treatment, protein_eigenPlot_cluster2$data_type, sep="_")
ggplot(protein_eigenPlot_cluster2, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,11) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_protein_cluster2.pdf")

protein_eigenPlot_cluster3$plot_label <- paste(protein_eigenPlot_cluster3$treatment, protein_eigenPlot_cluster3$data_type, sep="_")
ggplot(protein_eigenPlot_cluster3, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_protein_cluster3.pdf")

protein_eigenPlot_cluster4$plot_label <- paste(protein_eigenPlot_cluster4$treatment, protein_eigenPlot_cluster4$data_type, sep="_")
ggplot(protein_eigenPlot_cluster4, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,8) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_protein_cluster4.pdf")

protein_eigenPlot_cluster5$plot_label <- paste(protein_eigenPlot_cluster5$treatment, protein_eigenPlot_cluster5$data_type, sep="_")
ggplot(protein_eigenPlot_cluster5, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_protein_cluster5.pdf")

protein_eigenPlot_cluster6$plot_label <- paste(protein_eigenPlot_cluster6$treatment, protein_eigenPlot_cluster6$data_type, sep="_")
ggplot(protein_eigenPlot_cluster6, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,4) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_protein_cluster6.pdf")

protein_eigenPlot_cluster7$plot_label <- paste(protein_eigenPlot_cluster7$treatment, protein_eigenPlot_cluster7$data_type, sep="_")
ggplot(protein_eigenPlot_cluster7, aes(x=x, y=eigenvalue, group=plot_label)) + 
  geom_line(aes(color =plot_label, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_protein_cluster7.pdf")

####GO enrichment in each cluster####
#export list of all genes in all clusters to upload to batch entrez (https://www.ncbi.nlm.nih.gov/sites/batchentrez)
#v9 in the annotation file is the NCBI ID from the Mollusca database
#use all of the column names from gene_expMat_cluster_list and protein_expMat_cluster_list

GO_sort_function <-function(names_input) {
  #first we match names
  GO_names <- annotation[which(annotation$SeqName %in% names_input),]
  
  #then we pull out GO.IDs
  GO_names <- GO_names[,-c(4, 3)]
  
  #pull out BLAST
  BLAST <- GO_names[,c(1:2)]
  
  #now we get all of the unique GO IDs
  #there are 1344
  unique_GO <- sapply(GO_names[,c(4:40)], function(x) unique(x))
  unique_GO_vector <- unique(unlist(unique_GO))
  #remove the empty one
  unique_GO_vector <- unique_GO_vector[unique_GO_vector != ""]
  
  #get a count of how many times each one occurs
  #rows are the sequences, columns are the GO terms
  GO_ID_list <- matrix(nrow=length(GO_names$SeqName), ncol=length(unique_GO_vector))
  for (i in 1:length(unique_GO_vector)) {
    Values <- unique_GO_vector[i]
    GO_ID_list[,i] <- rowSums(GO_names == as.character(Values), na.rm = T)
  }
  
  all_genes <- as.character(GO_names$SeqName)
  
  gene_list <- data.frame()
  for (j in 1:length(unique_GO_vector)) {
    ind_term <- GO_ID_list[,j]
    row_genes <- all_genes[which(ind_term %in% 1)]
    row_BLAST <- as.character(BLAST[which(ind_term %in% 1),2])
    gene_list_ind <- character()
    BLAST_list_ind <- character()
    if (length(row_genes) > 0) {
      for (k in 1:length(row_genes)) {
        gene_list_ind <- as.character(paste(gene_list_ind, ", ", row_genes[k], sep=''))
        BLAST_list_ind <- as.character(paste(BLAST_list_ind, ", ", row_BLAST[k], sep=''))
      }
      gene_list_ind <- substring(gene_list_ind, 3)
      BLAST_list_ind <- substring(BLAST_list_ind, 3)
      gene_list[j,1] <- unique_GO_vector[j]
      gene_list[j,2] <- gene_list_ind
      gene_list[j,3] <- BLAST_list_ind
    }
  }
  
  #put row and column names back
  rownames(GO_ID_list) <- GO_names$SeqName
  colnames(GO_ID_list) <- unique_GO_vector

  
  #figure out which column has the most instances
  GO_ID_count <- colSums(GO_ID_list)
  sorted_GO_ID_count <- GO_ID_count[order(-GO_ID_count)] 
  #top 1% of terms: 13
  top1 <- round(length(unique_GO_vector)*.01)
  top_sorted_GO_ID_count <- sorted_GO_ID_count[1:top1]
  top10 <- round(length(unique_GO_vector)*.1)
  top10_sorted_GO_ID_count <- sorted_GO_ID_count[1:top10]
  list <- list(top_sorted_GO_ID_count, top10_sorted_GO_ID_count, sorted_GO_ID_count, gene_list)
  return(list)
}

#do it by cluster
#there are 7 gene clusters, 8 protein clusters
#this is the format: rownames(gene_expMat_cluster_list[[1]])
cluster1g_names <- rownames(gene_expMat_cluster_list[[1]])
gene_cluster1 <- GO_sort_function(cluster1g_names)
cluster2g_names <- rownames(gene_expMat_cluster_list[[2]])
gene_cluster2 <- GO_sort_function(cluster2g_names)
cluster3g_names <- rownames(gene_expMat_cluster_list[[3]])
gene_cluster3 <- GO_sort_function(cluster3g_names)
cluster4g_names <- rownames(gene_expMat_cluster_list[[4]])
gene_cluster4 <- GO_sort_function(cluster4g_names)
cluster5g_names <- rownames(gene_expMat_cluster_list[[5]])
gene_cluster5 <- GO_sort_function(cluster5g_names)
cluster6g_names <- rownames(gene_expMat_cluster_list[[6]])
gene_cluster6 <- GO_sort_function(cluster6g_names)


cluster1p_names <- rownames(protein_expMat_cluster_list[[1]])
protein_cluster1 <- GO_sort_function(cluster1p_names)
cluster2p_names <- rownames(protein_expMat_cluster_list[[2]])
protein_cluster2 <- GO_sort_function(cluster2p_names)
cluster3p_names <- rownames(protein_expMat_cluster_list[[3]])
protein_cluster3 <- GO_sort_function(cluster3p_names)
cluster4p_names <- rownames(protein_expMat_cluster_list[[4]])
protein_cluster4 <- GO_sort_function(cluster4p_names)
cluster5p_names <- rownames(protein_expMat_cluster_list[[5]])
protein_cluster5 <- GO_sort_function(cluster5p_names)
cluster6p_names <- rownames(protein_expMat_cluster_list[[6]])
protein_cluster6 <- GO_sort_function(cluster6p_names)
cluster7p_names <- rownames(protein_expMat_cluster_list[[7]])
protein_cluster7 <- GO_sort_function(cluster7p_names)


#export all of this in a comprehensive data table
#columns: cluster (+ gene/protein designation), GO term, count
export_GO <- data.frame(cluster=c(rep("gene1", length(gene_cluster1[[3]])), rep("gene2", length(gene_cluster2[[3]])),
                                  rep("gene3", length(gene_cluster3[[3]])), rep("gene4", length(gene_cluster4[[3]])),
                                  rep("gene5", length(gene_cluster5[[3]])), rep("gene6", length(gene_cluster6[[3]])),
                                  rep("protein1", length(protein_cluster1[[3]])),
                                  rep("protein2", length(protein_cluster2[[3]])), rep("protein3", length(protein_cluster3[[3]])),
                                  rep("protein4", length(protein_cluster4[[3]])), rep("protein5", length(protein_cluster5[[3]])),
                                  rep("protein6", length(protein_cluster6[[3]])), rep("protein7", length(protein_cluster7[[3]]))),
                        GO_term=c(names(gene_cluster1[[3]]), names(gene_cluster2[[3]]), names(gene_cluster3[[3]]),
                                  names(gene_cluster4[[3]]), names(gene_cluster5[[3]]), names(gene_cluster6[[3]]),
                                  names(protein_cluster1[[3]]), names(protein_cluster2[[3]]),
                                  names(protein_cluster3[[3]]), names(protein_cluster4[[3]]), names(protein_cluster5[[3]]),
                                  names(protein_cluster6[[3]]), names(protein_cluster7[[3]])),
                        count=c(unname(gene_cluster1[[3]]), unname(gene_cluster2[[3]]), unname(gene_cluster3[[3]]),
                                unname(gene_cluster4[[3]]), unname(gene_cluster5[[3]]), unname(gene_cluster6[[3]]),
                                unname(protein_cluster1[[3]]), unname(protein_cluster2[[3]]),
                                unname(protein_cluster3[[3]]), unname(protein_cluster4[[3]]), unname(protein_cluster5[[3]]),
                                unname(protein_cluster6[[3]]), unname(protein_cluster7[[3]])))

write.table(export_GO, "~/MCODE_GOfreqTable.txt", sep='\t', quote=F, row.names=F)

###############################
#use this: https://rdrr.io/cran/RVAideMemoire/man/fisher.multcomp.html
#multiple comparisons fishers exact test (uses p.adjust)

#get the "GO sort function" results for the 1521 transcripts to use as a null reference
all_names <- as.character(annotation$SeqName)
all_GO_counts <- GO_sort_function(all_names)
all_GO_counts <- all_GO_counts[[3]]
all_GO <- data.frame(GO_term=names(all_GO_counts), count=unname(all_GO_counts))

#contingency tables
#cluster 1 gene, first column is null second is cluster
c1g <- subset(export_GO, cluster == "gene1")
c1g$cluster <- NULL
c1g_cont <- merge(all_GO, c1g, by = "GO_term")
rownames(c1g_cont) <- c1g_cont$GO_term
c1g_cont$GO_term <- NULL
colnames(c1g_cont) <- c("null", "cluster")
c1g_cont <- as.matrix(t(c1g_cont))

c1g_fisher <- vector()
for (i in 1:length(c1g_cont[1,])) {
  test_mat <- matrix(c(c1g_cont[2,i], 410, c1g_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c1g_fisher <- c(c1g_fisher, p_val)
}

#adjust for multiple comparisons - this really kills our p values...
c1g_fisher_adj <- p.adjust(c1g_fisher, "BH")
names(c1g_fisher_adj) <- colnames(c1g_cont)

#cluster 2 gene, first column is null second is cluster
c2g <- subset(export_GO, cluster == "gene2")
c2g$cluster <- NULL
c2g_cont <- merge(all_GO, c2g, by = "GO_term")
rownames(c2g_cont) <- c2g_cont$GO_term
c2g_cont$GO_term <- NULL
colnames(c2g_cont) <- c("null", "cluster")
c2g_cont <- as.matrix(t(c2g_cont))

c2g_fisher <- vector()
for (i in 1:length(c2g_cont[1,])) {
  test_mat <- matrix(c(c2g_cont[2,i], 321, c2g_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c2g_fisher <- c(c2g_fisher, p_val)
}

#adjust for multiple comparisons - this really kills our p values...
c2g_fisher_adj <- p.adjust(c2g_fisher, "BH")
names(c2g_fisher_adj) <- colnames(c2g_cont)

#cluster 3 gene, first column is null second is cluster
c3g <- subset(export_GO, cluster == "gene3")
c3g$cluster <- NULL
c3g_cont <- merge(all_GO, c3g, by = "GO_term")
rownames(c3g_cont) <- c3g_cont$GO_term
c3g_cont$GO_term <- NULL
colnames(c3g_cont) <- c("null", "cluster")
c3g_cont <- as.matrix(t(c3g_cont))

c3g_fisher <- vector()
for (i in 1:length(c3g_cont[1,])) {
  test_mat <- matrix(c(c3g_cont[2,i], 126, c3g_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c3g_fisher <- c(c3g_fisher, p_val)
}

#adjust for multiple comparisons - this really kills our p values...
c3g_fisher_adj <- p.adjust(c3g_fisher, "BH")
names(c3g_fisher_adj) <- colnames(c3g_cont)

#cluster 4 gene, first column is null second is cluster
c4g <- subset(export_GO, cluster == "gene4")
c4g$cluster <- NULL
c4g_cont <- merge(all_GO, c4g, by = "GO_term")
rownames(c4g_cont) <- c4g_cont$GO_term
c4g_cont$GO_term <- NULL
colnames(c4g_cont) <- c("null", "cluster")
c4g_cont <- as.matrix(t(c4g_cont))

c4g_fisher <- vector()
for (i in 1:length(c4g_cont[1,])) {
  test_mat <- matrix(c(c4g_cont[2,i], 91, c4g_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c4g_fisher <- c(c4g_fisher, p_val)
}

#adjust for multiple comparisons - this really kills our p values...
c4g_fisher_adj <- p.adjust(c4g_fisher, "BH")
names(c4g_fisher_adj) <- colnames(c4g_cont)

#cluster 5 gene, first column is null second is cluster
c5g <- subset(export_GO, cluster == "gene5")
c5g$cluster <- NULL
c5g_cont <- merge(all_GO, c5g, by = "GO_term")
rownames(c5g_cont) <- c5g_cont$GO_term
c5g_cont$GO_term <- NULL
colnames(c5g_cont) <- c("null", "cluster")
c5g_cont <- as.matrix(t(c5g_cont))

c5g_fisher <- vector()
for (i in 1:length(c5g_cont[1,])) {
  test_mat <- matrix(c(c5g_cont[2,i], 6, c5g_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c5g_fisher <- c(c5g_fisher, p_val)
}

#adjust for multiple comparisons - this really kills our p values...
c5g_fisher_adj <- p.adjust(c5g_fisher, "BH")
names(c5g_fisher_adj) <- colnames(c5g_cont)

#cluster 6 gene, first column is null second is cluster
c6g <- subset(export_GO, cluster == "gene6")
c6g$cluster <- NULL
c6g_cont <- merge(all_GO, c6g, by = "GO_term")
rownames(c6g_cont) <- c6g_cont$GO_term
c6g_cont$GO_term <- NULL
colnames(c6g_cont) <- c("null", "cluster")
c6g_cont <- as.matrix(t(c6g_cont))

c6g_fisher <- vector()
for (i in 1:length(c6g_cont[1,])) {
  test_mat <- matrix(c(c6g_cont[2,i], 5, c6g_cont[1,i], 1521), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c6g_fisher <- c(c6g_fisher, p_val)
}

#adjust for multiple comparisons
c6g_fisher_adj <- p.adjust(c6g_fisher, "BH")
names(c6g_fisher_adj) <- colnames(c6g_cont)


#cluster 1 protein, first column is null second is cluster
c1p <- subset(export_GO, cluster == "protein1")
c1p$cluster <- NULL
c1p_cont <- merge(all_GO, c1p, by = "GO_term")
rownames(c1p_cont) <- c1p_cont$GO_term
c1p_cont$GO_term <- NULL
colnames(c1p_cont) <- c("null", "cluster")
c1p_cont <- as.matrix(t(c1p_cont))

c1p_fisher <- vector()
for (i in 1:length(c1p_cont[1,])) {
  test_mat <- matrix(c(c1p_cont[2,i], 52, c1p_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c1p_fisher <- c(c1p_fisher, p_val)
}

#adjust for multiple comparisons
c1p_fisher_adj <- p.adjust(c1p_fisher, "BH")
names(c1p_fisher_adj) <- colnames(c1p_cont)

#cluster 2 protein, first column is null second is cluster
c2p <- subset(export_GO, cluster == "protein2")
c2p$cluster <- NULL
c2p_cont <- merge(all_GO, c2p, by = "GO_term")
rownames(c2p_cont) <- c2p_cont$GO_term
c2p_cont$GO_term <- NULL
colnames(c2p_cont) <- c("null", "cluster")
c2p_cont <- as.matrix(t(c2p_cont))

c2p_fisher <- vector()
for (i in 1:length(c2p_cont[1,])) {
  test_mat <- matrix(c(c2p_cont[2,i], 11, c2p_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c2p_fisher <- c(c2p_fisher, p_val)
}

#adjust for multiple comparisons
c2p_fisher_adj <- p.adjust(c2p_fisher, "BH")
names(c2p_fisher_adj) <- colnames(c2p_cont)

#cluster 3 protein, first column is null second is cluster
c3p <- subset(export_GO, cluster == "protein3")
c3p$cluster <- NULL
c3p_cont <- merge(all_GO, c3p, by = "GO_term")
rownames(c3p_cont) <- c3p_cont$GO_term
c3p_cont$GO_term <- NULL
colnames(c3p_cont) <- c("null", "cluster")
c3p_cont <- as.matrix(t(c3p_cont))

c3p_fisher <- vector()
for (i in 1:length(c3p_cont[1,])) {
  test_mat <- matrix(c(c3p_cont[2,i], 16, c3p_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c3p_fisher <- c(c3p_fisher, p_val)
}

#adjust for multiple comparisons
c3p_fisher_adj <- p.adjust(c3p_fisher, "BH")
names(c3p_fisher_adj) <- colnames(c3p_cont)

#cluster 4 protein, first column is null second is cluster
c4p <- subset(export_GO, cluster == "protein4")
c4p$cluster <- NULL
c4p_cont <- merge(all_GO, c4p, by = "GO_term")
rownames(c4p_cont) <- c4p_cont$GO_term
c4p_cont$GO_term <- NULL
colnames(c4p_cont) <- c("null", "cluster")
c4p_cont <- as.matrix(t(c4p_cont))

c4p_fisher <- vector()
for (i in 1:length(c4p_cont[1,])) {
  test_mat <- matrix(c(c4p_cont[2,i], 8, c4p_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c4p_fisher <- c(c4p_fisher, p_val)
}

#adjust for multiple comparisons
c4p_fisher_adj <- p.adjust(c4p_fisher, "BH")
names(c4p_fisher_adj) <- colnames(c4p_cont)

#cluster 5 protein, first column is null second is cluster
c5p <- subset(export_GO, cluster == "protein5")
c5p$cluster <- NULL
c5p_cont <- merge(all_GO, c5p, by = "GO_term")
rownames(c5p_cont) <- c5p_cont$GO_term
c5p_cont$GO_term <- NULL
colnames(c5p_cont) <- c("null", "cluster")
c5p_cont <- as.matrix(t(c5p_cont))

c5p_fisher <- vector()
for (i in 1:length(c5p_cont[1,])) {
  test_mat <- matrix(c(c5p_cont[2,i], 55, c5p_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c5p_fisher <- c(c5p_fisher, p_val)
}

#adjust for multiple comparisons
c5p_fisher_adj <- p.adjust(c5p_fisher, "BH")
names(c5p_fisher_adj) <- colnames(c5p_cont)

#cluster 6 protein, first column is null second is cluster
c6p <- subset(export_GO, cluster == "protein6")
c6p$cluster <- NULL
c6p_cont <- merge(all_GO, c6p, by = "GO_term")
rownames(c6p_cont) <- c6p_cont$GO_term
c6p_cont$GO_term <- NULL
colnames(c6p_cont) <- c("null", "cluster")
c6p_cont <- as.matrix(t(c6p_cont))

c6p_fisher <- vector()
for (i in 1:length(c6p_cont[1,])) {
  test_mat <- matrix(c(c6p_cont[2,i], 4, c6p_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c6p_fisher <- c(c6p_fisher, p_val)
}

#adjust for multiple comparisons
c6p_fisher_adj <- p.adjust(c6p_fisher, "BH")
names(c6p_fisher_adj) <- colnames(c6p_cont)

#cluster 7 protein, first column is null second is cluster
c7p <- subset(export_GO, cluster == "protein7")
c7p$cluster <- NULL
c7p_cont <- merge(all_GO, c7p, by = "GO_term")
rownames(c7p_cont) <- c7p_cont$GO_term
c7p_cont$GO_term <- NULL
colnames(c7p_cont) <- c("null", "cluster")
c7p_cont <- as.matrix(t(c7p_cont))

c7p_fisher <- vector()
for (i in 1:length(c7p_cont[1,])) {
  test_mat <- matrix(c(c7p_cont[2,i], 30, c7p_cont[1,i], 1519), 2, 2)
  p_val <- fisher.test(test_mat, alternative = "greater")$p.value
  c7p_fisher <- c(c7p_fisher, p_val)
}

#adjust for multiple comparisons
c7p_fisher_adj <- p.adjust(c7p_fisher, "BH")
names(c7p_fisher_adj) <- colnames(c7p_cont)


#combine all of the p values
#also add cluster names
fisher_vals <- c(unname(c1g_fisher_adj), unname(c2g_fisher_adj), unname(c3g_fisher_adj), unname(c4g_fisher_adj), unname(c5g_fisher_adj), unname(c6g_fisher_adj),
                 unname(c1p_fisher_adj), unname(c2p_fisher_adj), unname(c3p_fisher_adj), unname(c4p_fisher_adj), unname(c5p_fisher_adj),
                 unname(c6p_fisher_adj), unname(c7p_fisher_adj))
fisher_names <- c(names(c1g_fisher_adj), names(c2g_fisher_adj), names(c3g_fisher_adj), names(c4g_fisher_adj), names(c5g_fisher_adj), names(c6g_fisher_adj),
                 names(c1p_fisher_adj), names(c2p_fisher_adj), names(c3p_fisher_adj), names(c4p_fisher_adj), names(c5p_fisher_adj),
                 names(c6p_fisher_adj), names(c7p_fisher_adj))
cluster_names <- c(rep("c1g", length(c1g_fisher_adj)), rep("c2g", length(c2g_fisher_adj)), rep("c3g", length(c3g_fisher_adj)), 
                   rep("c4g", length(c4g_fisher_adj)), rep("c5g", length(c5g_fisher_adj)), rep("c6g", length(c6g_fisher_adj)),
                   rep("c1p", length(c1p_fisher_adj)), rep("c2p", length(c2p_fisher_adj)), 
                   rep("c3p", length(c3p_fisher_adj)), rep("c4p", length(c4p_fisher_adj)), rep("c5p", length(c5p_fisher_adj)),
                   rep("c6p", length(c6p_fisher_adj)), rep("c7p", length(c7p_fisher_adj)))

fisher_adj_table <- data.frame(Cluster=cluster_names, GO_term=fisher_names, p_val=fisher_vals)

#pull out the significant ones
fisher_sig_all <- data.frame()
for (i in 1:length(fisher_adj_table$p_val)) {
  if (fisher_adj_table$p_val[i] < 0.05) {
    fisher_sig_all <- rbind(fisher_sig_all, fisher_adj_table[i,])
  }
}

#make the tables pretty for publication
enrich_all_clusterGO <- substring(fisher_adj_table$GO_term, 4)
fisher_adj_table$GO_names <- get_names(enrich_all_clusterGO)[2]$go_name
enrich_sig_clusterGO <- substring(fisher_sig_all$GO_term, 4)
fisher_sig_all$GO_names <- get_names(enrich_sig_clusterGO)[2]$go_name

fisher_adj_table$cluster_num <- substring(fisher_adj_table$Cluster, 2, 2)
fisher_adj_table$cluster_type <- substring(fisher_adj_table$Cluster, 3, 3)
fisher_adj_table$cluster_type <- gsub("g", "Gene", fisher_adj_table$cluster_type)
fisher_adj_table$cluster_type <- gsub("p", "Protein", fisher_adj_table$cluster_type)
fisher_adj_table$Cluster <- paste(fisher_adj_table$cluster_type, "_Cluster", fisher_adj_table$cluster_num, sep='')

all_cluster_enrichment <- fisher_adj_table[,c(1, 2, 4, 3)]
colnames(all_cluster_enrichment) <- c("Cluster", "GO ID", "GO name", "p value")

fisher_sig_all$cluster_num <- substring(fisher_sig_all$Cluster, 2, 2)
fisher_sig_all$cluster_type <- substring(fisher_sig_all$Cluster, 3, 3)
fisher_sig_all$cluster_type <- gsub("g", "Gene", fisher_sig_all$cluster_type)
fisher_sig_all$cluster_type <- gsub("p", "Protein", fisher_sig_all$cluster_type)
fisher_sig_all$Cluster <- paste(fisher_sig_all$cluster_type, "_Cluster", fisher_sig_all$cluster_num, sep='')

sig_cluster_enrichment <- fisher_sig_all[,c(1, 2, 4, 3)]
colnames(sig_cluster_enrichment) <- c("Cluster", "GO ID", "GO name", "p value")

#add in transcript names and BLAST in columns
all_cluster_transcripts_list <- list(gene_cluster1[[4]], gene_cluster2[[4]], gene_cluster3[[4]], gene_cluster4[[4]],
                                     gene_cluster5[[4]], gene_cluster6[[4]], protein_cluster1[[4]],
                                     protein_cluster2[[4]], protein_cluster3[[4]], protein_cluster4[[4]], protein_cluster5[[4]],
                                     protein_cluster6[[4]], protein_cluster7[[4]])

cluster_names <- c("Gene_Cluster1", "Gene_Cluster2", "Gene_Cluster3", "Gene_Cluster4", "Gene_Cluster5", "Gene_Cluster6",
                   "Protein_Cluster1", "Protein_Cluster2", "Protein_Cluster3", "Protein_Cluster4",
                   "Protein_Cluster5", "Protein_Cluster6", "Protein_Cluster7")

cluster_allEnrich_transcripts <- data.frame()
for (i in 1:length(cluster_names)) {
  cluster_allEnrich <- subset(all_cluster_enrichment, Cluster == cluster_names[i])
  cluster_transcripts <- all_cluster_transcripts_list[[i]]
  colnames(cluster_transcripts) <- c("GO ID", "transcripts", "BLAST")
  merged <- merge(cluster_allEnrich, cluster_transcripts, by = "GO ID")
  cluster_allEnrich_transcripts <- rbind(cluster_allEnrich_transcripts, merged)
}

#add in total number of genes for each cluster
cluster_totals <- data.frame(Cluster=c("Gene_Cluster1", "Gene_Cluster2", "Gene_Cluster3", "Gene_Cluster4", 
                             "Gene_Cluster5", "Gene_Cluster6", "Protein_Cluster1", "Protein_Cluster2",
                             "Protein_Cluster3", "Protein_Cluster4", "Protein_Cluster5", "Protein_Cluster6", "Protein_Cluster7"),
                             total = c(410, 321, 126, 91, 6, 5, 52, 11, 16, 8, 55, 4, 30))
all_cluster_enrich_totals <- merge(cluster_allEnrich_transcripts, cluster_totals, by = "Cluster")

colnames(all_cluster_enrich_totals) <- c("Cluster", "GO ID", "GO Name", "p value", "Transcripts List", "BLAST Descriptions List", "Total Number of Cluster Genes")

#export
write.table(all_cluster_enrich_totals, "~/cluster_enrichment_allPvals.txt", quote=F, sep="\t", row.names=F)


####transcripts that are in multiple clusters####
#there are (1249-1096=153) of them
unique_transcripts <- unique(all_names)
overlap_transcripts <- all_names[duplicated(all_names) == TRUE]

####transcripts that are in any of the pathways we already looked at####
#list of KEGG terms: 
KEGG_pathways <- c("path:map00010", "path:map03050", "path:map04120", "path:map04010", "path:map04210", "path:map00020")

#this gives all of the transcripts total in the selected pathways
all_pathways_list <- as.character(KEGG_annotation[which(KEGG_annotation$pathway_ID %in% KEGG_pathways),1])

gene_names <- vector()
for (i in 1:length(gene_expMat_cluster_list)) {
  gene_names <- c(gene_names, rownames(gene_expMat_cluster_list[[i]]))
}

protein_names <- vector()
for (i in 1:length(protein_expMat_cluster_list)) {
  protein_names <- c(protein_names, rownames(protein_expMat_cluster_list[[i]]))
}
#
all_names <- c(gene_names, protein_names)

#this is about 10% of the total
overlap <- all_names[which(all_names %in% all_pathways_list)]
length(unique(overlap))
length(unique(all_names))

####multivariate dispersion stats per cluster####
#code is from "mytilus variation bio levels.R"

library(ggplot2)
library(stringr)
library(plotly)
library(vegan)
library(compositions)
#library(pcaMethods)
library(multcompView)
library(colorspace)

#calculate multivariate dispersion per cluster
#use a matrix of expression values for each cluster - there are data frames in "gene_expMat_cluster_list" and "protein_expMat_cluster_list"
#the order of individuals is the same on everything, but we need a treatment name vector to use for all - same as "ind_treatments" above in cleaning section
treatment_vector <- c("FAE", "POL", "FAE", "FAP", "POH", "CGP", "POH", "FAE", "FAP", "FAE", "POL", "FAP", "FAP", "FAE",
                      "FAP", "POH", "CGP", "CGP", "POL", "FAE", "FAP", "POH", "POH", "POL", "CGP", "POL", "POL", "POH",
                      "POH", "POH", "POL", "POL", "FAE", "FAP", "POH", "CGP", "FAE", "FAP", "POH", "FAE", "FAP")


multivar_disp_function <- function(expr_dataframe) {
  mat <- t(expr_dataframe) #took out scale because it happens at the beginning
  dis <- dist(mat, method="euclidean")
  mod <- betadisper(dis, treatment_vector, type="median")
  #this gives MAD values
  #distances are from the centroid, so this is the distance from the median in multivariate space, you take the median of these absolute values to get MAD
  #10.17.19 use this as a reference from now on, much less complicated than older code
  distances <- mod$distances
  distances <- data.frame(distance=distances, Treatment=treatment_vector)
  CGP_dist <- subset(distances, Treatment == "CGP")
  FAE_dist <- subset(distances, Treatment == "FAE")
  FAP_dist <- subset(distances, Treatment == "FAP")
  POH_dist <- subset(distances, Treatment == "POH")
  POL_dist <- subset(distances, Treatment == "POL")
  MAD <- rep(0,5)
  MAD[1] <- median(abs(CGP_dist$distance))
  MAD[2] <- median(abs(FAE_dist$distance))
  MAD[3] <- median(abs(FAP_dist$distance))
  MAD[4] <- median(abs(POH_dist$distance))
  MAD[5] <- median(abs(POL_dist$distance))
  #this gives effect sizes
  effect_test1 <- cohen.d(CGP_dist$distance, FAE_dist$distance, hedges.correction = TRUE)
  effect_size1 <- effect_test1$estimate
  effect_test2 <- cohen.d(CGP_dist$distance, FAP_dist$distance, hedges.correction = TRUE)
  effect_size2 <- effect_test2$estimate
  effect_test3 <- cohen.d(CGP_dist$distance, POH_dist$distance, hedges.correction = TRUE)
  effect_size3 <- effect_test3$estimate
  effect_test4 <- cohen.d(CGP_dist$distance, POL_dist$distance, hedges.correction = TRUE)
  effect_size4 <- effect_test4$estimate
  effect_test5 <- cohen.d(FAE_dist$distance, FAP_dist$distance, hedges.correction = TRUE)
  effect_size5 <- effect_test5$estimate
  effect_test6 <- cohen.d(FAE_dist$distance, POH_dist$distance, hedges.correction = TRUE)
  effect_size6 <- effect_test6$estimate
  effect_test7 <- cohen.d(FAE_dist$distance, POL_dist$distance, hedges.correction = TRUE)
  effect_size7 <- effect_test7$estimate
  effect_test8 <- cohen.d(FAP_dist$distance, POH_dist$distance, hedges.correction = TRUE)
  effect_size8 <- effect_test8$estimate
  effect_test9 <- cohen.d(FAP_dist$distance, POL_dist$distance, hedges.correction = TRUE)
  effect_size9 <- effect_test9$estimate
  effect_test10 <- cohen.d(POL_dist$distance, POH_dist$distance, hedges.correction = TRUE)
  effect_size10 <- effect_test10$estimate
  effect_sizes <- data.frame(effect = c(effect_size1, effect_size2, effect_size3, effect_size4, effect_size5,
                                        effect_size6, effect_size7, effect_size8, effect_size9, effect_size10),
                             treat1 = c("CGP", "CGP", "CGP", "CGP", "FAE", "FAE", "FAE", "FAP", "FAP", "POH"),
                             treat2 = c("FAE", "FAP", "POH", "POL", "FAP", "POH", "POL", "POH", "POL", "POL"))
  #this gives p values
  table <-permutest(mod, pairwise = TRUE, permutations = 9999)
  pairwise_pVal <- table$pairwise$observed
  #adjust for multiple comparisons
  pairwise_pValadj <- p.adjust(pairwise_pVal, "BH")
  #this gives letters in case we want to plot
  letters <- multcompLetters(pairwise_pValadj, Letters=LETTERS)$Letters
  
  #return a list of these different elements to use
  #first one is MAD value and letters for plotting, second is the full pairwise comparisons adj p values, third is the effect sizes
  MAD_letters <- data.frame(MAD=MAD, sig_letters=letters)
  export <- list(MAD_letters, pairwise_pValadj, effect_sizes, pairwise_pVal)
  return(export)
}

#make giant data frames of this info
geneClusterMVD_MAD <- data.frame()
geneClusterMVD_multcompP <- vector()
geneClusterMVD_effSize <- data.frame()
geneClusterMVD_Pval <- vector()
for (i in 1:length(gene_expMat_cluster_list)) {
  MVD <- multivar_disp_function(gene_expMat_cluster_list[[i]])
  geneClusterMVD_MAD <- rbind(geneClusterMVD_MAD, MVD[[1]])
  geneClusterMVD_multcompP <- c(geneClusterMVD_multcompP, MVD[[2]])
  geneClusterMVD_effSize <- rbind(geneClusterMVD_effSize, MVD[[3]])
  geneClusterMVD_Pval <- c(geneClusterMVD_Pval, MVD[[4]])
}

#put in the cluster and data type labels
geneClusterMVD_MAD$cluster <- rep(rep(c(1:6), each=5), 2)
geneClusterMVD_MAD$type <- rep(c("gene", "protein"), each=30)
geneClusterMVD_MAD$treatment <- rep(c("CGP", "FAE", "FAP", "POH", "POL"), 12)
geneClusterMVD_pVal <- data.frame(comparison=names(geneClusterMVD_multcompP), effect_size=geneClusterMVD_effSize$effect, p_val=unname(geneClusterMVD_multcompP), cluster=rep(rep(c(1:6), each=10), 2),
                                  type=rep(c("gene", "protein"), each=60))
geneClusterMVD_pVal_noAdj <- data.frame(comparison=names(geneClusterMVD_Pval), p_val_adj=unname(geneClusterMVD_multcompP), p_val=unname(geneClusterMVD_Pval), cluster=rep(rep(c(1:6), each=10), 2),
                                  type=rep(c("gene", "protein"), each=60))
#reduce to only significant (before adjustment)
geneClusterMVD_pVal_noAdj <- geneClusterMVD_pVal_noAdj[(geneClusterMVD_pVal_noAdj$p_val<0.05),]

#export
write.table(geneClusterMVD_MAD, "~/geneClusterMVD_MAD_letters.txt", sep='\t', row.names=F, quote=F)
write.table(geneClusterMVD_pVal, "~/geneClusterMVD_multipleComp_pVal.txt", sep='\t', row.names=F, quote=F)
write.table(geneClusterMVD_pVal_noAdj, "~/geneCluster_MVD_multipleComp_pValOriginal.txt", sep='\t', row.names=F, quote=F)

#make giant data frames of this info
proteinClusterMVD_MAD <- data.frame()
proteinClusterMVD_multcompP <- vector()
proteinClusterMVD_effSize <- data.frame()
proteinClusterMVD_Pval <- vector()
for (i in 1:length(protein_expMat_cluster_list)) {
  MVD <- multivar_disp_function(protein_expMat_cluster_list[[i]])
  proteinClusterMVD_MAD <- rbind(proteinClusterMVD_MAD, MVD[[1]])
  proteinClusterMVD_multcompP <- c(proteinClusterMVD_multcompP, MVD[[2]])
  proteinClusterMVD_effSize <- rbind(proteinClusterMVD_effSize, MVD[[3]])
  proteinClusterMVD_Pval <- c(proteinClusterMVD_Pval, MVD[[4]])
}

#put in the cluster and data type labels
proteinClusterMVD_MAD$cluster <- rep(rep(c(1:7), each=5), 2)
proteinClusterMVD_MAD$type <- rep(c("gene", "protein"), each=35)
proteinClusterMVD_MAD$treatment <- rep(c("CGP", "FAE", "FAP", "POH", "POL"), 14)
proteinClusterMVD_pVal <- data.frame(comparison=names(proteinClusterMVD_multcompP), effect_size=proteinClusterMVD_effSize$effect, p_val=unname(proteinClusterMVD_multcompP), cluster=rep(rep(c(1:7), each=10), 2),
                                  type=rep(c("gene", "protein"), each=70))
proteinClusterMVD_pVal_noAdj <- data.frame(comparison=names(proteinClusterMVD_Pval), p_val_adj=unname(proteinClusterMVD_multcompP), p_val=unname(proteinClusterMVD_Pval), cluster=rep(rep(c(1:7), each=10), 2),
                                        type=rep(c("gene", "protein"), each=70))
#reduce to only significant (before adjustment)
proteinClusterMVD_pVal_noAdj <- proteinClusterMVD_pVal_noAdj[(proteinClusterMVD_pVal_noAdj$p_val<0.05),]

#export 
write.table(proteinClusterMVD_MAD, "~/proteinClusterMVD_MAD_letters.txt", sep='\t', row.names=F, quote=F)
write.table(proteinClusterMVD_pVal, "~/proteinClusterMVD_multipleComp_pVal.txt", sep='\t', row.names=F, quote=F)
write.table(proteinClusterMVD_pVal_noAdj, "~/proteinCluster_MVD_multipleComp_pValOriginal.txt", sep='\t', row.names=F, quote=F)


####plots for MAD within pathway between gene and protein####

#look at average over genes within a pathway
#across all treatments
MAD_avgAllGene_gC1 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[1]])))
MAD_avgAllProtein_gC1 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[7]])))
MAD_avgAllGene_gC2 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[2]])))
MAD_avgAllProtein_gC2 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[8]])))
MAD_avgAllGene_gC3 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[3]])))
MAD_avgAllProtein_gC3 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[9]])))
MAD_avgAllGene_gC4 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[4]])))
MAD_avgAllProtein_gC4 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[10]])))
MAD_avgAllGene_gC5 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[5]])))
MAD_avgAllProtein_gC5 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[11]])))
MAD_avgAllGene_gC6 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[6]])))
MAD_avgAllProtein_gC6 <- mean(rowMeans(na.omit(gene_cluster_MAD_list[[12]])))

MAD_avgAllGene_pC1 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[1]])))
MAD_avgAllProtein_pC1 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[8]])))
MAD_avgAllGene_pC2 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[2]])))
MAD_avgAllProtein_pC2 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[9]])))
MAD_avgAllGene_pC3 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[3]])))
MAD_avgAllProtein_pC3 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[10]])))
MAD_avgAllGene_pC4 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[4]])))
MAD_avgAllProtein_pC4 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[11]])))
MAD_avgAllGene_pC5 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[5]])))
MAD_avgAllProtein_pC5 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[12]])))
MAD_avgAllGene_pC6 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[6]])))
MAD_avgAllProtein_pC6 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[13]])))
MAD_avgAllGene_pC7 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[7]])))
MAD_avgAllProtein_pC7 <- mean(rowMeans(na.omit(protein_cluster_MAD_list[[14]])))

#within one treatment (POH?)
MAD_avgPOHgene_gC1 <- mean(na.omit(gene_cluster_MAD_list[[1]][,4]))
MAD_avgPOHprotein_gC1 <- mean(na.omit(gene_cluster_MAD_list[[7]][,4]))
MAD_avgPOHgene_gC2 <- mean(na.omit(gene_cluster_MAD_list[[2]][,4]))
MAD_avgPOHprotein_gC2 <- mean(na.omit(gene_cluster_MAD_list[[8]][,4]))
MAD_avgPOHgene_gC3 <- mean(na.omit(gene_cluster_MAD_list[[3]][,4]))
MAD_avgPOHprotein_gC3 <- mean(na.omit(gene_cluster_MAD_list[[9]][,4]))
MAD_avgPOHgene_gC4 <- mean(na.omit(gene_cluster_MAD_list[[4]][,4]))
MAD_avgPOHprotein_gC4 <- mean(na.omit(gene_cluster_MAD_list[[10]][,4]))
MAD_avgPOHgene_gC5 <- mean(na.omit(gene_cluster_MAD_list[[5]][,4]))
MAD_avgPOHprotein_gC5 <- mean(na.omit(gene_cluster_MAD_list[[11]][,4]))
MAD_avgPOHgene_gC6 <- mean(na.omit(gene_cluster_MAD_list[[6]][,4]))
MAD_avgPOHprotein_gC6 <- mean(na.omit(gene_cluster_MAD_list[[12]][,4]))


MAD_avgPOH_gene <- c(MAD_avgPOHgene_gC1, MAD_avgPOHgene_gC2, MAD_avgPOHgene_gC3, MAD_avgPOHgene_gC4, MAD_avgPOHgene_gC5, MAD_avgPOHgene_gC6,
                MAD_avgPOHprotein_gC1, MAD_avgPOHprotein_gC2, MAD_avgPOHprotein_gC3, MAD_avgPOHprotein_gC4, MAD_avgPOHprotein_gC5, MAD_avgPOHprotein_gC6)

MAD_avgPOHgene_pC1 <- mean(na.omit(protein_cluster_MAD_list[[1]][,4]))
MAD_avgPOHprotein_pC1 <- mean(na.omit(protein_cluster_MAD_list[[8]][,4]))
MAD_avgPOHgene_pC2 <- mean(na.omit(protein_cluster_MAD_list[[2]][,4]))
MAD_avgPOHprotein_pC2 <- mean(na.omit(protein_cluster_MAD_list[[9]][,4]))
MAD_avgPOHgene_pC3 <- mean(na.omit(protein_cluster_MAD_list[[3]][,4]))
MAD_avgPOHprotein_pC3 <- mean(na.omit(protein_cluster_MAD_list[[10]][,4]))
MAD_avgPOHgene_pC4 <- mean(na.omit(protein_cluster_MAD_list[[4]][,4]))
MAD_avgPOHprotein_pC4 <- mean(na.omit(protein_cluster_MAD_list[[11]][,4]))
MAD_avgPOHgene_pC5 <- mean(na.omit(protein_cluster_MAD_list[[5]][,4]))
MAD_avgPOHprotein_pC5 <- mean(na.omit(protein_cluster_MAD_list[[12]][,4]))
MAD_avgPOHgene_pC6 <- mean(na.omit(protein_cluster_MAD_list[[6]][,4]))
MAD_avgPOHprotein_pC6 <- mean(na.omit(protein_cluster_MAD_list[[13]][,4]))
MAD_avgPOHgene_pC7 <- mean(na.omit(protein_cluster_MAD_list[[7]][,4]))
MAD_avgPOHprotein_pC7 <- mean(na.omit(protein_cluster_MAD_list[[14]][,4]))

MAD_avgPOH_protein <- c(MAD_avgPOHgene_pC1, MAD_avgPOHgene_pC2, MAD_avgPOHgene_pC3, MAD_avgPOHgene_pC4, MAD_avgPOHgene_pC5, MAD_avgPOHgene_pC6, MAD_avgPOHgene_pC7,
                MAD_avgPOHprotein_pC1, MAD_avgPOHprotein_pC2, MAD_avgPOHprotein_pC3, MAD_avgPOHprotein_pC4, MAD_avgPOHprotein_pC5, MAD_avgPOHprotein_pC6, MAD_avgPOHprotein_pC7)

#within POL
MAD_avgPOLgene_gC1 <- mean(na.omit(gene_cluster_MAD_list[[1]][,5]))
MAD_avgPOLprotein_gC1 <- mean(na.omit(gene_cluster_MAD_list[[7]][,5]))
MAD_avgPOLgene_gC2 <- mean(na.omit(gene_cluster_MAD_list[[2]][,5]))
MAD_avgPOLprotein_gC2 <- mean(na.omit(gene_cluster_MAD_list[[8]][,5]))
MAD_avgPOLgene_gC3 <- mean(na.omit(gene_cluster_MAD_list[[3]][,5]))
MAD_avgPOLprotein_gC3 <- mean(na.omit(gene_cluster_MAD_list[[9]][,5]))
MAD_avgPOLgene_gC4 <- mean(na.omit(gene_cluster_MAD_list[[4]][,5]))
MAD_avgPOLprotein_gC4 <- mean(na.omit(gene_cluster_MAD_list[[10]][,5]))
MAD_avgPOLgene_gC5 <- mean(na.omit(gene_cluster_MAD_list[[5]][,5]))
MAD_avgPOLprotein_gC5 <- mean(na.omit(gene_cluster_MAD_list[[11]][,5]))
MAD_avgPOLgene_gC6 <- mean(na.omit(gene_cluster_MAD_list[[6]][,5]))
MAD_avgPOLprotein_gC6 <- mean(na.omit(gene_cluster_MAD_list[[12]][,5]))

MAD_avgPOL_gene <- c(MAD_avgPOLgene_gC1, MAD_avgPOLgene_gC2, MAD_avgPOLgene_gC3, MAD_avgPOLgene_gC4, MAD_avgPOLgene_gC5, MAD_avgPOLgene_gC6,
                     MAD_avgPOLprotein_gC1, MAD_avgPOLprotein_gC2, MAD_avgPOLprotein_gC3, MAD_avgPOLprotein_gC4, MAD_avgPOLprotein_gC5, MAD_avgPOLprotein_gC6)

MAD_avgPOLgene_pC1 <- mean(na.omit(protein_cluster_MAD_list[[1]][,5]))
MAD_avgPOLprotein_pC1 <- mean(na.omit(protein_cluster_MAD_list[[8]][,5]))
MAD_avgPOLgene_pC2 <- mean(na.omit(protein_cluster_MAD_list[[2]][,5]))
MAD_avgPOLprotein_pC2 <- mean(na.omit(protein_cluster_MAD_list[[9]][,5]))
MAD_avgPOLgene_pC3 <- mean(na.omit(protein_cluster_MAD_list[[3]][,5]))
MAD_avgPOLprotein_pC3 <- mean(na.omit(protein_cluster_MAD_list[[10]][,5]))
MAD_avgPOLgene_pC4 <- mean(na.omit(protein_cluster_MAD_list[[4]][,5]))
MAD_avgPOLprotein_pC4 <- mean(na.omit(protein_cluster_MAD_list[[11]][,5]))
MAD_avgPOLgene_pC5 <- mean(na.omit(protein_cluster_MAD_list[[5]][,5]))
MAD_avgPOLprotein_pC5 <- mean(na.omit(protein_cluster_MAD_list[[12]][,5]))
MAD_avgPOLgene_pC6 <- mean(na.omit(protein_cluster_MAD_list[[6]][,5]))
MAD_avgPOLprotein_pC6 <- mean(na.omit(protein_cluster_MAD_list[[13]][,5]))
MAD_avgPOLgene_pC7 <- mean(na.omit(protein_cluster_MAD_list[[7]][,5]))
MAD_avgPOLprotein_pC7 <- mean(na.omit(protein_cluster_MAD_list[[14]][,5]))

MAD_avgPOL_protein <- c(MAD_avgPOLgene_pC1, MAD_avgPOLgene_pC2, MAD_avgPOLgene_pC3, MAD_avgPOLgene_pC4, MAD_avgPOLgene_pC5, MAD_avgPOLgene_pC6, MAD_avgPOLgene_pC7,
                        MAD_avgPOLprotein_pC1, MAD_avgPOLprotein_pC2, MAD_avgPOLprotein_pC3, MAD_avgPOLprotein_pC4, MAD_avgPOLprotein_pC5, MAD_avgPOLprotein_pC6, MAD_avgPOLprotein_pC7)


#within CGP
MAD_avgCGPgene_gC1 <- mean(na.omit(gene_cluster_MAD_list[[1]][,1]))
MAD_avgCGPprotein_gC1 <- mean(na.omit(gene_cluster_MAD_list[[7]][,1]))
MAD_avgCGPgene_gC2 <- mean(na.omit(gene_cluster_MAD_list[[2]][,1]))
MAD_avgCGPprotein_gC2 <- mean(na.omit(gene_cluster_MAD_list[[8]][,1]))
MAD_avgCGPgene_gC3 <- mean(na.omit(gene_cluster_MAD_list[[3]][,1]))
MAD_avgCGPprotein_gC3 <- mean(na.omit(gene_cluster_MAD_list[[9]][,1]))
MAD_avgCGPgene_gC4 <- mean(na.omit(gene_cluster_MAD_list[[4]][,1]))
MAD_avgCGPprotein_gC4 <- mean(na.omit(gene_cluster_MAD_list[[10]][,1]))
MAD_avgCGPgene_gC5 <- mean(na.omit(gene_cluster_MAD_list[[5]][,1]))
MAD_avgCGPprotein_gC5 <- mean(na.omit(gene_cluster_MAD_list[[11]][,1]))
MAD_avgCGPgene_gC6 <- mean(na.omit(gene_cluster_MAD_list[[6]][,1]))
MAD_avgCGPprotein_gC6 <- mean(na.omit(gene_cluster_MAD_list[[12]][,1]))

MAD_avgCGP_gene <- c(MAD_avgCGPgene_gC1, MAD_avgCGPgene_gC2, MAD_avgCGPgene_gC3, MAD_avgCGPgene_gC4, MAD_avgCGPgene_gC5, MAD_avgCGPgene_gC6, 
                     MAD_avgCGPprotein_gC1, MAD_avgCGPprotein_gC2, MAD_avgCGPprotein_gC3, MAD_avgCGPprotein_gC4, MAD_avgCGPprotein_gC5, MAD_avgCGPprotein_gC6)

MAD_avgCGPgene_pC1 <- mean(na.omit(protein_cluster_MAD_list[[1]][,1]))
MAD_avgCGPprotein_pC1 <- mean(na.omit(protein_cluster_MAD_list[[8]][,1]))
MAD_avgCGPgene_pC2 <- mean(na.omit(protein_cluster_MAD_list[[2]][,1]))
MAD_avgCGPprotein_pC2 <- mean(na.omit(protein_cluster_MAD_list[[9]][,1]))
MAD_avgCGPgene_pC3 <- mean(na.omit(protein_cluster_MAD_list[[3]][,1]))
MAD_avgCGPprotein_pC3 <- mean(na.omit(protein_cluster_MAD_list[[10]][,1]))
MAD_avgCGPgene_pC4 <- mean(na.omit(protein_cluster_MAD_list[[4]][,1]))
MAD_avgCGPprotein_pC4 <- mean(na.omit(protein_cluster_MAD_list[[11]][,1]))
MAD_avgCGPgene_pC5 <- mean(na.omit(protein_cluster_MAD_list[[5]][,1]))
MAD_avgCGPprotein_pC5 <- mean(na.omit(protein_cluster_MAD_list[[12]][,1]))
MAD_avgCGPgene_pC6 <- mean(na.omit(protein_cluster_MAD_list[[6]][,1]))
MAD_avgCGPprotein_pC6 <- mean(na.omit(protein_cluster_MAD_list[[13]][,1]))
MAD_avgCGPgene_pC7 <- mean(na.omit(protein_cluster_MAD_list[[7]][,1]))
MAD_avgCGPprotein_pC7 <- mean(na.omit(protein_cluster_MAD_list[[14]][,1]))

MAD_avgCGP_protein <- c(MAD_avgCGPgene_pC1, MAD_avgCGPgene_pC2, MAD_avgCGPgene_pC3, MAD_avgCGPgene_pC4, MAD_avgCGPgene_pC5, MAD_avgCGPgene_pC6, MAD_avgCGPgene_pC7, 
                        MAD_avgCGPprotein_pC1, MAD_avgCGPprotein_pC2, MAD_avgCGPprotein_pC3, MAD_avgCGPprotein_pC4, MAD_avgCGPprotein_pC5, MAD_avgCGPprotein_pC6, MAD_avgCGPprotein_pC7)

#within FAE
MAD_avgFAEgene_gC1 <- mean(na.omit(gene_cluster_MAD_list[[1]][,2]))
MAD_avgFAEprotein_gC1 <- mean(na.omit(gene_cluster_MAD_list[[7]][,2]))
MAD_avgFAEgene_gC2 <- mean(na.omit(gene_cluster_MAD_list[[2]][,2]))
MAD_avgFAEprotein_gC2 <- mean(na.omit(gene_cluster_MAD_list[[8]][,2]))
MAD_avgFAEgene_gC3 <- mean(na.omit(gene_cluster_MAD_list[[3]][,2]))
MAD_avgFAEprotein_gC3 <- mean(na.omit(gene_cluster_MAD_list[[9]][,2]))
MAD_avgFAEgene_gC4 <- mean(na.omit(gene_cluster_MAD_list[[4]][,2]))
MAD_avgFAEprotein_gC4 <- mean(na.omit(gene_cluster_MAD_list[[10]][,2]))
MAD_avgFAEgene_gC5 <- mean(na.omit(gene_cluster_MAD_list[[5]][,2]))
MAD_avgFAEprotein_gC5 <- mean(na.omit(gene_cluster_MAD_list[[11]][,2]))
MAD_avgFAEgene_gC6 <- mean(na.omit(gene_cluster_MAD_list[[6]][,2]))
MAD_avgFAEprotein_gC6 <- mean(na.omit(gene_cluster_MAD_list[[12]][,2]))

MAD_avgFAE_gene <- c(MAD_avgFAEgene_gC1, MAD_avgFAEgene_gC2, MAD_avgFAEgene_gC3, MAD_avgFAEgene_gC4, MAD_avgFAEgene_gC5, MAD_avgFAEgene_gC6,
                     MAD_avgFAEprotein_gC1, MAD_avgFAEprotein_gC2, MAD_avgFAEprotein_gC3, MAD_avgFAEprotein_gC4, MAD_avgFAEprotein_gC5, MAD_avgFAEprotein_gC6)

MAD_avgFAEgene_pC1 <- mean(na.omit(protein_cluster_MAD_list[[1]][,2]))
MAD_avgFAEprotein_pC1 <- mean(na.omit(protein_cluster_MAD_list[[8]][,2]))
MAD_avgFAEgene_pC2 <- mean(na.omit(protein_cluster_MAD_list[[2]][,2]))
MAD_avgFAEprotein_pC2 <- mean(na.omit(protein_cluster_MAD_list[[9]][,2]))
MAD_avgFAEgene_pC3 <- mean(na.omit(protein_cluster_MAD_list[[3]][,2]))
MAD_avgFAEprotein_pC3 <- mean(na.omit(protein_cluster_MAD_list[[10]][,2]))
MAD_avgFAEgene_pC4 <- mean(na.omit(protein_cluster_MAD_list[[4]][,2]))
MAD_avgFAEprotein_pC4 <- mean(na.omit(protein_cluster_MAD_list[[11]][,2]))
MAD_avgFAEgene_pC5 <- mean(na.omit(protein_cluster_MAD_list[[5]][,2]))
MAD_avgFAEprotein_pC5 <- mean(na.omit(protein_cluster_MAD_list[[12]][,2]))
MAD_avgFAEgene_pC6 <- mean(na.omit(protein_cluster_MAD_list[[6]][,2]))
MAD_avgFAEprotein_pC6 <- mean(na.omit(protein_cluster_MAD_list[[13]][,2]))
MAD_avgFAEgene_pC7 <- mean(na.omit(protein_cluster_MAD_list[[7]][,2]))
MAD_avgFAEprotein_pC7 <- mean(na.omit(protein_cluster_MAD_list[[14]][,2]))

MAD_avgFAE_protein <- c(MAD_avgFAEgene_pC1, MAD_avgFAEgene_pC2, MAD_avgFAEgene_pC3, MAD_avgFAEgene_pC4, MAD_avgFAEgene_pC5, MAD_avgFAEgene_pC6, MAD_avgFAEgene_pC7,
                        MAD_avgFAEprotein_pC1, MAD_avgFAEprotein_pC2, MAD_avgFAEprotein_pC3, MAD_avgFAEprotein_pC4, MAD_avgFAEprotein_pC5, MAD_avgFAEprotein_pC6, MAD_avgFAEprotein_pC7)


#within FAP
MAD_avgFAPgene_gC1 <- mean(na.omit(gene_cluster_MAD_list[[1]][,3]))
MAD_avgFAPprotein_gC1 <- mean(na.omit(gene_cluster_MAD_list[[7]][,3]))
MAD_avgFAPgene_gC2 <- mean(na.omit(gene_cluster_MAD_list[[2]][,3]))
MAD_avgFAPprotein_gC2 <- mean(na.omit(gene_cluster_MAD_list[[8]][,3]))
MAD_avgFAPgene_gC3 <- mean(na.omit(gene_cluster_MAD_list[[3]][,3]))
MAD_avgFAPprotein_gC3 <- mean(na.omit(gene_cluster_MAD_list[[9]][,3]))
MAD_avgFAPgene_gC4 <- mean(na.omit(gene_cluster_MAD_list[[4]][,3]))
MAD_avgFAPprotein_gC4 <- mean(na.omit(gene_cluster_MAD_list[[10]][,3]))
MAD_avgFAPgene_gC5 <- mean(na.omit(gene_cluster_MAD_list[[5]][,3]))
MAD_avgFAPprotein_gC5 <- mean(na.omit(gene_cluster_MAD_list[[11]][,3]))
MAD_avgFAPgene_gC6 <- mean(na.omit(gene_cluster_MAD_list[[6]][,3]))
MAD_avgFAPprotein_gC6 <- mean(na.omit(gene_cluster_MAD_list[[12]][,3]))

MAD_avgFAP_gene <- c(MAD_avgFAPgene_gC1, MAD_avgFAPgene_gC2, MAD_avgFAPgene_gC3, MAD_avgFAPgene_gC4, MAD_avgFAPgene_gC5, MAD_avgFAPgene_gC6,
                     MAD_avgFAPprotein_gC1, MAD_avgFAPprotein_gC2, MAD_avgFAPprotein_gC3, MAD_avgFAPprotein_gC4, MAD_avgFAPprotein_gC5, MAD_avgFAPprotein_gC6)

MAD_avgFAPgene_pC1 <- mean(na.omit(protein_cluster_MAD_list[[1]][,3]))
MAD_avgFAPprotein_pC1 <- mean(na.omit(protein_cluster_MAD_list[[8]][,3]))
MAD_avgFAPgene_pC2 <- mean(na.omit(protein_cluster_MAD_list[[2]][,3]))
MAD_avgFAPprotein_pC2 <- mean(na.omit(protein_cluster_MAD_list[[9]][,3]))
MAD_avgFAPgene_pC3 <- mean(na.omit(protein_cluster_MAD_list[[3]][,3]))
MAD_avgFAPprotein_pC3 <- mean(na.omit(protein_cluster_MAD_list[[10]][,3]))
MAD_avgFAPgene_pC4 <- mean(na.omit(protein_cluster_MAD_list[[4]][,3]))
MAD_avgFAPprotein_pC4 <- mean(na.omit(protein_cluster_MAD_list[[11]][,3]))
MAD_avgFAPgene_pC5 <- mean(na.omit(protein_cluster_MAD_list[[5]][,3]))
MAD_avgFAPprotein_pC5 <- mean(na.omit(protein_cluster_MAD_list[[12]][,3]))
MAD_avgFAPgene_pC6 <- mean(na.omit(protein_cluster_MAD_list[[6]][,3]))
MAD_avgFAPprotein_pC6 <- mean(na.omit(protein_cluster_MAD_list[[13]][,3]))
MAD_avgFAPgene_pC7 <- mean(na.omit(protein_cluster_MAD_list[[7]][,3]))
MAD_avgFAPprotein_pC7 <- mean(na.omit(protein_cluster_MAD_list[[14]][,3]))

MAD_avgFAP_protein <- c(MAD_avgFAPgene_pC1, MAD_avgFAPgene_pC2, MAD_avgFAPgene_pC3, MAD_avgFAPgene_pC4, MAD_avgFAPgene_pC5, MAD_avgFAPgene_pC6, MAD_avgFAPgene_pC7,
                        MAD_avgFAPprotein_pC1, MAD_avgFAPprotein_pC2, MAD_avgFAPprotein_pC3, MAD_avgFAPprotein_pC4, MAD_avgFAPprotein_pC5, MAD_avgFAPprotein_pC6, MAD_avgFAPprotein_pC7)

#make a giant data frame for plotting
#rows = pathways x data type (12 total rows)
#columns = MVD, allMAD, POH_MAD, POL_MAD, FAE_MAD, FAP_MAD, CGP_MAD
plot_MAD_geneC <- data.frame(cluster=rep(c("gene C1", "gene C2", "gene C3", "gene C4", "gene C5", "gene C6"), 2),
                       type=rep(c("gene", "protein"), each=6),
                       MAD_all=c(MAD_avgAllGene_gC1, MAD_avgAllGene_gC2, MAD_avgAllGene_gC3, MAD_avgAllGene_gC4, MAD_avgAllGene_gC5, MAD_avgAllGene_gC6,
                                 MAD_avgAllProtein_gC1, MAD_avgAllProtein_gC2, MAD_avgAllProtein_gC3, MAD_avgAllProtein_gC4, MAD_avgAllProtein_gC5, MAD_avgAllProtein_gC6),
                       MAD_CGP=MAD_avgCGP_gene,
                       MAD_FAE=MAD_avgFAE_gene,
                       MAD_FAP=MAD_avgFAP_gene,
                       MAD_POL=MAD_avgPOL_gene,
                       MAD_POH=MAD_avgPOH_gene)

#eventual goal: 6 panel figure showing each of the treatments separately and one all together with all multivariate mapped on top
#bar plot: categorical pathway vs scaled MAD

#plot all together

ggplot(plot_MAD_geneC, aes(x=cluster, y=MAD_all, group=type)) +
  geom_col(aes(x=cluster, y=MAD_all, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#melt data frame so we can plot treatment x MAD by cluster 
library(reshape2)
plot_MAD_geneC_melt <- melt(plot_MAD_geneC)

gC1_meltMAD <- subset(plot_MAD_geneC_melt, cluster == "gene C1")
ggplot(gC1_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POH" = "POH", "MAD_POL" = "POL")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

gC2_meltMAD <- subset(plot_MAD_geneC_melt, cluster == "gene C2")
ggplot(gC2_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POH" = "POH", "MAD_POL" = "POL")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

gC3_meltMAD <- subset(plot_MAD_geneC_melt, cluster == "gene C3")
ggplot(gC3_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POH" = "POH", "MAD_POL" = "POL")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

gC4_meltMAD <- subset(plot_MAD_geneC_melt, cluster == "gene C4")
ggplot(gC4_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POH" = "POH", "MAD_POL" = "POL")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

gC5_meltMAD <- subset(plot_MAD_geneC_melt, cluster == "gene C5")
ggplot(gC5_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POH" = "POH", "MAD_POL" = "POL")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

gC6_meltMAD <- subset(plot_MAD_geneC_melt, cluster == "gene C6")
ggplot(gC6_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POH" = "POH", "MAD_POL" = "POL")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")


#POH
ggplot(plot_MAD_geneC, aes(x=cluster, y=MAD_POH, group=type)) +
  geom_col(aes(x=cluster, y=MAD_POH, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in POH") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot POL
ggplot(plot_MAD_geneC, aes(x=cluster, y=MAD_POL, group=type)) +
  geom_col(aes(x=cluster, y=MAD_POL, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in POL") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot CGP
ggplot(plot_MAD_geneC, aes(x=cluster, y=MAD_CGP, group=type)) +
  geom_col(aes(x=cluster, y=MAD_CGP, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in CGP") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot FAE
ggplot(plot_MAD_geneC, aes(x=cluster, y=MAD_FAE, group=type)) +
  geom_col(aes(x=cluster, y=MAD_FAE, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in FAE") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot FAP
ggplot(plot_MAD_geneC, aes(x=cluster, y=MAD_FAP, group=type)) +
  geom_col(aes(x=cluster, y=MAD_FAP, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in FAP") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

plot_MAD_proteinC <- data.frame(cluster=rep(c("protein C1", "protein C2", "protein C3", "protein C4", "protein C5", "protein C6", "protein C7"), 2),
                             type=rep(c("gene", "protein"), each=7),
                             MAD_all=c(MAD_avgAllGene_pC1, MAD_avgAllGene_pC2, MAD_avgAllGene_pC3, MAD_avgAllGene_pC4, MAD_avgAllGene_pC5, MAD_avgAllGene_pC6, MAD_avgAllGene_pC7, 
                                       MAD_avgAllProtein_pC1, MAD_avgAllProtein_pC2, MAD_avgAllProtein_pC3, MAD_avgAllProtein_pC4, MAD_avgAllProtein_pC5, MAD_avgAllProtein_pC6, MAD_avgAllProtein_pC7),
                             MAD_CGP=MAD_avgCGP_protein,
                             MAD_FAE=MAD_avgFAE_protein,
                             MAD_FAP=MAD_avgFAP_protein,
                             MAD_POL=MAD_avgPOL_protein,
                             MAD_POH=MAD_avgPOH_protein)

#eventual goal: 6 panel figure showing each of the treatments separately and one all together with all multivariate mapped on top
#bar plot: categorical pathway vs scaled MAD

#plot all together

ggplot(plot_MAD_proteinC, aes(x=cluster, y=MAD_all, group=type)) +
  geom_col(aes(x=cluster, y=MAD_all, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#melt data frame so we can plot treatment x MAD by cluster
library(reshape2)
plot_MAD_proteinC_melt <- melt(plot_MAD_proteinC)

pC1_meltMAD <- subset(plot_MAD_proteinC_melt, cluster == "protein C1")
ggplot(pC1_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POL" = "POL", "MAD_POH" = "POH")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

pC2_meltMAD <- subset(plot_MAD_proteinC_melt, cluster == "protein C2")
ggplot(pC2_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POL" = "POL", "MAD_POH" = "POH")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

pC3_meltMAD <- subset(plot_MAD_proteinC_melt, cluster == "protein C3")
ggplot(pC3_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POL" = "POL", "MAD_POH" = "POH")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

pC4_meltMAD <- subset(plot_MAD_proteinC_melt, cluster == "protein C4")
ggplot(pC4_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POL" = "POL", "MAD_POH" = "POH")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

pC5_meltMAD <- subset(plot_MAD_proteinC_melt, cluster == "protein C5")
ggplot(pC5_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POL" = "POL", "MAD_POH" = "POH")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

pC6_meltMAD <- subset(plot_MAD_proteinC_melt, cluster == "protein C6")
ggplot(pC6_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POL" = "POL", "MAD_POH" = "POH")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

pC7_meltMAD <- subset(plot_MAD_proteinC_melt, cluster == "protein C7")
ggplot(pC7_meltMAD, aes(x=variable, y=value, group=type)) +
  geom_col(aes(x=variable, y=value, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  scale_x_discrete(labels=c("MAD_all" = "All", "MAD_CGP" = "PCG", "MAD_FAE" = "EFA", "MAD_FAP" = "PFA", 
                            "MAD_POL" = "POL", "MAD_POH" = "POH")) +
  ylim(0, 0.85) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Average scaled MAD within Treatment") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.position = "none")

#plot POH
ggplot(plot_MAD_proteinC, aes(x=cluster, y=MAD_POH, group=type)) +
  geom_col(aes(x=cluster, y=MAD_POH, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in POH") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot POL
ggplot(plot_MAD_proteinC, aes(x=cluster, y=MAD_POL, group=type)) +
  geom_col(aes(x=cluster, y=MAD_POL, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in POL") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot CGP
ggplot(plot_MAD_proteinC, aes(x=cluster, y=MAD_CGP, group=type)) +
  geom_col(aes(x=cluster, y=MAD_CGP, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in CGP") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot FAE
ggplot(plot_MAD_proteinC, aes(x=cluster, y=MAD_FAE, group=type)) +
  geom_col(aes(x=cluster, y=MAD_FAE, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in FAE") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot FAP
ggplot(plot_MAD_proteinC, aes(x=cluster, y=MAD_FAP, group=type)) +
  geom_col(aes(x=cluster, y=MAD_FAP, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in FAP") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

####gene by gene analysis: GP, treatment, genes####
#same as done in the KEGG pathways
#needs data frames that are separated by cluster and by treatment - these are stored in "gene_expMAT_cluster_treat_list" and "protein_expMAT_cluster_treat_list"


residuals_function <- function(data_input) {
  #USE LOOPS
  residuals <- data.frame(matrix(NA, nrow=length(data_input[,1]), ncol=length(data_input[1,])))
  columns <- colnames(data_input)
  rows <- rownames(data_input)
  #calculate MAD per gene/protein
  for (i in 1:length(data_input[1,])) {
    sub <- na.omit(as.numeric(as.character(data_input[,i])))
    #find median for each gene
    median <- median(sub)
    #calculate deviations
    deviations <- abs(sub - median)
    residuals[,i] <- deviations
  }
  colnames(residuals) <- columns
  rownames(residuals) <- rows
  return(residuals)
}

t.test_treatmentGP_function <- function(gene_treat_res, protein_treat_res) { #input residual matrices for each treatment
  treat_comp_GP <- vector()
  for (i in 1:(length(gene_treat_res[1,]))) {
    gene_sub <- gene_treat_res[,i]
    protein_sub <- protein_treat_res[,i]
    if ((sum(!is.na(gene_sub)) > 1) && (sum(!is.na(protein_sub)) > 1)) {
      p_val <- t.test(gene_sub, protein_sub)$p.value
      stat <- t.test(gene_sub, protein_sub)$statistic
      effect_test <- cohen.d(gene_sub,protein_sub, hedges.correction = TRUE)
      effect_size <- effect_test$estimate
      add <- c(p_val, stat, effect_size, colnames(gene_treat_res)[i])
      #}
      treat_comp_GP <- rbind(treat_comp_GP, add)
      }
    }
  treat_comp_GP <- data.frame(treat_comp_GP)
  colnames(treat_comp_GP) <- c("p_val", "t","effect_size", "transcript")
  
  #multiple comparisons adjustment
  treat_comp_GP$p_val <- as.numeric(as.character(treat_comp_GP$p_val))
  treat_comp_GP$t <- as.numeric(as.character(treat_comp_GP$t))
  treat_comp_GP$effect_size <- as.numeric(as.character(treat_comp_GP$effect_size))
  treat_comp_GP$transcript <- as.character(treat_comp_GP$transcript)
  treat_comp_adj <- p.adjust(treat_comp_GP[,1], "BH")
  treat_comp_adj[is.na(treat_comp_adj)] <- 1
  
  #cut down to only significant ones
  treat_sigcomp_GP <- vector()
  for (i in 1:length(treat_comp_adj)) {
    if (treat_comp_adj[i] < 0.05) {
      sig <- c(treat_comp_adj[i], treat_comp_GP$t[i], treat_comp_GP$effect_size[i], treat_comp_GP$transcript[i])
      treat_sigcomp_GP <- rbind(treat_sigcomp_GP, sig)
    }
  }
  treat_sigcomp_GP <- data.frame(treat_sigcomp_GP)
  
  list_comp <- list(treat_comp_GP, treat_comp_adj, treat_sigcomp_GP)
  return(list_comp)
}

#TREAT FOR EFFECT SIZE REPORTING IS FLIPPED #fixed 5.28.20
anova_treatment_function <- function(treat_res_CGP, treat_res_FAE, treat_res_FAP, treat_res_POH, treat_res_POL) { #input residual matrices for each treatment
  treat_comp <- vector()
  for (i in 1:(length(treat_res_CGP[1,]))) {
    #contrasts(treat_res$treat) <- contr.sum #compares all to the grand mean, which is not what we want
    CGP_sub <- treat_res_CGP[,i]
    FAE_sub <- treat_res_FAE[,i]
    FAP_sub <- treat_res_FAP[,i]
    POH_sub <- treat_res_POH[,i]
    POL_sub <- treat_res_POL[,i]
    
    if ((sum(!is.na(CGP_sub)) > 1) && (sum(!is.na(FAE_sub)) > 1) && (sum(!is.na(FAP_sub)) > 1) && (sum(!is.na(POH_sub)) > 1) && (sum(!is.na(POL_sub)) > 1)) {
      effect_test1 <- cohen.d(CGP_sub, FAE_sub, hedges.correction = TRUE)
      effect_size1 <- effect_test1$estimate
      effect_test2 <- cohen.d(CGP_sub, FAP_sub, hedges.correction = TRUE)
      effect_size2 <- effect_test2$estimate
      effect_test3 <- cohen.d(CGP_sub, POH_sub, hedges.correction = TRUE)
      effect_size3 <- effect_test3$estimate
      effect_test4 <- cohen.d(CGP_sub, POL_sub, hedges.correction = TRUE)
      effect_size4 <- effect_test4$estimate
      effect_test5 <- cohen.d(FAE_sub, FAP_sub, hedges.correction = TRUE)
      effect_size5 <- effect_test5$estimate
      effect_test6 <- cohen.d(FAE_sub, POH_sub, hedges.correction = TRUE)
      effect_size6 <- effect_test6$estimate
      effect_test7 <- cohen.d(FAE_sub, POL_sub, hedges.correction = TRUE)
      effect_size7 <- effect_test7$estimate
      effect_test8 <- cohen.d(FAP_sub, POH_sub, hedges.correction = TRUE)
      effect_size8 <- effect_test8$estimate
      effect_test9 <- cohen.d(FAP_sub, POL_sub, hedges.correction = TRUE)
      effect_size9 <- effect_test9$estimate
      effect_test10 <- cohen.d(POL_sub, POH_sub, hedges.correction = TRUE)
      effect_size10 <- effect_test10$estimate
      effect_sizes <- data.frame(effect = c(effect_size1, effect_size2, effect_size3, effect_size4, effect_size5, 
                                            effect_size6, effect_size7, effect_size8, effect_size9, effect_size10),
                                 treat1 = c("CGP", "CGP", "CGP", "CGP", "FAE", "FAE", "FAE", "FAP", "FAP", "POH"),
                                 treat2 = c("FAE", "FAP", "POH", "POL", "FAP", "POH", "POL", "POH", "POL", "POL"))
      
      treat_list <- c(rep("CGP", length(CGP_sub)), rep("FAE", length(FAE_sub)), rep("FAP", length(FAP_sub)), 
                      rep("POH", length(POH_sub)), rep("POL", length(POL_sub)))
      treat_res <- data.frame(value=c(CGP_sub, FAE_sub, FAP_sub, POH_sub, POL_sub), treat=treat_list)
      #do an anova, if significant for whole dataset continue with multiple comparisons
      aov_sum <- summary(aov(value~treat, treat_res))
      aov_pVal <- aov_sum[[1]]$`Pr(>F)`[1]
      if (aov_pVal <0.05) { #if significant for all treatments
        #multiple comparisons now
        pair.comp <- pairwise.t.test(treat_res$value, treat_res$treat, p.adj="none")
        pair.compM <- melt(pair.comp$p.value)
        colnames(pair.compM) <- c("treat2", "treat1", "p_value")
        pair.comp_out <- data.frame(merge(pair.compM, effect_sizes), gene=colnames(treat_res_CGP)[i])
        treat_comp <- rbind(treat_comp, pair.comp_out)
        }
      }
  }
  if (length(treat_comp) > 0) {
    treat_comp <- data.frame(treat_comp)
    colnames(treat_comp) <- c("treat2", "treat1", "p_val", "effect_size", "gene")
    treat_comp <- na.omit(treat_comp)
    
    #multiple comparisons adjustment
    treat_comp$p_val <- as.numeric(as.character(treat_comp$p_val))
    treat_comp$effect_size <- as.numeric(as.character(treat_comp$effect_size))
    treat_comp$treat1 <- as.character(treat_comp$treat1)
    treat_comp$treat2 <- as.character(treat_comp$treat2)
    treat_comp$gene <- as.character(treat_comp$gene)
    treat_comp_adj <- p.adjust(treat_comp[,3], "BH")
    treat_comp_adj[is.na(treat_comp_adj)] <- 1
    
    #cut down to only significant ones
    treat_sigcomp <- vector()
    for (i in 1:length(treat_comp_adj)) {
      if (treat_comp_adj[i] < 0.05) {
        sig <- c(treat_comp_adj[i], treat_comp$effect_size[i], treat_comp$treat1[i], treat_comp$treat2[i], treat_comp$gene[i])
        treat_sigcomp <- rbind(treat_sigcomp, sig)
      }
    }
    treat_sigcomp <- data.frame(treat_sigcomp)
    
    list_comp <- list(treat_comp, treat_comp_adj, treat_sigcomp)
    return(list_comp)
  }
}

#TREAT FOR EFFECT SIZE IS NOT FLIPPED
t.test_genes_function <- function(treat_res) { #input residual matrices for each treatment
  gene_comp <- vector()
  for (i in 1:(length(treat_res[1,])-1)) {
    counter <- i+1
    gene1_sub <- treat_res[,i]
    for (j in counter:(length(treat_res[1,]))) {
      gene2_sub <- treat_res[,j]
      if ((sum(!is.na(gene1_sub)) > 0) && (sum(!is.na(gene2_sub)) > 0)) {
        p_val <- t.test(gene1_sub, gene2_sub)$p.value
        stat <- t.test(gene1_sub, gene2_sub)$statistic
        effect_test <- cohen.d(gene1_sub, gene2_sub, hedges.correction = TRUE)
        effect_size <- effect_test$estimate
        add <- c(p_val, stat, effect_size, colnames(treat_res)[i], colnames(treat_res)[j])
        gene_comp <- rbind(gene_comp, add)
      }
    }
  }
  gene_comp <- data.frame(gene_comp)
  gene_comp <- na.omit(gene_comp)
  colnames(gene_comp) <- c("p_val", "t", "effect_size", "gene1", "gene2")
  
  #multiple comparisons adjustment
  gene_comp$p_val <- as.numeric(as.character(gene_comp$p_val))
  gene_comp$t <- as.numeric(as.character(gene_comp$t))
  gene_comp$effect_size <- as.numeric(as.character(gene_comp$effect_size))
  gene_comp$gene1 <- as.character(gene_comp$gene1)
  gene_comp$gene2 <- as.character(gene_comp$gene2)
  gene_comp_adj <- p.adjust(gene_comp[,1], "BH")
  #replace NAs with 1 (non sig) - why are there NAs??
  gene_comp_adj[is.na(gene_comp_adj)] <- 1
  
  #cut down to only significant ones
  gene_sigcomp <- vector()
  for (k in 1:length(gene_comp_adj)) {
    if (gene_comp_adj[k] < 0.05) {
      sig <- c(gene_comp_adj[k], gene_comp$t[k], gene_comp$effect_size[k], gene_comp$gene1[k], gene_comp$gene2[k])
      gene_sigcomp <- rbind(gene_sigcomp, sig)
    }
  }
  gene_sigcomp <- data.frame(gene_sigcomp)
  
  list_comp <- list(gene_comp, gene_comp_adj, gene_sigcomp)
  return(list_comp)
}

#get residuals
#probably put this back into a list of data frames?
#use gene_treat_cluster: the format is CGP_Cluster 1_gene_gene (where the last one indicates it is a GENE CLUSTER and the second to last indicates it is a GENE DATASET)
gene_cluster_residuals_list <- setNames(replicate(length(gene_treat_cluster),data.frame()),paste("residuals_", gene_treat_cluster, sep=''))
for (i in 1:length(gene_cluster_residuals_list)){
  gene_cluster_residuals_list[[i]] <- residuals_function(t(gene_expMat_cluster_treat_list[[i]]))
}

#use protein_treat_cluster: the format is CGP_Cluster 1_protein_protein (where the last one indicates it is a protein CLUSTER and the second to last indicates it is a protein DATASET)
protein_cluster_residuals_list <- setNames(replicate(length(protein_treat_cluster),data.frame()),paste("residuals_", protein_treat_cluster, sep=''))
for (i in 1:length(protein_cluster_residuals_list)){
  protein_cluster_residuals_list[[i]] <- residuals_function(scale(t(protein_expMat_cluster_treat_list[[i]])))
}

#GP comparison
gene_cluster_names_treatment_df <- data.frame(treatment=rep(c("CGP", "FAE", "FAP", "POH", "POL"), length(gene_cluster_names)),
                                              cluster=rep(gene_cluster_names, each=5))
gene_cluster_names_treatment <- paste(gene_cluster_names_treatment_df$cluster, "_", gene_cluster_names_treatment_df$treatment, sep='')
geneCluster_GP_sigcomp_list <- setNames(replicate(length(gene_cluster_names_treatment),data.frame()),paste("GPsigcomp_", gene_cluster_names_treatment, sep=''))
for (i in 1:length(gene_cluster_names_treatment)) { #gene clusters first
  GP_list_temp <- t.test_treatmentGP_function(gene_cluster_residuals_list[[i]], gene_cluster_residuals_list[[i+length(gene_cluster_names_treatment)]])
  if (length(GP_list_temp[[3]]) > 0) {
    geneCluster_GP_sigcomp_list[[i]] <- GP_list_temp[[3]]
  }
}

#nothing significant here
protein_cluster_names_treatment_df <- data.frame(treatment=rep(c("CGP", "FAE", "FAP", "POH", "POL"), length(protein_cluster_names)),
                                              cluster=rep(protein_cluster_names, each=5))
protein_cluster_names_treatment <- paste(protein_cluster_names_treatment_df$cluster, "_", protein_cluster_names_treatment_df$treatment, sep='')
proteinCluster_GP_sigcomp_list <- setNames(replicate(length(protein_cluster_names_treatment),data.frame()),paste("GPsigcomp_", protein_cluster_names_treatment, sep=''))
for (i in 1:length(protein_cluster_names_treatment)) { #protein clusters first
  GP_list_temp <- t.test_treatmentGP_function(protein_cluster_residuals_list[[i]], protein_cluster_residuals_list[[i+length(protein_cluster_names_treatment)]])
  if (length(GP_list_temp[[3]]) > 0) {
    proteinCluster_GP_sigcomp_list[[i]] <- GP_list_temp[[3]]
  }
}

#treatment (anova) comparison
counter <- 0
geneCluster_anova_sigcomp_list <- setNames(replicate(length(gene_cluster_names_d),data.frame()),paste("anova_sigcomp_", gene_cluster_names_d, sep=''))
for (i in 1:length(gene_cluster_names_d)){
  anova_list_temp <- anova_treatment_function(gene_cluster_residuals_list[[i+counter]],
                                               gene_cluster_residuals_list[[i+counter+1]],
                                               gene_cluster_residuals_list[[i+counter+2]],
                                               gene_cluster_residuals_list[[i+counter+3]],
                                               gene_cluster_residuals_list[[i+counter+4]])
  counter <- counter+4
  if (length(anova_list_temp[[3]]) > 0) {
    geneCluster_anova_sigcomp_list[[i]] <- anova_list_temp[[3]]
  }
}

#no significant comparisons
counter <- 0
proteinCluster_anova_sigcomp_list <- setNames(replicate(length(protein_cluster_names_d),data.frame()),paste("anova_sigcomp_", protein_cluster_names_d, sep=''))
for (i in 1:length(protein_cluster_names_d)){
  anova_list_temp <- anova_treatment_function(protein_cluster_residuals_list[[i+counter]],
                                              protein_cluster_residuals_list[[i+counter+1]],
                                              protein_cluster_residuals_list[[i+counter+2]],
                                              protein_cluster_residuals_list[[i+counter+3]],
                                              protein_cluster_residuals_list[[i+counter+4]])
  counter <- counter+4
  if (length(anova_list_temp[[3]]) > 0) {
    proteinCluster_anova_sigcomp_list[[i]] <- anova_list_temp[[3]]
  }
}


#gene by gene within cluster comparison
#nothing significant here
geneCluster_gg_sigcomp_list <- setNames(replicate(length(gene_treat_cluster),data.frame()),paste("GGsigcomp_", gene_treat_cluster, sep=''))
for (i in 1:length(gene_treat_cluster)) { #gene clusters first
  gg_list_temp <- t.test_genes_function(gene_cluster_residuals_list[[i]])
  if (length(gg_list_temp[[3]]) > 0) {
    geneCluster_gg_sigcomp_list[[i]] <- GP_list_temp[[3]]
  }
}

#nothing significant here
proteinCluster_gg_sigcomp_list <- setNames(replicate(length(protein_treat_cluster),data.frame()),paste("GGsigcomp_", protein_treat_cluster, sep=''))
for (i in 1:length(protein_treat_cluster)) { #protein clusters first
  gg_list_temp <- t.test_genes_function(protein_cluster_residuals_list[[i]])
  if (length(gg_list_temp[[3]]) > 0) {
    proteinCluster_gg_sigcomp_list[[i]] <- GP_list_temp[[3]]
  }
}

#export final files
#only GP and ANOVA gene clusters had any significant comparisons
#geneCluster_anova_sigcomp_list
#geneCluster_GP_sigcomp_list
all_GP_sigcomp <- rbind(geneCluster_GP_sigcomp_list$`GPsigcomp_Cluster 1_POH`, geneCluster_GP_sigcomp_list$`GPsigcomp_Cluster 1_POL`, geneCluster_GP_sigcomp_list$`GPsigcomp_Cluster 2_POL`, geneCluster_GP_sigcomp_list$`GPsigcomp_Cluster 5_POH`)
all_GP_sigcomp$cluster_treatment <- c(rep("Cluster 1_POH", length(geneCluster_GP_sigcomp_list$`GPsigcomp_Cluster 1_POH`$X1)),
                                      rep("Cluster 1_POL", length(geneCluster_GP_sigcomp_list$`GPsigcomp_Cluster 1_POL`$X1)),
                                      rep("Cluster 2_POL", length(geneCluster_GP_sigcomp_list$`GPsigcomp_Cluster 2_POL`$X1)),
                                      rep("Cluster 5_POH", length(geneCluster_GP_sigcomp_list$`GPsigcomp_Cluster 5_POH`$X1)))
colnames(all_GP_sigcomp) <- c("p_val", "t", "effect_size", "gene", "cluster_treatment")

all_anova_sigcomp <- data.frame()
for (i in 1:(length(geneCluster_anova_sigcomp_list)-1)) {
  temp <- geneCluster_anova_sigcomp_list[[i]]
  if (length(temp) > 1) {
    temp$cluster <- names(geneCluster_anova_sigcomp_list[i])
    all_anova_sigcomp <- rbind(all_anova_sigcomp, temp)
  }
}

colnames(all_anova_sigcomp) <- c("p_val", "effect_size", "treat1", "treat2", "gene", "cluster")

#add in descriptions
BLAST <- read.delim("~/Full_Blast2GO_table.csv", sep=',')
colnames(BLAST)[2] <- "gene"
all_anova_sigcompBLAST <- merge(all_anova_sigcomp, BLAST, by = "gene")
all_anova_sigcompBLAST <- all_anova_sigcompBLAST[,c(1:6, 8)]

all_GP_sigcompBLAST <- merge(all_GP_sigcomp, BLAST, by = "gene")
all_GP_sigcompBLAST <- all_GP_sigcompBLAST[,c(1:5, 7)]

write.table(all_GP_sigcompBLAST, "~/all_geneCluster_GP_sigcomp.txt", quote = F, row.names = F, sep='\t')
write.table(all_anova_sigcompBLAST, "~/all_geneCluster_anova_sigcomp.txt", quote = F, row.names = F, sep='\t')

#this is irrespective of treatment for the end of the results section
gene_res <- residuals_function(gene_matrix)
protein_res <- residuals_function(protein_matrix)
protein_res <- protein_res[row.names(protein_res) %in% row.names(gene_res),]

gene_res <- t(gene_res)
protein_res <- t(protein_res)

global_GP <- t.test_treatmentGP_function(gene_res, protein_res) #nothing significant


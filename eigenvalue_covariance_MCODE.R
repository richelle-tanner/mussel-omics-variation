#does the same thing as eigenvalue_covariance_KEGG.R except for empirically derived clusters from Cytoscape

library("stringr")
library("naniar")
library("evolqg")
library("ggplot2")
library("reshape2")

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

#import data

ddr <- "~/Documents/R_Data/musselRNASeq/"

#matching table for gene names
names <- read.table("~/Documents/Cytoscape_projects/M_Cali/gene_names_labels.txt", header=TRUE)

#full matrices
gene <- read.table(paste(ddr, "final_files/reduced_GEmatrix_by_matching_proteins.txt", sep=''), header=TRUE)
protein <- read.table(paste(ddr, "final_files/reduced_PEmatrix_by_matching_genes.txt", sep=''), header=TRUE)

#annotation
#important! use "read.delim" to specify what the delimiter is, because some of the cells have multiple delimiter types within gene names
annotation <- read.delim(paste(ddr, "final_files/reduced_annotation_matching_genes_proteins.txt", sep=''), sep="\t", header=TRUE, fill=TRUE)

setwd("~/Documents/R_Data/musselRNASeq/eigenvalue_covariance/")

####prep gene/protein matrices####
#get rid of two rows that have no variance in the gene dataset
gene1 <- gene[-c(348,536),]
#get rid of those same ones in the protein dataset
#unfortunately they're not the same row numbers
eliminate <- c("Velvet_k35_Velvet_k35_Locus_1520_Transcript_1114_Confidence_0.589_Length_3521", "TRINITY_DN20953_c0_g2_i1")
elimRowNum <- which(rownames(protein) %in% eliminate) 
protein1 <- protein[-c(elimRowNum),]

#only use individuals that overlap between datasets
gene_sim <- gene1[,order(colnames(gene1))]
protein_sim <- protein1[,order(colnames(protein1))]
#protein columns to get rid of: 6, 12 (22b, 32p)
protein <- protein_sim[,-c(6,12)]
#gene columns to get rid of: 3, 11, 19, 34, 36, 40, 41, 42, 44, 49 (18b, 30b, 41p, 63b, 6e, 81b, 82b, 82o, 85o, 96b)
gene <- gene_sim[,-c(3,11,19,34,36,40,41,42,44,49)]

#row names
gene_names <- rownames(gene)
protein_names <- rownames(protein)

#make numeric
gene <- data.frame(sapply(gene, function(x) as.numeric(as.character(x))))
protein <- data.frame(sapply(protein, function(x) as.numeric(as.character(x))))

####separate out treatments####
gene <- t(gene)
protein <- t(protein)

#add treatment label column
treatments <- c("FAE", "POL", "FAE", "FAP", "POH", "CGP", "POH", "FAE", "FAP", "FAE", "POL", "FAP", "FAP", "FAE",
                "FAP", "POH", "CGP", "CGP", "POL", "FAE", "FAP", "POH", "POH", "POL", "CGP", "POL", "POL", "POH",
                "POH", "POH", "POL", "POL", "FAE", "FAP", "POH", "CGP", "FAE", "FAP", "POH", "FAE", "FAP")

gene <- cbind(treatments, gene)
protein <- cbind(treatments, protein)

gene <- data.frame(gene)
protein <- data.frame(protein)

colnames(gene) <- c("treatments", gene_names)
colnames(protein) <- c("treatments", protein_names)

CGP_gene <- gene[which(gene$treatments == "CGP"),]
CGP_gene$treatments <- NULL
CGP_protein <- protein[which(protein$treatments == "CGP"),]
CGP_protein$treatments <- NULL
FAE_gene <- gene[which(gene$treatments == "FAE"),]
FAE_gene$treatments <- NULL
FAE_protein <- protein[which(protein$treatments == "FAE"),]
FAE_protein$treatments <- NULL
FAP_gene <- gene[which(gene$treatments == "FAP"),]
FAP_gene$treatments <- NULL
FAP_protein <- protein[which(protein$treatments == "FAP"),]
FAP_protein$treatments <- NULL
POH_gene <- gene[which(gene$treatments == "POH"),]
POH_gene$treatments <- NULL
POH_protein <- protein[which(protein$treatments == "POH"),]
POH_protein$treatments <- NULL
POL_gene <- gene[which(gene$treatments == "POL"),]
POL_gene$treatments <- NULL
POL_protein <- protein[which(protein$treatments == "POL"),]
POL_protein$treatments <- NULL

####reduce datasets to only be from cluster####
cluster_names <- names[which(names$gene_numbers %in% cluster_example),]

cluster_gPOH <- POH_gene[,which(colnames(POH_gene) %in% cluster_names$gene_names)]
cluster_gPOL <- POL_gene[,which(colnames(POL_gene) %in% cluster_names$gene_names)]
cluster_gCGP <- CGP_gene[,which(colnames(CGP_gene) %in% cluster_names$gene_names)]
cluster_gFAE <- FAE_gene[,which(colnames(FAE_gene) %in% cluster_names$gene_names)]
cluster_gFAP <- FAP_gene[,which(colnames(FAP_gene) %in% cluster_names$gene_names)]

cluster_pPOH <- POH_protein[,which(colnames(POH_protein) %in% cluster_names$gene_names)]
cluster_pPOL <- POL_protein[,which(colnames(POL_protein) %in% cluster_names$gene_names)]
cluster_pCGP <- CGP_protein[,which(colnames(CGP_protein) %in% cluster_names$gene_names)]
cluster_pFAE <- FAE_protein[,which(colnames(FAE_protein) %in% cluster_names$gene_names)]
cluster_pFAP <- FAP_protein[,which(colnames(FAP_protein) %in% cluster_names$gene_names)]

#cluster_names <- colnames(cluster_gPOH)

#####calculate eigenvalues####
dfs_cluster <- list(cluster_gPOH, cluster_gPOL, cluster_gCGP, cluster_gFAE, cluster_gFAP, cluster_pPOH, cluster_pPOL, cluster_pCGP, cluster_pFAE, cluster_pFAP)
ev_cluster <- data.frame(eigenvalue=seq(1,length(cluster_names)), gPOH=rep(0, length(cluster_names)), gPOL=rep(0, length(cluster_names)),
                       gCGP=rep(0, length(cluster_names)), gFAE=rep(0, length(cluster_names)), gFAP=rep(0, length(cluster_names)),
                       pPOH=rep(0, length(cluster_names)), pPOL=rep(0, length(cluster_names)), pCGP=rep(0, length(cluster_names)),
                       pFAE=rep(0, length(cluster_names)), pFAP=rep(0, length(cluster_names)))
eigenStat_cluster <- data.frame(treatment=rep(c("POH", "POL", "CGP", "FAE", "FAP"), gene_protein=rep(c("gene", "protein"),each=5)),
                              ICV=rep(0, length(dfs_cluster)), SD=rep(0, length(dfs_cluster)))
for (i in 1:length(dfs_cluster)) {
  df <- sapply(dfs_cluster[[i]], function(x) as.numeric(as.character(x)))
  df <- scale(df) #this turns full columns of zeros into NaN
  #keep track of which columns those are, and take them out so the rest will run
  dropset <- colnames(df)[which(colSums(is.na(df)) != 0)]
  if (length(dropset) > 0) {
    dropset_colnum <- which(colnames(df) %in% dropset)
    dropset_vals <- data.frame(num=dropset_colnum, val=rep("NA", length(dropset_colnum)))
    keepset_colnum <- (1:44)[-dropset_colnum]
    df <- df[ , apply(df, 2, function(x) !any(is.na(x)))]
    #eigenvalu calcs
    covar <- cov(df)
    eigenvals <- eigen(covar)
    #insert NAs where we had to take columns out
    keepset_vals <- data.frame(num=keepset_colnum, val=eigenvals$values)
    eigenvals_full <- rbind(keepset_vals, dropset_vals)
    ev_cluster[,(i+1)] <- eigenvals_full$val
  }
  if (length(dropset) == 0) {
    covar <- cov(df)
    eigenvals <- eigen(covar)
    ev_cluster[,(i+1)] <- eigenvals$values
  }
  eigenStat_cluster[i,2] <- CalcICV(covar)
  eigenStat_cluster[i,3] <- CalcEigenSd(covar)
}

####plot####
ev_clusterplot <- melt(ev_cluster, id.vars = 'eigenvalue', variable.name = 'treatment')
ev_clusterplot$value <- as.numeric(ev_clusterplot$value)
#only show the first 15 values
ggplot(ev_clusterplot, aes(x=eigenvalue,y=value, group=treatment)) + 
  geom_line(aes(color = treatment), size=1.5) +
  theme_bw() +
  xlim(1,15)

####annotations in cluster####
#take "cluster_names" and pull out GO terms from annotation file
cluster_ann <- annotation[which(annotation$transcript_id %in% cluster_names), ]
cluster_GO <- as.character(cluster_ann$gene_ontology_blast)
cluster_GO <- strsplit(cluster_GO, '^', fixed=TRUE)
#first 10 characters is the first GO term, just look at those for now
cluster_GO1 <- as.character(cluster_ann$gene_ontology_blast)
cluster_GO1 <- as.character(sapply(cluster_GO1, function(x) substr(x, start=1, stop=10)))
table(cluster_GO1)
#three have more than one represented:
#GO0005737 (5) cytoplasm
#GO0005829 (6) cytosol
#GO0009986 (2) cell surface

#take "cluster_names" and pull out NCBI accession numbers
cluster_nr <- as.character(cluster_ann$nr_BLASTX)
#first 14 characters is the first XP term, just look at those for now
cluster_nr1 <- as.character(sapply(cluster_nr, function(x) substr(x, start=1, stop=14)))
table(cluster_nr1)
#none overlapping

#this does two things: calculates eigenvalues for selected pathways and makes pathway maps colored by MAD (p values included separately)
#run these processes separately! both require you to go through the end of "KEGG terms" section to run 
#KEGG pathways: glycolysis, proteasome, proteolysis, metabolic processes, TCA cycle, apoptosis, MAPK signaling pathway
#only significant p values will export!

library("KEGGREST") #no longer available in current version of R
library("stringr")
library("naniar")
#library("KEGGprofile")
library("evolqg")
library("ggplot2")
library("reshape2")
library("onewaytests")
library("effsize")

MAD_function <- function(data_input) {
  #USE LOOPS
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

####import data####
setwd("~/Documents/R_Data/musselRNASeq/")
ddr <- "~/Documents/R_Data/musselRNASeq/"
#this has all of the KEGG terms in different columns
red_GP_ann_KEGG <- read.table(paste(ddr, "reduced_annotation_matching_genes_proteins_KEGG.txt", sep=''), fill=TRUE, header=TRUE)
#import the pathways made in KEGG term extraction.R
pathways <- read.table(paste(ddr, "reduced_genes_proteins_pathways.txt", sep=''), header=TRUE)
enzymes <- read.table(paste(ddr, "enzymes_list_9.4.19.txt", sep=''), header=TRUE)
#full matrices
gene <- read.table("~/gene_matrixllsImpute_6.8.20.txt")
protein <- read.table("~/protein_matrixllsImpute_6.8.20.txt") #use lls impute method 

colnames(enzymes) <- c("transcript", "enzyme")

#scale the data here
gene_matrix_scaled <- scale(t(gene))
protein_matrix_scaled <- scale(t(protein))
#flip it back for now so that everything downstream works
gene <- t(gene_matrix_scaled)
protein <- t(protein_matrix_scaled)

####separate out treatments####
gene_names <- row.names(gene)
protein_names <- row.names(protein)

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
CGP_gene <- data.frame(sapply(CGP_gene, function(x) as.numeric(as.character(x))))
CGP_protein <- protein[which(protein$treatments == "CGP"),]
CGP_protein$treatments <- NULL
CGP_protein <- data.frame(sapply(CGP_protein, function(x) as.numeric(as.character(x))))
FAE_gene <- gene[which(gene$treatments == "FAE"),]
FAE_gene$treatments <- NULL
FAE_gene <- data.frame(sapply(FAE_gene, function(x) as.numeric(as.character(x))))
FAE_protein <- protein[which(protein$treatments == "FAE"),]
FAE_protein$treatments <- NULL
FAE_protein <- data.frame(sapply(FAE_protein, function(x) as.numeric(as.character(x))))
FAP_gene <- gene[which(gene$treatments == "FAP"),]
FAP_gene$treatments <- NULL
FAP_gene <- data.frame(sapply(FAP_gene, function(x) as.numeric(as.character(x))))
FAP_protein <- protein[which(protein$treatments == "FAP"),]
FAP_protein$treatments <- NULL
FAP_protein <- data.frame(sapply(FAP_protein, function(x) as.numeric(as.character(x))))
POH_gene <- gene[which(gene$treatments == "POH"),]
POH_gene$treatments <- NULL
POH_gene <- data.frame(sapply(POH_gene, function(x) as.numeric(as.character(x))))
POH_protein <- protein[which(protein$treatments == "POH"),]
POH_protein$treatments <- NULL
POH_protein <- data.frame(sapply(POH_protein, function(x) as.numeric(as.character(x))))
POL_gene <- gene[which(gene$treatments == "POL"),]
POL_gene$treatments <- NULL
POL_gene <- data.frame(sapply(POL_gene, function(x) as.numeric(as.character(x))))
POL_protein <- protein[which(protein$treatments == "POL"),]
POL_protein$treatments <- NULL
POL_protein <- data.frame(sapply(POL_protein, function(x) as.numeric(as.character(x))))

####KEGG terms####

#map00010 - glycolysis
m00010tran_list <- as.character(pathways[which(pathways$pathway_ID %in% "path:map00010"),1])
m00010_gPOH <- POH_gene[,which(colnames(POH_gene) %in% m00010tran_list)]
m00010_gPOL <- POL_gene[,which(colnames(POL_gene) %in% m00010tran_list)]
m00010_gCGP <- CGP_gene[,which(colnames(CGP_gene) %in% m00010tran_list)]
m00010_gFAE <- FAE_gene[,which(colnames(FAE_gene) %in% m00010tran_list)]
m00010_gFAP <- FAP_gene[,which(colnames(FAP_gene) %in% m00010tran_list)]
m00010_pPOH <- POH_protein[,which(colnames(POH_protein) %in% m00010tran_list)]
m00010_pPOL <- POL_protein[,which(colnames(POL_protein) %in% m00010tran_list)]
m00010_pCGP <- CGP_protein[,which(colnames(CGP_protein) %in% m00010tran_list)]
m00010_pFAE <- FAE_protein[,which(colnames(FAE_protein) %in% m00010tran_list)]
m00010_pFAP <- FAP_protein[,which(colnames(FAP_protein) %in% m00010tran_list)]

#map03050 = proteasome
m03050tran_list <- as.character(pathways[which(pathways$pathway_ID %in% "path:map03050"),1])
m03050_gPOH <- POH_gene[,which(colnames(POH_gene) %in% m03050tran_list)]
m03050_gPOL <- POL_gene[,which(colnames(POL_gene) %in% m03050tran_list)]
m03050_gCGP <- CGP_gene[,which(colnames(CGP_gene) %in% m03050tran_list)]
m03050_gFAE <- FAE_gene[,which(colnames(FAE_gene) %in% m03050tran_list)]
m03050_gFAP <- FAP_gene[,which(colnames(FAP_gene) %in% m03050tran_list)]
m03050_pPOH <- POH_protein[,which(colnames(POH_protein) %in% m03050tran_list)]
m03050_pPOL <- POL_protein[,which(colnames(POL_protein) %in% m03050tran_list)]
m03050_pCGP <- CGP_protein[,which(colnames(CGP_protein) %in% m03050tran_list)]
m03050_pFAE <- FAE_protein[,which(colnames(FAE_protein) %in% m03050tran_list)]
m03050_pFAP <- FAP_protein[,which(colnames(FAP_protein) %in% m03050tran_list)]

#map04120 - ubiquitin-mediated proteolysis
m04120tran_list <- as.character(pathways[which(pathways$pathway_ID %in% "path:map04120"),1])
m04120_gPOH <- POH_gene[,which(colnames(POH_gene) %in% m04120tran_list)]
m04120_gPOL <- POL_gene[,which(colnames(POL_gene) %in% m04120tran_list)]
m04120_gCGP <- CGP_gene[,which(colnames(CGP_gene) %in% m04120tran_list)]
m04120_gFAE <- FAE_gene[,which(colnames(FAE_gene) %in% m04120tran_list)]
m04120_gFAP <- FAP_gene[,which(colnames(FAP_gene) %in% m04120tran_list)]
m04120_pPOH <- POH_protein[,which(colnames(POH_protein) %in% m04120tran_list)]
m04120_pPOL <- POL_protein[,which(colnames(POL_protein) %in% m04120tran_list)]
m04120_pCGP <- CGP_protein[,which(colnames(CGP_protein) %in% m04120tran_list)]
m04120_pFAE <- FAE_protein[,which(colnames(FAE_protein) %in% m04120tran_list)]
m04120_pFAP <- FAP_protein[,which(colnames(FAP_protein) %in% m04120tran_list)]

#map04010/04013 - MAPK signaling pathway (2nd term is insect)
m04010tran_list <- as.character(pathways[which(pathways$pathway_ID %in% "path:map04010"),1])
m04010_gPOH <- POH_gene[,which(colnames(POH_gene) %in% m04010tran_list)]
m04010_gPOL <- POL_gene[,which(colnames(POL_gene) %in% m04010tran_list)]
m04010_gCGP <- CGP_gene[,which(colnames(CGP_gene) %in% m04010tran_list)]
m04010_gFAE <- FAE_gene[,which(colnames(FAE_gene) %in% m04010tran_list)]
m04010_gFAP <- FAP_gene[,which(colnames(FAP_gene) %in% m04010tran_list)]
m04010_pPOH <- POH_protein[,which(colnames(POH_protein) %in% m04010tran_list)]
m04010_pPOL <- POL_protein[,which(colnames(POL_protein) %in% m04010tran_list)]
m04010_pCGP <- CGP_protein[,which(colnames(CGP_protein) %in% m04010tran_list)]
m04010_pFAE <- FAE_protein[,which(colnames(FAE_protein) %in% m04010tran_list)]
m04010_pFAP <- FAP_protein[,which(colnames(FAP_protein) %in% m04010tran_list)]

#map04210 - apoptosis
m04210tran_list <- as.character(pathways[which(pathways$pathway_ID %in% "path:map04210"),1])
m04210_gPOH <- POH_gene[,which(colnames(POH_gene) %in% m04210tran_list)]
m04210_gPOL <- POL_gene[,which(colnames(POL_gene) %in% m04210tran_list)]
m04210_gCGP <- CGP_gene[,which(colnames(CGP_gene) %in% m04210tran_list)]
m04210_gFAE <- FAE_gene[,which(colnames(FAE_gene) %in% m04210tran_list)]
m04210_gFAP <- FAP_gene[,which(colnames(FAP_gene) %in% m04210tran_list)]
m04210_pPOH <- POH_protein[,which(colnames(POH_protein) %in% m04210tran_list)]
m04210_pPOL <- POL_protein[,which(colnames(POL_protein) %in% m04210tran_list)]
m04210_pCGP <- CGP_protein[,which(colnames(CGP_protein) %in% m04210tran_list)]
m04210_pFAE <- FAE_protein[,which(colnames(FAE_protein) %in% m04210tran_list)]
m04210_pFAP <- FAP_protein[,which(colnames(FAP_protein) %in% m04210tran_list)]

#map00020 - TCA cycle
m00020tran_list <- as.character(pathways[which(pathways$pathway_ID %in% "path:map00020"),1])
m00020_gPOH <- POH_gene[,which(colnames(POH_gene) %in% m00020tran_list)]
m00020_gPOL <- POL_gene[,which(colnames(POL_gene) %in% m00020tran_list)]
m00020_gCGP <- CGP_gene[,which(colnames(CGP_gene) %in% m00020tran_list)]
m00020_gFAE <- FAE_gene[,which(colnames(FAE_gene) %in% m00020tran_list)]
m00020_gFAP <- FAP_gene[,which(colnames(FAP_gene) %in% m00020tran_list)]
m00020_pPOH <- POH_protein[,which(colnames(POH_protein) %in% m00020tran_list)]
m00020_pPOL <- POL_protein[,which(colnames(POL_protein) %in% m00020tran_list)]
m00020_pCGP <- CGP_protein[,which(colnames(CGP_protein) %in% m00020tran_list)]
m00020_pFAE <- FAE_protein[,which(colnames(FAE_protein) %in% m00020tran_list)]
m00020_pFAP <- FAP_protein[,which(colnames(FAP_protein) %in% m00020tran_list)]


#get names to match later
#m01100_names <- colnames(m01100_gPOH)
m00010_names <- colnames(m00010_gPOH)
m03050_names <- colnames(m03050_gPOH)
m04120_names <- colnames(m04120_gPOH)
m04010_names <- colnames(m04010_gPOH)
m04210_names <- colnames(m04210_gPOH)
m00020_names <- colnames(m00020_gPOH)

#export names to compare with MCODE clusters
namesExport <- data.frame(transcript=c(m00010_names, m03050_names, m04120_names, m04010_names,
                                       m04210_names, m00020_names), 
                          KEGG_pathway=c(rep("m00010", length(m00010_names)),
                                         rep("m03050", length(m03050_names)), rep("m04120", length(m04120_names)),
                                         rep("m04010", length(m04010_names)), rep("m04210", length(m04210_names)),
                                         rep("m00020", length(m00020_names))))
write.table(namesExport, "~/KEGG_names_allPathways.txt", sep='\t', quote=F, row.names = F)

####make covariance matrices for each treatment for each pathway and calculate eigenvalues/stats####

#00010 - glycolysis
dfs_00010 <- list(m00010_gPOH, m00010_gPOL, m00010_gCGP, m00010_gFAE, m00010_gFAP, m00010_pPOH, m00010_pPOL, m00010_pCGP, m00010_pFAE, m00010_pFAP)
ev_00010 <- data.frame(eigenvalue=seq(1,length(m00010_names)), gPOH=rep(0, length(m00010_names)), gPOL=rep(0, length(m00010_names)),
                       gCGP=rep(0, length(m00010_names)), gFAE=rep(0, length(m00010_names)), gFAP=rep(0, length(m00010_names)),
                       pPOH=rep(0, length(m00010_names)), pPOL=rep(0, length(m00010_names)), pCGP=rep(0, length(m00010_names)),
                       pFAE=rep(0, length(m00010_names)), pFAP=rep(0, length(m00010_names)))
eigenStat_00010 <- data.frame(treatment=rep(c("POH", "POL", "CGP", "FAE", "FAP"), gene_protein=rep(c("gene", "protein"),each=5)),
                              ICV=rep(0, length(dfs_00010)), SD=rep(0, length(dfs_00010)))
for (i in 1:length(dfs_00010)) {
  df <- sapply(dfs_00010[[i]], function(x) as.numeric(as.character(x)))
  for (j in 1:length(df[1,])) {
    if (sum(df[,j]) == 0 | sum(df[,j]) == "NaN") {
      df[,j] <- NA
    }
  }
  dropset <- colnames(df)[which(colSums(is.na(df)) != 0)]
  if (length(dropset) > 0) {
    dropset_colnum <- which(colnames(df) %in% dropset)
    dropset_vals <- data.frame(num=dropset_colnum, val=rep("NA", length(dropset_colnum)))
    keepset_colnum <- (1:20)[-dropset_colnum]
    df <- df[ , apply(df, 2, function(x) !any(is.na(x)))]
    #eigenvalu calcs
    covar <- cov(df)
    eigenvals <- eigen(covar)
    #insert NAs where we had to take columns out
    keepset_vals <- data.frame(num=keepset_colnum, val=eigenvals$values)
    eigenvals_full <- rbind(keepset_vals, dropset_vals)
    ev_00010[,(i+1)] <- eigenvals_full$val
  }
  if (length(dropset) == 0) {
    covar <- cov(df)
    eigenvals <- eigen(covar)
    ev_00010[,(i+1)] <- eigenvals$values
  }
  eigenStat_00010[i,2] <- CalcICV(covar)
  eigenStat_00010[i,3] <- CalcEigenVar(covar, sd=TRUE, sample=length(eigenvals$values))
}

#03050 - proteasome
dfs_03050 <- list(m03050_gPOH, m03050_gPOL, m03050_gCGP, m03050_gFAE, m03050_gFAP, m03050_pPOH, m03050_pPOL, m03050_pCGP, m03050_pFAE, m03050_pFAP)
ev_03050 <- data.frame(eigenvalue=seq(1,length(m03050_names)), gPOH=rep(0, length(m03050_names)), gPOL=rep(0, length(m03050_names)),
                       gCGP=rep(0, length(m03050_names)), gFAE=rep(0, length(m03050_names)), gFAP=rep(0, length(m03050_names)),
                       pPOH=rep(0, length(m03050_names)), pPOL=rep(0, length(m03050_names)), pCGP=rep(0, length(m03050_names)),
                       pFAE=rep(0, length(m03050_names)), pFAP=rep(0, length(m03050_names)))
eigenStat_03050 <- data.frame(treatment=rep(c("POH", "POL", "CGP", "FAE", "FAP"), gene_protein=rep(c("gene", "protein"),each=5)),
                              ICV=rep(0, length(dfs_03050)), SD=rep(0, length(dfs_03050)))
for (i in 1:length(dfs_03050)) {
  df <- sapply(dfs_03050[[i]], function(x) as.numeric(as.character(x)))
  for (j in 1:length(df[1,])) {
    if (sum(df[,j]) == 0 | sum(df[,j]) == "NaN") {
      df[,j] <- NA
    }
  }
  dropset <- colnames(df)[which(colSums(is.na(df)) != 0)]
  if (length(dropset) > 0) {
    dropset_colnum <- which(colnames(df) %in% dropset)
    dropset_vals <- data.frame(num=dropset_colnum, val=rep("NA", length(dropset_colnum)))
    keepset_colnum <- (1:19)[-dropset_colnum]
    df <- df[ , apply(df, 2, function(x) !any(is.na(x)))]
    #eigenvalu calcs
    covar <- cov(df)
    eigenvals <- eigen(covar)
    #insert NAs where we had to take columns out
    keepset_vals <- data.frame(num=keepset_colnum, val=eigenvals$values)
    eigenvals_full <- rbind(keepset_vals, dropset_vals)
    ev_03050[,(i+1)] <- eigenvals_full$val
  }
  if (length(dropset) == 0) {
    covar <- cov(df)
    eigenvals <- eigen(covar)
    ev_03050[,(i+1)] <- eigenvals$values
  }
  eigenStat_03050[i,2] <- CalcICV(covar)
  eigenStat_03050[i,3] <- CalcEigenVar(covar, sd=TRUE, sample=length(eigenvals$values))
}

#04120 - ubiquitin-mediated proteolysis
dfs_04120 <- list(m04120_gPOH, m04120_gPOL, m04120_gCGP, m04120_gFAE, m04120_gFAP, m04120_pPOH, m04120_pPOL, m04120_pCGP, m04120_pFAE, m04120_pFAP)
ev_04120 <- data.frame(eigenvalue=seq(1,length(m04120_names)), gPOH=rep(0, length(m04120_names)), gPOL=rep(0, length(m04120_names)),
                       gCGP=rep(0, length(m04120_names)), gFAE=rep(0, length(m04120_names)), gFAP=rep(0, length(m04120_names)),
                       pPOH=rep(0, length(m04120_names)), pPOL=rep(0, length(m04120_names)), pCGP=rep(0, length(m04120_names)),
                       pFAE=rep(0, length(m04120_names)), pFAP=rep(0, length(m04120_names)))
eigenStat_04120 <- data.frame(treatment=rep(c("POH", "POL", "CGP", "FAE", "FAP"), gene_protein=rep(c("gene", "protein"),each=5)),
                              ICV=rep(0, length(dfs_04120)), SD=rep(0, length(dfs_04120)))
for (i in 1:length(dfs_04120)) {
  df <- sapply(dfs_04120[[i]], function(x) as.numeric(as.character(x)))
  for (j in 1:length(df[1,])) {
    if (sum(df[,j]) == 0 | sum(df[,j]) == "NaN") {
      df[,j] <- NA
    }
  }
  dropset <- colnames(df)[which(colSums(is.na(df)) != 0)]
  if (length(dropset) > 0) {
    dropset_colnum <- which(colnames(df) %in% dropset)
    dropset_vals <- data.frame(num=dropset_colnum, val=rep("NA", length(dropset_colnum)))
    keepset_colnum <- (1:4)[-dropset_colnum]
    df <- df[ , apply(df, 2, function(x) !any(is.na(x)))]
    #eigenvalu calcs
    covar <- cov(df)
    eigenvals <- eigen(covar)
    #insert NAs where we had to take columns out
    keepset_vals <- data.frame(num=keepset_colnum, val=eigenvals$values)
    eigenvals_full <- rbind(keepset_vals, dropset_vals)
    ev_04120[,(i+1)] <- eigenvals_full$val
  }
  if (length(dropset) == 0) {
    covar <- cov(df)
    eigenvals <- eigen(covar)
    ev_04120[,(i+1)] <- eigenvals$values
  }
  eigenStat_04120[i,2] <- CalcICV(covar)
  eigenStat_04120[i,3] <- CalcEigenVar(covar, sd=TRUE, sample=length(eigenvals$values))
}

#04010 - MAPK signaling pathway
dfs_04010 <- list(m04010_gPOH, m04010_gPOL, m04010_gCGP, m04010_gFAE, m04010_gFAP, m04010_pPOH, m04010_pPOL, m04010_pCGP, m04010_pFAE, m04010_pFAP)
ev_04010 <- data.frame(eigenvalue=seq(1,length(m04010_names)), gPOH=rep(0, length(m04010_names)), gPOL=rep(0, length(m04010_names)),
                       gCGP=rep(0, length(m04010_names)), gFAE=rep(0, length(m04010_names)), gFAP=rep(0, length(m04010_names)),
                       pPOH=rep(0, length(m04010_names)), pPOL=rep(0, length(m04010_names)), pCGP=rep(0, length(m04010_names)),
                       pFAE=rep(0, length(m04010_names)), pFAP=rep(0, length(m04010_names)))
eigenStat_04010 <- data.frame(treatment=rep(c("POH", "POL", "CGP", "FAE", "FAP"), gene_protein=rep(c("gene", "protein"),each=5)),
                              ICV=rep(0, length(dfs_04010)), SD=rep(0, length(dfs_04010)))
for (i in 1:length(dfs_04010)) {
  df <- sapply(dfs_04010[[i]], function(x) as.numeric(as.character(x)))
  for (j in 1:length(df[1,])) {
    if (sum(df[,j]) == 0 | sum(df[,j]) == "NaN") {
      df[,j] <- NA
    }
  }
  dropset <- colnames(df)[which(colSums(is.na(df)) != 0)]
  if (length(dropset) > 0) {
    dropset_colnum <- which(colnames(df) %in% dropset)
    dropset_vals <- data.frame(num=dropset_colnum, val=rep("NA", length(dropset_colnum)))
    keepset_colnum <- (1:18)[-dropset_colnum]
    df <- df[ , apply(df, 2, function(x) !any(is.na(x)))]
    #eigenvalu calcs
    covar <- cov(df)
    eigenvals <- eigen(covar)
    #insert NAs where we had to take columns out
    keepset_vals <- data.frame(num=keepset_colnum, val=eigenvals$values)
    eigenvals_full <- rbind(keepset_vals, dropset_vals)
    ev_04010[,(i+1)] <- eigenvals_full$val
  }
  if (length(dropset) == 0) {
    covar <- cov(df)
    eigenvals <- eigen(covar)
    ev_04010[,(i+1)] <- eigenvals$values
  }
  eigenStat_04010[i,2] <- CalcICV(covar)
  eigenStat_04010[i,3] <- CalcEigenVar(covar, sd=TRUE, sample=length(eigenvals$values))
}

#map04210 - apoptosis
dfs_04210 <- list(m04210_gPOH, m04210_gPOL, m04210_gCGP, m04210_gFAE, m04210_gFAP, m04210_pPOH, m04210_pPOL, m04210_pCGP, m04210_pFAE, m04210_pFAP)
ev_04210 <- data.frame(eigenvalue=seq(1,length(m04210_names)), gPOH=rep(0, length(m04210_names)), gPOL=rep(0, length(m04210_names)),
                       gCGP=rep(0, length(m04210_names)), gFAE=rep(0, length(m04210_names)), gFAP=rep(0, length(m04210_names)),
                       pPOH=rep(0, length(m04210_names)), pPOL=rep(0, length(m04210_names)), pCGP=rep(0, length(m04210_names)),
                       pFAE=rep(0, length(m04210_names)), pFAP=rep(0, length(m04210_names)))
eigenStat_04210 <- data.frame(treatment=rep(c("POH", "POL", "CGP", "FAE", "FAP"), gene_protein=rep(c("gene", "protein"),each=5)),
                              ICV=rep(0, length(dfs_04210)), SD=rep(0, length(dfs_04210)))
for (i in 1:length(dfs_04210)) {
  df <- sapply(dfs_04210[[i]], function(x) as.numeric(as.character(x)))
  for (j in 1:length(df[1,])) {
    if (sum(df[,j]) == 0 | sum(df[,j]) == "NaN") {
      df[,j] <- NA
    }
  }
  dropset <- colnames(df)[which(colSums(is.na(df)) != 0)]
  if (length(dropset) > 0) {
    dropset_colnum <- which(colnames(df) %in% dropset)
    dropset_vals <- data.frame(num=dropset_colnum, val=rep("NA", length(dropset_colnum)))
    keepset_colnum <- (1:20)[-dropset_colnum]
    df <- df[ , apply(df, 2, function(x) !any(is.na(x)))]
    #eigenvalu calcs
    covar <- cov(df)
    eigenvals <- eigen(covar)
    #insert NAs where we had to take columns out
    keepset_vals <- data.frame(num=keepset_colnum, val=eigenvals$values)
    eigenvals_full <- rbind(keepset_vals, dropset_vals)
    ev_04210[,(i+1)] <- eigenvals_full$val
  }
  if (length(dropset) == 0) {
    covar <- cov(df)
    eigenvals <- eigen(covar)
    ev_04210[,(i+1)] <- eigenvals$values
  }
  eigenStat_04210[i,2] <- CalcICV(covar)
  eigenStat_04210[i,3] <- CalcEigenVar(covar, sd=TRUE, sample=length(eigenvals$values))
}

#00020 - TCA cycle
dfs_00020 <- list(m00020_gPOH, m00020_gPOL, m00020_gCGP, m00020_gFAE, m00020_gFAP, m00020_pPOH, m00020_pPOL, m00020_pCGP, m00020_pFAE, m00020_pFAP)
ev_00020 <- data.frame(eigenvalue=seq(1,length(m00020_names)), gPOH=rep(0, length(m00020_names)), gPOL=rep(0, length(m00020_names)),
                       gCGP=rep(0, length(m00020_names)), gFAE=rep(0, length(m00020_names)), gFAP=rep(0, length(m00020_names)),
                       pPOH=rep(0, length(m00020_names)), pPOL=rep(0, length(m00020_names)), pCGP=rep(0, length(m00020_names)),
                       pFAE=rep(0, length(m00020_names)), pFAP=rep(0, length(m00020_names)))
eigenStat_00020 <- data.frame(treatment=rep(c("POH", "POL", "CGP", "FAE", "FAP"), gene_protein=rep(c("gene", "protein"),each=5)),
                              ICV=rep(0, length(dfs_00020)), SD=rep(0, length(dfs_00020)))
for (i in 1:length(dfs_00020)) {
  df <- sapply(dfs_00020[[i]], function(x) as.numeric(as.character(x)))
  for (j in 1:length(df[1,])) {
    if (sum(df[,j]) == 0 | sum(df[,j]) == "NaN") {
      df[,j] <- NA
    }
  }
  dropset <- colnames(df)[which(colSums(is.na(df)) != 0)]
  if (length(dropset) > 0) {
    dropset_colnum <- which(colnames(df) %in% dropset)
    dropset_vals <- data.frame(num=dropset_colnum, val=rep("NA", length(dropset_colnum)))
    keepset_colnum <- (1:22)[-dropset_colnum]
    df <- df[ , apply(df, 2, function(x) !any(is.na(x)))]
    #eigenvalu calcs
    covar <- cov(df)
    eigenvals <- eigen(covar)
    #insert NAs where we had to take columns out
    keepset_vals <- data.frame(num=keepset_colnum, val=eigenvals$values)
    eigenvals_full <- rbind(keepset_vals, dropset_vals)
    ev_00020[,(i+1)] <- eigenvals_full$val
  }
  if (length(dropset) == 0) {
    covar <- cov(df)
    eigenvals <- eigen(covar)
    ev_00020[,(i+1)] <- eigenvals$values
  }
  eigenStat_00020[i,2] <- CalcICV(covar)
  eigenStat_00020[i,3] <- CalcEigenVar(covar, sd=TRUE, sample=length(eigenvals$values))
}

####eigenvalue plots and exporting data####

ev_00010plot <- melt(ev_00010, id.vars = 'eigenvalue', variable.name = 'treatment')
ev_03050plot <- melt(ev_03050, id.vars = 'eigenvalue', variable.name = 'treatment')
ev_04120plot <- melt(ev_04120, id.vars = 'eigenvalue', variable.name = 'treatment')
ev_04010plot <- melt(ev_04010, id.vars = 'eigenvalue', variable.name = 'treatment')
ev_04210plot <- melt(ev_04210, id.vars = 'eigenvalue', variable.name = 'treatment')
ev_00020plot <- melt(ev_00020, id.vars = 'eigenvalue', variable.name = 'treatment')

ev_00010plot$value <- as.numeric(ev_00010plot$value)
ev_03050plot$value <- as.numeric(ev_03050plot$value)
ev_04120plot$value <- as.numeric(ev_04120plot$value)
ev_04010plot$value <- as.numeric(ev_04010plot$value)
ev_04210plot$value <- as.numeric(ev_04210plot$value)
ev_00020plot$value <- as.numeric(ev_00020plot$value)

ev_00010plot$data_type <- str_sub(ev_00010plot$treatment, 1, 1)
ev_03050plot$data_type <- str_sub(ev_03050plot$treatment, 1, 1)
ev_04120plot$data_type <- str_sub(ev_04120plot$treatment, 1, 1)
ev_04010plot$data_type <- str_sub(ev_04010plot$treatment, 1, 1)
ev_04210plot$data_type <- str_sub(ev_04210plot$treatment, 1, 1)
ev_00020plot$data_type <- str_sub(ev_00020plot$treatment, 1, 1)

#only show the first 15 values - except 04120 (only 4 eigenvalues)
#colors
eigenPlot_colorMan <- rep(c("red", "purple", "blue", "green", "orange"), 2)

ggplot(ev_00010plot, aes(x=eigenvalue,y=value, group=treatment)) + 
  geom_line(aes(color = treatment, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_00010.pdf")

ggplot(ev_03050plot, aes(x=eigenvalue,y=value, group=treatment)) + 
  geom_line(aes(color = treatment, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_03050.pdf")

ggplot(ev_04120plot, aes(x=eigenvalue,y=value, group=treatment)) + 
  geom_line(aes(color = treatment, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,4) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_04120.pdf")

ggplot(ev_04010plot, aes(x=eigenvalue,y=value, group=treatment)) + 
  geom_line(aes(color = treatment, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_04010.pdf")

ggplot(ev_04210plot, aes(x=eigenvalue,y=value, group=treatment)) + 
  geom_line(aes(color = treatment, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_04210.pdf")

ggplot(ev_00020plot, aes(x=eigenvalue,y=value, group=treatment)) + 
  geom_line(aes(color = treatment, linetype=data_type), size=1.5) +
  theme_bw() +
  scale_color_manual(values=eigenPlot_colorMan) + 
  xlim(1,15) +
  xlab("Eigenvector") +
  ylab("Eigenvalue") +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=25)) 
ggsave("~/Eigenvalue_gene_00020.pdf")

#export stats
write.table(eigenStat_00010, "eigen_variance_00010.txt", row.names=FALSE, sep="\t", quote=F)
write.table(eigenStat_03050, "eigen_variance_03050.txt", row.names=FALSE, sep="\t", quote=F)
write.table(eigenStat_04120, "eigen_variance_04120.txt", row.names=FALSE, sep="\t", quote=F)
write.table(eigenStat_04010, "eigen_variance_04010.txt", row.names=FALSE, sep="\t", quote=F)
write.table(eigenStat_04210, "eigen_variance_04210.txt", row.names=FALSE, sep="\t", quote=F)
write.table(eigenStat_00020, "eigen_variance_00020.txt", row.names=FALSE, sep="\t", quote=F)

#export eigenvalues
write.table(ev_00010, "eigenvalues_00010.txt", row.names=FALSE, sep="\t", quote=F)
write.table(ev_03050, "eigenvalues_03050.txt", row.names=FALSE, sep="\t", quote=F)
write.table(ev_04120, "eigenvalues_04120.txt", row.names=FALSE, sep="\t", quote=F)
write.table(ev_04010, "eigenvalues_04010.txt", row.names=FALSE, sep="\t", quote=F)
write.table(ev_04210, "eigenvalues_04210.txt", row.names=FALSE, sep="\t", quote=F)
write.table(ev_00020, "eigenvalues_00020.txt", row.names=FALSE, sep="\t", quote=F)

####calculate MAD for each treatment####

MAD_m00010_gPOH <- MAD_function(m00010_gPOH)
MAD_m00010_gPOL <- MAD_function(m00010_gPOL)
MAD_m00010_gCGP <- MAD_function(m00010_gCGP)
MAD_m00010_gFAE <- MAD_function(m00010_gFAE)
MAD_m00010_gFAP <- MAD_function(m00010_gFAP)
MAD_m00010_pPOH <- MAD_function(m00010_pPOH)
MAD_m00010_pPOL <- MAD_function(m00010_pPOL)
MAD_m00010_pCGP <- MAD_function(m00010_pCGP)
MAD_m00010_pFAE <- MAD_function(m00010_pFAE)
MAD_m00010_pFAP <- MAD_function(m00010_pFAP)

MAD_m03050_gPOH <- MAD_function(m03050_gPOH)
MAD_m03050_gPOL <- MAD_function(m03050_gPOL)
MAD_m03050_gCGP <- MAD_function(m03050_gCGP)
MAD_m03050_gFAE <- MAD_function(m03050_gFAE)
MAD_m03050_gFAP <- MAD_function(m03050_gFAP)
MAD_m03050_pPOH <- MAD_function(m03050_pPOH)
MAD_m03050_pPOL <- MAD_function(m03050_pPOL)
MAD_m03050_pCGP <- MAD_function(m03050_pCGP)
MAD_m03050_pFAE <- MAD_function(m03050_pFAE)
MAD_m03050_pFAP <- MAD_function(m03050_pFAP)

MAD_m04120_gPOH <- MAD_function(m04120_gPOH)
MAD_m04120_gPOL <- MAD_function(m04120_gPOL)
MAD_m04120_gCGP <- MAD_function(m04120_gCGP)
MAD_m04120_gFAE <- MAD_function(m04120_gFAE)
MAD_m04120_gFAP <- MAD_function(m04120_gFAP)
MAD_m04120_pPOH <- MAD_function(m04120_pPOH)
MAD_m04120_pPOL <- MAD_function(m04120_pPOL)
MAD_m04120_pCGP <- MAD_function(m04120_pCGP)
MAD_m04120_pFAE <- MAD_function(m04120_pFAE)
MAD_m04120_pFAP <- MAD_function(m04120_pFAP)

MAD_m04010_gPOH <- MAD_function(m04010_gPOH)
MAD_m04010_gPOL <- MAD_function(m04010_gPOL)
MAD_m04010_gCGP <- MAD_function(m04010_gCGP)
MAD_m04010_gFAE <- MAD_function(m04010_gFAE)
MAD_m04010_gFAP <- MAD_function(m04010_gFAP)
MAD_m04010_pPOH <- MAD_function(m04010_pPOH)
MAD_m04010_pPOL <- MAD_function(m04010_pPOL)
MAD_m04010_pCGP <- MAD_function(m04010_pCGP)
MAD_m04010_pFAE <- MAD_function(m04010_pFAE)
MAD_m04010_pFAP <- MAD_function(m04010_pFAP)

MAD_m04210_gPOH <- MAD_function(m04210_gPOH)
MAD_m04210_gPOL <- MAD_function(m04210_gPOL)
MAD_m04210_gCGP <- MAD_function(m04210_gCGP)
MAD_m04210_gFAE <- MAD_function(m04210_gFAE)
MAD_m04210_gFAP <- MAD_function(m04210_gFAP)
MAD_m04210_pPOH <- MAD_function(m04210_pPOH)
MAD_m04210_pPOL <- MAD_function(m04210_pPOL)
MAD_m04210_pCGP <- MAD_function(m04210_pCGP)
MAD_m04210_pFAE <- MAD_function(m04210_pFAE)
MAD_m04210_pFAP <- MAD_function(m04210_pFAP)

MAD_m00020_gPOH <- MAD_function(m00020_gPOH)
MAD_m00020_gPOL <- MAD_function(m00020_gPOL)
MAD_m00020_gCGP <- MAD_function(m00020_gCGP)
MAD_m00020_gFAE <- MAD_function(m00020_gFAE)
MAD_m00020_gFAP <- MAD_function(m00020_gFAP)
MAD_m00020_pPOH <- MAD_function(m00020_pPOH)
MAD_m00020_pPOL <- MAD_function(m00020_pPOL)
MAD_m00020_pCGP <- MAD_function(m00020_pCGP)
MAD_m00020_pFAE <- MAD_function(m00020_pFAE)
MAD_m00020_pFAP <- MAD_function(m00020_pFAP)


####KEGG pathway maps with MAD colors####
#export MAD values for reference
#each MAD value is for a different gene in the pathway

#make a color gradient based on min and max values across all treatments
#use the same palette for everything
pal <- colorRampPalette(c("yellow", "red")) #changed to have yellow (low) to red (high) 11.7.19!! it was reversed before
#since they're z scored, make ONE palette for both gene and protein datasets!
#make a new dataset that has a column for treatment, column for gene/protein, and column for gene name
treat_step <- c("POL", "POH", "CGP", "FAE", "FAP")

##00010----
#colors
colors_00010_MADdf <- data.frame(MAD=c(MAD_m00010_gPOL, MAD_m00010_gPOH, MAD_m00010_gCGP, MAD_m00010_gFAE, MAD_m00010_gFAP,
                                       MAD_m00010_pPOL, MAD_m00010_pPOH, MAD_m00010_pCGP, MAD_m00010_pFAE, MAD_m00010_pFAP),
                                 Treatment=rep(rep(treat_step, each=length(MAD_m00010_gPOH)),2), 
                                 Type=rep(c("gene", "protein"), each=(length(MAD_m00010_gPOH)*5)), 
                                 Name=rep(summed_names_00010, 10))
colors_00010_MADdf$order <- findInterval(colors_00010_MADdf$MAD, sort(colors_00010_MADdf$MAD))
#get the code for the vector
colors_00010_MADdf$colors <- pal(nrow(colors_00010_MADdf))[colors_00010_MADdf$order]

#make separate vectors for each treatment
#order: "POL" "POH" "CGP" "FAE" "FAP"
color_00010_MADg_POL <- colors_00010_MADdf[c(1:length(summed_names_00010)),6]
color_00010_MADg_POH <- colors_00010_MADdf[c((length(summed_names_00010)+1):(length(summed_names_00010)*2)),6]
color_00010_MADg_CGP <- colors_00010_MADdf[c((length(summed_names_00010)*2+1):(length(summed_names_00010)*3)),6]
color_00010_MADg_FAE <- colors_00010_MADdf[c((length(summed_names_00010)*3+1):(length(summed_names_00010)*4)),6]
color_00010_MADg_FAP <- colors_00010_MADdf[c((length(summed_names_00010)*4+1):(length(summed_names_00010)*5)),6]

color_00010_MADp_POL <- colors_00010_MADdf[c((length(summed_names_00010)*5+1):(length(summed_names_00010)*6)),6]
color_00010_MADp_POH <- colors_00010_MADdf[c((length(summed_names_00010)*6+1):(length(summed_names_00010)*7)),6]
color_00010_MADp_CGP <- colors_00010_MADdf[c((length(summed_names_00010)*7+1):(length(summed_names_00010)*8)),6]
color_00010_MADp_FAE <- colors_00010_MADdf[c((length(summed_names_00010)*8+1):(length(summed_names_00010)*9)),6]
color_00010_MADp_FAP <- colors_00010_MADdf[c((length(summed_names_00010)*9+1):(length(summed_names_00010)*10)),6]

#assign color vector
#black color for the writing - needs to be the length of the names list
black_color <- rep("#000000", length(summed_names_00010))
#make pathway map
url00010_MADg_POL <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADg_POL, black_color)
url00010_MADg_POH <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADg_POH, black_color)
url00010_MADg_CGP <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADg_CGP, black_color)
url00010_MADg_FAE <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADg_FAE, black_color)
url00010_MADg_FAP <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADg_FAP, black_color)

url00010_MADp_POL <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADp_POL, black_color)
url00010_MADp_POH <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADp_POH, black_color)
url00010_MADp_CGP <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADp_CGP, black_color)
url00010_MADp_FAE <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADp_FAE, black_color)
url00010_MADp_FAP <- color.pathway.by.objects("path:map00010", m00010_KO_terms, color_00010_MADp_FAP, black_color)

#make a list and export it
m00010_urls <- data.frame(URL=c(url00010_MADg_POL, url00010_MADg_POH, url00010_MADg_CGP, url00010_MADg_FAE, url00010_MADg_FAP, 
                                url00010_MADp_POL, url00010_MADp_POH, url00010_MADp_CGP, url00010_MADp_FAE, url00010_MADp_FAP), 
                          Treatment=rep(treat_step, 2), 
                          Type=rep(c("gene", "protein"),each=5))
write.table(m00010_urls, "pathway_maps/m00010_colored_expression_urls.txt", quote=F, sep="\t", row.names=F)
write.table(colors_00010_MADdf, "m00010_MAD_values.txt", quote=F, sep="\t", row.names=F)

##03050----
#names vector is not the summed one, it is the original one because there are no repeats
#colors
colors_03050_MADdf <- data.frame(MAD=c(MAD_m03050_gPOL, MAD_m03050_gPOH, MAD_m03050_gCGP, MAD_m03050_gFAE, MAD_m03050_gFAP,
                                       MAD_m03050_pPOL, MAD_m03050_pPOH, MAD_m03050_pCGP, MAD_m03050_pFAE, MAD_m03050_pFAP),
                                 Treatment=rep(rep(treat_step, each=length(MAD_m03050_gPOH)),2), 
                                 Type=rep(c("gene", "protein"), each=(length(MAD_m03050_gPOH)*5)), 
                                 Name=rep(m03050_names, 10))
colors_03050_MADdf$order <- findInterval(colors_03050_MADdf$MAD, sort(colors_03050_MADdf$MAD))
#get the code for the vector
colors_03050_MADdf$colors <- pal(nrow(colors_03050_MADdf))[colors_03050_MADdf$order]

#make separate vectors for each treatment
#order: "POL" "POH" "CGP" "FAE" "FAP"
color_03050_MADg_POL <- colors_03050_MADdf[c(1:length(m03050_names)),6]
color_03050_MADg_POH <- colors_03050_MADdf[c((length(m03050_names)+1):(length(m03050_names)*2)),6]
color_03050_MADg_CGP <- colors_03050_MADdf[c((length(m03050_names)*2+1):(length(m03050_names)*3)),6]
color_03050_MADg_FAE <- colors_03050_MADdf[c((length(m03050_names)*3+1):(length(m03050_names)*4)),6]
color_03050_MADg_FAP <- colors_03050_MADdf[c((length(m03050_names)*4+1):(length(m03050_names)*5)),6]

color_03050_MADp_POL <- colors_03050_MADdf[c((length(m03050_names)*5+1):(length(m03050_names)*6)),6]
color_03050_MADp_POH <- colors_03050_MADdf[c((length(m03050_names)*6+1):(length(m03050_names)*7)),6]
color_03050_MADp_CGP <- colors_03050_MADdf[c((length(m03050_names)*7+1):(length(m03050_names)*8)),6]
color_03050_MADp_FAE <- colors_03050_MADdf[c((length(m03050_names)*8+1):(length(m03050_names)*9)),6]
color_03050_MADp_FAP <- colors_03050_MADdf[c((length(m03050_names)*9+1):(length(m03050_names)*10)),6]

#assign color vector
#black color for the writing - needs to be the length of the names list
black_color <- rep("#000000", length(m03050_names))
#make pathway map
url03050_MADg_POL <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADg_POL, black_color)
url03050_MADg_POH <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADg_POH, black_color)
url03050_MADg_CGP <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADg_CGP, black_color)
url03050_MADg_FAE <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADg_FAE, black_color)
url03050_MADg_FAP <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADg_FAP, black_color)

url03050_MADp_POL <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADp_POL, black_color)
url03050_MADp_POH <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADp_POH, black_color)
url03050_MADp_CGP <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADp_CGP, black_color)
url03050_MADp_FAE <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADp_FAE, black_color)
url03050_MADp_FAP <- color.pathway.by.objects("path:map03050", m03050_KO_terms, color_03050_MADp_FAP, black_color)

#make a list and export it
m03050_urls <- data.frame(URL=c(url03050_MADg_POL, url03050_MADg_POH, url03050_MADg_CGP, url03050_MADg_FAE, url03050_MADg_FAP, 
                                url03050_MADp_POL, url03050_MADp_POH, url03050_MADp_CGP, url03050_MADp_FAE, url03050_MADp_FAP), 
                          Treatment=rep(treat_step, 2), 
                          Type=rep(c("gene", "protein"),each=5))
write.table(m03050_urls, "pathway_maps/m03050_colored_expression_urls.txt", quote=F, sep="\t", row.names=F)
write.table(colors_03050_MADdf, "m03050_MAD_values.txt", quote=F, sep="\t", row.names=F)

##04120----
#names vector is not the summed one, it is the original one because there are no repeats
#colors
colors_04120_MADdf <- data.frame(MAD=c(MAD_m04120_gPOL, MAD_m04120_gPOH, MAD_m04120_gCGP, MAD_m04120_gFAE, MAD_m04120_gFAP,
                                       MAD_m04120_pPOL, MAD_m04120_pPOH, MAD_m04120_pCGP, MAD_m04120_pFAE, MAD_m04120_pFAP),
                                 Treatment=rep(rep(treat_step, each=length(MAD_m04120_gPOH)),2), 
                                 Type=rep(c("gene", "protein"), each=(length(MAD_m04120_gPOH)*5)), 
                                 Name=rep(m04120_names, 10))
colors_04120_MADdf$order <- findInterval(colors_04120_MADdf$MAD, sort(colors_04120_MADdf$MAD))
#get the code for the vector
colors_04120_MADdf$colors <- pal(nrow(colors_04120_MADdf))[colors_04120_MADdf$order]

#make separate vectors for each treatment
#order: "POL" "POH" "CGP" "FAE" "FAP"
color_04120_MADg_POL <- colors_04120_MADdf[c(1:length(m04120_names)),6]
color_04120_MADg_POH <- colors_04120_MADdf[c((length(m04120_names)+1):(length(m04120_names)*2)),6]
color_04120_MADg_CGP <- colors_04120_MADdf[c((length(m04120_names)*2+1):(length(m04120_names)*3)),6]
color_04120_MADg_FAE <- colors_04120_MADdf[c((length(m04120_names)*3+1):(length(m04120_names)*4)),6]
color_04120_MADg_FAP <- colors_04120_MADdf[c((length(m04120_names)*4+1):(length(m04120_names)*5)),6]

color_04120_MADp_POL <- colors_04120_MADdf[c((length(m04120_names)*5+1):(length(m04120_names)*6)),6]
color_04120_MADp_POH <- colors_04120_MADdf[c((length(m04120_names)*6+1):(length(m04120_names)*7)),6]
color_04120_MADp_CGP <- colors_04120_MADdf[c((length(m04120_names)*7+1):(length(m04120_names)*8)),6]
color_04120_MADp_FAE <- colors_04120_MADdf[c((length(m04120_names)*8+1):(length(m04120_names)*9)),6]
color_04120_MADp_FAP <- colors_04120_MADdf[c((length(m04120_names)*9+1):(length(m04120_names)*10)),6]

#assign color vector
#black color for the writing - needs to be the length of the names list
black_color <- rep("#000000", length(m04120_names))
#make pathway map
url04120_MADg_POL <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADg_POL, black_color)
url04120_MADg_POH <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADg_POH, black_color)
url04120_MADg_CGP <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADg_CGP, black_color)
url04120_MADg_FAE <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADg_FAE, black_color)
url04120_MADg_FAP <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADg_FAP, black_color)

url04120_MADp_POL <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADp_POL, black_color)
url04120_MADp_POH <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADp_POH, black_color)
url04120_MADp_CGP <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADp_CGP, black_color)
url04120_MADp_FAE <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADp_FAE, black_color)
url04120_MADp_FAP <- color.pathway.by.objects("path:map04120", m04120_KO_terms, color_04120_MADp_FAP, black_color)

#make a list and export it
m04120_urls <- data.frame(URL=c(url04120_MADg_POL, url04120_MADg_POH, url04120_MADg_CGP, url04120_MADg_FAE, url04120_MADg_FAP, 
                                url04120_MADp_POL, url04120_MADp_POH, url04120_MADp_CGP, url04120_MADp_FAE, url04120_MADp_FAP), 
                          Treatment=rep(treat_step, 2), 
                          Type=rep(c("gene", "protein"),each=5))
write.table(m04120_urls, "pathway_maps/m04120_colored_expression_urls.txt", quote=F, sep="\t", row.names=F)
write.table(colors_04120_MADdf, "m04120_MAD_values.txt", quote=F, sep="\t", row.names=F)

##04010----
#names vector is not the summed one, it is the original one because there are no repeats
#colors
colors_04010_MADdf <- data.frame(MAD=c(MAD_m04010_gPOL, MAD_m04010_gPOH, MAD_m04010_gCGP, MAD_m04010_gFAE, MAD_m04010_gFAP,
                                       MAD_m04010_pPOL, MAD_m04010_pPOH, MAD_m04010_pCGP, MAD_m04010_pFAE, MAD_m04010_pFAP),
                                 Treatment=rep(rep(treat_step, each=length(MAD_m04010_gPOH)),2), 
                                 Type=rep(c("gene", "protein"), each=(length(MAD_m04010_gPOH)*5)), 
                                 Name=rep(summed_names_04010, 10))
colors_04010_MADdf$order <- findInterval(colors_04010_MADdf$MAD, sort(colors_04010_MADdf$MAD))
#get the code for the vector
colors_04010_MADdf$colors <- pal(nrow(colors_04010_MADdf))[colors_04010_MADdf$order]

#make separate vectors for each treatment
#order: "POL" "POH" "CGP" "FAE" "FAP"
color_04010_MADg_POL <- colors_04010_MADdf[c(1:length(summed_names_04010)),6]
color_04010_MADg_POH <- colors_04010_MADdf[c((length(summed_names_04010)+1):(length(summed_names_04010)*2)),6]
color_04010_MADg_CGP <- colors_04010_MADdf[c((length(summed_names_04010)*2+1):(length(summed_names_04010)*3)),6]
color_04010_MADg_FAE <- colors_04010_MADdf[c((length(summed_names_04010)*3+1):(length(summed_names_04010)*4)),6]
color_04010_MADg_FAP <- colors_04010_MADdf[c((length(summed_names_04010)*4+1):(length(summed_names_04010)*5)),6]

color_04010_MADp_POL <- colors_04010_MADdf[c((length(summed_names_04010)*5+1):(length(summed_names_04010)*6)),6]
color_04010_MADp_POH <- colors_04010_MADdf[c((length(summed_names_04010)*6+1):(length(summed_names_04010)*7)),6]
color_04010_MADp_CGP <- colors_04010_MADdf[c((length(summed_names_04010)*7+1):(length(summed_names_04010)*8)),6]
color_04010_MADp_FAE <- colors_04010_MADdf[c((length(summed_names_04010)*8+1):(length(summed_names_04010)*9)),6]
color_04010_MADp_FAP <- colors_04010_MADdf[c((length(summed_names_04010)*9+1):(length(summed_names_04010)*10)),6]

#assign color vector
#black color for the writing - needs to be the length of the names list
black_color <- rep("#000000", length(summed_names_04010))
#make pathway map
url04010_MADg_POL <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADg_POL, black_color)
url04010_MADg_POH <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADg_POH, black_color)
url04010_MADg_CGP <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADg_CGP, black_color)
url04010_MADg_FAE <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADg_FAE, black_color)
url04010_MADg_FAP <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADg_FAP, black_color)

url04010_MADp_POL <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADp_POL, black_color)
url04010_MADp_POH <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADp_POH, black_color)
url04010_MADp_CGP <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADp_CGP, black_color)
url04010_MADp_FAE <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADp_FAE, black_color)
url04010_MADp_FAP <- color.pathway.by.objects("path:map04010", m04010_KO_terms, color_04010_MADp_FAP, black_color)

#make a list and export it
m04010_urls <- data.frame(URL=c(url04010_MADg_POL, url04010_MADg_POH, url04010_MADg_CGP, url04010_MADg_FAE, url04010_MADg_FAP, 
                                url04010_MADp_POL, url04010_MADp_POH, url04010_MADp_CGP, url04010_MADp_FAE, url04010_MADp_FAP), 
                          Treatment=rep(treat_step, 2), 
                          Type=rep(c("gene", "protein"),each=5))
write.table(m04010_urls, "pathway_maps/m04010_colored_expression_urls.txt", quote=F, sep="\t", row.names=F)
write.table(colors_04010_MADdf, "m04010_MAD_values.txt", quote=F, sep="\t", row.names=F)

##04210----
#colors
colors_04210_MADdf <- data.frame(MAD=c(MAD_m04210_gPOL, MAD_m04210_gPOH, MAD_m04210_gCGP, MAD_m04210_gFAE, MAD_m04210_gFAP,
                                       MAD_m04210_pPOL, MAD_m04210_pPOH, MAD_m04210_pCGP, MAD_m04210_pFAE, MAD_m04210_pFAP),
                                 Treatment=rep(rep(treat_step, each=length(MAD_m04210_gPOH)),2), 
                                 Type=rep(c("gene", "protein"), each=(length(MAD_m04210_gPOH)*5)), 
                                 Name=rep(summed_names_04210, 10))
colors_04210_MADdf$order <- findInterval(colors_04210_MADdf$MAD, sort(colors_04210_MADdf$MAD))
#get the code for the vector
colors_04210_MADdf$colors <- pal(nrow(colors_04210_MADdf))[colors_04210_MADdf$order]

#make separate vectors for each treatment
#order: "POL" "POH" "CGP" "FAE" "FAP"
color_04210_MADg_POL <- colors_04210_MADdf[c(1:length(summed_names_04210)),6]
color_04210_MADg_POH <- colors_04210_MADdf[c((length(summed_names_04210)+1):(length(summed_names_04210)*2)),6]
color_04210_MADg_CGP <- colors_04210_MADdf[c((length(summed_names_04210)*2+1):(length(summed_names_04210)*3)),6]
color_04210_MADg_FAE <- colors_04210_MADdf[c((length(summed_names_04210)*3+1):(length(summed_names_04210)*4)),6]
color_04210_MADg_FAP <- colors_04210_MADdf[c((length(summed_names_04210)*4+1):(length(summed_names_04210)*5)),6]

color_04210_MADp_POL <- colors_04210_MADdf[c((length(summed_names_04210)*5+1):(length(summed_names_04210)*6)),6]
color_04210_MADp_POH <- colors_04210_MADdf[c((length(summed_names_04210)*6+1):(length(summed_names_04210)*7)),6]
color_04210_MADp_CGP <- colors_04210_MADdf[c((length(summed_names_04210)*7+1):(length(summed_names_04210)*8)),6]
color_04210_MADp_FAE <- colors_04210_MADdf[c((length(summed_names_04210)*8+1):(length(summed_names_04210)*9)),6]
color_04210_MADp_FAP <- colors_04210_MADdf[c((length(summed_names_04210)*9+1):(length(summed_names_04210)*10)),6]

#assign color vector
#black color for the writing - needs to be the length of the names list
black_color <- rep("#000000", length(summed_names_04210))
#make pathway map
url04210_MADg_POL <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADg_POL, black_color)
url04210_MADg_POH <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADg_POH, black_color)
url04210_MADg_CGP <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADg_CGP, black_color)
url04210_MADg_FAE <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADg_FAE, black_color)
url04210_MADg_FAP <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADg_FAP, black_color)

url04210_MADp_POL <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADp_POL, black_color)
url04210_MADp_POH <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADp_POH, black_color)
url04210_MADp_CGP <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADp_CGP, black_color)
url04210_MADp_FAE <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADp_FAE, black_color)
url04210_MADp_FAP <- color.pathway.by.objects("path:map04210", m04210_KO_terms, color_04210_MADp_FAP, black_color)

#make a list and export it
m04210_urls <- data.frame(URL=c(url04210_MADg_POL, url04210_MADg_POH, url04210_MADg_CGP, url04210_MADg_FAE, url04210_MADg_FAP, 
                                url04210_MADp_POL, url04210_MADp_POH, url04210_MADp_CGP, url04210_MADp_FAE, url04210_MADp_FAP), 
                          Treatment=rep(treat_step, 2), 
                          Type=rep(c("gene", "protein"),each=5))
write.table(m04210_urls, "pathway_maps/m04210_colored_expression_urls.txt", quote=F, sep="\t", row.names=F)
write.table(colors_04210_MADdf, "m04210_MAD_values.txt", quote=F, sep="\t", row.names=F)

##00020----
#colors
colors_00020_MADdf <- data.frame(MAD=c(MAD_m00020_gPOL, MAD_m00020_gPOH, MAD_m00020_gCGP, MAD_m00020_gFAE, MAD_m00020_gFAP,
                                       MAD_m00020_pPOL, MAD_m00020_pPOH, MAD_m00020_pCGP, MAD_m00020_pFAE, MAD_m00020_pFAP),
                                 Treatment=rep(rep(treat_step, each=length(MAD_m00020_gPOH)),2), 
                                 Type=rep(c("gene", "protein"), each=(length(MAD_m00020_gPOH)*5)), 
                                 Name=rep(summed_names_00020, 10))
colors_00020_MADdf$order <- findInterval(colors_00020_MADdf$MAD, sort(colors_00020_MADdf$MAD))
#get the code for the vector
colors_00020_MADdf$colors <- pal(nrow(colors_00020_MADdf))[colors_00020_MADdf$order]

#make separate vectors for each treatment
#order: "POL" "POH" "CGP" "FAE" "FAP"
color_00020_MADg_POL <- colors_00020_MADdf[c(1:length(summed_names_00020)),6]
color_00020_MADg_POH <- colors_00020_MADdf[c((length(summed_names_00020)+1):(length(summed_names_00020)*2)),6]
color_00020_MADg_CGP <- colors_00020_MADdf[c((length(summed_names_00020)*2+1):(length(summed_names_00020)*3)),6]
color_00020_MADg_FAE <- colors_00020_MADdf[c((length(summed_names_00020)*3+1):(length(summed_names_00020)*4)),6]
color_00020_MADg_FAP <- colors_00020_MADdf[c((length(summed_names_00020)*4+1):(length(summed_names_00020)*5)),6]

color_00020_MADp_POL <- colors_00020_MADdf[c((length(summed_names_00020)*5+1):(length(summed_names_00020)*6)),6]
color_00020_MADp_POH <- colors_00020_MADdf[c((length(summed_names_00020)*6+1):(length(summed_names_00020)*7)),6]
color_00020_MADp_CGP <- colors_00020_MADdf[c((length(summed_names_00020)*7+1):(length(summed_names_00020)*8)),6]
color_00020_MADp_FAE <- colors_00020_MADdf[c((length(summed_names_00020)*8+1):(length(summed_names_00020)*9)),6]
color_00020_MADp_FAP <- colors_00020_MADdf[c((length(summed_names_00020)*9+1):(length(summed_names_00020)*10)),6]

#assign color vector
#black color for the writing - needs to be the length of the names list
black_color <- rep("#000000", length(summed_names_00020))
#make pathway map
url00020_MADg_POL <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADg_POL, black_color)
url00020_MADg_POH <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADg_POH, black_color)
url00020_MADg_CGP <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADg_CGP, black_color)
url00020_MADg_FAE <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADg_FAE, black_color)
url00020_MADg_FAP <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADg_FAP, black_color)

url00020_MADp_POL <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADp_POL, black_color)
url00020_MADp_POH <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADp_POH, black_color)
url00020_MADp_CGP <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADp_CGP, black_color)
url00020_MADp_FAE <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADp_FAE, black_color)
url00020_MADp_FAP <- color.pathway.by.objects("path:map00020", m00020_KO_terms, color_00020_MADp_FAP, black_color)

#make a list and export it
m00020_urls <- data.frame(URL=c(url00020_MADg_POL, url00020_MADg_POH, url00020_MADg_CGP, url00020_MADg_FAE, url00020_MADg_FAP, 
                                url00020_MADp_POL, url00020_MADp_POH, url00020_MADp_CGP, url00020_MADp_FAE, url00020_MADp_FAP), 
                          Treatment=rep(treat_step, 2), 
                          Type=rep(c("gene", "protein"),each=5))
write.table(m00020_urls, "pathway_maps/m00020_colored_expression_urls.txt", quote=F, sep="\t", row.names=F)
write.table(colors_00020_MADdf, "m00020_MAD_values.txt", quote=F, sep="\t", row.names=F)


####multiple comparisons p values####
#use Brown-Forsythe test
#comparisons happen within each gene
#export p values

residuals_function <- function(data_input) {
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

#from MAD p values from residuals file
t.test_treatmentGP_function <- function(gene_treat_res, protein_treat_res) { #input residual matrices for each treatment
  treat_comp_GP <- vector()
  for (i in 1:(length(gene_treat_res[1,]))) {
    gene_sub <- gene_treat_res[,i]
    protein_sub <- protein_treat_res[,i]
    p_val <- t.test(gene_sub, protein_sub)$p.value
    stat <- t.test(gene_sub, protein_sub)$statistic
    effect_test <- cohen.d(gene_sub,protein_sub, hedges.correction = TRUE)
    effect_size <- effect_test$estimate
    add <- c(p_val, stat, effect_size, colnames(gene_treat_res[i]))
    treat_comp_GP <- rbind(treat_comp_GP, add)
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

anova_treatment_function <- function(treat_res_CGP, treat_res_FAE, treat_res_FAP, treat_res_POH, treat_res_POL) { #input residual matrices for each treatment
  treat_comp <- vector()
  for (i in 1:(length(treat_res_CGP[1,]))) {
    CGP_sub <- treat_res_CGP[,i]
    FAE_sub <- treat_res_FAE[,i]
    FAP_sub <- treat_res_FAP[,i]
    POH_sub <- treat_res_POH[,i]
    POL_sub <- treat_res_POL[,i]
    
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

t.test_genes_function <- function(treat_res) { #input residual matrices for each treatment
  gene_comp <- vector()
  for (i in 1:(length(treat_res[1,])-1)) {
    counter <- i+1
    gene1_sub <- treat_res[,i]
    for (j in counter:(length(treat_res[1,]))) {
      gene2_sub <- treat_res[,j]
      p_val <- t.test(gene1_sub, gene2_sub)$p.value
      stat <- t.test(gene1_sub, gene2_sub)$statistic
      effect_test <- cohen.d(gene1_sub, gene2_sub, hedges.correction = TRUE)
      effect_size <- effect_test$estimate
      add <- c(p_val, stat, effect_size, colnames(treat_res)[i], colnames(treat_res)[j])
    }
    gene_comp <- rbind(gene_comp, add)
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
  #replace NAs with 1 (non sig)
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

#00010 ----

residuals_m00010_gPOH <- residuals_function(m00010_gPOH)
residuals_m00010_gPOL <- residuals_function(m00010_gPOL)
residuals_m00010_gCGP <- residuals_function(m00010_gCGP)
residuals_m00010_gFAE <- residuals_function(m00010_gFAE)
residuals_m00010_gFAP <- residuals_function(m00010_gFAP)
residuals_m00010_pPOH <- residuals_function(m00010_pPOH)
residuals_m00010_pPOL <- residuals_function(m00010_pPOL)
residuals_m00010_pCGP <- residuals_function(m00010_pCGP)
residuals_m00010_pFAE <- residuals_function(m00010_pFAE)
residuals_m00010_pFAP <- residuals_function(m00010_pFAP)

residuals_m00010 <- rbind(residuals_m00010_gCGP, residuals_m00010_gFAE, residuals_m00010_gFAP, residuals_m00010_gPOH, residuals_m00010_gPOL,
                          residuals_m00010_pCGP, residuals_m00010_pFAE, residuals_m00010_pFAP, residuals_m00010_pPOH, residuals_m00010_pPOL)
residuals_m00010$treatment <- rep(c(rep("CGP", 5), rep("FAE", 9), rep("FAP", 9), rep("POH", 10), rep("POL", 8)),2)
residuals_m00010$type <- as.factor(rep(c("gene", "protein"), each=41))
residuals_m00010$treatment_type <- as.factor(paste(residuals_m00010$treatment, "_", residuals_m00010$type, sep=''))

#one significant comparison
m00010_t.comps_POL <- t.test_treatmentGP_function(residuals_m00010_gPOL, residuals_m00010_pPOL)
m00010_t.comps_POH <- t.test_treatmentGP_function(residuals_m00010_gPOH, residuals_m00010_pPOH) #NOT EMPTY
m00010_t.comps_CGP <- t.test_treatmentGP_function(residuals_m00010_gCGP, residuals_m00010_pCGP)
m00010_t.comps_FAE <- t.test_treatmentGP_function(residuals_m00010_gFAE, residuals_m00010_pFAE)
m00010_t.comps_FAP <- t.test_treatmentGP_function(residuals_m00010_gFAP, residuals_m00010_pFAP)

#compares treatments within data type
m00010_aov.comps_gene <- anova_treatment_function(residuals_m00010_gCGP, residuals_m00010_gFAE, residuals_m00010_gFAP, residuals_m00010_gPOH, residuals_m00010_gPOL) #NOT EMPTY 6.8
m00010_aov.comps_protein <- anova_treatment_function(residuals_m00010_pCGP, residuals_m00010_pFAE, residuals_m00010_pFAP, residuals_m00010_pPOH, residuals_m00010_pPOL) #NOT EMPTY 6.8

#compares across genes within treatment and data type
m00010_gene.comps_gCGP <- t.test_genes_function(residuals_m00010_gCGP)
m00010_gene.comps_gFAE <- t.test_genes_function(residuals_m00010_gFAE)
m00010_gene.comps_gFAP <- t.test_genes_function(residuals_m00010_gFAP)
m00010_gene.comps_gPOH <- t.test_genes_function(residuals_m00010_gPOH)
m00010_gene.comps_gPOL <- t.test_genes_function(residuals_m00010_gPOL)
m00010_gene.comps_pCGP <- t.test_genes_function(residuals_m00010_pCGP)
m00010_gene.comps_pFAE <- t.test_genes_function(residuals_m00010_pFAE)
m00010_gene.comps_pFAP <- t.test_genes_function(residuals_m00010_pFAP)
m00010_gene.comps_pPOH <- t.test_genes_function(residuals_m00010_pPOH)
m00010_gene.comps_pPOL <- t.test_genes_function(residuals_m00010_pPOL)


# #export significant treatment comparisons - m00010_sigcomp gets incorporated at the final export
m00010_sigcomp_gene <- m00010_aov.comps_gene[[3]]
m00010_sigcomp_gene$type <- "gene"
m00010_sigcomp_protein <- m00010_aov.comps_protein[[3]]
m00010_sigcomp_protein$type <- "protein"
m00010_sigcomp <- rbind(m00010_sigcomp_gene, m00010_sigcomp_protein)
colnames(m00010_sigcomp) <- c("p_val", "effect_size", "treat1", "treat2", "gene", "type")
rownames(m00010_sigcomp) <- NULL
m00010_sigcomp$pathway <- "m00010"
m00010_POH_sigcomp <- m00010_t.comps_POH[[3]]
colnames(m00010_POH_sigcomp) <- c("p_val", "t_val", "effect_size", "transcript")
rownames(m00010_POH_sigcomp) <- NULL
m00010_POH_sigcomp$treatment <- "POH"
m00010_POH_sigcomp$pathway <- "m00010"

#03050 ----

residuals_m03050_gPOH <- residuals_function(m03050_gPOH)
residuals_m03050_gPOL <- residuals_function(m03050_gPOL)
residuals_m03050_gCGP <- residuals_function(m03050_gCGP)
residuals_m03050_gFAE <- residuals_function(m03050_gFAE)
residuals_m03050_gFAP <- residuals_function(m03050_gFAP)
residuals_m03050_pPOH <- residuals_function(m03050_pPOH)
residuals_m03050_pPOL <- residuals_function(m03050_pPOL)
residuals_m03050_pCGP <- residuals_function(m03050_pCGP)
residuals_m03050_pFAE <- residuals_function(m03050_pFAE)
residuals_m03050_pFAP <- residuals_function(m03050_pFAP)

residuals_m03050 <- rbind(residuals_m03050_gCGP, residuals_m03050_gFAE, residuals_m03050_gFAP, residuals_m03050_gPOH, residuals_m03050_gPOL,
                          residuals_m03050_pCGP, residuals_m03050_pFAE, residuals_m03050_pFAP, residuals_m03050_pPOH, residuals_m03050_pPOL)
residuals_m03050$treatment <- rep(c(rep("CGP", 5), rep("FAE", 9), rep("FAP", 9), rep("POH", 10), rep("POL", 8)),2)
residuals_m03050$type <- as.factor(rep(c("gene", "protein"), each=41))
residuals_m03050$treatment_type <- as.factor(paste(residuals_m03050$treatment, "_", residuals_m03050$type, sep=''))

#no significant comparisons
m03050_t.comps_POL <- t.test_treatmentGP_function(residuals_m03050_gPOL, residuals_m03050_pPOL)
m03050_t.comps_POH <- t.test_treatmentGP_function(residuals_m03050_gPOH, residuals_m03050_pPOH)
m03050_t.comps_CGP <- t.test_treatmentGP_function(residuals_m03050_gCGP, residuals_m03050_pCGP)
m03050_t.comps_FAE <- t.test_treatmentGP_function(residuals_m03050_gFAE, residuals_m03050_pFAE)
m03050_t.comps_FAP <- t.test_treatmentGP_function(residuals_m03050_gFAP, residuals_m03050_pFAP)

#compares treatments within data type
m03050_aov.comps_gene <- anova_treatment_function(residuals_m03050_gCGP, residuals_m03050_gFAE, residuals_m03050_gFAP, residuals_m03050_gPOH, residuals_m03050_gPOL) #empty
m03050_aov.comps_protein <- anova_treatment_function(residuals_m03050_pCGP, residuals_m03050_pFAE, residuals_m03050_pFAP, residuals_m03050_pPOH, residuals_m03050_pPOL) #empty

#compares across genes within treatment and data type
m03050_gene.comps_gCGP <- t.test_genes_function(residuals_m03050_gCGP)
m03050_gene.comps_gFAE <- t.test_genes_function(residuals_m03050_gFAE)
m03050_gene.comps_gFAP <- t.test_genes_function(residuals_m03050_gFAP)
m03050_gene.comps_gPOH <- t.test_genes_function(residuals_m03050_gPOH)
m03050_gene.comps_gPOL <- t.test_genes_function(residuals_m03050_gPOL)
m03050_gene.comps_pCGP <- t.test_genes_function(residuals_m03050_pCGP)
m03050_gene.comps_pFAE <- t.test_genes_function(residuals_m03050_pFAE)
m03050_gene.comps_pFAP <- t.test_genes_function(residuals_m03050_pFAP) 
m03050_gene.comps_pPOH <- t.test_genes_function(residuals_m03050_pPOH) 
m03050_gene.comps_pPOL <- t.test_genes_function(residuals_m03050_pPOL)

#NO EXPORT

#04120 ----

residuals_m04120_gPOH <- residuals_function(m04120_gPOH)
residuals_m04120_gPOL <- residuals_function(m04120_gPOL)
residuals_m04120_gCGP <- residuals_function(m04120_gCGP)
residuals_m04120_gFAE <- residuals_function(m04120_gFAE)
residuals_m04120_gFAP <- residuals_function(m04120_gFAP)
residuals_m04120_pPOH <- residuals_function(m04120_pPOH)
residuals_m04120_pPOL <- residuals_function(m04120_pPOL)
residuals_m04120_pCGP <- residuals_function(m04120_pCGP)
residuals_m04120_pFAE <- residuals_function(m04120_pFAE)
residuals_m04120_pFAP <- residuals_function(m04120_pFAP)

residuals_m04120 <- rbind(residuals_m04120_gCGP, residuals_m04120_gFAE, residuals_m04120_gFAP, residuals_m04120_gPOH, residuals_m04120_gPOL,
                          residuals_m04120_pCGP, residuals_m04120_pFAE, residuals_m04120_pFAP, residuals_m04120_pPOH, residuals_m04120_pPOL)
residuals_m04120$treatment <- rep(c(rep("CGP", 5), rep("FAE", 9), rep("FAP", 9), rep("POH", 10), rep("POL", 8)),2)
residuals_m04120$type <- as.factor(rep(c("gene", "protein"), each=41))
residuals_m04120$treatment_type <- as.factor(paste(residuals_m04120$treatment, "_", residuals_m04120$type, sep=''))

#12.12.19: POL and FAP have sig ones
m04120_t.comps_POL <- t.test_treatmentGP_function(residuals_m04120_gPOL, residuals_m04120_pPOL) 
m04120_t.comps_POH <- t.test_treatmentGP_function(residuals_m04120_gPOH, residuals_m04120_pPOH)
m04120_t.comps_CGP <- t.test_treatmentGP_function(residuals_m04120_gCGP, residuals_m04120_pCGP) 
m04120_t.comps_FAE <- t.test_treatmentGP_function(residuals_m04120_gFAE, residuals_m04120_pFAE)
m04120_t.comps_FAP <- t.test_treatmentGP_function(residuals_m04120_gFAP, residuals_m04120_pFAP) 

#compares treatments within data type
m04120_aov.comps_gene <- anova_treatment_function(residuals_m04120_gCGP, residuals_m04120_gFAE, residuals_m04120_gFAP, residuals_m04120_gPOH, residuals_m04120_gPOL) #NOT EMPTY 6.8
m04120_aov.comps_protein <- anova_treatment_function(residuals_m04120_pCGP, residuals_m04120_pFAE, residuals_m04120_pFAP, residuals_m04120_pPOH, residuals_m04120_pPOL)

#compares across genes within treatment and data type
m04120_gene.comps_gCGP <- t.test_genes_function(residuals_m04120_gCGP)
m04120_gene.comps_gFAE <- t.test_genes_function(residuals_m04120_gFAE)
m04120_gene.comps_gFAP <- t.test_genes_function(residuals_m04120_gFAP)
m04120_gene.comps_gPOH <- t.test_genes_function(residuals_m04120_gPOH)
m04120_gene.comps_gPOL <- t.test_genes_function(residuals_m04120_gPOL)
m04120_gene.comps_pCGP <- t.test_genes_function(residuals_m04120_pCGP)
m04120_gene.comps_pFAE <- t.test_genes_function(residuals_m04120_pFAE)
m04120_gene.comps_pFAP <- t.test_genes_function(residuals_m04120_pFAP) 
m04120_gene.comps_pPOH <- t.test_genes_function(residuals_m04120_pPOH)
m04120_gene.comps_pPOL <- t.test_genes_function(residuals_m04120_pPOL) 

#export significant treatment comparisons - m04120_sigcomp gets incorporated at the final export
m04120_sigcomp_gene <- m04120_aov.comps_gene[[3]]
m04120_sigcomp_gene$type <- "gene"
m04120_sigcomp <- m04120_sigcomp_gene
colnames(m04120_sigcomp) <- c("p_val", "effect_size", "treat1", "treat2", "gene", "type")
rownames(m04120_sigcomp) <- NULL
m04120_sigcomp$pathway <- "m04120"

#04010 ----

residuals_m04010_gPOH <- residuals_function(m04010_gPOH)
residuals_m04010_gPOL <- residuals_function(m04010_gPOL)
residuals_m04010_gCGP <- residuals_function(m04010_gCGP)
residuals_m04010_gFAE <- residuals_function(m04010_gFAE)
residuals_m04010_gFAP <- residuals_function(m04010_gFAP)
residuals_m04010_pPOH <- residuals_function(m04010_pPOH)
residuals_m04010_pPOL <- residuals_function(m04010_pPOL)
residuals_m04010_pCGP <- residuals_function(m04010_pCGP)
residuals_m04010_pFAE <- residuals_function(m04010_pFAE)
residuals_m04010_pFAP <- residuals_function(m04010_pFAP)

residuals_m04010 <- rbind(residuals_m04010_gCGP, residuals_m04010_gFAE, residuals_m04010_gFAP, residuals_m04010_gPOH, residuals_m04010_gPOL,
                          residuals_m04010_pCGP, residuals_m04010_pFAE, residuals_m04010_pFAP, residuals_m04010_pPOH, residuals_m04010_pPOL)
residuals_m04010$treatment <- rep(c(rep("CGP", 5), rep("FAE", 9), rep("FAP", 9), rep("POH", 10), rep("POL", 8)),2)
residuals_m04010$type <- as.factor(rep(c("gene", "protein"), each=41))
residuals_m04010$treatment_type <- as.factor(paste(residuals_m04010$treatment, "_", residuals_m04010$type, sep=''))

#only POH has sig ones
m04010_t.comps_POL <- t.test_treatmentGP_function(residuals_m04010_gPOL, residuals_m04010_pPOL)
m04010_t.comps_POH <- t.test_treatmentGP_function(residuals_m04010_gPOH, residuals_m04010_pPOH)
m04010_t.comps_CGP <- t.test_treatmentGP_function(residuals_m04010_gCGP, residuals_m04010_pCGP)
m04010_t.comps_FAE <- t.test_treatmentGP_function(residuals_m04010_gFAE, residuals_m04010_pFAE) #NOT EMPTY
m04010_t.comps_FAP <- t.test_treatmentGP_function(residuals_m04010_gFAP, residuals_m04010_pFAP) #NOT EMPTY

#compares treatments within data type
m04010_aov.comps_gene <- anova_treatment_function(residuals_m04010_gCGP, residuals_m04010_gFAE, residuals_m04010_gFAP, residuals_m04010_gPOH, residuals_m04010_gPOL) #NOT EMPTY 6.8
m04010_aov.comps_protein <- anova_treatment_function(residuals_m04010_pCGP, residuals_m04010_pFAE, residuals_m04010_pFAP, residuals_m04010_pPOH, residuals_m04010_pPOL)

#compares across genes within treatment and data type
m04010_gene.comps_gCGP <- t.test_genes_function(residuals_m04010_gCGP)
m04010_gene.comps_gFAE <- t.test_genes_function(residuals_m04010_gFAE)
m04010_gene.comps_gFAP <- t.test_genes_function(residuals_m04010_gFAP)
m04010_gene.comps_gPOH <- t.test_genes_function(residuals_m04010_gPOH)
m04010_gene.comps_gPOL <- t.test_genes_function(residuals_m04010_gPOL) 
m04010_gene.comps_pCGP <- t.test_genes_function(residuals_m04010_pCGP)
m04010_gene.comps_pFAE <- t.test_genes_function(residuals_m04010_pFAE)
m04010_gene.comps_pFAP <- t.test_genes_function(residuals_m04010_pFAP)
m04010_gene.comps_pPOH <- t.test_genes_function(residuals_m04010_pPOH) 
m04010_gene.comps_pPOL <- t.test_genes_function(residuals_m04010_pPOL) 

#export significant treatment comparisons - m04010_sigcomp gets incorporated at the final export 2.19.20
m04010_sigcomp_gene <- m04010_aov.comps_gene[[3]]
m04010_sigcomp_gene$type <- "gene"
m04010_sigcomp <- m04010_sigcomp_gene
colnames(m04010_sigcomp) <- c("p_val", "effect_size", "treat1", "treat2", "gene", "type")
rownames(m04010_sigcomp) <- NULL
m04010_sigcomp$pathway <- "m04010"
#do we want to get rid of the intercepts? leave for now

#export FAE significant p vals file
m04010_FAE_sigcomp <- m04010_t.comps_FAE[[3]]
colnames(m04010_FAE_sigcomp) <- c("p_val", "t_val", "effect_size", "transcript")
rownames(m04010_FAE_sigcomp) <- NULL
m04010_FAE_sigcomp$treatment <- "FAE"
m04010_FAE_sigcomp$pathway <- "m04010"
m04010_FAP_sigcomp <- m04010_t.comps_FAP[[3]]
colnames(m04010_FAP_sigcomp) <- c("p_val", "t_val", "effect_size", "transcript")
rownames(m04010_FAP_sigcomp) <- NULL
m04010_FAP_sigcomp$treatment <- "FAP"
m04010_FAP_sigcomp$pathway <- "m04010"

m04010_all_sigcomp <- rbind(m04010_FAE_sigcomp, m04010_FAP_sigcomp)

#04210 ----

residuals_m04210_gPOH <- residuals_function(m04210_gPOH)
residuals_m04210_gPOL <- residuals_function(m04210_gPOL)
residuals_m04210_gCGP <- residuals_function(m04210_gCGP)
residuals_m04210_gFAE <- residuals_function(m04210_gFAE)
residuals_m04210_gFAP <- residuals_function(m04210_gFAP)
residuals_m04210_pPOH <- residuals_function(m04210_pPOH)
residuals_m04210_pPOL <- residuals_function(m04210_pPOL)
residuals_m04210_pCGP <- residuals_function(m04210_pCGP)
residuals_m04210_pFAE <- residuals_function(m04210_pFAE)
residuals_m04210_pFAP <- residuals_function(m04210_pFAP)

residuals_m04210 <- rbind(residuals_m04210_gCGP, residuals_m04210_gFAE, residuals_m04210_gFAP, residuals_m04210_gPOH, residuals_m04210_gPOL,
                          residuals_m04210_pCGP, residuals_m04210_pFAE, residuals_m04210_pFAP, residuals_m04210_pPOH, residuals_m04210_pPOL)
residuals_m04210$treatment <- rep(c(rep("CGP", 5), rep("FAE", 9), rep("FAP", 9), rep("POH", 10), rep("POL", 8)),2)
residuals_m04210$type <- as.factor(rep(c("gene", "protein"), each=41))
residuals_m04210$treatment_type <- as.factor(paste(residuals_m04210$treatment, "_", residuals_m04210$type, sep=''))

#only POH has sig ones
m04210_t.comps_POL <- t.test_treatmentGP_function(residuals_m04210_gPOL, residuals_m04210_pPOL)
m04210_t.comps_POH <- t.test_treatmentGP_function(residuals_m04210_gPOH, residuals_m04210_pPOH)
m04210_t.comps_CGP <- t.test_treatmentGP_function(residuals_m04210_gCGP, residuals_m04210_pCGP)
m04210_t.comps_FAE <- t.test_treatmentGP_function(residuals_m04210_gFAE, residuals_m04210_pFAE)
m04210_t.comps_FAP <- t.test_treatmentGP_function(residuals_m04210_gFAP, residuals_m04210_pFAP)

#compares treatments within data type
m04210_aov.comps_gene <- anova_treatment_function(residuals_m04210_gCGP, residuals_m04210_gFAE, residuals_m04210_gFAP, residuals_m04210_gPOH, residuals_m04210_gPOL) #not empty 6.8
m04210_aov.comps_protein <- anova_treatment_function(residuals_m04210_pCGP, residuals_m04210_pFAE, residuals_m04210_pFAP, residuals_m04210_pPOH, residuals_m04210_pPOL)

#compares across genes within treatment and data type
m04210_gene.comps_gCGP <- t.test_genes_function(residuals_m04210_gCGP)
m04210_gene.comps_gFAE <- t.test_genes_function(residuals_m04210_gFAE)
m04210_gene.comps_gFAP <- t.test_genes_function(residuals_m04210_gFAP)
m04210_gene.comps_gPOH <- t.test_genes_function(residuals_m04210_gPOH)
m04210_gene.comps_gPOL <- t.test_genes_function(residuals_m04210_gPOL)
m04210_gene.comps_pCGP <- t.test_genes_function(residuals_m04210_pCGP)
m04210_gene.comps_pFAE <- t.test_genes_function(residuals_m04210_pFAE)
m04210_gene.comps_pFAP <- t.test_genes_function(residuals_m04210_pFAP)
m04210_gene.comps_pPOH <- t.test_genes_function(residuals_m04210_pPOH)
m04210_gene.comps_pPOL <- t.test_genes_function(residuals_m04210_pPOL)

# #export significant treatment comparisons - m04210_sigcomp gets incorporated at the final export
m04210_sigcomp_gene <- m04210_aov.comps_gene[[3]]
m04210_sigcomp_gene$type <- "gene"
m04210_sigcomp <- m04210_sigcomp_gene
colnames(m04210_sigcomp) <- c("p_val", "effect_size", "treat1", "treat2", "gene", "type")
rownames(m04210_sigcomp) <- NULL
m04210_sigcomp$pathway <- "m04210"

#00020 ----

residuals_m00020_gPOH <- residuals_function(m00020_gPOH)
residuals_m00020_gPOL <- residuals_function(m00020_gPOL)
residuals_m00020_gCGP <- residuals_function(m00020_gCGP)
residuals_m00020_gFAE <- residuals_function(m00020_gFAE)
residuals_m00020_gFAP <- residuals_function(m00020_gFAP)
residuals_m00020_pPOH <- residuals_function(m00020_pPOH)
residuals_m00020_pPOL <- residuals_function(m00020_pPOL)
residuals_m00020_pCGP <- residuals_function(m00020_pCGP)
residuals_m00020_pFAE <- residuals_function(m00020_pFAE)
residuals_m00020_pFAP <- residuals_function(m00020_pFAP)

residuals_m00020 <- rbind(residuals_m00020_gCGP, residuals_m00020_gFAE, residuals_m00020_gFAP, residuals_m00020_gPOH, residuals_m00020_gPOL,
                          residuals_m00020_pCGP, residuals_m00020_pFAE, residuals_m00020_pFAP, residuals_m00020_pPOH, residuals_m00020_pPOL)
residuals_m00020$treatment <- rep(c(rep("CGP", 5), rep("FAE", 9), rep("FAP", 9), rep("POH", 10), rep("POL", 8)),2)
residuals_m00020$type <- as.factor(rep(c("gene", "protein"), each=41))
residuals_m00020$treatment_type <- as.factor(paste(residuals_m00020$treatment, "_", residuals_m00020$type, sep=''))


#only FAP and POH sig ones
m00020_t.comps_POL <- t.test_treatmentGP_function(residuals_m00020_gPOL, residuals_m00020_pPOL) #NOT EMPTY 6.8
m00020_t.comps_POH <- t.test_treatmentGP_function(residuals_m00020_gPOH, residuals_m00020_pPOH)
m00020_t.comps_CGP <- t.test_treatmentGP_function(residuals_m00020_gCGP, residuals_m00020_pCGP)
m00020_t.comps_FAE <- t.test_treatmentGP_function(residuals_m00020_gFAE, residuals_m00020_pFAE)
m00020_t.comps_FAP <- t.test_treatmentGP_function(residuals_m00020_gFAP, residuals_m00020_pFAP) #NOT EMPTY 6.8

#compares treatments within data type
m00020_aov.comps_gene <- anova_treatment_function(residuals_m00020_gCGP, residuals_m00020_gFAE, residuals_m00020_gFAP, residuals_m00020_gPOH, residuals_m00020_gPOL) #NOT EMPTY 6.8
m00020_aov.comps_protein <- anova_treatment_function(residuals_m00020_pCGP, residuals_m00020_pFAE, residuals_m00020_pFAP, residuals_m00020_pPOH, residuals_m00020_pPOL)

#compares across genes within treatment and data type
m00020_gene.comps_gCGP <- t.test_genes_function(residuals_m00020_gCGP)
m00020_gene.comps_gFAE <- t.test_genes_function(residuals_m00020_gFAE) 
m00020_gene.comps_gFAP <- t.test_genes_function(residuals_m00020_gFAP) 
m00020_gene.comps_gPOH <- t.test_genes_function(residuals_m00020_gPOH)
m00020_gene.comps_gPOL <- t.test_genes_function(residuals_m00020_gPOL)
m00020_gene.comps_pCGP <- t.test_genes_function(residuals_m00020_pCGP)
m00020_gene.comps_pFAE <- t.test_genes_function(residuals_m00020_pFAE) 
m00020_gene.comps_pFAP <- t.test_genes_function(residuals_m00020_pFAP)
m00020_gene.comps_pPOH <- t.test_genes_function(residuals_m00020_pPOH)
m00020_gene.comps_pPOL <- t.test_genes_function(residuals_m00020_pPOL)

#export significant treatment comparisons - m00020_sigcomp gets incorporated at the final export
m00020_sigcomp_gene <- m00020_aov.comps_gene[[3]]
m00020_sigcomp_gene$type <- "gene"
m00020_sigcomp <- m00020_sigcomp_gene
colnames(m00020_sigcomp) <- c("p_val", "effect_size", "treat1", "treat2", "gene", "type")
rownames(m00020_sigcomp) <- NULL
m00020_sigcomp$pathway <- "m00020"

#export FAP significant p vals file
m00020_FAP_sigcomp <- m00020_t.comps_FAP[[3]]
colnames(m00020_FAP_sigcomp) <- c("p_val", "t_val", "effect_size", "transcript")
rownames(m00020_FAP_sigcomp) <- NULL
m00020_FAP_sigcomp$treatment <- "FAP"
m00020_FAP_sigcomp$pathway <- "m00020"
m00020_POL_sigcomp <- m00020_t.comps_POL[[3]]
colnames(m00020_POL_sigcomp) <- c("p_val", "t_val", "effect_size", "transcript")
rownames(m00020_POL_sigcomp) <- NULL
m00020_POL_sigcomp$treatment <- "POL"
m00020_POL_sigcomp$pathway <- "m00020"

m00020_all_sigcomp <- rbind(m00020_FAP_sigcomp, m00020_POL_sigcomp)

####put all significant residual comparisons together and export####
#add in annotation
BLAST <- read.delim("~/Documents/R_Data/musselRNASeq/final_files/Full_Blast2GO_table.csv", sep=',')
colnames(BLAST)[2] <- "transcript"

all_sigcomp <- rbind(m00010_POH_sigcomp, m04010_all_sigcomp, m00020_all_sigcomp)
all_sigcompBLAST <- merge(all_sigcomp, BLAST, by = "transcript")
all_sigcompBLAST <- all_sigcompBLAST[,c(1:6, 8)]

write.table(all_sigcompBLAST, "~/Documents/R_Data/musselRNASeq/compare_KEGG_MCODE/KEGG_allPath_sigcomp_GPtreat_6.8.20.txt", sep='\t', quote=F, row.names = F, col.names = T)
# 
# #they all have something except 03050
all_treat_sigcomp <- rbind(m00010_sigcomp, m04120_sigcomp, m04010_sigcomp, m04210_sigcomp, m00020_sigcomp)
colnames(all_treat_sigcomp)[5] <- "transcript"
all_treat_sigcompBLAST <- merge(all_treat_sigcomp, BLAST, by = "transcript")
all_treat_sigcompBLAST <- all_treat_sigcompBLAST[,c(1:7, 9)]


write.table(all_treat_sigcompBLAST, "~/Documents/R_Data/musselRNASeq/compare_KEGG_MCODE/KEGG_allPath_sigcomp_treat_6.8.20.txt", sep='\t', quote=F, row.names = F, col.names = T)


####multivariate MAD####
#this is copied from MCODE output analysis
#original code is from "mytilus variation bio levels.R"

#recreate scaled datasets by pathway here
m00010_allGz <- rbind(m00010_gCGP, m00010_gFAE, m00010_gFAP, m00010_gPOH, m00010_gPOL)
m00020_allGz <- rbind(m00020_gCGP, m00020_gFAE, m00020_gFAP, m00020_gPOH, m00020_gPOL)
m03050_allGz <- rbind(m03050_gCGP, m03050_gFAE, m03050_gFAP, m03050_gPOH, m03050_gPOL)
m04010_allGz <- rbind(m04010_gCGP, m04010_gFAE, m04010_gFAP, m04010_gPOH, m04010_gPOL)
m04120_allGz <- rbind(m04120_gCGP, m04120_gFAE, m04120_gFAP, m04120_gPOH, m04120_gPOL)
m04210_allGz <- rbind(m04210_gCGP, m04210_gFAE, m04210_gFAP, m04210_gPOH, m04210_gPOL)

m00010_allPz <- rbind(m00010_pCGP, m00010_pFAE, m00010_pFAP, m00010_pPOH, m00010_pPOL)
m00020_allPz <- rbind(m00020_pCGP, m00020_pFAE, m00020_pFAP, m00020_pPOH, m00020_pPOL)
m03050_allPz <- rbind(m03050_pCGP, m03050_pFAE, m03050_pFAP, m03050_pPOH, m03050_pPOL)
m04010_allPz <- rbind(m04010_pCGP, m04010_pFAE, m04010_pFAP, m04010_pPOH, m04010_pPOL)
m04120_allPz <- rbind(m04120_pCGP, m04120_pFAE, m04120_pFAP, m04120_pPOH, m04120_pPOL)
m04210_allPz <- rbind(m04210_pCGP, m04210_pFAE, m04210_pFAP, m04210_pPOH, m04210_pPOL)

#input is a list of data frames - change from other script - gets put in already transposed and scaled
#use the scaled datasets that look like this: m00010_allGz
KEGG_listG <- setNames(replicate(6,data.frame()),c("m00010", "m00020", "m03050", "m04010", "m04120", "m04210"))
KEGG_listG[[1]] <- m00010_allGz
KEGG_listG[[2]] <- m00020_allGz
KEGG_listG[[3]] <- m03050_allGz
KEGG_listG[[4]] <- m04010_allGz
KEGG_listG[[5]] <- m04120_allGz
KEGG_listG[[6]] <- m04210_allGz

KEGG_listP <- setNames(replicate(6,data.frame()),c("m00010", "m00020", "m03050", "m04010", "m04120", "m04210"))
KEGG_listP[[1]] <- m00010_allPz
KEGG_listP[[2]] <- m00020_allPz
KEGG_listP[[3]] <- m03050_allPz
KEGG_listP[[4]] <- m04010_allPz
KEGG_listP[[5]] <- m04120_allPz
KEGG_listP[[6]] <- m04210_allPz


library(ggplot2)
library(stringr)
library(plotly)
library(vegan)
library(compositions)
#library(pcaMethods)
library(multcompView)
library(colorspace)

#calculate multivariate dispersion per cluster
#use a matrix of expression values for each cluster - there are data frames in "KEGG_listG" and "KEGG_listP"
#the order is CGP (5), FAE (9), FAP (9), POH (10), and POL (8)
treatment_vector <- c(rep("CGP", 5), rep("FAE", 9), rep("FAP", 9), rep("POH", 10), rep("POL", 8))

multivar_disp_function <- function(expr_dataframe) {
  mat <- expr_dataframe
  dis <- dist(mat, method="euclidean")
  mod <- betadisper(dis, treatment_vector, type="median")
  #this gives MAD values
  #distances are from the centroid, so this is the distance from the median in multivariate space, you take the median of these absolute values to get MAD
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
  #adjust for multiple comparisons - EDIT 2.19.20
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
geneKEGG_MVD_MAD <- data.frame()
geneKEGG_MVD_multcompP <- vector()
geneKEGGMVD_effSize <- data.frame()
geneKEGGMVD_Pval <- vector()
for (i in 1:length(KEGG_listG)) {
  MVD <- multivar_disp_function(KEGG_listG[[i]])
  geneKEGG_MVD_MAD <- rbind(geneKEGG_MVD_MAD, MVD[[1]])
  geneKEGG_MVD_multcompP <- c(geneKEGG_MVD_multcompP, MVD[[2]])
  geneKEGGMVD_effSize <- rbind(geneKEGGMVD_effSize, MVD[[3]])
  geneKEGGMVD_Pval <- c(geneKEGGMVD_Pval, MVD[[4]])
}

#put in the cluster and data type labels
geneKEGG_MVD_MAD$pathway <- rep(rep(c("m00010", "m00020", "m03050", "m04010", "m04120", "m04210"), each=5))
geneKEGG_MVD_MAD$treatment <- rep(c("CGP", "FAE", "FAP", "POH", "POL"), 6)
geneKEGG_MVD_pVal <- data.frame(comparison=names(geneKEGG_MVD_multcompP), p_val=unname(geneKEGG_MVD_multcompP), 
                                pathway=rep(rep(c("m00010", "m00020", "m03050", "m04010", "m04120", "m04210"), each=10)))
geneKEGG_MVD_pVal_noAdj <- data.frame(comparison=names(geneKEGG_MVD_multcompP), p_val_adj=unname(geneKEGG_MVD_multcompP), p_val=unname(geneKEGGMVD_Pval),
                                pathway=rep(rep(c("m00010", "m00020", "m03050", "m04010", "m04120", "m04210"), each=10)))
#reduce to only significant (before adjustment)
geneKEGG_MVD_pVal_noAdj <- geneKEGG_MVD_pVal_noAdj[(geneKEGG_MVD_pVal_noAdj$p_val<0.05),]

#export
write.table(geneKEGG_MVD_MAD, "~/geneKEGG_MVD_MAD_letters.txt", sep='\t', row.names=F, quote=F)
write.table(geneKEGG_MVD_pVal, "~/geneKEGG_MVD_multipleComp_pVal.txt", sep='\t', row.names=F, quote=F)
write.table(geneKEGG_MVD_pVal_noAdj, "~/geneKEGG_MVD_multipleComp_pValOriginal.txt", sep='\t', row.names=F, quote=F)

#make giant data frames of this info
proteinKEGG_MVD_MAD <- data.frame()
proteinKEGG_MVD_multcompP <- vector()
proteinKEGGMVD_effSize <- data.frame()
proteinKEGGMVD_Pval <- vector()
for (i in 1:length(KEGG_listP)) {
  MVD <- multivar_disp_function(KEGG_listP[[i]])
  proteinKEGG_MVD_MAD <- rbind(proteinKEGG_MVD_MAD, MVD[[1]])
  proteinKEGG_MVD_multcompP <- c(proteinKEGG_MVD_multcompP, MVD[[2]])
  proteinKEGGMVD_effSize <- rbind(proteinKEGGMVD_effSize, MVD[[3]])
  proteinKEGGMVD_Pval <- c(proteinKEGGMVD_Pval, MVD[[4]])
}

#put in the cluster and data type labels
proteinKEGG_MVD_MAD$pathway <- rep(rep(c("m00010", "m00020", "m03050", "m04010", "m04120", "m04210"), each=5))
proteinKEGG_MVD_MAD$treatment <- rep(c("CGP", "FAE", "FAP", "POH", "POL"), 6)
proteinKEGG_MVD_pVal <- data.frame(comparison=names(proteinKEGG_MVD_multcompP), p_val=unname(proteinKEGG_MVD_multcompP), 
                                pathway=rep(rep(c("m00010", "m00020", "m03050", "m04010", "m04120", "m04210"), each=10)))
proteinKEGG_MVD_pVal_noAdj <- data.frame(comparison=names(proteinKEGG_MVD_multcompP), p_val_adj=unname(proteinKEGG_MVD_multcompP), p_val=unname(proteinKEGGMVD_Pval),
                                      pathway=rep(rep(c("m00010", "m00020", "m03050", "m04010", "m04120", "m04210"), each=10)))
#reduce to only significant (before adjustment)
proteinKEGG_MVD_pVal_noAdj <- proteinKEGG_MVD_pVal_noAdj[(proteinKEGG_MVD_pVal_noAdj$p_val<0.05),]


#export 
write.table(proteinKEGG_MVD_MAD, "~/proteinKEGG_MVD_MAD_letters.txt", sep='\t', row.names=F, quote=F)
write.table(proteinKEGG_MVD_pVal, "~/proteinKEGG_MVD_multipleComp_pVal.txt", sep='\t', row.names=F, quote=F)

#column 2 is regular p vals and column 5 is adjusted p vals
#pull out only significant pairwise comparisons
geneKEGG_MVD_sigPVal <- data.frame()
for (i in 1:length(geneKEGG_MVD_pVal[,1])) {
  if (geneKEGG_MVD_pVal[i,2] < 0.05) {
    geneKEGG_MVD_sigPVal <- rbind(geneKEGG_MVD_sigPVal, geneKEGG_MVD_pVal[i,])
  }
}

proteinKEGG_MVD_sigPVal <- data.frame()
for (i in 1:length(proteinKEGG_MVD_pVal[,1])) {
  if (proteinKEGG_MVD_pVal[i,2] < 0.05) {
    proteinKEGG_MVD_sigPVal <- rbind(proteinKEGG_MVD_sigPVal, proteinKEGG_MVD_pVal[i,])
  }
}

#export 
write.table(geneKEGG_MVD_sigPVal, "~/geneKEGG_MVD_multipleComp_sig_pVal.txt", sep='\t', row.names=F, quote=F)

####plots for MAD within pathway between gene and protein####
#look at multivariate
#use geneKEGG_MVD_MAD and proteinKEGG_MVD_MAD
MVD_avgGene_00010 <- mean(geneKEGG_MVD_MAD[c(1:5),1])
MVD_avgGene_00020 <- mean(geneKEGG_MVD_MAD[c(6:10),1])
MVD_avgGene_03050 <- mean(geneKEGG_MVD_MAD[c(11:15),1])
MVD_avgGene_04010 <- mean(geneKEGG_MVD_MAD[c(16:20),1])
MVD_avgGene_04120 <- mean(geneKEGG_MVD_MAD[c(21:25),1])
MVD_avgGene_04210 <- mean(geneKEGG_MVD_MAD[c(26:30),1])
MVD_avgProtein_00010 <- mean(proteinKEGG_MVD_MAD[c(1:5),1])
MVD_avgProtein_00020 <- mean(proteinKEGG_MVD_MAD[c(6:10),1])
MVD_avgProtein_03050 <- mean(proteinKEGG_MVD_MAD[c(11:15),1])
MVD_avgProtein_04010 <- mean(proteinKEGG_MVD_MAD[c(16:20),1])
MVD_avgProtein_04120 <- mean(proteinKEGG_MVD_MAD[c(21:25),1])
MVD_avgProtein_04210 <- mean(proteinKEGG_MVD_MAD[c(26:30),1])

#reorganize the treatment-specific values so they can be easily plotted
MVD_CGPgene <- geneKEGG_MVD_MAD[c(1,6,11,16,21,26),1]
MVD_FAEgene <- geneKEGG_MVD_MAD[c(2,7,12,17,22,27),1]
MVD_FAPgene <- geneKEGG_MVD_MAD[c(3,8,13,18,23,28),1]
MVD_POHgene <- geneKEGG_MVD_MAD[c(4,9,14,19,24,29),1]
MVD_POLgene <- geneKEGG_MVD_MAD[c(5,10,15,20,25,30),1]

MVD_CGPprotein <- proteinKEGG_MVD_MAD[c(1,6,11,16,21,26),1]
MVD_FAEprotein <- proteinKEGG_MVD_MAD[c(2,7,12,17,22,27),1]
MVD_FAPprotein <- proteinKEGG_MVD_MAD[c(3,8,13,18,23,28),1]
MVD_POHprotein <- proteinKEGG_MVD_MAD[c(4,9,14,19,24,29),1]
MVD_POLprotein <- proteinKEGG_MVD_MAD[c(5,10,15,20,25,30),1]

#look at average over genes within a pathway
#across all treatments
MAD_avgAllGene_00010 <- mean(rowMeans(cbind(MAD_m00010_gCGP, MAD_m00010_gFAE, MAD_m00010_gFAP, MAD_m00010_gPOH, MAD_m00010_gPOL)))
MAD_avgAllProtein_00010 <- mean(rowMeans(cbind(MAD_m00010_pCGP, MAD_m00010_pFAE, MAD_m00010_pFAP, MAD_m00010_pPOH, MAD_m00010_pPOL)))
MAD_avgAllGene_03050 <- mean(rowMeans(cbind(MAD_m03050_gCGP, MAD_m03050_gFAE, MAD_m03050_gFAP, MAD_m03050_gPOH, MAD_m03050_gPOL)))
MAD_avgAllProtein_03050 <- mean(rowMeans(cbind(MAD_m03050_pCGP, MAD_m03050_pFAE, MAD_m03050_pFAP, MAD_m03050_pPOH, MAD_m03050_pPOL)))
MAD_avgAllGene_04010 <- mean(rowMeans(cbind(MAD_m04010_gCGP, MAD_m04010_gFAE, MAD_m04010_gFAP, MAD_m04010_gPOH, MAD_m04010_gPOL)))
MAD_avgAllProtein_04010 <- mean(rowMeans(cbind(MAD_m04010_pCGP, MAD_m04010_pFAE, MAD_m04010_pFAP, MAD_m04010_pPOH, MAD_m04010_pPOL)))
MAD_avgAllGene_04120 <- mean(rowMeans(cbind(MAD_m04120_gCGP, MAD_m04120_gFAE, MAD_m04120_gFAP, MAD_m04120_gPOH, MAD_m04120_gPOL)))
MAD_avgAllProtein_04120 <- mean(rowMeans(cbind(MAD_m04120_pCGP, MAD_m04120_pFAE, MAD_m04120_pFAP, MAD_m04120_pPOH, MAD_m04120_pPOL)))
MAD_avgAllGene_04210 <- mean(rowMeans(cbind(MAD_m04210_gCGP, MAD_m04210_gFAE, MAD_m04210_gFAP, MAD_m04210_gPOH, MAD_m04210_gPOL)))
MAD_avgAllProtein_04210 <- mean(rowMeans(cbind(MAD_m04210_pCGP, MAD_m04210_pFAE, MAD_m04210_pFAP, MAD_m04210_pPOH, MAD_m04210_pPOL)))
MAD_avgAllGene_00020 <- mean(rowMeans(cbind(MAD_m00020_gCGP, MAD_m00020_gFAE, MAD_m00020_gFAP, MAD_m00020_gPOH, MAD_m00020_gPOL)))
MAD_avgAllProtein_00020 <- mean(rowMeans(cbind(MAD_m00020_pCGP, MAD_m00020_pFAE, MAD_m00020_pFAP, MAD_m00020_pPOH, MAD_m00020_pPOL)))

#within one treatment (POH?)
MAD_avgPOHgene_00010 <- mean(MAD_m00010_gPOH)
MAD_avgPOHprotein_00010 <- mean(MAD_m00010_pPOH)
MAD_avgPOHgene_03050 <- mean(MAD_m03050_gPOH)
MAD_avgPOHprotein_03050 <- mean(MAD_m03050_pPOH)
MAD_avgPOHgene_04010 <- mean(MAD_m04010_gPOH)
MAD_avgPOHprotein_04010 <- mean(MAD_m04010_pPOH)
MAD_avgPOHgene_04120 <- mean(MAD_m04120_gPOH)
MAD_avgPOHprotein_04120 <- mean(MAD_m04120_pPOH)
MAD_avgPOHgene_04210 <- mean(MAD_m04210_gPOH)
MAD_avgPOHprotein_04210 <- mean(MAD_m04210_pPOH)
MAD_avgPOHgene_00020 <- mean(MAD_m00020_gPOH)
MAD_avgPOHprotein_00020 <- mean(MAD_m00020_pPOH)

MAD_avgPOH <- c(MAD_avgPOHgene_00010, MAD_avgPOHgene_00020, MAD_avgPOHgene_03050, MAD_avgPOHgene_04010, MAD_avgPOHgene_04120, MAD_avgPOHgene_04210,
                      MAD_avgPOHprotein_00010, MAD_avgPOHprotein_00020, MAD_avgPOHprotein_03050, MAD_avgPOHprotein_04010, MAD_avgPOHprotein_04120, MAD_avgPOHprotein_04210)


#within POL
MAD_avgPOLgene_00010 <- mean(MAD_m00010_gPOL)
MAD_avgPOLprotein_00010 <- mean(MAD_m00010_pPOL)
MAD_avgPOLgene_03050 <- mean(MAD_m03050_gPOL)
MAD_avgPOLprotein_03050 <- mean(MAD_m03050_pPOL)
MAD_avgPOLgene_04010 <- mean(MAD_m04010_gPOL)
MAD_avgPOLprotein_04010 <- mean(MAD_m04010_pPOL)
MAD_avgPOLgene_04120 <- mean(MAD_m04120_gPOL)
MAD_avgPOLprotein_04120 <- mean(MAD_m04120_pPOL)
MAD_avgPOLgene_04210 <- mean(MAD_m04210_gPOL)
MAD_avgPOLprotein_04210 <- mean(MAD_m04210_pPOL)
MAD_avgPOLgene_00020 <- mean(MAD_m00020_gPOL)
MAD_avgPOLprotein_00020 <- mean(MAD_m00020_pPOL)

MAD_avgPOL <- c(MAD_avgPOLgene_00010, MAD_avgPOLgene_00020, MAD_avgPOLgene_03050, MAD_avgPOLgene_04010, MAD_avgPOLgene_04120, MAD_avgPOLgene_04210,
                      MAD_avgPOLprotein_00010, MAD_avgPOLprotein_00020, MAD_avgPOLprotein_03050, MAD_avgPOLprotein_04010, MAD_avgPOLprotein_04120, MAD_avgPOLprotein_04210)


#within CGP
MAD_avgCGPgene_00010 <- mean(MAD_m00010_gCGP)
MAD_avgCGPprotein_00010 <- mean(MAD_m00010_pCGP)
MAD_avgCGPgene_03050 <- mean(MAD_m03050_gCGP)
MAD_avgCGPprotein_03050 <- mean(MAD_m03050_pCGP)
MAD_avgCGPgene_04010 <- mean(MAD_m04010_gCGP)
MAD_avgCGPprotein_04010 <- mean(MAD_m04010_pCGP)
MAD_avgCGPgene_04120 <- mean(MAD_m04120_gCGP)
MAD_avgCGPprotein_04120 <- mean(MAD_m04120_pCGP)
MAD_avgCGPgene_04210 <- mean(MAD_m04210_gCGP)
MAD_avgCGPprotein_04210 <- mean(MAD_m04210_pCGP)
MAD_avgCGPgene_00020 <- mean(MAD_m00020_gCGP)
MAD_avgCGPprotein_00020 <- mean(MAD_m00020_pCGP)

MAD_avgCGP <- c(MAD_avgCGPgene_00010, MAD_avgCGPgene_00020, MAD_avgCGPgene_03050, MAD_avgCGPgene_04010, MAD_avgCGPgene_04120, MAD_avgCGPgene_04210,
                      MAD_avgCGPprotein_00010, MAD_avgCGPprotein_00020, MAD_avgCGPprotein_03050, MAD_avgCGPprotein_04010, MAD_avgCGPprotein_04120, MAD_avgCGPprotein_04210)

#within FAE
MAD_avgFAEgene_00010 <- mean(MAD_m00010_gFAE)
MAD_avgFAEprotein_00010 <- mean(MAD_m00010_pFAE)
MAD_avgFAEgene_03050 <- mean(MAD_m03050_gFAE)
MAD_avgFAEprotein_03050 <- mean(MAD_m03050_pFAE)
MAD_avgFAEgene_04010 <- mean(MAD_m04010_gFAE)
MAD_avgFAEprotein_04010 <- mean(MAD_m04010_pFAE)
MAD_avgFAEgene_04120 <- mean(MAD_m04120_gFAE)
MAD_avgFAEprotein_04120 <- mean(MAD_m04120_pFAE)
MAD_avgFAEgene_04210 <- mean(MAD_m04210_gFAE)
MAD_avgFAEprotein_04210 <- mean(MAD_m04210_pFAE)
MAD_avgFAEgene_00020 <- mean(MAD_m00020_gFAE)
MAD_avgFAEprotein_00020 <- mean(MAD_m00020_pFAE)

MAD_avgFAE <- c(MAD_avgFAEgene_00010, MAD_avgFAEgene_00020, MAD_avgFAEgene_03050, MAD_avgFAEgene_04010, MAD_avgFAEgene_04120, MAD_avgFAEgene_04210,
                      MAD_avgFAEprotein_00010, MAD_avgFAEprotein_00020, MAD_avgFAEprotein_03050, MAD_avgFAEprotein_04010, MAD_avgFAEprotein_04120, MAD_avgFAEprotein_04210)


#within FAP
MAD_avgFAPgene_00010 <- mean(MAD_m00010_gFAP)
MAD_avgFAPprotein_00010 <- mean(MAD_m00010_pFAP)
MAD_avgFAPgene_03050 <- mean(MAD_m03050_gFAP)
MAD_avgFAPprotein_03050 <- mean(MAD_m03050_pFAP)
MAD_avgFAPgene_04010 <- mean(MAD_m04010_gFAP)
MAD_avgFAPprotein_04010 <- mean(MAD_m04010_pFAP)
MAD_avgFAPgene_04120 <- mean(MAD_m04120_gFAP)
MAD_avgFAPprotein_04120 <- mean(MAD_m04120_pFAP)
MAD_avgFAPgene_04210 <- mean(MAD_m04210_gFAP)
MAD_avgFAPprotein_04210 <- mean(MAD_m04210_pFAP)
MAD_avgFAPgene_00020 <- mean(MAD_m00020_gFAP)
MAD_avgFAPprotein_00020 <- mean(MAD_m00020_pFAP)

MAD_avgFAP <- c(MAD_avgFAPgene_00010, MAD_avgFAPgene_00020, MAD_avgFAPgene_03050, MAD_avgFAPgene_04010, MAD_avgFAPgene_04120, MAD_avgFAPgene_04210,
                      MAD_avgFAPprotein_00010, MAD_avgFAPprotein_00020, MAD_avgFAPprotein_03050, MAD_avgFAPprotein_04010, MAD_avgFAPprotein_04120, MAD_avgFAPprotein_04210)

#make a giant data frame for plotting
#rows = pathways x data type (12 total rows)
#columns = MVD, allMAD, POH_MAD, POL_MAD, FAE_MAD, FAP_MAD, CGP_MAD
plot_MAD <- data.frame(pathway=rep(c("m00010", "m00020", "m03050", "m04010", "m4120", "m04210"), 2),
                       type=rep(c("gene", "protein"), each=6),
                       MAD_all=c(MAD_avgAllGene_00010, MAD_avgAllGene_00020, MAD_avgAllGene_03050, MAD_avgAllGene_04010, MAD_avgAllGene_04120, MAD_avgAllGene_04210,
                                 MAD_avgAllProtein_00010, MAD_avgAllProtein_00020, MAD_avgAllProtein_03050, MAD_avgAllProtein_04010, MAD_avgAllProtein_04120, MAD_avgAllProtein_04210),
                       MAD_CGP=MAD_avgCGP,
                       MAD_FAE=MAD_avgFAE,
                       MAD_FAP=MAD_avgFAP,
                       MAD_POL=MAD_avgPOL,
                       MAD_POH=MAD_avgPOH)

#plot all together
ggplot(plot_MAD, aes(x=pathway, y=MAD_all, group=type)) +
  geom_col(aes(x=pathway, y=MAD_all, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD across transcripts") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

library(reshape2)
plot_MAD_melt <- melt(plot_MAD)

m00010_meltMAD <- subset(plot_MAD_melt, pathway == "m00010")
ggplot(m00010_meltMAD, aes(x=variable, y=value, group=type)) +
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

m00020_meltMAD <- subset(plot_MAD_melt, pathway == "m00020")
ggplot(m00020_meltMAD, aes(x=variable, y=value, group=type)) +
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

m03050_meltMAD <- subset(plot_MAD_melt, pathway == "m03050")
ggplot(m03050_meltMAD, aes(x=variable, y=value, group=type)) +
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

m04010_meltMAD <- subset(plot_MAD_melt, pathway == "m04010")
ggplot(m04010_meltMAD, aes(x=variable, y=value, group=type)) +
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

m04120_meltMAD <- subset(plot_MAD_melt, pathway == "m4120")
ggplot(m04120_meltMAD, aes(x=variable, y=value, group=type)) +
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

m04210_meltMAD <- subset(plot_MAD_melt, pathway == "m04210")
ggplot(m04210_meltMAD, aes(x=variable, y=value, group=type)) +
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
ggplot(plot_MAD, aes(x=pathway, y=MAD_POH, group=type)) +
  geom_col(aes(x=pathway, y=MAD_POH, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in POH") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot POL
ggplot(plot_MAD, aes(x=pathway, y=MAD_POL, group=type)) +
  geom_col(aes(x=pathway, y=MAD_POL, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in POL") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot CGP
ggplot(plot_MAD, aes(x=pathway, y=MAD_CGP, group=type)) +
  geom_col(aes(x=pathway, y=MAD_CGP, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in CGP") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot FAE
ggplot(plot_MAD, aes(x=pathway, y=MAD_FAE, group=type)) +
  geom_col(aes(x=pathway, y=MAD_FAE, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in FAE") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))

#plot FAP
ggplot(plot_MAD, aes(x=pathway, y=MAD_FAP, group=type)) +
  geom_col(aes(x=pathway, y=MAD_FAP, group=type, fill=type), color="black", position = "dodge") +
  scale_fill_manual(values=c(gene="white", protein="black")) +
  theme_bw() +
  ylab("Average scaled MAD within cluster in FAP") +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))



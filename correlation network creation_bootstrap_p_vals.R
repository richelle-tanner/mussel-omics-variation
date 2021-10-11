#create a correlation network to import into Cytoscape
#multiple types: genes, proteins, by treatment, overall
#nodes/edges need to be transcripts, not individuals!

#install packages
library(igraph)
library(RCy3)

#import data
ddr <- "~/Documents/Cytoscape_projects/M_Cali/"
gene <- read.table("~/Documents/R_Data/musselRNASeq/final_files/gene_matrixllsImpute_6.8.20.txt")
protein <- read.table("~/Documents/R_Data/musselRNASeq/final_files/protein_matrixllsImpute_6.8.20.txt") #use lls impute method 


#row names
gene_names <- rownames(gene)
protein_names <- rownames(protein)

#make numeric
gene <- data.frame(sapply(gene, function(x) as.numeric(as.character(x))))
#sometimes this makes NAs? Hmmm
protein <- data.frame(sapply(protein, function(x) as.numeric(as.character(x))))

####cull using p values####
#import p value files
boot_folder <- "~/Documents/R_Data/musselRNASeq/bootstrapping/"
p_val_gene <- read.table(paste(boot_folder, "bootstrap_p_vals_allMAD_gene_6.15.20.txt", sep=''), header=TRUE)
p_val_protein <- read.table(paste(boot_folder, "bootstrap_p_vals_allMAD_protein_6.15.20.txt", sep=''), header=TRUE)
#import name file for cytoscape
cytoscape_folder <- "~/Documents/Cytoscape_projects/M_Cali/"
gene_labels <- read.table(paste(cytoscape_folder, "gene_names_labels.txt", sep=''), header=TRUE)

#match names and cytoscape name IDs
colnames(p_val_gene)[2] <- "gene_names"
p_val_gene_labeled <- merge(p_val_gene, gene_labels, by = "gene_names")
colnames(p_val_gene_labeled)[4] <- "shared.name"

colnames(p_val_protein)[2] <- "gene_names"
p_val_protein_labeled <- merge(p_val_protein, gene_labels, by = "gene_names")
colnames(p_val_protein_labeled)[4] <- "shared.name"

#get list of non-sig p values
non_sig_gene <- p_val_gene_labeled[ which(p_val_gene_labeled$p_vals_gene > 0.051), ]
non_sig_protein <- p_val_protein_labeled[ which(p_val_protein_labeled$p_vals_protein > 0.051), ]

#make separate data frames for each treatment
gene_z <- scale(t(gene))
protein_z <- scale(t(protein))

#add treatment label column
#thankfully they're in the same order
treatments <- c("FAE", "POL", "FAE", "FAP", "POH", "CGP", "POH", "FAE", "FAP", "FAE", "POL", "FAP", "FAP", "FAE",
                "FAP", "POH", "CGP", "CGP", "POL", "FAE", "FAP", "POH", "POH", "POL", "CGP", "POL", "POL", "POH",
                "POH", "POH", "POL", "POL", "FAE", "FAP", "POH", "CGP", "FAE", "FAP", "POH", "FAE", "FAP")

gene_z <- cbind(treatments, gene_z)
protein_z <- cbind(treatments, protein_z)

gene_z <- data.frame(gene_z)
protein_z <- data.frame(protein_z)

colnames(gene_z) <- c("treatments", gene_names)
colnames(protein_z) <- c("treatments", protein_names)

CGP_gene <- gene_z[which(gene_z$treatments == "CGP"),]
CGP_gene$treatments <- NULL
CGP_protein <- protein_z[which(protein_z$treatments == "CGP"),]
CGP_protein$treatments <- NULL
FAE_gene <- gene_z[which(gene_z$treatments == "FAE"),]
FAE_gene$treatments <- NULL
FAE_protein <- protein_z[which(protein_z$treatments == "FAE"),]
FAE_protein$treatments <- NULL
FAP_gene <- gene_z[which(gene_z$treatments == "FAP"),]
FAP_gene$treatments <- NULL
FAP_protein <- protein_z[which(protein_z$treatments == "FAP"),]
FAP_protein$treatments <- NULL
POH_gene <- gene_z[which(gene_z$treatments == "POH"),]
POH_gene$treatments <- NULL
POH_protein <- protein_z[which(protein_z$treatments == "POH"),]
POH_protein$treatments <- NULL
POL_gene <- gene_z[which(gene_z$treatments == "POL"),]
POL_gene$treatments <- NULL
POL_protein <- protein_z[which(protein_z$treatments == "POL"),]
POL_protein$treatments <- NULL

#subset by taking out non-sig p value genes
CGP_non_sig_gene <- non_sig_gene[which(non_sig_gene$treatment == "CGP"),]
CGP_gene_sig_list <- which(colnames(CGP_gene) %in% CGP_non_sig_gene$gene_names)
CGP_gene_all_sig <- CGP_gene[,-CGP_gene_sig_list]

FAE_non_sig_gene <- non_sig_gene[which(non_sig_gene$treatment == "FAE"),]
FAE_gene_sig_list <- which(colnames(FAE_gene) %in% FAE_non_sig_gene$gene_names)
FAE_gene_all_sig <- FAE_gene[,-FAE_gene_sig_list]

FAP_non_sig_gene <- non_sig_gene[which(non_sig_gene$treatment == "FAP"),]
FAP_gene_sig_list <- which(colnames(FAP_gene) %in% FAP_non_sig_gene$gene_names)
FAP_gene_all_sig <- FAP_gene[,-FAP_gene_sig_list]

POH_non_sig_gene <- non_sig_gene[which(non_sig_gene$treatment == "POH"),]
POH_gene_sig_list <- which(colnames(POH_gene) %in% POH_non_sig_gene$gene_names)
POH_gene_all_sig <- POH_gene[,-POH_gene_sig_list]

POL_non_sig_gene <- non_sig_gene[which(non_sig_gene$treatment == "POL"),]
POL_gene_sig_list <- which(colnames(POL_gene) %in% POL_non_sig_gene$gene_names)
POL_gene_all_sig <- POL_gene[,-POL_gene_sig_list]

CGP_non_sig_protein <- non_sig_protein[which(non_sig_protein$treatment == "CGP"),]
CGP_protein_sig_list <- which(colnames(CGP_protein) %in% CGP_non_sig_protein$gene_names)
CGP_protein_all_sig <- CGP_protein[,-CGP_protein_sig_list]

FAE_non_sig_protein <- non_sig_protein[which(non_sig_protein$treatment == "FAE"),]
FAE_protein_sig_list <- which(colnames(FAE_protein) %in% FAE_non_sig_protein$gene_names)
FAE_protein_all_sig <- FAE_protein[,-FAE_protein_sig_list]

FAP_non_sig_protein <- non_sig_protein[which(non_sig_protein$treatment == "FAP"),]
FAP_protein_sig_list <- which(colnames(FAP_protein) %in% FAP_non_sig_protein$gene_names)
FAP_protein_all_sig <- FAP_protein[,-FAP_protein_sig_list]

POH_non_sig_protein <- non_sig_protein[which(non_sig_protein$treatment == "POH"),]
POH_protein_sig_list <- which(colnames(POH_protein) %in% POH_non_sig_protein$gene_names)
POH_protein_all_sig <- POH_protein[,-POH_protein_sig_list]

POL_non_sig_protein <- non_sig_protein[which(non_sig_protein$treatment == "POL"),]
POL_protein_sig_list <- which(colnames(POL_protein) %in% POL_non_sig_protein$gene_names)
POL_protein_all_sig <- POL_protein[,-POL_protein_sig_list]

####predetermine high and low MAD (1st and last quintiles) genes by treatments and build networks separately####

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

CGPg_MAD <- MAD_function(CGP_gene_all_sig)
FAEg_MAD <- MAD_function(FAE_gene_all_sig)
FAPg_MAD <- MAD_function(FAP_gene_all_sig)
POHg_MAD <- MAD_function(POH_gene_all_sig)
POLg_MAD <- MAD_function(POL_gene_all_sig)

CGPp_MAD <- MAD_function(CGP_protein_all_sig)
FAEp_MAD <- MAD_function(FAE_protein_all_sig)
FAPp_MAD <- MAD_function(FAP_protein_all_sig)
POHp_MAD <- MAD_function(POH_protein_all_sig)
POLp_MAD <- MAD_function(POL_protein_all_sig)

#put the gene names back in
CGPg_MAD <- data.frame(transcript=colnames(CGP_gene_all_sig), MAD=CGPg_MAD)
FAEg_MAD <- data.frame(transcript=colnames(FAE_gene_all_sig), MAD=FAEg_MAD)
FAPg_MAD <- data.frame(transcript=colnames(FAP_gene_all_sig), MAD=FAPg_MAD)
POHg_MAD <- data.frame(transcript=colnames(POH_gene_all_sig), MAD=POHg_MAD)
POLg_MAD <- data.frame(transcript=colnames(POL_gene_all_sig), MAD=POLg_MAD)

CGPp_MAD <- data.frame(transcript=colnames(CGP_protein_all_sig), MAD=CGPp_MAD)
FAEp_MAD <- data.frame(transcript=colnames(FAE_protein_all_sig), MAD=FAEp_MAD)
FAPp_MAD <- data.frame(transcript=colnames(FAP_protein_all_sig), MAD=FAPp_MAD)
POHp_MAD <- data.frame(transcript=colnames(POH_protein_all_sig), MAD=POHp_MAD)
POLp_MAD <- data.frame(transcript=colnames(POL_protein_all_sig), MAD=POLp_MAD)

#put them in order from low to high MAD
#remove NAs, if any
CGPg_MAD <- na.omit(CGPg_MAD[rev(order(CGPg_MAD$MAD)),])
FAEg_MAD <- na.omit(FAEg_MAD[rev(order(FAEg_MAD$MAD)),])
FAPg_MAD <- na.omit(FAPg_MAD[rev(order(FAPg_MAD$MAD)),])
POHg_MAD <- na.omit(POHg_MAD[rev(order(POHg_MAD$MAD)),])
POLg_MAD <- na.omit(POLg_MAD[rev(order(POLg_MAD$MAD)),])

CGPp_MAD <- na.omit(CGPp_MAD[rev(order(CGPp_MAD$MAD)),])
FAEp_MAD <- na.omit(FAEp_MAD[rev(order(FAEp_MAD$MAD)),])
FAPp_MAD <- na.omit(FAPp_MAD[rev(order(FAPp_MAD$MAD)),])
POHp_MAD <- na.omit(POHp_MAD[rev(order(POHp_MAD$MAD)),])
POLp_MAD <- na.omit(POLp_MAD[rev(order(POLp_MAD$MAD)),])

#get top 20% of genes/proteins organized by MAD (extract row numbers)
#20% of 1521 is 304.2 (=304) -- we still want this number because we want top 20% of all genes, regardless of sig.
#we want first 304 rows because it is already in reverse order
CGPg_MAD20 <- CGPg_MAD[c(1:304),]
FAEg_MAD20 <- FAEg_MAD[c(1:304),]
FAPg_MAD20 <- FAPg_MAD[c(1:304),]
POHg_MAD20 <- POHg_MAD[c(1:304),]
POLg_MAD20 <- POLg_MAD[c(1:304),]

CGPp_MAD20 <- CGPp_MAD[c(1:304),]
FAEp_MAD20 <- FAEp_MAD[c(1:304),]
FAPp_MAD20 <- FAPp_MAD[c(1:304),]
POHp_MAD20 <- POHp_MAD[c(1:304),]
POLp_MAD20 <- POLp_MAD[c(1:304),]

#create new DFs with new, reduced row numbers
#use gene and protein data frames to pull from
#only want individuals from that treatment
#needs to be genes x individuals
gene_CGPtop <- CGP_gene[,which(colnames(CGP_gene) %in% CGPg_MAD20$transcript)]
gene_FAEtop <- FAE_gene[,which(colnames(FAE_gene) %in% FAEg_MAD20$transcript)]
gene_FAPtop <- FAP_gene[,which(colnames(FAP_gene) %in% FAPg_MAD20$transcript)]
gene_POHtop <- POH_gene[,which(colnames(POH_gene) %in% POHg_MAD20$transcript)]
gene_POLtop <- POL_gene[,which(colnames(POL_gene) %in% POLg_MAD20$transcript)]

protein_CGPtop <- CGP_protein[,which(colnames(CGP_protein) %in% CGPp_MAD20$transcript)]
protein_FAEtop <- FAE_protein[,which(colnames(FAE_protein) %in% FAEp_MAD20$transcript)]
protein_FAPtop <- FAP_protein[,which(colnames(FAP_protein) %in% FAPp_MAD20$transcript)]
protein_POHtop <- POH_protein[,which(colnames(POH_protein) %in% POHp_MAD20$transcript)]
protein_POLtop <- POL_protein[,which(colnames(POL_protein) %in% POLp_MAD20$transcript)]

#make numeric again
gene_CGPtop <- sapply(gene_CGPtop, function(x) as.numeric(as.character(x)))
gene_FAEtop <- sapply(gene_FAEtop, function(x) as.numeric(as.character(x)))
gene_FAPtop <- sapply(gene_FAPtop, function(x) as.numeric(as.character(x)))
gene_POHtop <- sapply(gene_POHtop, function(x) as.numeric(as.character(x)))
gene_POLtop <- sapply(gene_POLtop, function(x) as.numeric(as.character(x)))

protein_CGPtop <- sapply(protein_CGPtop, function(x) as.numeric(as.character(x)))
protein_FAEtop <- sapply(protein_FAEtop, function(x) as.numeric(as.character(x)))
protein_FAPtop <- sapply(protein_FAPtop, function(x) as.numeric(as.character(x)))
protein_POHtop <- sapply(protein_POHtop, function(x) as.numeric(as.character(x)))
protein_POLtop <- sapply(protein_POLtop, function(x) as.numeric(as.character(x)))

#make them genes x individuals
gene_CGPtop <- t(gene_CGPtop)
gene_FAEtop <- t(gene_FAEtop)
gene_FAPtop <- t(gene_FAPtop)
gene_POHtop <- t(gene_POHtop)
gene_POLtop <- t(gene_POLtop)

protein_CGPtop <- t(protein_CGPtop)
protein_FAEtop <- t(protein_FAEtop)
protein_FAPtop <- t(protein_FAPtop)
protein_POHtop <- t(protein_POHtop)
protein_POLtop <- t(protein_POLtop)


#make correlation matrices with new DFs
geneCorr_POLtop <- cor(t(gene_POLtop))
geneCorr_POHtop <- cor(t(gene_POHtop))
geneCorr_CGPtop <- cor(t(gene_CGPtop))
geneCorr_FAEtop <- cor(t(gene_FAEtop))
geneCorr_FAPtop <- cor(t(gene_FAPtop))

proteinCorr_POLtop <- cor(t(protein_POLtop))
proteinCorr_POHtop <- cor(t(protein_POHtop))
proteinCorr_CGPtop <- cor(t(protein_CGPtop))
proteinCorr_FAEtop <- cor(t(protein_FAEtop))
proteinCorr_FAPtop <- cor(t(protein_FAPtop))

#create network
net_genePOLtop <- graph_from_adjacency_matrix(geneCorr_POLtop, weighted=T, mode="directed", diag=F)
net_genePOHtop <- graph_from_adjacency_matrix(geneCorr_POHtop, weighted=T, mode="directed", diag=F)

net_proteinPOHtop <- graph_from_adjacency_matrix(proteinCorr_POHtop, weighted=T, mode="directed", diag=F)
net_proteinPOLtop <- graph_from_adjacency_matrix(proteinCorr_POLtop, weighted=T, mode="directed", diag=F)

net_geneCGPtop <- graph_from_adjacency_matrix(geneCorr_CGPtop, weighted=T, mode="directed", diag=F)
net_proteinCGPtop <- graph_from_adjacency_matrix(proteinCorr_CGPtop, weighted=T, mode="directed", diag=F)

net_geneFAEtop <- graph_from_adjacency_matrix(geneCorr_FAEtop, weighted=T, mode="directed", diag=F)
net_proteinFAEtop <- graph_from_adjacency_matrix(proteinCorr_FAEtop, weighted=T, mode="directed", diag=F)

net_geneFAPtop <- graph_from_adjacency_matrix(geneCorr_FAPtop, weighted=T, mode="directed", diag=F)
net_proteinFAPtop <- graph_from_adjacency_matrix(proteinCorr_FAPtop, weighted=T, mode="directed", diag=F)

#export network to Cytoscape
createNetworkFromIgraph(net_genePOLtop,"POL_gene_top20MAD")
createNetworkFromIgraph(net_genePOHtop,"POH_gene_top20MAD")
createNetworkFromIgraph(net_proteinPOHtop,"POH_protein_top20MAD")
createNetworkFromIgraph(net_proteinPOLtop,"POL_protein_top20MAD")

createNetworkFromIgraph(net_geneCGPtop,"CGP_gene_top20MAD")
createNetworkFromIgraph(net_proteinCGPtop,"CGP_protein_top20MAD")
createNetworkFromIgraph(net_geneFAEtop,"FAE_gene_top20MAD")
createNetworkFromIgraph(net_proteinFAEtop,"FAE_protein_top20MAD")
createNetworkFromIgraph(net_geneFAPtop,"FAP_gene_top20MAD")
createNetworkFromIgraph(net_proteinFAPtop,"FAP_protein_top20MAD")

####export COVARIANCE matrices for evolQG analysis####
geneCovar_POLtop <- cov(t(gene_POLtop))
geneCovar_POHtop <- cov(t(gene_POHtop))
geneCovar_CGPtop <- cov(t(gene_CGPtop))
geneCovar_FAEtop <- cov(t(gene_FAEtop))
geneCovar_FAPtop <- cov(t(gene_FAPtop))

proteinCovar_POLtop <- cov(t(protein_POLtop))
proteinCovar_POHtop <- cov(t(protein_POHtop))
proteinCovar_CGPtop <- cov(t(protein_CGPtop))
proteinCovar_FAEtop <- cov(t(protein_FAEtop))
proteinCovar_FAPtop <- cov(t(protein_FAPtop))

setwd("~/Documents/Cytoscape_projects/M_Cali/covariance_matrices/")
write.table(geneCovar_POLtop, "covarMat_gene_POLtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(geneCovar_POHtop, "covarMat_gene_POHtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(geneCovar_CGPtop, "covarMat_gene_CGPtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(geneCovar_FAEtop, "covarMat_gene_FAEtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(geneCovar_FAPtop, "covarMat_gene_FAPtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")

write.table(proteinCovar_POLtop, "covarMat_protein_POLtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(proteinCovar_POHtop, "covarMat_protein_POHtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(proteinCovar_CGPtop, "covarMat_protein_CGPtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(proteinCovar_FAEtop, "covarMat_protein_FAEtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(proteinCovar_FAPtop, "covarMat_protein_FAPtop20MAD.txt", row.names = T, col.names = T, quote = F, sep = "\t")


####predetermine high DE genes overall (top quintile)####
#first do this only with two treatments: outplant low and high
#z score transform data
gene_z <- scale(t(gene))
protein_z <- scale(t(protein))

#add treatment label column
#thankfully they're in the same order
treatments <- c("FAE", "POL", "FAE", "FAP", "POH", "CGP", "POH", "FAE", "FAP", "FAE", "POL", "FAP", "FAP", "FAE",
                "FAP", "POH", "CGP", "CGP", "POL", "FAE", "FAP", "POH", "POH", "POL", "CGP", "POL", "POL", "POH",
                "POH", "POH", "POL", "POL", "FAE", "FAP", "POH", "CGP", "FAE", "FAP", "POH", "FAE", "FAP")
#gene_z <- data.frame(cbind(treatments, gene_z))
#protein_z <- data.frame(cbind(treatments, protein_z))

gene_z <- cbind(treatments, gene_z)
protein_z <- cbind(treatments, protein_z)

#take out only POL and POH
gene_z_POL <- subset(gene_z, treatments == "POL")
gene_z_POH <- subset(gene_z, treatments == "POH")
protein_z_POL <- subset(protein_z, treatments == "POL")
protein_z_POH <- subset(protein_z, treatments == "POH")

#do an anova and sort by p value
#empty p value vectors - only one comparison, so just a vector for now
DEsig_g <- rep((length(gene_z[1,])-1), 0)
DEsig_p <- rep((length(gene_z[1,])-1), 0)
for (i in 1:(length(gene_z[1,])-1)) {
  #separate by treatment
  subset_g <- subset(gene_z, treatments == treat_step[i])
  subset_p <- subset(protein_z, treatments == treat_step[i])
  subset_g <- subset_g[,-1]
  subset_p <- subset_p[,-1]
  #calculate MAD per gene/protein
  for (j in 1:(length(gene_z[1,])-1)) {
    sub_gG <- as.numeric(as.character(subset_g[,j]))
    sub_pP <- as.numeric(as.character(subset_p[,j]))
    #find median for each gene
    median_g <- median(sub_gG)
    median_p <- median(sub_pP)
    #calculate deviations
    deviations_g <- sub_gG - median_g
    deviations_p <- sub_pP - median_p
    #calculate and store MAD
    MADg[j,i] <- median(abs(deviations_g))
    MADp[j,i] <- median(abs(deviations_p))
  }
}


#####using all samples####
#make correlation matrices - all genes/proteins, not separated by treatment
#this is Pearson's!

gene_allMAT <- cor(t(gene))
protein_allMAT <- cor(t(protein))

#add rownames back in
rownames(gene_allMAT) <- gene_names
rownames(protein_allMAT) <- protein_names

write.table(gene_allMAT, "~/full_gene_network_adjacency_matrix.txt", quote=F, sep="\t", row.names=T)
write.table(protein_allMAT, "~/full_protein_network_adjacency_matrix.txt", quote=F, sep="\t", row.names=T)

#get the transcript IDs for the rows and export separately
gene_names_transcripts <- data.frame(gene_names=gene_names, gene_numbers=seq(1, length(gene_names)))
write.table(gene_names_transcripts, "~/gene_names_labels.txt", quote=F, sep="\t", row.names=T)

#get rid of anything that's under a correlation of 0.6 - threshold is from Ecker et al. 2017
#if you don't do this, it will not export (too big)
gene_allMAT[gene_allMAT<0.6]=0
protein_allMAT[protein_allMAT<0.6]=0

networkG <- graph_from_adjacency_matrix(gene_allMAT, weighted=T, mode="directed", diag=F, add.rownames=T)
#plot(networkG)

networkP <- graph_from_adjacency_matrix(protein_allMAT, weighted=T, mode="directed", diag=F, add.rownames=T)
#plot(networkP)

createNetworkFromIgraph(networkG,"allTreat_geneIG")
createNetworkFromIgraph(networkP, "allTreat_proteinIG")

####export the transcripts that are high MAD for every treatment####
#the column is POH_MADg$name
#we want rows 1218:1521 (304 of them)
#add in a label for treatment

top20MAD_protein_transcripts <- data.frame(treatment=rep(c("POH", "POL", "CGP", "FAE", "FAP"),each=304),
                                   transcript=c(POH_MADp$name[1218:1521], POL_MADp$name[1218:1521], 
                                                CGP_MADp$name[1218:1521], FAE_MADp$name[1218:1521],
                                                FAP_MADp$name[1218:1521]))

top20MAD_gene_transcripts <- data.frame(treatment=rep(c("POH", "POL", "CGP", "FAE", "FAP"),each=304),
                                        transcript=c(POH_MADg$name[1218:1521], POL_MADg$name[1218:1521], 
                                                     CGP_MADg$name[1218:1521], FAE_MADg$name[1218:1521],
                                                     FAP_MADg$name[1218:1521]))

setwd("~/Documents/R_Data/musselRNASeq/")
write.table(top20MAD_protein_transcripts, "top20MAD_protein_transcripts.txt", quote=F, sep="\t", row.names=F)
write.table(top20MAD_gene_transcripts, "top20MAD_gene_transcripts.txt", quote=F, sep="\t", row.names=F)

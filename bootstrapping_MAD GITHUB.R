#we want p values so we can select the genes/proteins to highlight in network analyses
library(boot)

#import data
ddr <- "~/Documents/R_Data/musselRNASeq/final_files/"

gene <- read.table("~/gene_matrixllsImpute_6.8.20.txt")
protein <- read.table("~/protein_matrixllsImpute_6.8.20.txt") #use lls impute method 

#row names
gene_names <- rownames(gene)
protein_names <- rownames(protein)

#make numeric
gene <- data.frame(sapply(gene, function(x) as.numeric(as.character(x))))
#sometimes this makes NAs? Hmmm
protein <- data.frame(sapply(protein, function(x) as.numeric(as.character(x))))

gene_z <- scale(t(gene))
protein_z <- scale(t(protein))

#add treatment label column
#thankfully they're in the same order
treatments <- c("FAE", "POL", "FAE", "FAP", "POH", "CGP", "POH", "FAE", "FAP", "FAE", "POL", "FAP", "FAP", "FAE",
                "FAP", "POH", "CGP", "CGP", "POL", "FAE", "FAP", "POH", "POH", "POL", "CGP", "POL", "POL", "POH",
                "POH", "POH", "POL", "POL", "FAE", "FAP", "POH", "CGP", "FAE", "FAP", "POH", "FAE", "FAP")

gene_z <- cbind(treatments, gene_z)
protein_z <- cbind(treatments, protein_z)

####make treatment subsets####
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

MAD_function <- function(data_input,indices) {
  data_input <- data_input[indices,]
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

####run bootstrapping####
results_geneCGP <- boot(data=CGP_gene, statistic=MAD_function, R=10000, sim="ordinary")
results_genePOL <- boot(data=POL_gene, statistic=MAD_function, R=10000, sim="ordinary")
results_genePOH <- boot(data=POH_gene, statistic=MAD_function, R=10000, sim="ordinary")
results_geneFAE <- boot(data=FAE_gene, statistic=MAD_function, R=10000, sim="ordinary")
results_geneFAP <- boot(data=FAP_gene, statistic=MAD_function, R=10000, sim="ordinary")

results_proteinCGP <- boot(data=CGP_protein, statistic=MAD_function, R=10000, sim="ordinary")
results_proteinPOL <- boot(data=POL_protein, statistic=MAD_function, R=10000, sim="ordinary")
results_proteinPOH <- boot(data=POH_protein, statistic=MAD_function, R=10000, sim="ordinary")
results_proteinFAE <- boot(data=FAE_protein, statistic=MAD_function, R=10000, sim="ordinary")
results_proteinFAP <- boot(data=FAP_protein, statistic=MAD_function, R=10000, sim="ordinary")


#generate p values
p_vals_geneCGP <- rep(0, length(results_geneCGP$t[1,]))
p_vals_genePOL <- rep(0, length(results_genePOL$t[1,]))
p_vals_genePOH <- rep(0, length(results_genePOH$t[1,]))
p_vals_geneFAE <- rep(0, length(results_geneFAE$t[1,]))
p_vals_geneFAP <- rep(0, length(results_geneFAP$t[1,]))

for (i in 1:length(p_vals_geneCGP)) {
  results.under.H0_gene <- results_geneCGP$t[,i] - mean(results_geneCGP$t[,i])
  p_vals_geneCGP[i] <- mean(abs(results.under.H0_gene) > abs(results_geneCGP$t0[i]))
}

for (i in 1:length(p_vals_genePOL)) {
  results.under.H0_gene <- results_genePOL$t[,i] - mean(results_genePOL$t[,i])
  p_vals_genePOL[i] <- mean(abs(results.under.H0_gene) > abs(results_genePOL$t0[i]))
}

for (i in 1:length(p_vals_genePOH)) {
  results.under.H0_gene <- results_genePOH$t[,i] - mean(results_genePOH$t[,i])
  p_vals_genePOH[i] <- mean(abs(results.under.H0_gene) > abs(results_genePOH$t0[i]))
}

for (i in 1:length(p_vals_geneFAE)) {
  results.under.H0_gene <- results_geneFAE$t[,i] - mean(results_geneFAE$t[,i])
  p_vals_geneFAE[i] <- mean(abs(results.under.H0_gene) > abs(results_geneFAE$t0[i]))
}

for (i in 1:length(p_vals_geneFAP)) {
  results.under.H0_gene <- results_geneFAP$t[,i] - mean(results_geneFAP$t[,i])
  p_vals_geneFAP[i] <- mean(abs(results.under.H0_gene) > abs(results_geneFAP$t0[i]))
}

p_vals_proteinCGP <- rep(0, length(results_proteinCGP$t[1,]))
p_vals_proteinPOL <- rep(0, length(results_proteinPOL$t[1,]))
p_vals_proteinPOH <- rep(0, length(results_proteinPOH$t[1,]))
p_vals_proteinFAE <- rep(0, length(results_proteinFAE$t[1,]))
p_vals_proteinFAP <- rep(0, length(results_proteinFAP$t[1,]))

for (i in 1:length(p_vals_proteinCGP)) {
  results.under.H0_protein <- results_proteinCGP$t[,i] - mean(results_proteinCGP$t[,i])
  p_vals_proteinCGP[i] <- mean(abs(results.under.H0_protein) > abs(results_proteinCGP$t0[i]))
}

for (i in 1:length(p_vals_proteinPOL)) {
  results.under.H0_protein <- results_proteinPOL$t[,i] - mean(results_proteinPOL$t[,i])
  p_vals_proteinPOL[i] <- mean(abs(results.under.H0_protein) > abs(results_proteinPOL$t0[i]))
}

for (i in 1:length(p_vals_proteinPOH)) {
  results.under.H0_protein <- results_proteinPOH$t[,i] - mean(results_proteinPOH$t[,i])
  p_vals_proteinPOH[i] <- mean(abs(results.under.H0_protein) > abs(results_proteinPOH$t0[i]))
}

for (i in 1:length(p_vals_proteinFAE)) {
  results.under.H0_protein <- results_proteinFAE$t[,i] - mean(results_proteinFAE$t[,i])
  p_vals_proteinFAE[i] <- mean(abs(results.under.H0_protein) > abs(results_proteinFAE$t0[i]))
}

for (i in 1:length(p_vals_proteinFAP)) {
  results.under.H0_protein <- results_proteinFAP$t[,i] - mean(results_proteinFAP$t[,i])
  p_vals_proteinFAP[i] <- mean(abs(results.under.H0_protein) > abs(results_proteinFAP$t0[i]))
}


#export
#put the gene and protein names back in
#every 304 is something out of the treatment vector POL, POH, CGP, FAE, FAP
treatment <- rep(c("CGP", "POL", "POH", "FAE", "FAP"),each=1519)
p_vals_gene_names <- c(colnames(CGP_gene), colnames(POL_gene), colnames(POH_gene), colnames(FAE_gene), colnames(FAP_gene))
p_vals_gene <- c(p_vals_geneCGP, p_vals_genePOL, p_vals_genePOH, p_vals_geneFAE, p_vals_geneFAP)
p_vals_gene_table <- cbind(treatment,p_vals_gene_names, p_vals_gene)
p_vals_protein_names <- c(colnames(CGP_protein), colnames(POL_protein), colnames(POH_protein), colnames(FAE_protein), colnames(FAP_protein))
p_vals_protein <- c(p_vals_proteinCGP, p_vals_proteinPOL, p_vals_proteinPOH, p_vals_proteinFAE, p_vals_proteinFAP)
p_vals_protein_table <- cbind(treatment,p_vals_protein_names, p_vals_protein)

write.table(p_vals_gene_table, file="~/bootstrap_p_vals_allMAD_gene.txt", quote=F, sep="\t", row.names=F)
write.table(p_vals_protein_table, file="~/bootstrapping/bootstrap_p_vals_allMAD_protein.txt", quote=F, sep="\t", row.names=F)

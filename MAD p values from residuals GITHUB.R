#this only does p values for every single gene, every MAD comparison
#use this for the MCODE cluster and KEGG pathway comparisons
#calculates absolute residuals for every point and then does an ANOVA or t test
#uses names from the "compare_KEGG_MCODE" folder

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

#import data
gene_matrix <- read.table("~/Documents/R_Data/musselRNASeq/final_files/gene_exp_noScale_noTrans_12.11.19.txt")
protein_matrix <- read.table("~/Documents/R_Data/musselRNASeq/final_files/protein_exp_noScale_noTrans_12.11.19.txt")

####prep data####
gene_matrix_org <- gene_matrix[order(row.names(gene_matrix)),]
protein_matrix_org <- protein_matrix[order(row.names(protein_matrix)),]
gene_names <- row.names(gene_matrix_org)
protein_names <- row.names(protein_matrix_org)

#transpose
gene_matrix <- t(gene_matrix_org)
protein_matrix <- t(protein_matrix_org)

#scale the data so we can compare gene and protein
gene <- scale(gene_matrix, center = TRUE, scale = TRUE)
protein <- scale(protein_matrix, center = TRUE, scale = TRUE)

#add treatment label column
#thankfully they're in the same order
treatments <- c("FAE", "POL", "FAE", "FAP", "POH", "CGP", "POH", "FAE", "FAP", "FAE", "POL", "FAP", "FAP", "FAE",
                "FAP", "POH", "CGP", "CGP", "POL", "FAE", "FAP", "POH", "POH", "POL", "CGP", "POL", "POL", "POH",
                "POH", "POH", "POL", "POL", "FAE", "FAP", "POH", "CGP", "FAE", "FAP", "POH", "FAE", "FAP")

gene <- cbind(treatments, gene)
protein <- cbind(treatments, protein)

gene <- data.frame(gene)
protein <- data.frame(protein)

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

####looking for significant residual comparisons####
#do the residuals calculation
residuals_CGPg <- residuals_function(CGP_gene)
residuals_FAEg <- residuals_function(FAE_gene)
residuals_FAPg <- residuals_function(FAP_gene)
residuals_POHg <- residuals_function(POH_gene)
residuals_POLg <- residuals_function(POL_gene)

residuals_CGPp <- residuals_function(CGP_protein)
residuals_FAEp <- residuals_function(FAE_protein)
residuals_FAPp <- residuals_function(FAP_protein)
residuals_POHp <- residuals_function(POH_protein)
residuals_POLp <- residuals_function(POL_protein)

#put them back together in a data frame
residuals_all <- rbind(residuals_CGPg, residuals_FAEg, residuals_FAPg, residuals_POHg, residuals_POLg,
                          residuals_CGPp, residuals_FAEp, residuals_FAPp, residuals_POHp, residuals_POLp)
residuals_all$treatment <- rep(c(rep("CGP", 5), rep("FAE", 9), rep("FAP", 9), rep("POH", 10), rep("POL", 8)),2)
residuals_all$type <- as.factor(rep(c("gene", "protein"), each=41))
residuals_all$treatment_type <- as.factor(paste(residuals_all$treatment, "_", residuals_all$type, sep=''))

#loop through to print p values for all multiple comparisons
allTukey_comp <- vector()
for (i in 1:(length(residuals_all[1,])-3)) {
  last_3 <- ncol(residuals_all)
  sub <- residuals_all[,c(i, (last_3-2):last_3)]
  lev_aov <- aov(sub[,1] ~ sub[,4])
  tukey <- TukeyHSD(lev_aov)$`sub[, 4]`
  for (j in 1:length(tukey[,4])) {
    add <- c(tukey[j,], rownames(tukey)[j], colnames(residuals_all[i]))
    allTukey_comp <- rbind(allTukey_comp, add)
  }
  
}

#multiple comparisons adjustment
allTukey_compAdj <- p.adjust(allTukey_comp[,4], "BH")

#cut down to only significant ones
allTukey_sigcomp <- vector()
for (i in 1:length(allTukey_compAdj)) {
  if (allTukey_compAdj[i] < 0.05) {
    sig <- c(allTukey_compAdj[i], allTukey_comp[i,5], allTukey_comp[i,6])
    allTukey_sigcomp <- rbind(allTukey_sigcomp, sig)
  }
}
colnames(allTukey_sigcomp) <- c("p_val", "comparison", "transcript")

#compare among treatments in separate gene and protein datasets
residuals_all_G <- residuals_all[which(residuals_all$type %in% "gene"),]
residuals_all_P <- residuals_all[which(residuals_all$type %in% "protein"),]
residuals_all_G$treatment <- as.factor(residuals_all_G$treatment)
residuals_all_P$treatment <- as.factor(residuals_all_P$treatment)

treat_comp_G <- vector()
contrasts(residuals_all_G$treatment) <- contr.sum
for (i in 1:(length(residuals_all_G[1,])-3)) {
  last_3 <- ncol(residuals_all_G)
  sub <- residuals_all_G[,c(i, (last_3-2):last_3)]
  lev_aov <- aov(sub[,1] ~ sub[,2])
  tukey <- TukeyHSD(lev_aov)$`sub[, 2]`
  for (j in 1:length(tukey[,4])) {
    # do this afterward so we can correct for multiple comparisons first 9.30.19
      #if (tukey[j,4] < 0.05) {
      #  sig <- c(tukey[j,], rownames(tukey)[j], colnames(residuals_all_G[i]))
      #  treat_comp_G <- rbind(treat_comp_G, sig)
      #}
    add <- c(tukey[j,], rownames(tukey)[j], colnames(residuals_all_G[i]))
    treat_comp_G <- rbind(treat_comp_G, add)
  }
}

#multiple comparisons adjustment
treat_comp_G <- data.frame(treat_comp_G)
treat_comp_G$p.adj <- as.numeric(as.character(treat_comp_G$p.adj))
treatG_compAdj <- p.adjust(treat_comp_G$p.adj, "BH")

#cut down to only significant ones
treat_sigcomp_G <- vector()
for (i in 1:length(treatG_compAdj)) {
  if (treatG_compAdj[i] < 0.05) {
    sig <- c(treatG_compAdj[i], treat_comp_G[i,5], treat_comp_G[i,6])
    treat_sigcomp_G <- rbind(treat_sigcomp_G, sig)
  }
}
colnames(treat_sigcomp_G) <- c("p_val", "comparison", "transcript")



treat_comp_P <- vector()
contrasts(residuals_all_P$treatment) <- contr.sum
for (i in 1:(length(residuals_all_P[1,])-3)) {
  last_3 <- ncol(residuals_all_P)
  sub <- residuals_all_P[,c(i, (last_3-2):last_3)]
  lev_aov <- aov(sub[,1] ~ sub[,2])
  tukey <- TukeyHSD(lev_aov)$`sub[, 2]`
  for (j in 1:length(tukey[,4])) {
    add <- c(tukey[j,], rownames(tukey)[j], colnames(residuals_all_P[i]))
    treat_comp_P <- rbind(treat_comp_P, add)
  }
}

#multiple comparisons adjustment
treat_comp_P <- data.frame(treat_comp_P)
treat_comp_P$p.adj <- as.numeric(as.character(treat_comp_P$p.adj))
treatP_compAdj <- p.adjust(treat_comp_P$p.adj, "BH")

#cut down to only significant ones
treat_sigcomp_P <- vector()
for (i in 1:length(treatP_compAdj)) {
  if (treatP_compAdj[i] < 0.05) {
    sig <- c(treatP_compAdj[i], treat_comp_P[i,5], treat_comp_P[i,6])
    treat_sigcomp_P <- rbind(treat_sigcomp_P, sig)
  }
}
colnames(treat_sigcomp_G) <- c("p_val", "comparison", "transcript")

#do this to compare gene vs protein
#use sum to zero for this comparison
all_comp_GP <- vector()
contrasts(residuals_all$type) <- contr.sum
for (i in 1:(length(residuals_all[1,])-3)) {
  last_3 <- ncol(residuals_all)
  sub <- residuals_all[,c(i, (last_3-2):last_3)]
  p_val <- summary(lm(sub[,1] ~ sub[,3] + sub[,2]))$coefficients[1,4]
  add <- c(p_val, colnames(residuals_all[i]))
  all_comp_GP <- rbind(all_comp_GP, add)
}

#multiple comparisons adjustment
allGP_compAdj <- p.adjust(all_comp_GP[,1], "BH")

#cut down to only significant ones
all_sigcomp_GP <- vector()
for (i in 1:length(allGP_compAdj)) {
  if (allGP_compAdj[i] < 0.05) {
    sig <- c(allGP_compAdj[i], all_comp_GP[i,2])
    all_sigcomp_GP <- rbind(all_sigcomp_GP, sig)
  }
}
colnames(all_sigcomp_GP) <- c("p_val", "transcript")

#do this to compare gene and protein within each treatment separately
#make separate data frames for each treatment and then do a welch's t-test
residuals_CGP <- rbind(residuals_CGPg, residuals_CGPp)
residuals_CGP$type <- as.factor(rep(c("gene", "protein"), each=5))
residuals_FAE <- rbind(residuals_FAEg, residuals_FAEp)
residuals_FAE$type <- as.factor(rep(c("gene", "protein"), each=9))
residuals_FAP <- rbind(residuals_FAPg, residuals_FAPp)
residuals_FAP$type <- as.factor(rep(c("gene", "protein"), each=9))
residuals_POH <- rbind(residuals_POHg, residuals_POHp)
residuals_POH$type <- as.factor(rep(c("gene", "protein"), each=10))
residuals_POL <- rbind(residuals_POLg, residuals_POLp)
residuals_POL$type <- as.factor(rep(c("gene", "protein"), each=8))

CGP_comp_GP <- vector()
for (i in 1:(length(residuals_CGP[1,])-1)) {
  sub <- residuals_CGP[,c(i, length(residuals_CGP[1,]))]
  gene_sub <- subset(sub, type == "gene")
  protein_sub <- subset(sub, type == "protein")
  p_val <- t.test(gene_sub[,1], protein_sub[,1])$p.value
  stat <- t.test(gene_sub[,1], protein_sub[,1])$statistic
  add <- c(p_val, stat, colnames(residuals_CGP[i]))
  CGP_comp_GP <- rbind(CGP_comp_GP, add)
}
CGP_comp_GP <- data.frame(CGP_comp_GP)
colnames(CGP_comp_GP) <- c("p_val", "t", "transcript")

#multiple comparisons adjustment
CGP_comp_GP$p_val <- as.numeric(as.character(CGP_comp_GP$p_val))
CGP_comp_GP$t <- as.numeric(as.character(CGP_comp_GP$t))
CGP_comp_GP$transcript <- as.character(CGP_comp_GP$transcript)
CGP_compAdj <- p.adjust(CGP_comp_GP[,1], "BH")

#cut down to only significant ones
CGP_sigcomp_GP <- vector()
for (i in 1:length(CGP_compAdj)) {
  if (CGP_compAdj[i] < 0.05) {
    sig <- c(CGP_compAdj[i], CGP_comp_GP$t[i], CGP_comp_GP$transcript[i])
    CGP_sigcomp_GP <- rbind(CGP_sigcomp_GP, sig)
  }
}

FAE_comp_GP <- vector()
for (i in 1:(length(residuals_FAE[1,])-1)) {
  sub <- residuals_FAE[,c(i, length(residuals_FAE[1,]))]
  gene_sub <- subset(sub, type == "gene")
  protein_sub <- subset(sub, type == "protein")
  p_val <- t.test(gene_sub[,1], protein_sub[,1])$p.value
  stat <- t.test(gene_sub[,1], protein_sub[,1])$statistic
  add <- c(p_val, stat, colnames(residuals_FAE[i]))
  FAE_comp_GP <- rbind(FAE_comp_GP, add)
}
FAE_comp_GP <- data.frame(FAE_comp_GP)
colnames(FAE_comp_GP) <- c("p_val", "t", "transcript")

#multiple comparisons adjustment
FAE_comp_GP$p_val <- as.numeric(as.character(FAE_comp_GP$p_val))
FAE_comp_GP$t <- as.numeric(as.character(FAE_comp_GP$t))
FAE_comp_GP$transcript <- as.character(FAE_comp_GP$transcript)
FAE_compAdj <- p.adjust(FAE_comp_GP[,1], "BH")

#cut down to only significant ones
FAE_sigcomp_GP <- vector()
for (i in 1:length(FAE_compAdj)) {
  if (FAE_compAdj[i] < 0.05) {
    sig <- c(FAE_compAdj[i], FAE_comp_GP$t[i], FAE_comp_GP$transcript[i])
    FAE_sigcomp_GP <- rbind(FAE_sigcomp_GP, sig)
  }
}

FAP_comp_GP <- vector()
for (i in 1:(length(residuals_FAP[1,])-1)) {
  sub <- residuals_FAP[,c(i, length(residuals_FAP[1,]))]
  gene_sub <- subset(sub, type == "gene")
  protein_sub <- subset(sub, type == "protein")
  p_val <- t.test(gene_sub[,1], protein_sub[,1])$p.value
  stat <- t.test(gene_sub[,1], protein_sub[,1])$statistic
  add <- c(p_val, stat, colnames(residuals_FAP[i]))
  FAP_comp_GP <- rbind(FAP_comp_GP, add)
}
FAP_comp_GP <- data.frame(FAP_comp_GP)
colnames(FAP_comp_GP) <- c("p_val", "t", "transcript")

#multiple comparisons adjustment
FAP_comp_GP$p_val <- as.numeric(as.character(FAP_comp_GP$p_val))
FAP_comp_GP$t <- as.numeric(as.character(FAP_comp_GP$t))
FAP_comp_GP$transcript <- as.character(FAP_comp_GP$transcript)
FAP_compAdj <- p.adjust(FAP_comp_GP[,1], "BH")

#cut down to only significant ones
FAP_sigcomp_GP <- vector()
for (i in 1:length(FAP_compAdj)) {
  if (FAP_compAdj[i] < 0.05) {
    sig <- c(FAP_compAdj[i], FAP_comp_GP$t[i], FAP_comp_GP$transcript[i])
    FAP_sigcomp_GP <- rbind(FAP_sigcomp_GP, sig)
  }
}

POH_comp_GP <- vector()
for (i in 1:(length(residuals_POH[1,])-1)) {
  sub <- residuals_POH[,c(i, length(residuals_POH[1,]))]
  gene_sub <- subset(sub, type == "gene")
  protein_sub <- subset(sub, type == "protein")
  p_val <- t.test(gene_sub[,1], protein_sub[,1])$p.value
  stat <- t.test(gene_sub[,1], protein_sub[,1])$statistic
  add <- c(p_val, stat, colnames(residuals_POH[i]))
  POH_comp_GP <- rbind(POH_comp_GP, add)
}
POH_comp_GP <- data.frame(POH_comp_GP)
colnames(POH_comp_GP) <- c("p_val", "t", "transcript")

#multiple comparisons adjustment
POH_comp_GP$p_val <- as.numeric(as.character(POH_comp_GP$p_val))
POH_comp_GP$t <- as.numeric(as.character(POH_comp_GP$t))
POH_comp_GP$transcript <- as.character(POH_comp_GP$transcript)
POH_compAdj <- p.adjust(POH_comp_GP[,1], "BH")

#cut down to only significant ones
POH_sigcomp_GP <- vector()
for (i in 1:length(POH_compAdj)) {
  if (POH_compAdj[i] < 0.05) {
    sig <- c(POH_compAdj[i], POH_comp_GP$t[i], POH_comp_GP$transcript[i])
    POH_sigcomp_GP <- rbind(POH_sigcomp_GP, sig)
  }
}
POH_sigcomp_GP <- data.frame(POH_sigcomp_GP)

POL_comp_GP <- vector()
for (i in 1:(length(residuals_POL[1,])-1)) {
  sub <- residuals_POL[,c(i, length(residuals_POL[1,]))]
  gene_sub <- subset(sub, type == "gene")
  protein_sub <- subset(sub, type == "protein")
  p_val <- t.test(gene_sub[,1], protein_sub[,1])$p.value
  stat <- t.test(gene_sub[,1], protein_sub[,1])$statistic
  add <- c(p_val, stat, colnames(residuals_POL[i]))
  POL_comp_GP <- rbind(POL_comp_GP, add)
}
POL_comp_GP <- data.frame(POL_comp_GP)
colnames(POL_comp_GP) <- c("p_val", "t", "transcript")

#multiple comparisons adjustment
POL_comp_GP$p_val <- as.numeric(as.character(POL_comp_GP$p_val))
POL_comp_GP$t <- as.numeric(as.character(POL_comp_GP$t))
POL_comp_GP$transcript <- as.character(POL_comp_GP$transcript)
POL_compAdj <- p.adjust(POL_comp_GP[,1], "BH")

#cut down to only significant ones
POL_sigcomp_GP <- vector()
for (i in 1:length(POL_compAdj)) {
  if (POL_compAdj[i] < 0.05) {
    sig <- c(POL_compAdj[i], POL_comp_GP$t[i], POL_comp_GP$transcript[i])
    POL_sigcomp_GP <- rbind(POL_sigcomp_GP, sig)
  }
}
POL_sigcomp_GP <- data.frame(POL_sigcomp_GP)

#rename all of the columns so they will match with the KEGG and MCODE files
#colnames(allTukey_sigcomp) <- c("diff", "lwr", "upr", "p_val", "comparison", "transcript")
#colnames(all_sigcomp_GP) <- c("p_val", "transcript")
#colnames(CGP_sigcomp_GP) <- c("p_val", "t_val", "transcript")
#colnames(FAE_sigcomp_GP) <- c("p_val", "t_val", "transcript")
#colnames(FAP_sigcomp_GP) <- c("p_val", "t_val", "transcript")
colnames(POH_sigcomp_GP) <- c("p_val", "t_val", "transcript")
colnames(POL_sigcomp_GP) <- c("p_val", "t_val", "transcript")
#CGP_sigcomp_GP$treatment <- "CGP"
#FAE_sigcomp_GP$treatment <- "FAE"
#FAP_sigcomp_GP$treatment <- "FAP"
POH_sigcomp_GP$treatment <- "POH"
POL_sigcomp_GP$treatment <- "POL"
allInd_sigcomp_GP <- rbind(POH_sigcomp_GP, POL_sigcomp_GP)

####look for significant terms in KEGG and MCODE pathways####
ddr <- "~/Documents/R_Data/musselRNASeq/compare_KEGG_MCODE/"
KEGG <- read.table(paste(ddr, "KEGG_names_allPathways.txt", sep=''), header=TRUE)
MCODE <- read.table(paste(ddr, "all transcripts_MCODE clusters.txt", sep=''),row.names=NULL)
colnames(MCODE)[6] <- "transcript"

#pull out everything from each sigcomp DF for each dataset
KEGG_allTukey <- merge(allTukey_sigcomp, KEGG, by="transcript")
KEGG_allGP <- merge(all_sigcomp_GP, KEGG, by="transcript")
KEGG_indGP <- merge(allInd_sigcomp_GP, KEGG, by="transcript")

KEGG_allTukey$sig_test <- "Tukey"
KEGG_allGP$sig_test <- "GP_sumZero"
KEGG_indGP$sig_test <- "GP_treat_t_test"

MCODE_allTukey <- merge(allTukey_sigcomp, MCODE, by="transcript")
MCODE_allGP <- merge(all_sigcomp_GP, MCODE, by="transcript")
MCODE_indGP <- merge(allInd_sigcomp_GP, MCODE, by="transcript")

MCODE_allTukey$sig_test <- "Tukey"
MCODE_allGP$sig_test <- "GP_sumZero"
MCODE_indGP$sig_test <- "GP_treat_t_test"

#export all of this
write.table(KEGG_allTukey, paste(ddr, "KEGG_allTukey_10.4.19.txt", sep=''), sep='\t', quote=F, row.names = F)
write.table(KEGG_allGP, paste(ddr, "KEGG_allGP_12.11.19.txt", sep=''), sep='\t', quote=F, row.names = F)
write.table(KEGG_indGP, paste(ddr, "KEGG_indGP_12.11.19.txt", sep=''), sep='\t', quote=F, row.names = F)


write.table(MCODE_allTukey, paste(ddr, "MCODE_allTukey_10.4.19.txt", sep=''), sep='\t', quote=F, row.names = F)
write.table(MCODE_allGP, paste(ddr, "MCODE_allGP_12.11.19.txt", sep=''), sep='\t', quote=F, row.names = F)
write.table(MCODE_indGP, paste(ddr, "MCODE_indGP_12.11.19.txt", sep=''), sep='\t', quote=F, row.names = F)

####explore some genes of interest####
#THIS IS FOR KEGG
#match enzyme names
enzymes <- read.table("~/Documents/R_Data/musselRNASeq/enzymes_list_9.4.19.txt", header=TRUE)
colnames(enzymes) <- c("transcript", "enzyme")
#apparently not all of the pathways are in enzyme names, so we can look at KO terms?
KO <- read.table("~/Documents/R_Data/musselRNASeq/KO_list_9.24.19.txt", header=TRUE)
colnames(KO) <- c("transcript", "KO")

#this reduces the KEGG_indGP dataset, so it's not a 1:1 comparison 
KEGG_indGP_enz <- merge(KEGG_indGP, enzymes, by="transcript")
KEGG_indGP_enz <- unique(KEGG_indGP_enz)

KEGG_indGP_KO <- merge(KEGG_indGP, KO, by="transcript")

#TCA cycle
m00020 <- KEGG_indGP_KO[which(KEGG_indGP_KO$KEGG_pathway %in% "m00020"),]
length(unique(m00020$KO))
#glycolysis
m00010 <- KEGG_indGP_KO[which(KEGG_indGP_KO$KEGG_pathway %in% "m00010"),]
length(unique(m00010$KO))
#metabolic processes
m01100 <- KEGG_indGP_KO[which(KEGG_indGP_KO$KEGG_pathway %in% "m01100"),]
length(unique(m01100$KO))
#proteasome
m03050 <- KEGG_indGP_KO[which(KEGG_indGP_KO$KEGG_pathway %in% "m03050"),]
length(unique(m03050$KO))
#ubiquitin mediated proteolysis
#this only has one significantly different KO term
m04120 <- KEGG_indGP_KO[which(KEGG_indGP_KO$KEGG_pathway %in% "m04120"),]
length(unique(m04120$KO))
#MAP-K signaling
m04010 <- KEGG_indGP_KO[which(KEGG_indGP_KO$KEGG_pathway %in% "m04010"),]
length(unique(m04010$KO))
#apoptosis
m04210 <- KEGG_indGP_KO[which(KEGG_indGP_KO$KEGG_pathway %in% "m04210"),]
length(unique(m04210$KO))

#MCODE
MCODE_indGP$type_cluster <- paste(MCODE_indGP$type, MCODE_indGP$MCODE_Cluster, sep="_")
gene1 <- subset(MCODE_indGP, type_cluster == "gene_1")
gene2 <- subset(MCODE_indGP, type_cluster == "gene_2")
gene3 <- subset(MCODE_indGP, type_cluster == "gene_3")
gene4 <- subset(MCODE_indGP, type_cluster == "gene_4")
gene5 <- subset(MCODE_indGP, type_cluster == "gene_5")
gene6 <- subset(MCODE_indGP, type_cluster == "gene_6")
gene7 <- subset(MCODE_indGP, type_cluster == "gene_7")
protein1 <- subset(MCODE_indGP, type_cluster == "protein_1")
protein2 <- subset(MCODE_indGP, type_cluster == "protein_2")
protein3 <- subset(MCODE_indGP, type_cluster == "protein_3")
protein4 <- subset(MCODE_indGP, type_cluster == "protein_4")
protein5 <- subset(MCODE_indGP, type_cluster == "protein_5")
protein6 <- subset(MCODE_indGP, type_cluster == "protein_6")
protein7 <- subset(MCODE_indGP, type_cluster == "protein_7")
protein8 <- subset(MCODE_indGP, type_cluster == "protein_8")


####look at directionality of MAD effects in significant genes####
MAD_CGPg <- MAD_function(CGP_gene)
MAD_FAEg <- MAD_function(FAE_gene)
MAD_FAPg <- MAD_function(FAP_gene)
MAD_POHg <- MAD_function(POH_gene)
MAD_POLg <- MAD_function(POL_gene)

MAD_CGPp <- MAD_function(CGP_protein)
MAD_FAEp <- MAD_function(FAE_protein)
MAD_FAPp <- MAD_function(FAP_protein)
MAD_POHp <- MAD_function(POH_protein)
MAD_POLp <- MAD_function(POL_protein)

MAD_all <- data.frame(transcript=gene_names, MAD_CGPg, MAD_FAEg, MAD_FAPg, MAD_POHg, MAD_POLg, MAD_CGPp, MAD_FAEp, MAD_FAPp, MAD_POHp, MAD_POLp)

m00020_MAD <- merge(MAD_all, m00020, by="transcript")
m00010_MAD <- merge(MAD_all, m00010, by="transcript")
m01100_MAD <- merge(MAD_all, m01100, by="transcript")
m03050_MAD <- merge(MAD_all, m03050, by="transcript")
m04010_MAD <- merge(MAD_all, m04010, by="transcript")
m04120_MAD <- merge(MAD_all, m04120, by="transcript")
m04210_MAD <- merge(MAD_all, m04210, by="transcript")

####test by doing the exact CBP analysis####

dis_GE <-dist(gene, method= "euclidean")   #appears to work and give same distance matrix as next line
#dis_GE <- vegdist(adGEt, method="euclidean", binary=FALSE)
mod_GE <- betadisper(dis_GE, treatments, type="median")
distance_GE <- mod_GE$distance

distance_GE <- data.frame(distance=distance_GE, treatment=treatments)

CGP_dev <- subset(distance_GE, treatment == "CGP")
FAE_dev <- subset(distance_GE, treatment == "FAE")
FAP_dev <- subset(distance_GE, treatment == "FAP")
POH_dev <- subset(distance_GE, treatment == "POH")
POL_dev <- subset(distance_GE, treatment == "POL")

#calculate MAD for each treatment and put into a DF
GE_MAD <- rep(0,5)
GE_MAD[1] <- median(abs(CGP_dev$distance))
GE_MAD[2] <- median(abs(FAE_dev$distance))
GE_MAD[3] <- median(abs(FAP_dev$distance))
GE_MAD[4] <- median(abs(POH_dev$distance))
GE_MAD[5] <- median(abs(POL_dev$distance))


#make a df of just the absolute deviations to put into the levene test
GE_MAD_distance <- data.frame(Outplant=c(rep("CG", length(CGP_dev$distance)), rep("field", length(FAE_dev$distance)),
                                         rep("field", length(FAP_dev$distance)), rep("outH", length(POH_dev$distance)),
                                         rep("outL", length(POL_dev$distance))),
                              Origin=c(rep("pro", length(CGP_dev$distance)), rep("exp", length(FAE_dev$distance)),
                                       rep("pro", length(FAP_dev$distance)), rep("pro", length(POH_dev$distance)),
                                       rep("pro", length(POL_dev$distance))),
                              Deviation=c(abs(CGP_dev$distance), abs(FAE_dev$distance), abs(FAP_dev$distance), 
                                          abs(POH_dev$distance), abs(POL_dev$distance)))


anova(mod_GE) #this says NOT SIGNIFICANT, BUT PAIRWISE COMPARISONS ARE SIGNIFICANT...
table_GE<-permutest(mod_GE, pairwise = TRUE, permutations = 9999)
ad_GE_pVal <- table_GE$pairwise$observed
ad_GE_letters <- multcompLetters(ad_GE_pVal, Letters=LETTERS)$Letters

#bar graph
#colors:
#only exposed and protected!
#this means that dark colors need to be exposed/protected low
cols <- c("gray68", "black")
ggplot(GE_MAD_distance, aes(x=Outplant, y=Deviation, color=Origin)) +
  geom_boxplot() +
  theme_bw() +
  #make the fill correspond to the two groupings within origin
  scale_color_manual(values=cols) +
  scale_x_discrete(limits=c("field", "CG", "outL", "outH")) +
  #add in the lettering
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

dis_PE <-dist(protein, method= "euclidean") 
mod_PE <- betadisper(dis_PE, treatments, type="median")
distance_PE <- mod_PE$distance

distance_PE <- data.frame(distance=distance_PE, treatment=treatments)

CGP_dev <- subset(distance_PE, treatment == "CGP")
FAE_dev <- subset(distance_PE, treatment == "FAE")
FAP_dev <- subset(distance_PE, treatment == "FAP")
POH_dev <- subset(distance_PE, treatment == "POH")
POL_dev <- subset(distance_PE, treatment == "POL")

#calculate MAD for each treatment and put into a DF
PE_MAD <- rep(0,5)
PE_MAD[1] <- median(abs(CGP_dev$distance))
PE_MAD[2] <- median(abs(FAE_dev$distance))
PE_MAD[3] <- median(abs(FAP_dev$distance))
PE_MAD[4] <- median(abs(POH_dev$distance))
PE_MAD[5] <- median(abs(POL_dev$distance))

#make a df of just the absolute deviations to put into the levene test
PE_MAD_distances <- data.frame(Outplant=c(rep("CG", length(CGP_dev$distance)), rep("field", length(FAE_dev$distance)),
                                          rep("field", length(FAP_dev$distance)), rep("outH", length(POH_dev$distance)),
                                          rep("outL", length(POL_dev$distance))),
                               Origin=c(rep("pro", length(CGP_dev$distance)), rep("exp", length(FAE_dev$distance)),
                                        rep("pro", length(FAP_dev$distance)), rep("pro", length(POH_dev$distance)),
                                        rep("pro", length(POL_dev$distance))),
                               Deviation=c(abs(CGP_dev$distance), abs(FAE_dev$distance), abs(FAP_dev$distance), 
                                           abs(POH_dev$distance), abs(POL_dev$distance)))

anova(mod_PE) #this says SIGNIFICANT, BUT PAIRWISE COMPARISONS ARE SIGNIFICANT...
table_PE<-permutest(mod_PE, pairwise = TRUE, permutations = 9999)
ad_PE_pVal <- table_PE$pairwise$observed
ad_PE_letters <- multcompLetters(ad_PE_pVal, Letters=LETTERS)$Letters


#bar graph
#colors:
#only exposed and protected!
#this means that dark colors need to be exposed/protected low (?)
cols <- c("gray68", "black")
ggplot(PE_MAD_distances, aes(x=Outplant, y=Deviation, color=Origin)) +
  geom_boxplot() +
  theme_bw() +
  #make the fill correspond to the two groupings within origin
  scale_color_manual(values=cols) +
  scale_x_discrete(limits=c("field", "CG", "outL", "outH")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

####look at relative expression in gene vs protein for KEGG components####
#since it is gene minus protein, a negative value means that protein has higher expression

CGP_gene <- sapply(CGP_gene, function(x) as.numeric(as.character(x)))
CGP_meanG <- colMeans(CGP_gene)
CGP_protein <- sapply(CGP_protein, function(x) as.numeric(as.character(x)))
CGP_meanP <- colMeans(CGP_protein)
CGP_diff <- CGP_meanG - CGP_meanP

FAE_gene <- sapply(FAE_gene, function(x) as.numeric(as.character(x)))
FAE_meanG <- colMeans(FAE_gene)
FAE_protein <- sapply(FAE_protein, function(x) as.numeric(as.character(x)))
FAE_meanP <- colMeans(FAE_protein)
FAE_diff <- FAE_meanG - FAE_meanP

FAP_gene <- sapply(FAP_gene, function(x) as.numeric(as.character(x)))
FAP_meanG <- colMeans(FAP_gene)
FAP_protein <- sapply(FAP_protein, function(x) as.numeric(as.character(x)))
FAP_meanP <- colMeans(FAP_protein)
FAP_diff <- FAP_meanG - FAP_meanP

POH_gene <- sapply(POH_gene, function(x) as.numeric(as.character(x)))
POH_meanG <- colMeans(POH_gene)
POH_protein <- sapply(POH_protein, function(x) as.numeric(as.character(x)))
POH_meanP <- colMeans(POH_protein)
POH_diff <- POH_meanG - POH_meanP

POL_gene <- sapply(POL_gene, function(x) as.numeric(as.character(x)))
POL_meanG <- colMeans(POL_gene)
POL_protein <- sapply(POL_protein, function(x) as.numeric(as.character(x)))
POL_meanP <- colMeans(POL_protein)
POL_diff <- POL_meanG - POL_meanP

m00010_diff_POH <- POH_diff[which(names(POH_diff) %in% m00010$transcript)]
m00010_diff_POL <- POL_diff[which(names(POL_diff) %in% m00010$transcript)]

m00020_diff_POH <- POH_diff[which(names(POH_diff) %in% m00020$transcript)]
m00020_diff_POL <- POL_diff[which(names(POL_diff) %in% m00020$transcript)]

m03050_diff_POH <- POH_diff[which(names(POH_diff) %in% m03050$transcript)]
m03050_diff_POL <- POL_diff[which(names(POL_diff) %in% m03050$transcript)]

m04010_diff_POH <- POH_diff[which(names(POH_diff) %in% m04010$transcript)]
m04010_diff_POL <- POL_diff[which(names(POL_diff) %in% m04010$transcript)]

m04120_diff_POH <- POH_diff[which(names(POH_diff) %in% m04120$transcript)]
m04120_diff_POL <- POL_diff[which(names(POL_diff) %in% m04120$transcript)]

m04210_diff_POH <- POH_diff[which(names(POH_diff) %in% m04210$transcript)]
m04210_diff_POL <- POL_diff[which(names(POL_diff) %in% m04210$transcript)]

#export for cytoscape here!
#export MAD_all because the clusters are formed within the whole network
#put in the matching transcript to cytoscape number for transcripts
names_cytoscape <- read.table("~/match_cytoscapeNum_transcriptID.txt", header=TRUE)
colnames(names_cytoscape) <- c("Cytoscape_ID", "transcript")
GO_terms <- 
MAD_all_cytoscape <- merge(MAD_all, names_cytoscape, by = "transcript")
write.table(MAD_all_cytoscape, paste(ddr, "MAD_all1521.txt", sep=''), sep='\t', quote=F, row.names = F)
#write.table(MAD_gene1, paste(ddr, "MAD_Cytoscape/MAD_gene1.txt", sep=''), sep='\t', quote=F, row.names = F)

MAD_gene1$transcript <- NULL
MAD_gene2$transcript <- NULL
MAD_gene3$transcript <- NULL
MAD_gene4$transcript <- NULL
MAD_gene5$transcript <- NULL
MAD_gene6$transcript <- NULL
MAD_gene7$transcript <- NULL
MAD_protein1$transcript <- NULL
MAD_protein2$transcript <- NULL
MAD_protein3$transcript <- NULL
MAD_protein4$transcript <- NULL
MAD_protein5$transcript <- NULL
MAD_protein6$transcript <- NULL
MAD_protein7$transcript <- NULL
MAD_protein8$transcript <- NULL

MAD_gene1_avg <- colMeans(MAD_gene1)
MAD_gene2_avg <- colMeans(MAD_gene2)
MAD_gene3_avg <- colMeans(MAD_gene3)
MAD_gene4_avg <- colMeans(MAD_gene4)
MAD_gene5_avg <- colMeans(MAD_gene5)
MAD_gene6_avg <- colMeans(MAD_gene6)
MAD_gene7_avg <- colMeans(MAD_gene7)
MAD_protein1_avg <- colMeans(MAD_protein1)
MAD_protein2_avg <- colMeans(MAD_protein2)
MAD_protein3_avg <- colMeans(MAD_protein3)
MAD_protein4_avg <- colMeans(MAD_protein4)
MAD_protein5_avg <- colMeans(MAD_protein5)
MAD_protein6_avg <- colMeans(MAD_protein6)
MAD_protein7_avg <- colMeans(MAD_protein7)
MAD_protein8_avg <- colMeans(MAD_protein8)

MCODE_avgMADoverall <- c(mean(MAD_gene1_avg[1:5]), mean(MAD_gene1_avg[6:10]), 
                         mean(MAD_gene2_avg[1:5]), mean(MAD_gene2_avg[6:10]),
                         mean(MAD_gene3_avg[1:5]), mean(MAD_gene3_avg[6:10]),
                         mean(MAD_gene4_avg[1:5]), mean(MAD_gene4_avg[6:10]),
                         mean(MAD_gene5_avg[1:5]), mean(MAD_gene5_avg[6:10]),
                         mean(MAD_gene6_avg[1:5]), mean(MAD_gene6_avg[6:10]),
                         mean(MAD_gene7_avg[1:5]), mean(MAD_gene7_avg[6:10]),
                         mean(MAD_protein1_avg[1:5]), mean(MAD_protein1_avg[6:10]), 
                         mean(MAD_protein2_avg[1:5]), mean(MAD_protein2_avg[6:10]),
                         mean(MAD_protein3_avg[1:5]), mean(MAD_protein3_avg[6:10]),
                         mean(MAD_protein4_avg[1:5]), mean(MAD_protein4_avg[6:10]),
                         mean(MAD_protein5_avg[1:5]), mean(MAD_protein5_avg[6:10]),
                         mean(MAD_protein6_avg[1:5]), mean(MAD_protein6_avg[6:10]),
                         mean(MAD_protein7_avg[1:5]), mean(MAD_protein7_avg[6:10]),
                         mean(MAD_protein8_avg[1:5]), mean(MAD_protein8_avg[6:10]))

MCODE_MAD_plot <- data.frame(MAD=MCODE_avgMADoverall, pathway=rep(c("gene1", "gene2", "gene3", "gene4", 
                                                                    "gene5", "gene6", "gene7", "protein1", 
                                                                    "protein 2", "protein 3", "protein 4", 
                                                                    "protein 5", "protein 6", "protein 7",
                                                                    "protein 8"), each=2),
                            type=rep(rep(c("gene", "protein"), 2), 15))

#MAD POH by pathway
MCODE_avg_MAD_POH <- c(MAD_gene1_avg[4], MAD_gene1_avg[9],
                       MAD_gene2_avg[4], MAD_gene2_avg[9],
                       MAD_gene3_avg[4], MAD_gene3_avg[9],
                       MAD_gene4_avg[4], MAD_gene4_avg[9],
                       MAD_gene5_avg[4], MAD_gene5_avg[9],
                       MAD_gene6_avg[4], MAD_gene6_avg[9],
                       MAD_gene7_avg[4], MAD_gene7_avg[9],
                       MAD_protein1_avg[4], MAD_protein1_avg[9],
                       MAD_protein2_avg[4], MAD_protein2_avg[9],
                       MAD_protein3_avg[4], MAD_protein3_avg[9],
                       MAD_protein4_avg[4], MAD_protein4_avg[9],
                       MAD_protein5_avg[4], MAD_protein5_avg[9],
                       MAD_protein6_avg[4], MAD_protein6_avg[9],
                       MAD_protein7_avg[4], MAD_protein7_avg[9],
                       MAD_protein8_avg[4], MAD_protein8_avg[9])

MCODE_MAD_POH_plot <- data.frame(MAD=MCODE_avg_MAD_POH, pathway=rep(c("gene1", "gene2", "gene3", "gene4", 
                                                                     "gene5", "gene6", "gene7", "protein1", 
                                                                     "protein 2", "protein 3", "protein 4", 
                                                                     "protein 5", "protein 6", "protein 7",
                                                                     "protein 8"), each=2),
                                 type=rep(rep(c("gene", "protein"), 2), 15))

ggplot(MCODE_MAD_plot, aes(x=pathway, y=MAD, group=type, fill=type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw() +
  ylim(0, 0.25) +
  xlab("Pathway/Cluster") +
  ylab("MAD (variance measure)") +
  theme(axis.title = element_text(size=15))

ggplot(MCODE_MAD_POH_plot, aes(x=pathway, y=MAD, group=type, fill=type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw()


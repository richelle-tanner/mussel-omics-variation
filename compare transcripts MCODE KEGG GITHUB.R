#compare transcripts in the MCODE clusters and KEGG pathways
#files come from "eigenvalue_covariance_KEGG.R" and "MCODE output analysis.R"

#import data
ddr <- "~/Documents/R_Data/musselRNASeq/compare_KEGG_MCODE/"
KEGG <- read.table(paste(ddr, "KEGG_names_allPathways.txt", sep=''), header=TRUE)
MCODE <- read.table(paste(ddr, "all transcripts_MCODE clusters.txt", sep=''),row.names=NULL)

MCODE$type_cluster <- paste(MCODE$type, MCODE$MCODE_Cluster, sep='_')
MCODE <- MCODE[,c(6,8)]

#are there any duplicates among pathways/clusters? - yes, about a few hundred
KEGG_unique <- data.frame(table(KEGG$transcript))
MCODE_unique <- data.frame(table(MCODE$TRUE.))

#what is the overlap between KEGG and MCODE
KEGG_names <- as.character(KEGG_unique$Var1)
MCODE_names <- as.character(MCODE_unique$Var1)

overlap <- intersect(KEGG_names, MCODE_names)

#put the frequency columns for KEGG and MCODE to see whether any of these show up in multiple clusters/pathways
overlap_freq <- data.frame(Names = overlap,
                           KEGG_freq = KEGG_unique[which(KEGG_unique$Var1 %in% overlap),2],
                           MCODE_freq = MCODE_unique[which(MCODE_unique$Var1 %in% overlap),2])

colnames(MCODE) <- c("Names", "type_cluster")
colnames(KEGG) <- c("Names", "pathway")

overlap_freq <- merge(overlap_freq, MCODE, by="Names")
overlap_freq <- merge(overlap_freq, KEGG, by="Names")

####anova-specific analyses####
MCODE_anova <- read.delim("~/all_geneCluster_anova_sigcomp_2.25.20.txt", header=TRUE, sep='\t')
KEGG_anova <- read.delim("~/KEGG_allPath_sigcomp_treat_2.25.20.txt", header=TRUE, sep='\t')

#look specifically for transcripts that are higher in stressful treatments (just take POH for now)
table(MCODE_anova$treat1)
table(MCODE_anova$treat2)
#treat 1 > 2 if effect size is negative
#treat 1 < 2 if effect size is positive
#remember the labels got switched in the function code so its the opposite of what you would expect
MCODE_anova_POH1 <- subset(MCODE_anova, treat1 == "POH")
MCODE_anova_POH1h <- subset(MCODE_anova_POH1, effect_size < 0)
MCODE_anova_POH2 <- subset(MCODE_anova, treat2 == "POH")
MCODE_anova_POH2h <- subset(MCODE_anova_POH2, effect_size > 0)
MCODE_anova_POH <- rbind(MCODE_anova_POH1, MCODE_anova_POH2)
MCODE_anova_POHh <- rbind(MCODE_anova_POH1h, MCODE_anova_POH2h) #395 of 534 observations involving POH show higher variation in POH

#it isn't POL for MCODE
MCODE_anova_POL1 <- subset(MCODE_anova_POH, treat1 == "POL")
MCODE_anova_POL1h <- subset(MCODE_anova_POL1, effect_size < 0) #something here
MCODE_anova_POL2 <- subset(MCODE_anova_POH, treat2 == "POL")
MCODE_anova_POL2h <- subset(MCODE_anova_POL2, effect_size > 0)

MCODE_anova_FAP1 <- subset(MCODE_anova_POH, treat1 == "FAP")
MCODE_anova_FAP1h <- subset(MCODE_anova_FAP1, effect_size < 0)
MCODE_anova_FAP2 <- subset(MCODE_anova_POH, treat2 == "FAP")
MCODE_anova_FAP2h <- subset(MCODE_anova_FAP2, effect_size > 0) #something here

MCODE_anova_FAE1 <- subset(MCODE_anova_POH, treat1 == "FAE")
MCODE_anova_FAE1h <- subset(MCODE_anova_FAE1, effect_size < 0)
MCODE_anova_FAE2 <- subset(MCODE_anova_POH, treat2 == "FAE")
MCODE_anova_FAE2h <- subset(MCODE_anova_FAE2, effect_size > 0) #something here

table(KEGG_anova$treat1)
table(KEGG_anova$treat2)
#treat 1 > 2 if effect size is negative
#treat 1 < 2 if effect size is positive
#remember the labels got switched in the function code so its the opposite of what you would expect
KEGG_anova_POH1 <- subset(KEGG_anova, treat1 == "POH")
KEGG_anova_POH1h <- subset(KEGG_anova_POH1, effect_size < 0)
KEGG_anova_POH2 <- subset(KEGG_anova, treat2 == "POH")
KEGG_anova_POH2h <- subset(KEGG_anova_POH2, effect_size > 0)
KEGG_anova_POHh <- rbind(KEGG_anova_POH1h, KEGG_anova_POH2h) #33 of 43 observation involving POH show higher variation in POH (all of the other ones have higher variation in POL)

KEGG_anova_POL1 <- subset(KEGG_anova_POH, treat1 == "POL")
KEGG_anova_POL1h <- subset(KEGG_anova_POL1, effect_size < 0)
KEGG_anova_POL2 <- subset(KEGG_anova_POH, treat2 == "POL")
KEGG_anova_POL2h <- subset(KEGG_anova_POL2, effect_size > 0) #the rest are higher in POL (of the ones from POH)



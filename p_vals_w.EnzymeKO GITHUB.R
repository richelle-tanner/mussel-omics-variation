#puts in the enzyme and KO codes so you can look up gene names online
#end section also has some superficial analyses needed for the text of the paper
#this end section includes looking at gene pairs

#Author: R Tanner
#Last edit: 2.26.20

enzymes <- read.delim("~/Documents/R_Data/musselRNASeq/enzymes_list_9.4.19.txt", header=TRUE)
KO <- read.delim("~/Documents/R_Data/musselRNASeq/KO_list_9.24.19.txt", header=TRUE)
BLAST <- read.delim("~/Documents/R_Data/musselRNASeq/final_files/Full_Blast2GO_table_tab.txt", header=TRUE, sep='\t', fill=TRUE)

colnames(BLAST)[1] <- "transcript_ID"

combined <- merge(KO, enzymes, by = "transcript_ID", all.x=TRUE)
combined_blastKO <- merge(BLAST, KO, by = "transcript_ID", all.x=TRUE)

p_vals <- read.table("~/Documents/R_Data/musselRNASeq/compare_KEGG_MCODE/all_geneCluster_GP_sigcomp_2.25.20.txt", header=TRUE, sep='\t')
colnames(p_vals)[4] <- "transcript_ID"

combined_p <- merge(p_vals, combined_blastKO, by = "transcript_ID", all.x=TRUE)

write.table(combined_p, "~/Documents/R_Data/musselRNASeq/annotated_sigGeneCluster_GP_2.26.20.txt", quote=F, row.names = F, sep='\t')

#for the gene by gene comparison
p_vals <- read.table("~/Documents/R_Data/musselRNASeq/compare_KEGG_MCODE/KEGG_allPath_sigcomp_genes_2.25.20.txt", header=TRUE)
#separate out the second gene column
gene2_sep <- data.frame(transcript_ID=as.character(p_vals$gene2),p_val=p_vals$p_val)

colnames(p_vals)[4] <- "transcript_ID"

combined_p <- merge(p_vals, combined_blastKO, by = "transcript_ID", all.x=TRUE)
gene2_p <- merge(gene2_sep,combined_blastKO, by = "transcript_ID", all.x=TRUE)

write.table(combined_p, "~/Documents/R_Data/musselRNASeq/annotated_sigKEGGgene1_2.26.20.txt", quote=F, row.names = F, sep='\t')


####looking at gene pairs analysis####

#import data
KEGGgenePairs <- read.delim("~/Documents/R_Data/musselRNASeq/annotated_sigKEGGgenes_2.26.20.txt", header=TRUE)

#find gene pairs that are in gene and protein datasets
KEGGgenePairs_sort <- KEGGgenePairs[order(KEGGgenePairs$gene1, KEGGgenePairs$gene2),]

KEGGgenePairs_m00010 <- subset(KEGGgenePairs_noEnz, pathway == "m00010")
KEGGgenePairs_m00010pro <- subset(KEGGgenePairs_m00010, type == "protein")
KEGGgenePairs_m00010gen <- subset(KEGGgenePairs_m00010, type == "gene")
mean(KEGGgenePairs_m00010pro$effect_size)
sd(KEGGgenePairs_m00010pro$effect_size)/sqrt(length(KEGGgenePairs_m00010pro$effect_size))
mean(KEGGgenePairs_m00010gen$effect_size)
sd(KEGGgenePairs_m00010gen$effect_size)/sqrt(length(KEGGgenePairs_m00010gen$effect_size))

#how many times a single gene appears
genes_list <- c(as.character(KEGGgenePairs$gene1), as.character(KEGGgenePairs$gene2))
genes_list_table <- data.frame(table(genes_list))
genes_list_table <- genes_list_table[order(genes_list_table$Freq),]

#by treatment
table(KEGGgenePairs$treatment)
table(KEGGgenePairs$pathway)
table(KEGGgenePairs$type)

####looking at GP analysis####
#import data
KEGG_GP <- read.delim("~/Documents/R_Data/musselRNASeq/annotated_sigKEGG_GP_2.26.20.txt", header=TRUE)
MCODE_GP <- read.delim("~/Documents/R_Data/musselRNASeq/annotated_sigGeneCluster_GP_2.26.20.txt", header=TRUE, sep='\t') #there is no protein one

MCODE_GP$cluster <- substr(MCODE_GP$cluster_treatment, 1, 9)
MCODE_GP$treatment <- substr(MCODE_GP$cluster_treatment, 11, 13)

#mean effect size by transcript/protein
table(KEGG_GP$pathway)
mean(KEGG_GP[(KEGG_GP$effect_size < 0), 4])
sd(KEGG_GP[(KEGG_GP$effect_size < 0), 4])/sqrt(length(KEGG_GP[(KEGG_GP$effect_size < 0), 4]))
mean(KEGG_GP[(KEGG_GP$effect_size > 0), 4])
sd(KEGG_GP[(KEGG_GP$effect_size > 0), 4])/sqrt(length(KEGG_GP[(KEGG_GP$effect_size > 0), 4]))

#just for TCA cycle
TCA_GP <- subset(KEGG_GP, pathway == "m00020")
mean(TCA_GP[(TCA_GP$effect_size < 0), 4])
sd(TCA_GP[(TCA_GP$effect_size < 0), 4])/sqrt(length(TCA_GP[(TCA_GP$effect_size < 0), 4]))

#mean effect size by transcript/protein
table(MCODE_GP$cluster)
table(MCODE_GP$treatment)
mean(MCODE_GP[(MCODE_GP$effect_size < 0), 4])
mean(MCODE_GP[(MCODE_GP$effect_size > 0), 4])

MCODE_GP_POH <- subset(MCODE_GP, treatment == "POH")


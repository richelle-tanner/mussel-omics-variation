#file that uses imputation method to replace zeros in transcript and protein expression
#updated 6.3.20 to also replace zeros for transcript expression in addition to protein expression

####protein expression####
ddr <- "~/Documents/Cytoscape_projects/M_Cali/"
protein <- read.delim(paste(ddr, "reduced_PEmatrix_by_matching_genes.txt", sep=''))

protein <- data.frame(sapply(protein, function(x) as.numeric(as.character(x))))

protein[protein == 0] <- NA


#use impute function
#expression matrix with genes in the rows, samples in the columns

#get rid of two individuals and two genes that aren't in the final set
protein_old <- read.delim("~/protein_exp_noScale_noTrans.txt")

protein_new <- protein[which(row.names(protein) %in% row.names(gene_llsImpute_export)),]
protein_new <- protein_new[,which(colnames(protein_new) %in% colnames(protein_old))]

protein_genenames <- row.names(protein_new)
protein_indnames <- colnames(protein_new)

#get rid of two transcripts that had all zeros in the gene expression dataset
protein_new <- data.frame(sapply(protein_new, function(x) as.numeric(as.character(x))))

protein_new[protein_new == 0] <- NA

library(pcaMethods)
#local least squares impute method (llsImpute)
#requires samples (row) x genes (col)
protein_newT <- t(protein_new)

#recommended to use the nni() wrapper function, but we want to see what is going on
protein_llsImpute <- llsImpute(protein_newT, k=10, center=FALSE, completeObs=TRUE, correlation="pearson", allVariables=TRUE)
protein_llsImpute_export <- data.frame(protein_llsImpute@completeObs)
protein_llsImpute_export <- data.frame(t(protein_llsImpute_export))
row.names(protein_llsImpute_export) <- protein_genenames
colnames(protein_llsImpute_export) <- protein_indnames


####gene expression####
gene <- read.delim("~/Documents/R_Data/musselRNASeq/final_files/gene_exp_noScale_noTrans_12.11.19.txt")

#remove columns that have ALL zeros
gene[gene == 0] <- NA

#local least squares impute method (llsImpute)
#requires samples (row) x genes (col)
geneT <- t(gene)
geneF <- geneT[,colSums(is.na(geneT))<nrow(geneT)] #there are two genes we need to get rid of in the protein dataset!

#recommended to use the nni() wrapper function, but we want to see what is going on
gene_llsImpute <- llsImpute(geneF, k=10, center=FALSE, completeObs=TRUE, correlation="pearson", allVariables=TRUE)
gene_llsImpute_export <- data.frame(gene_llsImpute@completeObs)
gene_llsImpute_export <- data.frame(t(gene_llsImpute_export))
row.names(gene_llsImpute_export) <- colnames(geneF)
colnames(gene_llsImpute_export) <- row.names(geneF)


write.table(gene_llsImpute_export, "~/gene_matrixllsImpute.txt", quote=F, sep='\t')
write.table(protein_llsImpute_export, "~/protein_matrixllsImpute.txt", quote=F, sep='\t')

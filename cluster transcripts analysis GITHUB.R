#use the results from ClusterViz to get a better look at what is in each cluster, annotation-wise
#only look at most complex clusters

####import the data####
ddr <- "~/Documents/R_Data/musselRNASeq/"
full_gene <- read.csv(paste(ddr, "full_gene_clusters_MCODE.csv", sep=''))
transc_ann <- read.table(paste(ddr, "final_files/reduced_annotation_full_transcriptome.txt", sep=''), fill=TRUE, stringsAsFactors = FALSE)
gene_adj_matrix <- read.table(paste(ddr, "full_gene_network_adjacency_matrix.txt", sep=''), header=TRUE)
gene_matrix <- read.table(paste(ddr, "final_files/reduced_GEmatrix_by_matching_proteins.txt", sep=''), fill=TRUE)
protein_matrix <- read.table(paste(ddr, "final_files/reduced_PEmatrix_by_matching_genes.txt", sep=''), fill=TRUE)

####cut down gene and protein matrices####
#this will match the matrices that were used to create the networks (in correlation network creation.R)
#input is same as output: gene_matrix and protein_matrix

#get rid of two rows that have no variance in the gene dataset
gene1 <- gene_matrix[-c(348,536),]
#get rid of those same ones in the protein dataset
#unfortunately they're not the same row numbers
eliminate <- c("Velvet_k35_Velvet_k35_Locus_1520_Transcript_1114_Confidence_0.589_Length_3521", "TRINITY_DN20953_c0_g2_i1")
elimRowNum <- which(rownames(protein_matrix) %in% eliminate) 
protein1 <- protein_matrix[-c(elimRowNum),]

#only use individuals that overlap between datasets
gene_sim <- gene1[,order(colnames(gene1))]
protein_sim <- protein1[,order(colnames(protein1))]
protein_matrix <- protein_sim[,-c(6,12)]
gene_matrix <- gene_sim[,-c(3,11,19,34,36,40,41,42,44,49)]


####match row numbers and names####
#need to match the row numbers from the "correlation network creation.r" file to the transcript IDs
#extract row names and put them with numbers in a new data frame
row_name_gene <- data.frame(names=rownames(gene_adj_matrix), numbers=seq(1,length(gene_adj_matrix$X1)))


####first cluster####
cluster1 <- full_gene[1,]
cluster1 <- cluster1[,-c(1:3)]
rows_selected <- which(row_name_gene$numbers %in% cluster1)
cluster1_genes <- row_name_gene[rows_selected,]

#put in nr database return (V9), KEGG term (V15), and GO term (V16)
rows_selected_ann <- which(transc_ann$V3 %in% cluster1_genes$names)
cluster1_ann <- transc_ann[rows_selected_ann, c(3, 9, 15, 16)]

#get the reduced gene matrix so we can summarize DV and DE
rows_selected_matrix <- which(rownames(gene_matrix) %in% cluster1_ann$V3)
cluster1_matrix <- gene_matrix[rows_selected_matrix,]

#alphabetical order 
cluster1_matrix <- cluster1_matrix[order(rownames(cluster1_matrix)),]

#calculate DV and DE
gene_names_cluster1 <- rownames(cluster1_matrix)
cluster1_matrix <- data.frame(lapply(cluster1_matrix, function(x) as.numeric(as.character(x))))
cluster1_matrixT <- scale(t(cluster1_matrix))

treatments <- c("FAE", "POL", "FAE", "FAP", "POH", "CGP", "POH", "FAE", "FAP", "FAE", "POL", "FAP", "FAP", "FAE",
                "FAP", "POH", "CGP", "CGP", "POL", "FAE", "FAP", "POH", "POH", "POL", "CGP", "POL", "POL", "POH",
                "POH", "POH", "POL", "POL", "FAE", "FAP", "POH", "CGP", "FAE", "FAP", "POH", "FAE", "FAP")

cluster1_matrixT <- cbind(treatments, cluster1_matrixT)

#empty MAD vectors, columns are treatment, rows are genes
cluster1_MADg <- matrix(ncol=5, nrow=(length(cluster1_matrixT[1,])-1))
cluster1_meang <- matrix(ncol=5, nrow=(length(cluster1_matrixT[1,])-1))
#treatment vector to step through
treat_step <- c("POL", "POH", "CGP", "FAE", "FAP")
for (i in 1:length(treat_step)) {
  #separate by treatment
  subset_g <- subset(cluster1_matrixT, treatments == treat_step[i])
  subset_g <- subset_g[,-1]
  #calculate MAD per gene/protein
  for (j in 1:(length(cluster1_matrixT[1,])-1)) {
    sub_gG <- as.numeric(as.character(subset_g[,j]))
    #find median for each gene
    median_g <- median(sub_gG)
    #calculate deviations
    deviations_g <- sub_gG - median_g
    #calculate and store MAD
    cluster1_MADg[j,i] <- median(abs(deviations_g))
    cluster1_meang[j,i] <- mean(sub_gG)
  }
}

#put names back in
cluster1_MADg <- data.frame(cluster1_MADg)
cluster1_MADg$name <- gene_names_cluster1

cluster1_meang <- data.frame(cluster1_meang)
cluster1_meang$name <- gene_names_cluster1

#put the annotations in
colnames(cluster1_ann) <- c("name", "BLASTX_nr", "KEGG", "GO")
colnames(cluster1_MADg) <- c("POL", "POH", "CGP", "FAE", "FAP", "name")
cluster1_MADg_ann <- merge(cluster1_MADg, cluster1_ann, by = "name")

#cut down GO column to the first GO term
cluster1_MADg_ann$BLASTX_nr <- substring(cluster1_MADg_ann$BLASTX_nr, 1, 14)

#cut down BLAST column to the first term
cluster1_MADg_ann$GO <- substring(cluster1_MADg_ann$GO, 1, 10)

#KEGG terms
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

cluster1_MADg_ann$KEGG <- substrRight(cluster1_MADg_ann$KEGG, 9)


table(cluster1_MADg_ann$GO)
table(cluster1_MADg_ann$BLASTX_nr)

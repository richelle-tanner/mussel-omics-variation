# mussel-omics-variation
R scripts and Cytoscape files associated with the Tanner, Gleason, &amp; Dowd manuscript exploring the role of inter-individual variation in gene and protein networks of Mytilus mussels under heat stress.

#expression replace zeros.R ---- data cleaning of raw gene and protein expression datasets using impute function

####Scripts and files associated with subnetwork creation in Cytoscape application:
#correlation network creation_bootstrap_p_vals.R ---- takes in raw gene and protein expression datasets (all treatments) to generate sif files for Cytoscape import
#allTreat_geneIG_6.16.20.sif ---- created by above R file for gene expression
#allTreat_proteinIG_6.16.20.sif ---- created by above R file for protein expression
#MCODE_6.16.20.cys ---- Cytoscape working file with the above sif files imported and MCODE analysis run

####Subnetwork analyses and comparisons between KEGG and MCODE
#MCODE output analysis.R ---- post-Cytoscape identification of gene functions and MCODE cluster (subnetwork) characteristics 
#compare transcripts MCODE KEGG.R ---- compares the transcripts present in MCODE clusters and selected KEGG pathways in gene and protein expression datasets

####Median Absolute Deviation calculations:
#bootstrapping_MAD.R ---- bootstraps MAD values for all pathways and clusters
#MAD p values from residuals.R ---- calculates p values for all MAD values for all pathways and clusters

####Eigenvalue covariance calculations:
#eigenvalue_covariance_KEGG.R ---- calculates eigenvalues SD and ICV for all KEGG pathways
#eigenvalue_covariance_MCODE.R ---- calculates eigenvalues SD and ICV for all MCODE clusters
#pathway rank analyses.R ---- conducts a rank analysis of eigenvalues SD and ICV for all pathways and clusters

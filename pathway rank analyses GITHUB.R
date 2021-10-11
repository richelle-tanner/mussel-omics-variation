#script for doing mean rank for eigenvalue statistics

library(crank)

#import tables S4-6 in one file
allSDs <- read.delim("~/EigenSDsALL.csv", sep=',')
pathways <- as.character(unique(allSDs$Pathway))

####mean rank across pathways ####
#SDs first
#transform SDs into ranks by pathway
SDranks <- data.frame()
for (i in 1:length(pathways)) {
  pathway <- subset(allSDs, Pathway == pathways[i])
  ranks <- rank(-pathway$SD, na.last = TRUE) #puts NAs last, wanted to keep as NA but meanrank isn't performing
  add <- data.frame(pathway$Pathway, pathway$Treatment, pathway$SD, ranks)
  SDranks <- rbind(SDranks, add)
}
colnames(SDranks) <- c("Pathway", "Treatment", "SD", "Rank")

#matrix should be pathway x treatment (calculates means by column)
PCG_transcript <- subset(SDranks, Treatment == "PCG-transcript exp.")
EFA_transcript <- subset(SDranks, Treatment == "EFA-transcript exp.")
PFA_transcript <- subset(SDranks, Treatment == "PFA-transcript exp.")
POL_transcript <- subset(SDranks, Treatment == "POL-transcript exp.")
POH_transcript <- subset(SDranks, Treatment == "POH-transcript exp.")
PCG_protein <- subset(SDranks, Treatment == "PCG-protein exp.")
EFA_protein <- subset(SDranks, Treatment == "EFA-protein exp.")
PFA_protein <- subset(SDranks, Treatment == "PFA-protein exp.")
POL_protein <- subset(SDranks, Treatment == "POL-protein exp.")
POH_protein <- subset(SDranks, Treatment == "POH-protein exp.")

rank_mat <- cbind(PCG_transcript$Rank, EFA_transcript$Rank, PFA_transcript$Rank, 
                  POL_transcript$Rank, POH_transcript$Rank,
                  PCG_protein$Rank, EFA_protein$Rank, PFA_protein$Rank, 
                  POL_protein$Rank, POH_protein$Rank)

mean_rank <- meanranks(rank_mat)

#Also do ICV because there are no NAs
ICVranks <- data.frame()
for (i in 1:length(pathways)) {
  pathway <- subset(allSDs, Pathway == pathways[i])
  ranks <- rank(-pathway$ICV, na.last = TRUE) #puts NAs last, wanted to keep as NA but meanrank isn't performing
  add <- data.frame(pathway$Pathway, pathway$Treatment, pathway$ICV, ranks)
  ICVranks <- rbind(ICVranks, add)
}
colnames(ICVranks) <- c("Pathway", "Treatment", "ICV", "Rank")

#matrix should be pathway x treatment (calculates means by column)
PCG_transcriptICV <- subset(ICVranks, Treatment == "PCG-transcript exp.")
EFA_transcriptICV <- subset(ICVranks, Treatment == "EFA-transcript exp.")
PFA_transcriptICV <- subset(ICVranks, Treatment == "PFA-transcript exp.")
POL_transcriptICV <- subset(ICVranks, Treatment == "POL-transcript exp.")
POH_transcriptICV <- subset(ICVranks, Treatment == "POH-transcript exp.")
PCG_proteinICV <- subset(ICVranks, Treatment == "PCG-protein exp.")
EFA_proteinICV <- subset(ICVranks, Treatment == "EFA-protein exp.")
PFA_proteinICV <- subset(ICVranks, Treatment == "PFA-protein exp.")
POL_proteinICV <- subset(ICVranks, Treatment == "POL-protein exp.")
POH_proteinICV <- subset(ICVranks, Treatment == "POH-protein exp.")

rank_matICV <- cbind(PCG_transcriptICV$Rank, EFA_transcriptICV$Rank, PFA_transcriptICV$Rank, 
                  POL_transcriptICV$Rank, POH_transcriptICV$Rank,
                  PCG_proteinICV$Rank, EFA_proteinICV$Rank, PFA_proteinICV$Rank, 
                  POL_proteinICV$Rank, POH_proteinICV$Rank)

mean_rankICV <- meanranks(rank_matICV)

allMeanRanks <- data.frame(Treatment = c("PCG-transcript exp.", "EFA-transcript exp.", "PFA-transcript exp.",
                                         "POL-transcript exp.", "POH-transcript exp.",
                                         "PCG-protein exp.", "EFA-protein exp.", "PFA-protein exp.",
                                         "POL-protein exp.", "POH-protein exp."),
                           SD_meanrank = mean_rank$mean.ranks, ICV_meanrank = mean_rankICV$mean.ranks)

write.table(allMeanRanks, "~/allMeanRanks.csv", sep=',', quote = F, row.names = F)



#######################################################
#      ***Filter taxa using OTU_ID in Phyloseq***     #
#                                                     #
#                  Gian MN Benucci, PhD               #
#               Michigan State University             #
#                   benucci@msu.edu                   #      
#                                                     #
# Usage:                                              #
# badTaxa = c("OTU_1", "OTU_123", "OTU_456", ...)     #
#	biom_prun = filter_taxa(biom, badTaxa)              #
#######################################################


filter_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}




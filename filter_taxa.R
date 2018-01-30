###########################################################
# Filter taxa using OTUID across samples
# in a phyloseq object
#
# Usage:
# badTaxa = c("OTU_1", "OTU_123", "OTU_456", ...) 
#	biom_prun = filter_taxa(biom, badTaxa)
############################################################


filter_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}




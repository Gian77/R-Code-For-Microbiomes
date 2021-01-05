# Simple funciton to remove unwanted OTUs form a phyloseq object

# Example: badTaxa = c("OTU_1", "OTU_123", "OTU_456", ...) 
#	   filterTaxa(phyloseq_object, badTaxa)

filterTaxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}


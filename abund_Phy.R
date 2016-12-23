# This is a test function!
# This function calculate the relative abundance at Phylum level of your community table
# to load: source("abund_Phy.R")
# usage: abund_phy(biom, nsample)
# biom = biom_table imported in Phyloseq()
# nsample = number of samples in your map_file

abund_Phy <- function(biom, nsample){

require(phyloseq)
require(plyr)

biom_phyla = tax_glom(biom, "Phylum")
biom_phyla_prop = transform_sample_counts(biom_phyla, function(x) x/sum(x))

if (sum(taxa_sums(biom_phyla_prop)/nsample)==1) {
result <- taxa_sums(biom_phyla_prop)/nsample
}

# you can print a list instead of rearranging the table
#abund <- result
#Phy <- tax_table(biom_phyla_prop)
#newList <- list("Phylum abundances" = abund, "Phylum names" = Phy)
#return(newList)}

tax_table <- data.frame(tax_table(biom_phyla))
tax_table <- tax_table[c(2)]
tax_table$abundance <- as.vector(taxa_sums(biom_phyla_prop)/nsample)
tax_table <- arrange(tax_table, abundance, decreasing=TRUE)

return(tax_table)}



#############################################################
# To filter the tax_table assignment perfomed using QIIME 
#############################################################

tax_table(biom)[, "Kingdom"] <- gsub("k__", "", tax_table(biom)[, "Kingdom"])
tax_table(biom)[, "Phylum"] <- gsub("p__", "", tax_table(biom)[, "Phylum"])
tax_table(biom)[, "Class"] <- gsub("c__", "", tax_table(biom)[, "Class"])
tax_table(biom)[, "Order"] <- gsub("o__", "", tax_table(biom)[, "Order"])
tax_table(biom)[, "Family"] <- gsub("f__", "", tax_table(biom)[, "Family"])
tax_table(biom)[, "Genus"] <- gsub("g__", "", tax_table(biom)[, "Genus"])
tax_table(biom)[, "Species"] <- gsub("s__", "", tax_table(biom)[, "Species"])

# for 16s reads also consider
biom_16s <- subset_taxa(biom_16s, Phylum!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Class!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Order!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Family!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Genus!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Phylum!="chloroplast")
biom_16s <- subset_taxa(biom_16s, Class!="chloroplast")
biom_16s <- subset_taxa(biom_16s, Order!="chloroplast")
biom_16s <- subset_taxa(biom_16s, Family!="chloroplast")
biom_16s <- subset_taxa(biom_16s, Genus!="chloroplast")

biom_16s <- subset_taxa(biom_16s, Phylum!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Class!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Order!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Family!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Genus!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Phylum!="mitochondria")
biom_16s <- subset_taxa(biom_16s, Class!="mitochondria")
biom_16s <- subset_taxa(biom_16s, Order!="mitochondria")
biom_16s <- subset_taxa(biom_16s, Family!="mitochondria")
biom_16s <- subset_taxa(biom_16s, Genus!="mitochondria")






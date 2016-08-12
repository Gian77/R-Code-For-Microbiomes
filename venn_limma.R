### Plotting Venn Diagram form Phyloseq object using limma package

library(phyloseq)
library(limma)

## import biom_table in R as phyloseq objects
biom = import_biom("biom_table.biom")
map = import_qiime_sample_data("mapping.txt")
sample_data(biom) <- map
colnames(tax_table(biom)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")

## merge samples according a metadata variable
biom_soil = merge_samples(biom, "Soil")
sample_data(biom_soil)

## create the object for intersections
otu_table_soil <- t(otu_table(subset_samples(biom_ITS_soil, Soil%in%c(1,2,4,5,6))))
venn_counts <- vennCounts(otu_table_ITS_soil)
venn_counts

## plotting diagrams
vennDiagram(venn_counts_ITS)
vennDiagram(venn_counts_ITS, cex=c(1,1.2,0.8),
            names = c("BC", "CA", "NC","NWA", "OSU"),
            circle.col = c("red","blue","green","grey","yellow"))

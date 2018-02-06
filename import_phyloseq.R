##########################################
##  ***Importing data in phyloseq ***   ##
##        Gian MN Benucci, PhD          ##
##     Michigan State University        ##
##         benucci@msu.edu              ##
##       February 2, 2018               ##
##########################################


# removes all variables in the global environment for a fresh start!
rm(list = ls(all=TRUE)) 

library("phyloseq")
library("Biostrings")

# apprach 1)
biom = import_biom("otu_table.biom")
map = import_qiime_sample_data("mapping.txt")
sample_data(biom) <- map
colnames(tax_table(biom)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_set <- readDNAStringSet("otus.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
physeq_obj <- merge_phyloseq(biom, otus_rep_set)

# apprach 2)
# import OTU table
otu_tab <- read.delim("otutab_UPARSE_ITS.txt", row.names=1, header=TRUE, sep="\t") 
head(otu_tab)
dim(otu_tab)






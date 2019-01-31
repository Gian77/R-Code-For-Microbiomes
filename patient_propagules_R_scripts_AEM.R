# Patient propagules: do long-term soil archives preserve legacy fungal and bacterial communities?
# Benucci GMN, Rennick B, Bonito G.
# Michigan State University
# Applied anbd Environmental Microbiology 

# Benucci et al. - manuscripot R scripts ---------------

options(scipen = 999) 
options(max.print=100000000)

library(phyloseq)
library(Biostrings)
library(ggplot2)

# color palettes -----------------------------------------
palette_CB6 <-c("#332288","#88CCEE","#117733","#DDCC77","#FF899d","#AA4499")


# importing datasets -------------------------------------
biom_ITS_uparse = import_biom("otutab_UPARSE_ITS_tax_json.biom")
map_ITS = import_qiime_sample_data("metadata_its.txt")
sample_data(biom_ITS_uparse) <- map_ITS
colnames(tax_table(biom_ITS_uparse)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_uparse <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_uparse <- merge_phyloseq(biom_ITS_uparse, otus_rep_ITS_uparse)
head(tax_table(biom_ITS_uparse))
refseq(biom_ITS_uparse)
biom_ITS_uparse


biom_16s_uparse = import_biom("otutab_UPARSE_16s_tax_json.biom")
map_16s = import_qiime_sample_data("metadata_16s.txt")
sample_data(biom_16s_uparse) <- map_16s
colnames(tax_table(biom_16s_uparse)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_16s_uparse <- readDNAStringSet("otus_16s.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_16s_uparse <- merge_phyloseq(biom_16s_uparse, otus_rep_16s_uparse)
sample_data(biom_16s_uparse)
head(tax_table(biom_16s_uparse))
biom_16s_uparse


# cleaning and filtering ----------------------------------
tax_table(biom_ITS_uparse)[, "Kingdom"] <- gsub("PMI", "", tax_table(biom_ITS_uparse)[, "Kingdom"]) # based on the reference taxonomy
tax_table(biom_ITS_uparse)[, "Kingdom"] <- gsub("NVP", "", tax_table(biom_ITS_uparse)[, "Kingdom"])

tax_table(biom_ITS_uparse)[, "Kingdom"] <- gsub("k__", "", tax_table(biom_ITS_uparse)[, "Kingdom"])
tax_table(biom_ITS_uparse)[, "Phylum"] <- gsub("p__", "", tax_table(biom_ITS_uparse)[, "Phylum"])
tax_table(biom_ITS_uparse)[, "Class"] <- gsub("c__", "", tax_table(biom_ITS_uparse)[, "Class"])
tax_table(biom_ITS_uparse)[, "Order"] <- gsub("o__", "", tax_table(biom_ITS_uparse)[, "Order"])
tax_table(biom_ITS_uparse)[, "Family"] <- gsub("f__", "", tax_table(biom_ITS_uparse)[, "Family"])
tax_table(biom_ITS_uparse)[, "Genus"] <- gsub("g__", "", tax_table(biom_ITS_uparse)[, "Genus"])
tax_table(biom_ITS_uparse)[, "Species"] <- gsub("s__", "", tax_table(biom_ITS_uparse)[, "Species"])

biom_ITS_uparse <- subset_taxa(biom_ITS_uparse, Kingdom == "Fungi")
tax_table(biom_ITS_uparse)[tax_table(biom_ITS_uparse)=="unidentified"]<- NA
tax_table(biom_ITS_uparse)[is.na(tax_table(biom_ITS_uparse))]<-"Unclassified"

head(otu_table(biom_ITS_uparse))
head(tax_table(biom_ITS_uparse))
head(sample_data(biom_ITS_uparse))
biom_ITS_uparse


tax_table(biom_16s_uparse)[, "Kingdom"] <- gsub("k__", "", tax_table(biom_16s_uparse)[, "Kingdom"])
tax_table(biom_16s_uparse)[, "Phylum"] <- gsub("p__", "", tax_table(biom_16s_uparse)[, "Phylum"])
tax_table(biom_16s_uparse)[, "Class"] <- gsub("c__", "", tax_table(biom_16s_uparse)[, "Class"])
tax_table(biom_16s_uparse)[, "Order"] <- gsub("o__", "", tax_table(biom_16s_uparse)[, "Order"])
tax_table(biom_16s_uparse)[, "Family"] <- gsub("f__", "", tax_table(biom_16s_uparse)[, "Family"])
tax_table(biom_16s_uparse)[, "Genus"] <- gsub("g__", "", tax_table(biom_16s_uparse)[, "Genus"])
tax_table(biom_16s_uparse)[, "Species"] <- gsub("s__", "", tax_table(biom_16s_uparse)[, "Species"])

biom_16s_uparse <- subset_taxa(biom_16s_uparse, Phylum!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Class!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Order!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Family!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Genus!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Phylum!="chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Class!="chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Order!="chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Family!="chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Genus!="chloroplast")

biom_16s_uparse <- subset_taxa(biom_16s_uparse, Phylum!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Class!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Order!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Family!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Genus!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Phylum!="mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Class!="mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Order!="mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Family!="mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Genus!="mitochondria")


tax_table(biom_16s_uparse)[tax_table(biom_16s_uparse)==""]<- NA
tax_table(biom_16s_uparse)[is.na(tax_table(biom_16s_uparse))]<-"Unclassified"
head(otu_table(biom_16s_uparse))
head(tax_table(biom_16s_uparse))
head(sample_data(biom_16s_uparse))
biom_16s_uparse



# filtering contaminants ----------------------------------
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

badTaxa_uparse_ITS =c("OTU_8","OTU_39","OTU_269")

biom_ITS_uparse_qf_prun = pop_taxa(biom_ITS_uparse, badTaxa_uparse_ITS)
biom_ITS_uparse_qf_prun

badTaxa_uparse_16s =c("OTU_4","OTU_9","OTU_8176","OTU_5073","OTU_37","OTU_47","OTU_109","OTU_135","OTU_133",
                      "OTU_405","OTU_603","OTU_397","OTU_865","OTU_1303","OTU_970","OTU_2480",
                      "OTU_2015","OTU_2488","OTU_2642","OTU_2815","OTU_2390","OTU_12514",
                      "OTU_9394","OTU_7279","OTU_8176")

biom_16s_uparse_qf_prun = pop_taxa(biom_16s_uparse, badTaxa_uparse_16s)
biom_16s_uparse_qf_prun


# Filtering out OTUs -----------------------------------
# Oliver et al. 2015, PCR errors and tag switching
# Lindhal et al. 2013, tag switching - that's a good  one!

biom_ITS_uparse_qf_prun -> biom_ITS_uparse_filt
otu_table(biom_ITS_uparse_filt) <- otu_table(biom_ITS_uparse_filt)[which(rowSums(otu_table(biom_ITS_uparse_filt)) >= 10),]
biom_ITS_uparse_filt

tax_table(biom_ITS_uparse_filt)
sample_data(biom_ITS_uparse_filt)


biom_16s_uparse_qf_prun -> biom_16s_uparse_filt
otu_table(biom_16s_uparse_filt) <- otu_table(biom_16s_uparse_filt)[which(rowSums(otu_table(biom_16s_uparse_filt)) >= 10),] 
biom_16s_uparse_filt

tax_table(biom_16s_uparse_filt)
sample_data(biom_16s_uparse_filt)



# filtering out PCR sequencing control samples ------------------------------

sample_data(biom_ITS_uparse_filt)
otu_table(biom_ITS_uparse_filt) <- subset(otu_table(biom_ITS_uparse_filt), select = -c(DF33bis, control2))
otu_table(biom_ITS_uparse_filt) <- otu_table(biom_ITS_uparse_filt)[which(rowSums(otu_table(biom_ITS_uparse_filt)) >= 1),] 
biom_ITS_uparse_filt

sum(rowSums(otu_table(biom_ITS_uparse_filt)))
any(taxa_sums(biom_ITS_uparse_filt) == 0)

sample_data(biom_16s_uparse_filt)
otu_table(biom_16s_uparse_filt) <- subset(otu_table(biom_16s_uparse_filt),
                          select = -c(DF33bis, FJ34B, P615A, C4FC, p589, empty1, empty3, control3, control2, control1))
otu_table(biom_16s_uparse_filt) <- otu_table(biom_16s_uparse_filt)[which(rowSums(otu_table(biom_16s_uparse_filt)) >= 1),] 
biom_16s_uparse_filt

sum(rowSums(otu_table(biom_16s_uparse_filt)))
any(taxa_sums(biom_16s_uparse_filt) == 0)

# model species richness vs. read number  -----------------------------------
library("dplyr")
library("ggpubr")
source("../lm_eqn.R")

otu_ITS_uparse <- as.data.frame(otu_table(biom_ITS_uparse_filt))

sums_ITS_uparse <- data.frame(colSums(otu_table(biom_ITS_uparse_filt)))
colnames(sums_ITS_uparse) <- "readNO"
sums_ITS_uparse$richness <- specnumber(otu_ITS_uparse, MARGIN = 2)
sums_ITS_uparse$rarefied <- rarefy(otu_ITS_uparse, min(colSums(otu_ITS_uparse)), se = FALSE, MARGIN = 2)
sums_ITS_uparse$sample <- row.names(sums_ITS_uparse)
sums_ITS_uparse <- arrange(sums_ITS_uparse, sample)
sums_ITS_uparse$group <- c(rep("Forest", times=32),rep("Poplar", times=27))
sums_ITS_uparse
str(sums_ITS_uparse)

plot(sums_ITS_uparse$readNO, sums_ITS_uparse$richness)
cor(sums_ITS_uparse$readNO, sums_ITS_uparse$richness, method = "pearson")
cor.test(sums_ITS_uparse$readNO, sums_ITS_uparse$richness)$p.value
cor(sums_ITS_uparse$readNO, sums_ITS_uparse$rarefied)
cor.test(sums_ITS_uparse$readNO, sums_ITS_uparse$rarefied)$p.value
fit_ob_ITS <- lm(richness ~ readNO, data = sums_ITS_uparse)
fit_ra_ITS <- lm(rarefied ~ readNO, data = sums_ITS_uparse)
summary(fit_ob_ITS)
summary(fit_ra_ITS)
plot(fit_ob_ITS)
plot(fit_ra_ITS)

# Additional exploratory figure ------------
par(mfrow=c(2,2)) 
plot(fit_ob_ITS) 
par(mfrow=c(1,1)) 


mod_ITS = ggplot(sums_ITS_uparse, aes(x=sums_ITS_uparse$readNO, y=as.numeric(sums_ITS_uparse$richness))) +
  geom_point(aes(size=readNO, shape=group, color=group)) +
  labs(x="Read number", y="Observed richness", title="") +
  geom_smooth(method="lm") + 
  annotate("text", x=30000, y=1000, label=lm_eqn(lm(richness ~ readNO, data = sums_ITS_uparse)), size = 3, parse=TRUE) +
  theme(title = element_text(size = 8)) +
  theme_bw() +
  theme(legend.position="none")

mod_ITS


otu_16s_uparse <- as.data.frame(otu_table(biom_16s_uparse_filt))

sums_16s_uparse <- data.frame(colSums(otu_table(biom_16s_uparse_filt)))
colnames(sums_16s_uparse) <- "readNO"
sums_16s_uparse$richness <- specnumber(otu_16s_uparse, MARGIN = 2)
sums_16s_uparse$rarefied <- rarefy(otu_16s_uparse, min(colSums(otu_16s_uparse)), se = FALSE, MARGIN = 2)
sums_16s_uparse$sample <- row.names(sums_16s_uparse)
sums_16s_uparse <- arrange(sums_16s_uparse, sample)
sums_16s_uparse$group <- c(rep("Forest", times=32),rep("Poplar", times=27))
sums_16s_uparse
str(sums_16s_uparse)

plot(sums_16s_uparse$readNO, sums_16s_uparse$richness)
cor(sums_16s_uparse$readNO, sums_16s_uparse$richness)
cor.test(sums_16s_uparse$readNO, sums_16s_uparse$richness)$p.value
cor(sums_16s_uparse$readNO, sums_16s_uparse$rarefied)
cor.test(sums_16s_uparse$readNO, sums_16s_uparse$rarefied)$p.value
fit_ob_16s <- lm(richness ~ readNO, data = sums_16s_uparse)
fit_ra_16s <- lm(rarefied ~ readNO, data = sums_16s_uparse)
summary(fit_ob_16s)
summary(fit_ra_16s)
plot(fit_ob_16s)
plot(fit_ra_16s)

# Additional exploratory figure ------------
par(mfrow=c(2,2)) 
plot(fit_ob_16s) 
par(mfrow=c(1,1)) 

mod_16s = ggplot(sums_16s_uparse, aes(x=sums_16s_uparse$readNO, y=as.numeric(sums_16s_uparse$richness))) +
  geom_point(aes(size=readNO, shape=group, color=group)) +
  labs(x="Read number", y="Observed richness", title="") +
  geom_smooth(method="lm") + 
  annotate("text", x=45000, y=5000, label=lm_eqn(lm(richness ~ readNO, data = sums_16s_uparse)),size = 3, parse=TRUE) +
  theme(title = element_text(size = 8)) +
  expand_limits(y = 0) +
  theme_bw() +
  theme(legend.position="none")

mod_16s

# FIGURE S2 -----------------------------------------------
ggarrange(mod_ITS, mod_16s,
          labels = c("A", "B"),
          widths = c(1, 1),
          align = "none", ncol = 2, nrow = 1)



# rarefy the whole dataset ----------------------------------

min(sample_sums(biom_ITS_uparse_filt))
biom_ITS_uparse_filt_ev = rarefy_even_depth(biom_ITS_uparse_filt,
              sample.size = min(sample_sums(biom_ITS_uparse_filt)),
              rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_ITS_uparse_filt_ev
colSums(otu_table(biom_ITS_uparse_filt_ev))
any(taxa_sums(biom_ITS_uparse_filt_ev) == 0)

min(sample_sums(biom_16s_uparse_filt))
biom_16s_uparse_filt_ev = rarefy_even_depth(biom_16s_uparse_filt,
                sample.size = min(sample_sums(biom_16s_uparse_filt)),
                rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_16s_uparse_filt_ev
colSums(otu_table(biom_16s_uparse_filt_ev))
any(taxa_sums(biom_16s_uparse_filt_ev) == 0)


# Extracting data of different experiments -----------------------------


# this is EXPERIMENT 1 ITS dataset ---------------
biom_ITS_exp1 <- subset_samples(biom_ITS_uparse_filt, exp1%in%c("exp1"))
otu_table(biom_ITS_exp1) <- otu_table(biom_ITS_exp1)[which(rowSums(otu_table(biom_ITS_exp1)) >= 1),] 
biom_ITS_exp1
sample_data(biom_ITS_exp1)
tax_table(biom_ITS_exp1)

biom_ITS_exp1_ev <- subset_samples(biom_ITS_uparse_filt_ev, exp1%in%c("exp1"))
otu_table(biom_ITS_exp1_ev) <- otu_table(biom_ITS_exp1_ev)[which(rowSums(otu_table(biom_ITS_exp1_ev)) >= 1),] 
biom_ITS_exp1_ev
sample_data(biom_ITS_exp1_ev)
tax_table(biom_ITS_exp1_ev)

# this is EXPERIMENT 2 ITS dataset --------------
biom_ITS_exp2 <- subset_samples(biom_ITS_uparse_filt, exp2%in%c("exp2"))
otu_table(biom_ITS_exp2) <- otu_table(biom_ITS_exp2)[which(rowSums(otu_table(biom_ITS_exp2)) >=1),] 
biom_ITS_exp2
sample_data(biom_ITS_exp2)

biom_ITS_exp2_ev <- subset_samples(biom_ITS_uparse_filt_ev, exp2%in%c("exp2"))
otu_table(biom_ITS_exp2_ev) <- otu_table(biom_ITS_exp2_ev)[which(rowSums(otu_table(biom_ITS_exp2_ev)) >=1),] 
biom_ITS_exp2_ev
sample_data(biom_ITS_exp2_ev)

# this is EXPERIMENT 3 ITS dataset --------------
biom_ITS_exp3 <- subset_samples(biom_ITS_uparse_filt, exp3%in%c("exp3"))
otu_table(biom_ITS_exp3) <- otu_table(biom_ITS_exp3)[which(rowSums(otu_table(biom_ITS_exp3)) >=1),] 
biom_ITS_exp3
sample_data(biom_ITS_exp3)

biom_ITS_exp3_ev <- subset_samples(biom_ITS_uparse_filt_ev, exp3%in%c("exp3"))
otu_table(biom_ITS_exp3_ev) <- otu_table(biom_ITS_exp3_ev)[which(rowSums(otu_table(biom_ITS_exp3_ev)) >=1),] 
biom_ITS_exp3_ev
sample_data(biom_ITS_exp3_ev)




# this is EXPERIMENT 1 16s dataset ---------------
biom_16s_exp1 <- subset_samples(biom_16s_uparse_filt, exp1%in%c("exp1"))
otu_table(biom_16s_exp1) <- otu_table(biom_16s_exp1)[which(rowSums(otu_table(biom_16s_exp1)) >= 1),] 
biom_16s_exp1
sample_data(biom_16s_exp1)
tax_table(biom_16s_exp1)

biom_16s_exp1_ev <- subset_samples(biom_16s_uparse_filt_ev, exp1%in%c("exp1"))
otu_table(biom_16s_exp1_ev) <- otu_table(biom_16s_exp1_ev)[which(rowSums(otu_table(biom_16s_exp1_ev)) >= 1),] 
biom_16s_exp1_ev
sample_data(biom_16s_exp1_ev)
tax_table(biom_16s_exp1_ev)

# this is EXPERIMENT 2 16s dataset ---------------
biom_16s_exp2 <- subset_samples(biom_16s_uparse_filt, exp2%in%c("exp2"))
otu_table(biom_16s_exp2) <- otu_table(biom_16s_exp2)[which(rowSums(otu_table(biom_16s_exp2)) >=1),] 
biom_16s_exp2
sample_data(biom_16s_exp2)

biom_16s_exp2_ev <- subset_samples(biom_16s_uparse_filt_ev, exp2%in%c("exp2"))
otu_table(biom_16s_exp2_ev) <- otu_table(biom_16s_exp2_ev)[which(rowSums(otu_table(biom_16s_exp2_ev)) >=1),] 
biom_16s_exp2_ev
sample_data(biom_16s_exp2_ev)


# this is EXPERIMENT 3 16s dataset ---------------
biom_16s_exp3 <- subset_samples(biom_16s_uparse_filt, exp3%in%c("exp3"))
otu_table(biom_16s_exp3) <- otu_table(biom_16s_exp3)[which(rowSums(otu_table(biom_16s_exp3)) >=1),] 
biom_16s_exp3
sample_data(biom_16s_exp3)

biom_16s_exp3_ev <- subset_samples(biom_16s_uparse_filt_ev, exp3%in%c("exp3"))
otu_table(biom_16s_exp3_ev) <- otu_table(biom_16s_exp3_ev)[which(rowSums(otu_table(biom_16s_exp3_ev)) >=1),] 
biom_16s_exp3_ev
sample_data(biom_16s_exp3_ev)




# alpha-diversity -----------------------------------

sample_data(biom_ITS_exp1)$year <- as.factor(sample_data(biom_ITS_exp1)$year)
sample_data(biom_ITS_exp1)$year <- factor(sample_data(biom_ITS_exp1)$year, levels=c("2015","2014","2010","2005","2000","1995"))

alpha_ITS_exp1 = plot_richness(biom_ITS_exp1, x="year", color="year", shape="habitat", measures = c("Observed")) + # nrow=3, sortby = "Chao1"
  geom_point(size=0.1, alpha = 0.9) + 
  geom_point(position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(17, 16))  +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  facet_grid(~habitat, scales = "free_x", space="free_x") + 
  theme(strip.text.x = element_text(size = 8)) +
  ylim(0, 700) +
  theme_pubclean() +
  scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000")) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, size=7)) +
  theme(axis.text.y = element_text(hjust=0.5, size=7)) +
  labs(title="", x="Year", y="Observed richness") +
  theme(axis.title = element_text(size = 8)) +  
  theme(legend.position="none")

alpha_ITS_exp1

sample_data(biom_16s_exp1)$year <- as.factor(sample_data(biom_16s_exp1)$year)
sample_data(biom_16s_exp1)$year <- factor(sample_data(biom_16s_exp1)$year, levels=c("2015","2014","2010","2005","2000","1995"))

alpha_16s_exp1 = plot_richness(biom_16s_exp1, x="year", color="year", shape="habitat", measures = c("Observed")) + # nrow=3, sortby = "Chao1"
  geom_point(size=0.1, alpha = 0.9) + 
  geom_point(position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(17, 16))  +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  facet_grid(~habitat, scales = "free_x", space="free_x") + 
  theme(strip.text.x = element_text(size = 8)) +
  ylim(0, 4500) +
  theme_pubclean() +
  scale_colour_manual(values=c("#0000FF","#0000FF","#0000FF","#0000FF","#0000FF","#0000FF")) + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, size=7)) +
  theme(axis.text.y = element_text(hjust=0.5, size=7)) +
  labs(title="", x="Year", y="Observed richness") +
  theme(axis.title = element_text(size = 8)) +  
  theme(legend.position="none")

alpha_16s_exp1


sample_data(biom_ITS_exp2)$year <- as.factor(sample_data(biom_ITS_exp2)$year)
sample_data(biom_ITS_exp2)$year <- factor(sample_data(biom_ITS_exp2)$year, levels=c("2015","2014","2010","2005","2000"))

alpha_ITS_exp2 = plot_richness(biom_ITS_exp2, x="year", color="habitat", shape="habitat", measures = c("Observed")) + # nrow=3, sortby = "Chao1"
  geom_point(size=1, alpha = 0.9) + 
  scale_shape_manual(values = c(17, 16))  +
  ylim(0, 700) +
  facet_grid(~habitat, scales = "free_x") + 
  theme(strip.text.x = element_text(size = 8)) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  scale_colour_manual(values=c("#000000")) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, size=7)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + # element_text(hjust=0.5, size=7)
  theme(axis.title = element_text(size = 8)) + 
  labs(title="", x="Year", y = "") +
  expand_limits(x = 0, y = 0) +
  theme(legend.position="none")

alpha_ITS_exp2


sample_data(biom_16s_exp2)$year <- as.factor(sample_data(biom_16s_exp2)$year)
sample_data(biom_16s_exp2)$year <- factor(sample_data(biom_16s_exp2)$year, levels=c("2015","2014","2010","2005","2000"))

alpha_16s_exp2 = plot_richness(biom_16s_exp2, x="year", color="habitat", shape="habitat", measures = c("Observed")) + # nrow=3, sortby = "Chao1"
  geom_point(size=1, alpha = 0.9) + 
  scale_shape_manual(values = c(17, 16))  +
  ylim(0, 4500) +
  facet_grid(~habitat, scales = "free_x") + 
  theme(strip.text.x = element_text(size = 8)) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  scale_colour_manual(values=c("#0000FF")) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, size=7)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(axis.title = element_text(size = 8)) + 
  labs(title="", x="Year", y = "") +
  expand_limits(x = 0, y = 0) +
  theme(legend.position="none")

alpha_16s_exp2


alpha_ITS_exp3 = plot_richness(biom_ITS_exp3, x="site", color="year", shape="habitat", measures = c("Observed")) + # nrow=3, sortby = "Chao1"
  geom_point(size=1, alpha = 0.9) + 
  scale_shape_manual(values = c(17, 16))  +
  ylim(0, 700) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  facet_grid(~habitat, scales = "free_x") + 
  theme(strip.text.x = element_text(size = 8)) +
  scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000")) + 
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, size=7)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(axis.title = element_text(size = 8)) +  
  labs(title="", x="Site", y = "") +
  expand_limits(x = 0, y = 0) +
  theme(legend.position="none")

alpha_ITS_exp3


alpha_16s_exp3 = plot_richness(biom_16s_exp3, x="site", color="year", shape="habitat", measures = c("Observed")) + # nrow=3, sortby = "Chao1"
  geom_point(size=1, alpha = 0.9) + 
  scale_shape_manual(values = c(17, 16))  +
  ylim(0, 4500) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  facet_grid(~habitat, scales = "free_x") + 
  theme(strip.text.x = element_text(size = 8)) +
  scale_colour_manual(values=c("#0000FF","#0000FF","#0000FF","#0000FF","#0000FF","#0000FF")) + 
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, size=7)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(axis.title = element_text(size = 8)) +  
  labs(title="", x="Site", y = "") +
  expand_limits(x = 0, y = 0) +
  theme(legend.position="none")

alpha_16s_exp3

# > FIGURE 1 -----------------------------------------------
library("ggpubr")

ggarrange(alpha_ITS_exp1, alpha_ITS_exp2, alpha_ITS_exp3,
          alpha_16s_exp1, alpha_16s_exp2, alpha_16s_exp3,
          labels = c("A","B","C","D","E","F"),
          widths = c(1,0.45,0.7,1,0.45,0.7),
          align = "h", ncol = 6, nrow = 1)








# linear models: richness - year -------------------------
library("vegan")

otus_ITS_exp1 <- as.data.frame(otu_table(biom_ITS_exp1))
metadata_ITS_exp1 <- as.data.frame(as.matrix(sample_data(biom_ITS_exp1)))

div_ITS_exp1 <- as.data.frame(specnumber(otus_ITS_exp1, MARGIN = 2))
colnames(div_ITS_exp1) <- "richness"
div_ITS_exp1$rarefied <- rarefy(otus_ITS_exp1, min(colSums(otus_ITS_exp1)), se = FALSE, MARGIN = 2)
div_ITS_exp1$habitat <- metadata_ITS_exp1$habitat
div_ITS_exp1$site <- metadata_ITS_exp1$site
div_ITS_exp1$year <- as.numeric(as.character(metadata_ITS_exp1$year))
div_ITS_exp1$time <- as.vector(2015-div_ITS_exp1$year)
div_ITS_exp1

otus_16s_exp1 <- as.data.frame(otu_table(biom_16s_exp1))
metadata_16s_exp1 <- as.data.frame(as.matrix(sample_data(biom_16s_exp1)))

div_16s_exp1 <- as.data.frame(specnumber(otus_16s_exp1, MARGIN = 2))
colnames(div_16s_exp1) <- "richness"
div_16s_exp1$rarefied <- rarefy(otus_16s_exp1, min(colSums(otus_16s_exp1)), se = FALSE, MARGIN = 2)
div_16s_exp1$habitat <- metadata_16s_exp1$habitat
div_16s_exp1$site <- metadata_16s_exp1$site
div_16s_exp1$year <- as.numeric(as.character(metadata_16s_exp1$year))
div_16s_exp1$time <- as.vector(2015-div_16s_exp1$year)
div_16s_exp1

otus_ITS_exp1_PS <- subset(div_ITS_exp1, habitat%in%c("PS"))
otus_ITS_exp1_DF <- subset(div_ITS_exp1, habitat%in%c("DF"))
otus_16s_exp1_PS <- subset(div_16s_exp1, habitat%in%c("PS"))
otus_16s_exp1_DF <- subset(div_16s_exp1, habitat%in%c("DF"))


cor(otus_ITS_exp1_PS$richness, div_pop$time)
cor.test(otus_ITS_exp1_PS$richness, otus_ITS_exp1_PS$time)$p.value
fit_ITS_exp1_PS <- lm(richness ~ time, data = otus_ITS_exp1_PS)
summary(fit_ITS_exp1_PS)

cor(otus_ITS_exp1_DF$richness, otus_ITS_exp1_DF$time)
cor.test(otus_ITS_exp1_DF$richness, otus_ITS_exp1_DF$time)$p.value
fit_ITS_exp1_DF <- lm(richness ~ time, data = otus_ITS_exp1_DF)
summary(fit_ITS_exp1_DF)

cor(otus_16s_exp1_PS$richness, otus_16s_exp1_PS$time)
cor.test(otus_16s_exp1_PS$richness, otus_16s_exp1_PS$time)$p.value
fit_16s_exp1_PS <- lm(richness ~ time, data = otus_16s_exp1_PS)
summary(fit_16s_exp1_PS)

cor(otus_16s_exp1_DF$richness, otus_16s_exp1_DF$time)
cor.test(otus_16s_exp1_DF$richness, otus_16s_exp1_DF$time)$p.value
fit_16s_exp1_DF <- lm(richness ~ time, data = otus_16s_exp1_DF)
summary(fit_16s_exp1_DF)



lm_ITS_exp1_PS = ggplot(otus_ITS_exp1_PS, aes(x=otus_ITS_exp1_PS$time, y=otus_ITS_exp1_PS$richness)) +
  geom_point(alpha = 0.8, shape=16, size=3, color="black") +
  labs(x="Time", y="Observed richness", title="") +
  geom_smooth(method="lm", col="red") + 
  annotate("text", x=10, y=700, label=lm_eqn(lm(richness ~ time, otus_ITS_exp1_PS)), size = 2.5, parse=TRUE) +
  theme(title = element_text(size = 8)) +
  expand_limits(y = 0) +
  scale_x_continuous(breaks=c(0,1,5,10,15,20)) +
  theme_bw() +
  theme(legend.position="none")

lm_ITS_exp1_PS

lm_ITS_exp1_DF = ggplot(otus_ITS_exp1_DF, aes(x=otus_ITS_exp1_DF$time, y=otus_ITS_exp1_DF$richness)) +
  geom_point(alpha = 0.8, shape=17, size=3, color="black") +
  labs(x="Time", y="Observed richness", title="") +
  geom_smooth(method="lm", col="red") + 
  annotate("text", x=7.5, y=700, label=lm_eqn(lm(richness ~ time, otus_ITS_exp1_DF)), size = 2.5, parse=TRUE) +
  theme(title = element_text(size = 8)) +
  expand_limits(y = 0) +
  scale_x_continuous(breaks=c(0,1,5,10,15)) +
  theme_bw() +
  theme(legend.position="none")

lm_ITS_exp1_DF

lm_16s_exp1_PS = ggplot(otus_16s_exp1_PS, aes(x=otus_16s_exp1_PS$time, y=otus_16s_exp1_PS$richness)) +
  geom_point(alpha = 0.8, shape=16, size=3, color="blue") +
  labs(x="Time", y="Observed richness", title="") +
  geom_smooth(method="lm", col="red") + 
  annotate("text", x=10, y=3700, label=lm_eqn(lm(richness ~ time, otus_16s_exp1_PS)), size = 2.5, parse=TRUE) +
  theme(title = element_text(size = 8)) +
  expand_limits(y = 0) +
  scale_x_continuous(breaks=c(0,1,5,10,15,20)) +
  theme_bw() +
  theme(legend.position="none")

lm_16s_exp1_PS

lm_16s_exp1_DF = ggplot(otus_16s_exp1_DF, aes(x=otus_16s_exp1_DF$time, y=otus_16s_exp1_DF$richness)) +
  geom_point(alpha = 0.8, shape=17, size=3, color="blue") +
  labs(x="Time", y="Observed richness", title="") +
  geom_smooth(method="lm", col="red") + 
  annotate("text", x=7.5, y=3700, label=lm_eqn(lm(richness ~ time, otus_16s_exp1_DF)), size = 2.5, parse=TRUE) +
  theme(title = element_text(size = 8)) +
  expand_limits(y = 0) +
  scale_x_continuous(breaks=c(0,1,5,10,15)) +
  theme_bw() +
  theme(legend.position="none")

lm_16s_exp1_DF

# > FIGURE 2 -----------------------------------------------
ggarrange(lm_ITS_exp1_PS, lm_16s_exp1_PS,
          lm_ITS_exp1_DF, lm_16s_exp1_DF,
          labels = c("A", "B", "C", "D"),
          widths = c(1, 1, 1, 1),
          align = "h", ncol = 2, nrow = 2)



# INDICATOR SPECIES ANALYSIS -----------------------------------------------------------------------------------------------------------------------------

#indicator species analysis (isa) Fungi
library("indicspecies")

isa_ITS_exp1 <- multipatt(as.data.frame(t(otus_ITS_exp1)), metadata_ITS_exp1$year, control=how(nperm=9999))
summary(isa_ITS_exp1, indvalcomp=TRUE)
head(isa_ITS_exp1$sign)
isa_ITS_exp1 -> isa_ITS_exp1_fdr
isa_ITS_exp1_fdr$sign$p.value<-p.adjust(isa_ITS_exp1$sign$p.value, "fdr")
summary(isa_ITS_exp1_fdr)

isa_16s_exp1 <- multipatt(as.data.frame(t(otus_16s_exp1)), metadata_16s_exp1$year, control=how(nperm=9999))
summary(isa_16s_exp1)
head(isa_16s_exp1$sign)
isa_16s_exp1 -> isa_16s_exp1_fdr
isa_16s_exp1_fdr$sign$p.value<-p.adjust(isa_16s_exp1$sign$p.value, "fdr")
summary(isa_uparse_16s_df_fdr)


# beta diversity: non-metric multidimensional scaling --------------------------

nmds_ITS_exp1_ev <- ordinate(biom_ITS_exp1_ev, method ="NMDS", distance="bray")
nmds_ITS_exp1_ev
nmds_16s_exp1_ev <- ordinate(biom_16s_exp1_ev, method ="NMDS", distance="bray")
nmds_16s_exp1_ev

plot_nmds_ITS_exp1_ev = plot_ordination(biom_ITS_exp1_ev, nmds_ITS_exp1_ev, color="year", shape="habitat") + 
  geom_point(size=2.5, alpha=0.9, aes(shape=habitat)) +
  scale_shape_manual(values = c(17, 16), labels=c("DF","PS")) +
  theme_bw() +
  scale_colour_manual(values=palette_CB6) +
  theme(legend.position="right")

plot_nmds_ITS_exp1_ev

plot_nmds_16s_exp1_ev = plot_ordination(biom_16s_exp1_ev, nmds_16s_exp1_ev, color="year", shape="habitat") + 
  geom_point(size=2.5, alpha=0.9, aes(shape=habitat)) +
  scale_shape_manual(values = c(17, 16), labels=c("DF","PS")) +
  theme_bw() +
  scale_colour_manual(values=palette_CB6) +
  theme(legend.position="right")

plot_nmds_16s_exp1_ev

otus_ITS_exp1_ev <- as.data.frame(otu_table(biom_ITS_exp1_ev))
nmds_ITS_exp1_ev_vegan = metaMDS(t(otus_ITS_exp1_ev), k=2, distance="bray")
nmds_ITS_exp1_ev_vegan
stressplot(nmds_ITS_exp1_ev_vegan)

otus_16s_exp1_ev <- as.data.frame(otu_table(biom_16s_exp1_ev))
nmds_16s_exp1_ev_vegan = metaMDS(t(otus_16s_exp1_ev), k=2, distance="bray")
nmds_16s_exp1_ev_vegan
stressplot(nmds_16s_exp1_ev_vegan)


# FIGURE S9 ----------------------------------------------------------
par(mfrow=c(1,2)) 
stressplot(nmds_ITS_exp1_ev_vegan, main="A", adj="0")
stressplot(nmds_16s_exp1_ev_vegan, main="B", adj="0")
par(mfrow=c(1,1)) 



# CAP models -------------------------

otus_ITS_exp1_ev <- as.data.frame(otu_table(biom_ITS_exp1_ev))
metadata_ITS_exp1$year <- as.factor(metadata_ITS_exp1$year)

fullModel_ITS_exp1_ev <- capscale(t(otus_ITS_exp1_ev) ~ year * habitat * sample, metadata_ITS_exp1, distance = "bray")
fullModel_ITS_exp1_ev <- capscale(t(otus_ITS_exp1_ev) ~ year * site, metadata_ITS_exp1, distance = "bray")
fullModel_ITS_exp1_ev <- capscale(t(otus_ITS_exp1_ev) ~ year * habitat, metadata_ITS_exp1, distance = "bray")
fullModel_ITS_exp1_ev <- capscale(t(otus_ITS_exp1_ev) ~ year + habitat + year:habitat, metadata_ITS_exp1, distance = "bray")
anova(fullModel_ITS_exp1_ev, permutations = 9999)

# backward selection
parsimoniousModel_ITS_exp1_ev <- ordistep(fullModel_ITS_exp1_ev, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
parsimoniousModel_ITS_exp1_ev

# extract the parsimonious model formula for later use
FormulaParsimoniousModel_ITS_exp1_ev <- formula(parsimoniousModel_ITS_exp1_ev)[1:3]
FormulaParsimoniousModel_ITS_exp1_ev

# performing capscale analysis
capscale(FormulaParsimoniousModel_ITS_exp1_ev, metadata_ITS_exp1, distance = "bray") -> CAP
head(summary(CAP))
anova(CAP, permutations = 9999)

RsquareAdj(CAP)$adj.r.squared



otus_16s_exp1_ev <- as.data.frame(otu_table(biom_16s_exp1_ev))
metadata_16s_exp1$year <- as.factor(metadata_16s_exp1$year)

fullModel_16s_exp1_ev <- capscale(t(otus_16s_exp1_ev) ~ year * habitat * sample, metadata_16s_exp1, distance = "bray")
fullModel_16s_exp1_ev <- capscale(t(otus_16s_exp1_ev) ~ year * site, metadata_16s_exp1, distance = "bray")
fullModel_16s_exp1_ev <- capscale(t(otus_16s_exp1_ev) ~ year * habitat, metadata_16s_exp1, distance = "bray")
fullModel_16s_exp1_ev <- capscale(t(otus_16s_exp1_ev) ~ year + habitat + year:habitat, metadata_16s_exp1, distance = "bray")
anova(fullModel_16s_exp1_ev, permutations = 9999)

# backward selection
parsimoniousModel_16s_exp1_ev <- ordistep(fullModel_16s_exp1_ev, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
parsimoniousModel_16s_exp1_ev

# extract the parsimonious model formula for later use
FormulaParsimoniousModel_16s_exp1_ev <- formula(parsimoniousModel_16s_exp1_ev)[1:3]
FormulaParsimoniousModel_16s_exp1_ev

# performing capscale analysis
capscale(FormulaParsimoniousModel_16s_exp1_ev, metadata_16s_exp1, distance = "bray") -> CAP
head(summary(CAP))
anova(CAP, permutations = 9999)

RsquareAdj(CAP)$adj.r.squared



cap_ITS_exp1_ev = ordinate(biom_ITS_exp1_ev, "CAP", "bray", ~ year + habitat + year:habitat)
cap_16s_exp1_ev = ordinate(biom_16s_exp1_ev, "CAP", "bray", ~ year + habitat + year:habitat)

plot_cap_ITS_exp1_ev = plot_ordination(biom_ITS_exp1_ev, cap_ITS_exp1_ev, color="year", shape="habitat") + 
  geom_point(size=2.5, alpha=0.9, aes(shape=habitat)) +
  scale_shape_manual(values = c(17, 16), labels=c("DF","PS")) +
  theme_bw() +
  scale_colour_manual(values=palette_CB6) +
  theme(legend.position="right")

plot_cap_ITS_exp1_ev


plot_cap_16s_exp1_ev = plot_ordination(biom_16s_exp1_ev, cap_16s_exp1_ev, color="year", shape="habitat") + 
  geom_point(size=2.5, alpha=0.9, aes(shape=habitat)) +
  scale_shape_manual(values = c(17, 16), labels=c("DF","PS")) +
  theme_bw() +
  scale_colour_manual(values=palette_CB6) +
  scale_x_continuous() +
  theme(legend.position="right")

plot_cap_16s_exp1_ev


# > FIGURE 3 -------------------------------------------------
library("ggpubr")
ggarrange(plot_nmds_16s_exp1_ev, plot_cap_ITS_exp1_ev,
          plot_nmds_16s_exp1_ev, plot_cap_16s_exp1_ev,
          labels = c("A", "B","C","D"),
          widths = c(1,1,1,1),
          align = "hv", ncol = 2, nrow = 2,
          legend="right", common.legend = TRUE)



#  PERMANOVA  -------------------------------------------------------

adonis_ITS_exp1_ev <- adonis(t(otus_ITS_exp1_ev) ~ year * habitat * sample, data=metadata_ITS_exp1, permutations=9999)
adonis_ITS_exp1_ev <- adonis(t(otus_ITS_exp1_ev) ~ year * site, data=metadata_ITS_exp1, permutations=9999)
adonis_ITS_exp1_ev <- adonis(t(otus_ITS_exp1_ev) ~ year * habitat, data=metadata_ITS_exp1, permutations=9999)
adonis_ITS_exp1_ev

adonis_16s_exp1_ev <- adonis(t(otus_16s_exp1_ev) ~ year * habitat * sample, data=metadata_16s_exp1, permutations=9999)
adonis_16s_exp1_ev <- adonis(t(otus_16s_exp1_ev) ~ year * site, data=metadata_16s_exp1, permutations=9999)
adonis_16s_exp1_ev <- adonis(t(otus_16s_exp1_ev) ~ year * habitat, data=metadata_16s_exp1, permutations=9999)
adonis_16s_exp1_ev


# Test the homogenity of group variances ------------------------------
vegan::vegdist(t(otus_ITS_exp1_ev), method="bray") -> dist_ITS_exp1_ev

permdisp_dist_ITS_exp1_ev_year <- betadisper(dist_ITS_exp1_ev, metadata_ITS_exp1$year)
permdisp_dist_ITS_exp1_ev_habitat <- betadisper(dist_ITS_exp1_ev, metadata_ITS_exp1$habitat)

permdisp_dist_ITS_exp1_ev_year
anova(permdisp_dist_ITS_exp1_ev_year, permutations = 9999)
permutest(permdisp_dist_ITS_exp1_ev_year, permutations = 9999, pairwise = T)

permdisp_dist_ITS_exp1_ev_habitat
anova(permdisp_dist_ITS_exp1_ev_habitat, permutations = 9999)
permutest(permdisp_dist_ITS_exp1_ev_habitat, permutations = 9999, pairwise = T)

# FIGURE S10 -----------------------------------------------------
par(mfrow=c(2,3)) # Change the panel layout to 2 x 3
plot(permdisp_dist_ITS_exp1_ev_year, main="PCoA (year)", las=1)
boxplot(permdisp_dist_ITS_exp1_ev_year, main="Distance from centroid\n distribution (year)", xlab="year", las=2)
plot(TukeyHSD(permdisp_ITS_pop_df_ev_year), las=1)
plot(permdisp_dist_ITS_exp1_ev_habitat, main="PCoA (habitat)", las=1)
boxplot(permdisp_dist_ITS_exp1_ev_habitat, main="Distance from centroid\n distribution (habitat)", las=2)
plot(TukeyHSD(permdisp_ITS_pop_df_ev_treat), las=1)
dev.off()
par(mfrow=c(1,1)) # Change back to 1 x 1



vegan::vegdist(t(otus_16s_exp1_ev), method="bray") -> dist_16s_exp1_ev

permdisp_dist_16s_exp1_ev_year <- betadisper(dist_16s_exp1_ev, metadata_16s_exp1$year)
permdisp_dist_16s_exp1_ev_habitat <- betadisper(dist_16s_exp1_ev, metadata_16s_exp1$habitat)

permdisp_dist_16s_exp1_ev_year
anova(permdisp_dist_16s_exp1_ev_year, permutations = 9999)
permutest(permdisp_dist_16s_exp1_ev_year, permutations = 9999, pairwise = T)

permdisp_dist_16s_exp1_ev_habitat
anova(permdisp_dist_16s_exp1_ev_habitat, permutations = 9999)
permutest(permdisp_dist_16s_exp1_ev_habitat, permutations = 9999, pairwise = T)

# FIGURE S11 ----------------------------------------------------------
par(mfrow=c(2,3)) 
plot(permdisp_dist_16s_exp1_ev_year, main="PCoA (year)", las=1)
boxplot(permdisp_dist_16s_exp1_ev_year, main="Distance from centroid\n distribution (year)", xlab="year", las=2)
plot(TukeyHSD(permdisp_16s_pop_df_ev_year), las=1)
plot(permdisp_dist_16s_exp1_ev_habitat, main="PCoA (habitat)", las=1)
boxplot(permdisp_dist_16s_exp1_ev_habitat, main="Distance from centroid\n distribution (habitat)", las=2)
plot(TukeyHSD(permdisp_16s_pop_df_ev_treat), las=1)
dev.off()
par(mfrow=c(1,1)) 



# plotting HEATMAPS -------------------------------------------

# order level
biom_ITS_exp1_order = tax_glom(biom_ITS_exp1, "Order") 
biom_ITS_exp1_order

otu_table(biom_ITS_exp1_order) <- otu_table(biom_ITS_exp1_order)[which
                                  (rowSums(otu_table(biom_ITS_exp1_order)) > 50),] # ~>1% of the most abundant OTU in a single sample 
biom_ITS_exp1_order

# manual correction of the tax_table to include fungi in the culture
#write.csv(tax_table(biom_ITS_exp1_order) , file = "tax_ITS_exp1_order.csv") 
#write.csv(otu_table(biom_ITS_exp1_order) , file = "otus_ITS_exp1_order.csv")

tax_ITS_exp1_order <- read.csv("tax_ITS_exp1_order.csv", header=T, row.names=1)
tax_table(biom_ITS_exp1_order) <- tax_table(as.matrix(tax_ITS_exp1_order))
biom_ITS_exp1_order

otu_ITS_exp1_ord = as.data.frame(rowSums(otu_table(biom_ITS_exp1_order)))
colnames(otu_ITS_exp1_ord) <- "counts"
otu_ITS_exp1_ord$OTU_ID <- row.names(otu_ITS_exp1_ord)
otu_ITS_exp1_ord$Order <- as.data.frame(tax_table(biom_ITS_exp1_order))$Order
otu_ITS_exp1_ord <- arrange(otu_ITS_exp1_ord, counts) # to reorder by sample desc(counts)
head(otu_ITS_exp1_ord)
tail(otu_ITS_exp1_ord)
otu_ITS_exp1_ord

# > FUGURE 4 ---------------------------------
heat_ITS_ord <- phyloseq::plot_heatmap(biom_ITS_exp1_order, sample.order = "year", sample.label = "year", 
                                       taxa.label = "Order", low="#000033", high="#66CCFF", na.value="light grey",
                                       taxa.order=otu_ITS_uparse_pop_df_ord$OTU_ID, trans=log_trans(2)) +
  facet_grid(~habitat+site, labeller = label_bquote(cols = .(site)), scales = "free_x", space="free_x") +
  theme(axis.text.y = element_text(size = 6)) +
  theme(axis.text.x = element_text(size = 6)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(legend.position="bottom") 

heat_ITS_ord


# genus level - dropping out unclassified taxa
tax_table(biom_ITS_exp1)[tax_table(biom_ITS_exp1)=="Unclassified"]<- NA

biom_ITS_exp1_genus = tax_glom(biom_ITS_exp1, "Genus")
otu_table(biom_ITS_exp1_genus) <- otu_table(biom_ITS_exp1_genus)[which
                                  (rowSums(otu_table(biom_ITS_exp1_genus)) > 50),]

write.csv(tax_table(biom_ITS_exp1_genus), file = "tax_ITS_exp1_genus.csv")
write.csv(otu_table(biom_ITS_exp1_genus), file = "otus_ITS_exp1_genus.csv")

otu_ITS_exp1_genus <- read.csv("otus_ITS_exp1_genus.csv", header=T, row.names=1)
otu_table(biom_ITS_exp1_genus) <- otu_table(as.matrix(otu_ITS_exp1_genus), taxa_are_rows=TRUE)

otu_ITS_exp1_genus = as.data.frame(rowSums(otu_table(biom_ITS_exp1_genus)))
colnames(otu_ITS_exp1_genus) <- "counts"
otu_ITS_exp1_genus$OTU_ID <- row.names(otu_ITS_exp1_genus)
otu_ITS_exp1_genus$Order <- as.data.frame(tax_table(biom_ITS_exp1_genus))$Order
otu_ITS_exp1_genus <- arrange(otu_ITS_exp1_genus, counts) # to reorder by sample desc(counts)
head(otu_ITS_exp1_genus)
tail(otu_ITS_exp1_genus)
otu_ITS_exp1_genus

# FIGURE S13 ----------------------------
heat_ITS_GEN <- phyloseq::plot_heatmap(biom_ITS_exp1_genus, sample.order = "year", sample.label = "year", 
                                       taxa.label = "Genus", low="#000033", high="#66CCFF", na.value="light grey",
                                       taxa.order=otu_ITS_exp1_genus$OTU_ID, trans=log_trans(2)) +
  facet_grid(~habitat+site, labeller = label_bquote(cols = .(site)), scales = "free_x", space="free_x") +
  theme(axis.text.y = element_text(face = "italic", size = 3)) +
  theme(axis.text.x = element_text(size = 6)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(legend.position="bottom") 

heat_ITS_GEN


tax_table(biom_16s_exp1)[tax_table(biom_16s_exp1)=="Unclassified"]<- NA

biom_16s_exp1_genus = tax_glom(biom_16s_exp1, "Genus")
otu_table(biom_16s_exp1_genus) <- otu_table(biom_16s_exp1_genus)[which
                                  (rowSums(otu_table(biom_16s_exp1_genus)) > 50),]

write.csv(tax_table(biom_16s_exp1_genus), file = "tax_16s_exp1_genus.csv")
write.csv(otu_table(biom_16s_exp1_genus), file = "otus_16s_exp1_genus.csv")

otu_16s_exp1_genus <- read.csv("FunGuild/otus_16s_uparse_pop_df_genus_new.csv", header=T, row.names=1)
otu_table(biom_16s_exp1_genus) <- otu_table(as.matrix(otu_16s_exp1_genus), taxa_are_rows=TRUE)

otu_16s_exp1_genus = as.data.frame(rowSums(otu_table(biom_16s_exp1_genus)))
colnames(otu_16s_exp1_genus) <- "counts"
otu_16s_exp1_genus$OTU_ID <- row.names(otu_16s_exp1_genus)
otu_16s_exp1_genus$Order <- as.data.frame(tax_table(biom_16s_exp1_genus))$Order
otu_16s_exp1_genus <- arrange(otu_16s_exp1_genus, counts) # to reorder by sample desc(counts)
head(otu_16s_exp1_genus)
tail(otu_16s_exp1_genus)
otu_16s_exp1_genus

# FIGURE S15 -----------------------------------------
heat_16s_GEN <- phyloseq::plot_heatmap(biom_16s_exp1_genus, sample.order = "year", sample.label = "year", 
                                       taxa.label = "Genus", low="#000033", high = "#66FF66", na.value="light grey",
                                       taxa.order=otu_16s_exp1_genus$OTU_ID, trans=log_trans(2)) +
  facet_grid(~habitat+site, labeller = label_bquote(cols = .(site)), scales = "free_x", space="free_x") +
  theme(axis.text.y = element_text(face = "italic", size = 4)) +
  theme(axis.text.x = element_text(size = 6)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(legend.position="bottom") 

heat_16s_GEN


# reimporting guilds form FUNGuilds ---------------------------------
guilds_ITS <- read.table("guilds.csv", row.names=1, sep=",", header=T)
guilds_ITS
dim(guilds_ITS)

identical(rownames(tax_table(biom_ITS_exp1_genus)), rownames(guilds_ITS))

tax_table(biom_ITS_exp1_genus) <- cbind(tax_table(biom_ITS_exp1_genus),as.matrix(guilds_ITS))
tax_table(biom_ITS_exp1_genus)

## Add OTU labels for retriving the identified Genus
tax_ITS_exp1_genus <- as.data.frame(as.matrix(tax_table(biom_ITS_exp1_genus)))
head(tax_ITS_exp1_genus)

tax_table(biom_ITS_exp1_genus) <- cbind(tax_table(biom_ITS_exp1_genus),
                                  as.matrix(paste(tax_ITS_exp1_genus$Genus," (",rownames(tax_ITS_exp1_genus),")",sep="")))

colnames(tax_table(biom_ITS_exp1_genus)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                                                    "guilds", "OTU_ID", "OTU_labels")
head(tax_table(biom_ITS_exp1_genus))

# extracting only Myrrhizal Fungi -----------------
biom_ITS_exp1_MF <- subset_taxa(biom_ITS_exp1_genus, guilds%in%c("ECM","VAM"))
tax_table(biom_ITS_exp1_MF)

# ranking abundances to reorder -------------------
library(dplyr)

otu_ITS_exp1_reordered_MF = as.data.frame(rowSums(otu_table(biom_ITS_exp1_MF)))
colnames(otu_ITS_exp1_reordered_MF) <- "counts"
otu_ITS_exp1_reordered_MF$OTU_ID <- row.names(otu_ITS_exp1_reordered_MF)
otu_ITS_exp1_reordered_MF$Genus <- as.data.frame(tax_table(biom_ITS_exp1_MF))$Genus
otu_ITS_exp1_reordered_MF <- arrange(otu_ITS_exp1_reordered_MF, counts) # to reorder by sample desc(counts)
head(otu_ITS_exp1_reordered_MF)
tail(otu_ITS_exp1_reordered_MF)
otu_ITS_exp1_reordered_MF

#biom_ITS_uparse_pop_df_gn_test_ECM -> biom_ITS_exp1_MF
#colnames(sample_data(biom_ITS_exp1_MF)) <- c("X.SampleID","site","year","habitat")

# > FIGURE 5 ---------------------------------------------
library(scales)

heat_ITS_myco = phyloseq::plot_heatmap(biom_ITS_exp1_MF, sample.order = "year", sample.label = "year", 
                taxa.label = "Genus", low="#000033", high="#FF3300", na.value="light grey",
                taxa.order=otu_ITS_exp1_reordered_MF$OTU_ID, trans=log_trans(2)) + 
  facet_grid(~habitat+site, labeller = label_bquote(cols = .(site)), scales = "free_x", space="free_x") +
  theme(axis.text.y = element_text(face = "italic", size = 6)) +
  theme(axis.text.x = element_text(size = 6)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(legend.position="bottom") 

heat_ITS_myco

ggsave("heat_myco_test.tiff", heat_ITS_myco, dpi=900, width=6)

tax_table(biom_ITS_exp1_MF)
otu_table(biom_ITS_exp1_MF)
max(otu_table(biom_ITS_exp1_MF))


# heatmap for prokaryotes 
tax_table(biom_16s_exp1)
tax_table(biom_16s_exp1)[tax_table(biom_16s_exp1)=="Unclassified"]<- NA

biom_16s_exp1_class = tax_glom(biom_16s_exp1, "Class")
otu_table(biom_16s_exp1_class) <- otu_table(biom_16s_exp1_class)[which
                                  (rowSums(otu_table(biom_16s_exp1_class)) > 50),] 
biom_16s_exp1_class
tax_table(biom_16s_exp1_class)

otu_16s_exp1_class = as.data.frame(rowSums(otu_table(biom_16s_exp1_class)))
colnames(otu_16s_exp1_class) <- "counts"
otu_16s_exp1_class$OTU_ID <- row.names(otu_16s_exp1_class)
otu_16s_exp1_class$Class <- as.data.frame(tax_table(biom_16s_exp1_class))$Class
otu_16s_exp1_class <- arrange(otu_16s_exp1_class, counts) # to reorder by sample desc(counts)
head(otu_16s_exp1_class)
otu_16s_exp1_class


# FIGURE S14 ------------------------------------------------------
heat_16s_Class = phyloseq::plot_heatmap(biom_16s_exp1_class, sample.order = "year", sample.label = "year", 
                                        taxa.label = "Class", low="#000033", high = "#66FF66", na.value="light grey",
                                        taxa.order=otu_16s_exp1_class$OTU_ID, trans=log_trans(2)) +
  facet_grid(~habitat+site, labeller = label_bquote(cols = .(site)), scales = "free_x", space="free_x") +
  theme(axis.text.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(size = 6)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(legend.position="bottom") 

heat_16s_Class



# OTHER FIGURES -----------------------------------------

#colnames(sample_data(biom_ITS_exp1_PAT)) <- c("X.SampleID","site","year","habitat")
biom_ITS_exp1_PAT_PAT <- subset_taxa(biom_ITS_exp1_genus, guilds%in%c("PAT"))
tax_table(biom_ITS_exp1_PAT)
sample_data(biom_ITS_exp1_PAT)

otu_ITS_exp1_reordered_PAT = as.data.frame(rowSums(otu_table(biom_ITS_exp1_PAT)))
colnames(otu_ITS_exp1_reordered_PAT) <- "counts"
otu_ITS_exp1_reordered_PAT$OTU_ID <- row.names(otu_ITS_exp1_reordered_PAT)
otu_ITS_exp1_reordered_PAT$Genus <- as.data.frame(tax_table(biom_ITS_exp1_PAT))$Genus
otu_ITS_exp1_reordered_PAT <- arrange(otu_ITS_exp1_reordered_PAT, counts) # to reorder by sample desc(counts)
head(otu_ITS_exp1_reordered_PAT)
tail(otu_ITS_exp1_reordered_PAT)
otu_ITS_exp1_reordered_PAT

heat_ITS_patho = phyloseq::plot_heatmap(biom_ITS_exp1_PAT, sample.order = "year", sample.label = "year", 
                                        taxa.label = "Genus", low = "#ff0000", high = "#15ff00", na.value="grey",
                                        taxa.order=otu_ITS_exp1_reordered_PAT$OTU_ID, trans=log_trans(2)) +
  facet_grid(~habitat+site, labeller = label_bquote(cols = .(site)), scales = "free_x", space="free_x") +
  theme(axis.text.y = element_text(face = "italic", size = 6)) +
  theme(axis.text.x = element_text(size = 6)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(legend.position="bottom") 

heat_ITS_patho

biom_ITS_exp1_SAP <- subset_taxa(biom_ITS_exp1_genus, guilds%in%c("SAP"))
tax_table(biom_ITS_exp1_SAP)
sample_data(biom_ITS_exp1_SAP)

otu_ITS_exp1_reordered_SAP = as.data.frame(rowSums(otu_table(biom_ITS_exp1_SAP)))
colnames(otu_ITS_exp1_reordered_SAP) <- "counts"
otu_ITS_exp1_reordered_SAP$OTU_ID <- row.names(otu_ITS_exp1_reordered_SAP)
otu_ITS_exp1_reordered_SAP$Genus <- as.data.frame(tax_table(biom_ITS_exp1_SAP))$Genus
otu_ITS_exp1_reordered_SAP <- arrange(otu_ITS_exp1_reordered_SAP, counts)
head(otu_ITS_exp1_reordered_SAP)
tail(otu_ITS_exp1_reordered_SAP)
otu_ITS_exp1_reordered_SAP

heat_ITS_sapro = phyloseq::plot_heatmap(biom_ITS_exp1_SAP, sample.order = "year", sample.label = "year", 
                                        taxa.label = "Genus", low = "#ff0000", high = "#15ff00", na.value="grey",
                                        taxa.order=otu_ITS_exp1_reordered_SAP$OTU_ID, trans=log_trans(2)) +
  facet_grid(~habitat+site, labeller = label_bquote(cols = .(site)), scales = "free_x", space="free_x") +
  theme(axis.text.y = element_text(face = "italic", size = 6)) +
  theme(axis.text.x = element_text(size = 6)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(legend.position="bottom") 

heat_ITS_sapro



# Culturable fungi -------------

read.csv("fungi_cultures.csv", header = TRUE, row.names=1, sep = ",") -> fungi_count
fungi_count

fungi_count$year <- as.factor(fungi_count$year)

# FIGURE S16 -------------------------
ggplot(fungi_count, aes(x=year, y=counts, fill= habitat)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=3) +
  theme_classic() +
  theme(legend.position = "bottom")
  









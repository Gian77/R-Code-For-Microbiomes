# *******Script for data analyses ******************************** -----
# Manuscript:   Microbial Communities Associated with Cultivated Morel (Morchella spp.) Mushrooms
# Authors:      Longley R, Benucci GMN, ... Bonito G.
# Affiliation:  Michigan State University
# Journal:      FEMS Microbiology Ecology 
# Date:         January 23, 2019
# ******************************************************************** -----

# WORKING ENVIRONMENT SETUP ------------------------------------------------------------
options(scipen = 999) #to use decimals
options(max.print=100000000) # to print more lines on the display
options(verbose=TRUE)

# >>> importing datasets --------------------------------------------------------------
library(ggplot2)
library(phyloseq)
library(vegan)
library(dplyr)
library(Biostrings)
library(ape)


# Import data ------------
otu_table_fungi <- read.csv("Datasets//otu_table_fungi.csv", header=T, row.names =1)
otu_table_fungi_phy <-otu_table(otu_table_fungi,taxa_are_rows = TRUE)

tax_table_fungi <- read.csv("Datasets//tax_table_fungi.csv", header=T, row.names =1)
tax_table_fungi_phy <- tax_table(as.matrix(tax_table_fungi))

sample_data_fungi <- read.csv("Datasets/sample_data_fungi.csv", header=T, row.names =1)
sample_data_fungi_phy <-sample_data(sample_data_fungi)

otus_fungi <- readDNAStringSet("Datasets/otus_fungi.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

# create a phyloseq object
physeq_fungi <- merge_phyloseq(otu_table_fungi_phy, 
                              tax_table_fungi_phy,
                              sample_data_fungi_phy,
                              otus_fungi)

physeq_fungi
sample_data(physeq_fungi)
otu_table(physeq_fungi)
tax_table(physeq_fungi)


otu_table_prokaryote <- read.csv("Datasets//otu_table_prokaryote.csv", header=T, row.names =1)
otu_table_prokaryote_phy <-otu_table(otu_table_prokaryote,taxa_are_rows = TRUE)

tax_table_prokaryote <- read.csv("Datasets//tax_table_prokaryote.csv", header=T, row.names =1)
tax_table_prokaryote_phy <- tax_table(as.matrix(tax_table_prokaryote))

sample_data_prokaryote <- read.csv("Datasets/sample_data_prokaryote.csv", header=T, row.names =1)
sample_data_prokaryote_phy <-sample_data(sample_data_prokaryote)

otus_prokaryote <- readDNAStringSet("Datasets/otus_prokaryote.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

# create a phyloseq object
physeq_prokaryote <- merge_phyloseq(otu_table_prokaryote_phy, 
                                    tax_table_prokaryote_phy,
                                    sample_data_prokaryote_phy,
                                    otus_prokaryote)

physeq_prokaryote
sample_data(physeq_prokaryote)
otu_table(physeq_prokaryote)
tax_table(physeq_prokaryote)


# extracting MOCK positive samples ---------------
physeq_fungi -> biom_ITS_uparse

biom_ITS_mock <- subset_samples(biom_ITS_uparse, Description%in%c("MOCK4"))
otu_table(biom_ITS_mock) <- otu_table(biom_ITS_mock)[which(rowSums(otu_table(biom_ITS_mock)) >= 1),] 
biom_ITS_mock
otu_table(biom_ITS_mock)
tax_table(biom_ITS_mock)
write.table(refseq(biom_ITS_mock), "MOCK1_Morel_ITS.fasta")

physeq_prokaryote -> biom_16s_uparse

biom_16s_mock <- subset_samples(biom_16s_uparse, Description%in%c("MOCK1"))
otu_table(biom_16s_mock) <- otu_table(biom_16s_mock)[which(rowSums(otu_table(biom_16s_mock)) >= 1),] 
biom_16s_mock
otu_table(biom_16s_mock)
tax_table(biom_16s_mock)
write.table(refseq(biom_16s_mock), "MOCK1_Morel_16s.fasta")


# >>> Filtering out OTUs ------------------------------
# Oliver et al. 2015, PCR errors and tag switching
# Lindhal et al. 2013, tag switching - that's a good  one!'''

# Barberan et al. 2012, removing OTUs that appear in less than x samples


biom_ITS_uparse -> biom_ITS_uparse_filt
#otu_table(biom_ITS_uparse_filt)[otu_table(biom_ITS_uparse_filt) <= 4] <- 0 ### tag switching
#otu_table(biom_16s_qc) <- otu_table(biom_16s_qc)[rowSums(otu_table(biom_16s_qc) > 0) >= 3, ] ### PCR errors  
otu_table(biom_ITS_uparse_filt) <- otu_table(biom_ITS_uparse_filt)[which(rowSums(otu_table(biom_ITS_uparse_filt)) >= 10),] ### PCR Errors 
biom_ITS_uparse_filt

biom_16s_uparse -> biom_16s_uparse_filt
#otu_table(biom_16s_uparse_filt)[otu_table(biom_16s_uparse_filt) <= 4 ] <- 0 ### tag switching
#otu_table(biom_16s_qc) <- otu_table(biom_16s_qc)[rowSums(otu_table(biom_16s_qc) > 0) >= 3, ] ### PCR errors  
otu_table(biom_16s_uparse_filt) <- otu_table(biom_16s_uparse_filt)[which(rowSums(otu_table(biom_16s_uparse_filt)) >= 10),] ### PCR Errors 
biom_16s_uparse_filt


# function to remove dab taxa -----------
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

badTaxa_ITS =c("OTU_17", "OTU_51")

biom_ITS_uparse_filt = remove_taxa(biom_ITS_uparse_filt, badTaxa_ITS)
biom_ITS_uparse_filt
sample_data(biom_ITS_uparse_filt)


badTaxa_16s =c("OTU_260","OTU_1933", "OTU_245","OTU_1155","OTU_2789","OTU_2693","OTU_9938","OTU_2645",
               "OTU_348","OTU_1790","OTU_414","OTU_2021","OTU_11","OTU_3102","OTU_5590","OTU_1530",
               "OTU_111","OTU_811","OTU_78","OTU_9041","OTU_104","OTU_8262","OTU_2701","OTU_5","OTU_150",
               "OTU_112","OTU_68","OTU_460","OTU_161","OTU_786","OTU_883","OTU_532","OTU_1149","OTU_889",
               "OTU_1244","OTU_482","OTU_1001","OTU_1172","OTU_692","OTU_957","OTU_1127","OTU_4226",
               "OTU_601","OTU_977")

biom_16s_uparse_filt = remove_taxa(biom_16s_uparse_filt, badTaxa_16s)
biom_16s_uparse_filt
sample_data(biom_16s_uparse_filt)

# > subsetting datasets -------------------------------------

# Morel ITS outdoor ------
biom_ITS_uparse_filt
sample_data(biom_ITS_uparse_filt)

biom_ITS_uparse_out <- subset_samples(biom_ITS_uparse_filt, Study%in%c("Outdoor"))
otu_table(biom_ITS_uparse_out) <- otu_table(biom_ITS_uparse_out)[which(rowSums(otu_table(biom_ITS_uparse_out)) >= 1),] 
biom_ITS_uparse_out
sample_data(biom_ITS_uparse_out)

count(sample_data(biom_ITS_uparse_out), vars = Growhouse)

# Morel ITS indoor -------
biom_ITS_uparse_in <- subset_samples(biom_ITS_uparse_filt, Study%in%c("Indoor"))
otu_table(biom_ITS_uparse_in) <- otu_table(biom_ITS_uparse_in)[which(rowSums(otu_table(biom_ITS_uparse_in)) >= 1),] 
biom_ITS_uparse_in
sample_data(biom_ITS_uparse_in)


# Morel 16S outdoor ------
sample_data(biom_16s_uparse_filt)

biom_16s_uparse_out <- subset_samples(biom_16s_uparse_filt, Study%in%c("Outdoor"))
otu_table(biom_16s_uparse_out) <- otu_table(biom_16s_uparse_out)[which(rowSums(otu_table(biom_16s_uparse_out)) >= 1),] 
biom_16s_uparse_out
sample_data(biom_16s_uparse_out)

count(sample_data(biom_16s_uparse_out), vars = Stage)
count(sample_data(biom_16s_uparse_out), vars = Origin)
count(sample_data(biom_16s_uparse_out), vars = Growhouse)

# Morel 16S indoor ------
biom_16s_uparse_in <- subset_samples(biom_16s_uparse_filt, Study%in%c("Indoor"))
otu_table(biom_16s_uparse_in) <- otu_table(biom_16s_uparse_in)[which(rowSums(otu_table(biom_16s_uparse_in)) >= 1),] 
biom_16s_uparse_in
sample_data(biom_16s_uparse_in)


# rarefy the whole dataset -------------------------------
set.seed(2018)
min(sample_sums(biom_ITS_uparse_out))
data.frame(colSums(otu_table(biom_ITS_uparse_out)))

biom_ITS_uparse_out_ev = rarefy_even_depth(biom_ITS_uparse_out, 
                                           sample.size = min(sample_sums(biom_ITS_uparse_out)), 
                                                 rngseed = FALSE, replace = TRUE, 
                                                 trimOTUs = TRUE, verbose = TRUE) 
otu_table(biom_ITS_uparse_out_ev) <- otu_table(biom_ITS_uparse_out_ev)[which
                                        (rowSums(otu_table(biom_ITS_uparse_out_ev)) >= 1),]
biom_ITS_uparse_out_ev
colSums(otu_table(biom_ITS_uparse_out_ev))
any(taxa_sums(biom_ITS_uparse_out_ev) == 0)


data.frame(colSums(otu_table(biom_ITS_uparse_in)))
biom_ITS_uparse_in_ev = rarefy_even_depth(biom_ITS_uparse_in,
                          sample.size = min(sample_sums(biom_ITS_uparse_in)),
                          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_ITS_uparse_in_ev
colSums(otu_table(biom_ITS_uparse_in_ev))
any(taxa_sums(biom_ITS_uparse_in_ev) == 0)


min(sample_sums(biom_16s_uparse_out))
data.frame(colSums(otu_table(biom_16s_uparse_out)))

biom_16s_uparse_out_ev = rarefy_even_depth(biom_16s_uparse_out, sample.size = min(sample_sums(biom_16s_uparse_out)), 
                                           rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) 

otu_table(biom_16s_uparse_out_ev) <- otu_table(biom_16s_uparse_out_ev)[which
                                      (rowSums(otu_table(biom_16s_uparse_out_ev)) >= 1),]
biom_16s_uparse_out_ev
colSums(otu_table(biom_16s_uparse_out_ev))

data.frame(colSums(otu_table(biom_16s_uparse_in)))
biom_16s_uparse_in_ev = rarefy_even_depth(biom_16s_uparse_in,
                        sample.size = min(sample_sums(biom_16s_uparse_in)),
                        rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_16s_uparse_in_ev
colSums(otu_table(biom_16s_uparse_in_ev))
any(taxa_sums(biom_16s_uparse_in_ev) == 0)



# >>> STACKED BARPLOTS  - example with with the fungal dataset -----------------------------------
source("Rscript_2_plot_ordered_bar.R")

# Fungi Outdoor G1 ---------------

biom_ITS_uparse_out_ev_g1 <- subset_samples(biom_ITS_uparse_out_ev, Growhouse%in%c("G1"))
otu_table(biom_ITS_uparse_out_ev_g1) <- otu_table(biom_ITS_uparse_out_ev_g1)[which(rowSums(otu_table(biom_ITS_uparse_out_ev_g1)) >= 1),] 
biom_ITS_uparse_out_ev_g1
sample_data(biom_ITS_uparse_out_ev_g1)

biom_ITS_uparse_out_ev_g1_genus = tax_glom(biom_ITS_uparse_out_ev_g1, "Genus")
biom_ITS_uparse_out_ev_g1_genus_stage = merge_samples(biom_ITS_uparse_out_ev_g1_genus, "Stage")
sample_data(biom_ITS_uparse_out_ev_g1_genus_stage)$Stage <- levels(sample_data(biom_ITS_uparse_out_ev_g1)$Stage)
biom_ITS_uparse_out_ev_g1_genus_stage_ab = transform_sample_counts(biom_ITS_uparse_out_ev_g1_genus_stage, function(x) 100 * x/sum(x))
sample_data(biom_ITS_uparse_out_ev_g1_genus_stage_ab)

top_15_ITS_out_g1 <- names(sort(taxa_sums(biom_ITS_uparse_out_ev_g1_genus_stage_ab), TRUE)[1:15])
top_15_ITS_out_g1
biom_ITS_uparse_out_ev_g1_genus_stage_ab_15 <- prune_taxa(top_15_ITS_out_g1, biom_ITS_uparse_out_ev_g1_genus_stage_ab)
biom_ITS_uparse_out_ev_g1_genus_stage_ab_15
sample_data(biom_ITS_uparse_out_ev_g1_genus_stage_ab_15)

sample_data(biom_ITS_uparse_out_ev_g1_genus_stage_ab_15)$Stage <- factor(sample_data(biom_ITS_uparse_out_ev_g1_genus_stage_ab_15)$Stage, 
                                                                      levels=c("Young", "Mature"))

tax_table(biom_ITS_uparse_out_ev_g1_genus_stage_ab_15)


# creating plot --------

ITS_Outdoor_Stage_g1 = plot_ordered_bar(biom_ITS_uparse_out_ev_g1_genus_stage_ab_15, x = "Stage",fill="Genus", 
                                        leg_size = 0.1) +
  guides(fill=guide_legend(ncol=1)) +
  labs(title = "", x = "G1", y ="Relative Abundance")+
  theme_classic()+
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.text=element_text(size=8),
        legend.key.size = unit(1,"line"),
        legend.position = c("right"))

ITS_Outdoor_Stage_g1


# Fungi Outdoor G2 ---------------
biom_ITS_uparse_out_ev_g2 <- subset_samples(biom_ITS_uparse_out_ev, Growhouse%in%c("G2"))
otu_table(biom_ITS_uparse_out_ev_g2) <- otu_table(biom_ITS_uparse_out_ev_g2)[which(rowSums(otu_table(biom_ITS_uparse_out_ev_g2)) >= 1),] 
biom_ITS_uparse_out_ev_g2
sample_data(biom_ITS_uparse_out_ev_g2)

biom_ITS_uparse_out_ev_g2_genus = tax_glom(biom_ITS_uparse_out_ev_g2, "Genus")
biom_ITS_uparse_out_ev_g2_genus_stage = merge_samples(biom_ITS_uparse_out_ev_g2_genus, "Stage")
sample_data(biom_ITS_uparse_out_ev_g2_genus_stage)$Stage <- levels(sample_data(biom_ITS_uparse_out_ev_g2)$Stage)
biom_ITS_uparse_out_ev_g2_genus_stage_ab = transform_sample_counts(biom_ITS_uparse_out_ev_g2_genus_stage, function(x) 100 * x/sum(x))
sample_data(biom_ITS_uparse_out_ev_g2_genus_stage_ab)

top_15_ITS_out_g2 <- names(sort(taxa_sums(biom_ITS_uparse_out_ev_g2_genus_stage_ab), TRUE)[1:15])
top_15_ITS_out_g2
biom_ITS_uparse_out_ev_g2_genus_stage_ab_15 <- prune_taxa(top_15_ITS_out_g2, biom_ITS_uparse_out_ev_g2_genus_stage_ab)
biom_ITS_uparse_out_ev_g2_genus_stage_ab_15
sample_data(biom_ITS_uparse_out_ev_g2_genus_stage_ab_15)

sample_data(biom_ITS_uparse_out_ev_g2_genus_stage_ab_15)$Stage <- factor(sample_data(biom_ITS_uparse_out_ev_g2_genus_stage_ab_15)$Stage, 
                                                                         levels=c("Young", "Mature"))

tax_table(biom_ITS_uparse_out_ev_g2_genus_stage_ab_15)


# creating plot --------

ITS_Outdoor_Stage_g2 = plot_ordered_bar(biom_ITS_uparse_out_ev_g2_genus_stage_ab_15, x = "Stage",fill="Genus", 
                                        leg_size = 0.1) +
  guides(fill=guide_legend(ncol=1)) +
  labs(title = "", x = "G2", y ="Relative Abundance")+
  theme_classic()+
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.text=element_text(size=8),
        legend.key.size = unit(1,"line"),
        legend.position = c("right"))

ITS_Outdoor_Stage_g2


# Fungi Indoor by Stage variable ---------------

biom_ITS_uparse_in_genus = tax_glom(biom_ITS_uparse_in_ev, "Genus")
biom_ITS_uparse_in_genus= merge_samples(biom_ITS_uparse_in_genus, "Stage")
sample_data(biom_ITS_uparse_in_genus)$Stage <- levels(sample_data(biom_ITS_uparse_in_ev)$Stage)
biom_ITS_uparse_in_genus_ab = transform_sample_counts(biom_ITS_uparse_in_genus, 
                                                      function(x) 100 * x/sum(x))

top_15_ITS_in <- names(sort(taxa_sums(biom_ITS_uparse_in_genus_ab), TRUE)[1:15]) 
biom_ITS_uparse_in_genus_ab <- prune_taxa(top_15_ITS_in, biom_ITS_uparse_in_genus_ab)
biom_ITS_uparse_in_genus_ab
sample_data(biom_ITS_uparse_in_genus_ab)


sample_data(biom_ITS_uparse_in_genus_ab)$Stage <- factor(sample_data(biom_ITS_uparse_in_genus_ab)$Stage, levels=c("Primordia", "Differentiation", "Fundamental", "Late fruiting", "Late non fruiting"))

ITS_Indoor_Stage = plot_ordered_bar(biom_ITS_uparse_in_genus_ab, x = "Stage", fill="Genus", 
                                    leg_size = 0.1) +
  guides(fill=guide_legend(ncol=1)) +
  labs(title = "", x = "Indoor", y ="")+
  theme_classic()+
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.text=element_text(size=8),
        legend.key.size = unit(1,"line"),
        legend.position = c("right"))

ITS_Indoor_Stage


# *** FIGURE 1 ---------------------------------------
# NOTE. We corrected the taxonomy before plotting by BLASING OTUs rep seq in NCBI

library("ggpubr")
ggarrange(ITS_Outdoor_Stage_g1,
          ITS_Outdoor_Stage_g2,
          ITS_Indoor_Stage,
          labels = c("A", "B" ,"C"),
          widths = c(1.6, 1.6,2.0),
          align = "h", ncol = 3, nrow = 1)





# >>> ALPHA DIVERSITY -------------------------------------------

# OUTDOOR-INDOOR alpha Fungi ------------------------------------
library("ggpubr")

library("gridExtra")
library("grid")
library("cowplot")

label_names <- c(Observed="Richness", Shannon="Shannon")
label_names

alpha_fungi_out = plot_richness(biom_ITS_uparse_out_ev, x="Growhouse", 
                                color="Stage", measures = c("Shannon", "Observed")) +
  theme_pubclean() +
  expand_limits(x = 0, y = 0) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 0, shape = 16) +
  labs(title="", x="", y = "Alpha Diversity Measure") +
  scale_x_discrete("Growhouse", labels = c("Growhouse 1" = "G1","Growhouse 2" = "G2")) +
  scale_colour_manual("Substrate",breaks = c("Mature", "Young"),
                      values = c("Mature"="green", "Young"="red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())

# alpha_fungi_out ---------
alpha_fungi_out
alpha_fungi_out$layers
alpha_fungi_out$layers <- alpha_fungi_out$layers[-1]
alpha_fungi_out$layers <- alpha_fungi_out$layers[-3]

# rename facets on the plot 
levels(alpha_fungi_out$data$variable)
levels(alpha_fungi_out$data$variable)[levels
      (alpha_fungi_out$data$variable)=="Observed"] <- "Richness"


# alpha_fungi_in ------------
levels(sample_data(biom_ITS_uparse_in_ev)$Stage)[levels
      (sample_data(biom_ITS_uparse_in_ev)$Stage)=="Late_fruiting"] <- "Late fruiting"

levels(sample_data(biom_ITS_uparse_in_ev)$Stage)[levels
      (sample_data(biom_ITS_uparse_in_ev)$Stage)=="Late_non_fruiting"] <- "Late non fruiting"


sample_data(biom_ITS_uparse_in_ev)$Stage <- factor(sample_data(biom_ITS_uparse_in_ev)$Stage,
      levels=c("Primordia","Fundamental","Differentiation","Late fruiting","Late non fruiting"))


alpha_fungi_in = plot_richness(biom_ITS_uparse_in_ev, x="Stage", color = "Stage",
                               measures = c("Shannon", "Observed")) +
  theme_pubclean() +
  geom_point(size = 0, shape = 16) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  expand_limits(x = 0, y = 0) +
  labs(title="", x="Stage", y = "") +
  scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 10)) +
  theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="none") +
  theme(legend.title=element_blank())

alpha_fungi_in
alpha_fungi_in$layers
alpha_fungi_in$layers <- alpha_fungi_in$layers[-1]
alpha_fungi_in$layers <- alpha_fungi_in$layers[-1]


# rename facets on the plot 
levels(alpha_fungi_in$data$variable)
levels(alpha_fungi_in$data$variable)[levels
    (alpha_fungi_in$data$variable)=="Observed"] <- "Richness"


alpha_fungi_out_A = ggarrange(alpha_fungi_out,
                              labels = "A",
                              widths = 1,
                              ncol = 1,
                              nrow = 1,
                              legend = "bottom")
alpha_fungi_out_A

alpha_fungi_out_B = ggarrange(alpha_fungi_in,
                              labels = "B",
                              widths = 1,
                              ncol = 1,
                              nrow = 1,
                              legend = "none")
alpha_fungi_out_B

# *** FIGURE 3 ---------------------------------
library("ggpubr")

ggarrange(alpha_fungi_out_A,
          alpha_fungi_out_B,
          widths = c(1, 1),
          align = "v",
          ncol = 2, nrow = 1, 
          legend = c("bottom"))


# OUTDOOR- INDOOR alpha Prokaryotes -------------

# alpha_16s_out1 -----------
alpha_16s_out1 = plot_richness(biom_16s_uparse_out_ev, x="Growhouse",
                               color="Stage", measures = c("Shannon", "Observed")) + 
  theme_pubclean() +
  expand_limits(x = 0, y = 0) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  labs(title="", x="Growhouse", y = "Alpha Diversity Measure") +
  scale_x_discrete("Growhouse", labels = c("Growhouse 1" = "G1","Growhouse 2" = "G2")) +
  scale_colour_manual("Mushroom",breaks = c("Mature", "Young"),
                      values = c("Mature"="cyan", "Young"="orange")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  theme(legend.key = element_blank()) +
  theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())

alpha_16s_out1
alpha_16s_out1$layers
alpha_16s_out1$layers <- alpha_16s_out1$layers[-1]

# rename facets on the plot 
levels(alpha_16s_out1$data$variable)
levels(alpha_16s_out1$data$variable)[levels
       (alpha_16s_out1$data$variable)=="Observed"] <- "Richness"


# alpha_16s_out2 ---------------
sample_data(biom_16s_uparse_out_ev)$Origin <- factor(sample_data(biom_16s_uparse_out_ev)$Origin,
                                                     levels=c("Cap","Stem","Soil"))

alpha_16s_out2 = plot_richness(biom_16s_uparse_out_ev, x="Growhouse", 
                               color="Origin", measures = c("Shannon", "Observed")) + 
  theme_pubclean() +
  expand_limits(x = 0, y = 0) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  labs(title="", x="Growhouse", y = "") +
  scale_x_discrete("Growhouse", labels = c("Growhouse 1" = "G1","Growhouse 2" = "G2")) +
  scale_colour_manual("Mushroom",breaks = c("Cap", "Stem", "Soil"),
                      values = c("Cap"="purple", "Stem"="grey", Soil="pink")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +  
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(legend.key = element_blank()) +
  theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())

alpha_16s_out2
alpha_16s_out2$layers
alpha_16s_out2$layers <- alpha_16s_out2$layers[-1]

# rename facets on the plot 
levels(alpha_16s_out2$data$variable)
levels(alpha_16s_out2$data$variable)[levels
      (alpha_16s_out2$data$variable)=="Observed"] <- "Richness"

# alpha_prokaryote_in --------------
levels(sample_data(biom_16s_uparse_in_ev)$Stage)[levels
       (sample_data(biom_16s_uparse_in_ev)$Stage)=="Late_fruiting"] <- "Late fruiting"

levels(sample_data(biom_16s_uparse_in_ev)$Stage)[levels
       (sample_data(biom_16s_uparse_in_ev)$Stage)=="Late_non_fruiting"] <- "Late non fruiting"

sample_data(biom_16s_uparse_in_ev)$Stage <- factor(sample_data(biom_16s_uparse_in_ev)$Stage,
           levels=c("Primordia","Fundamental","Differentiation","Late fruiting","Late non fruiting"))

alpha_prokaryote_in = plot_richness(biom_16s_uparse_in_ev, x="Stage", color ="Stage", 
                               measures = c("Shannon", "Observed")) +
  theme_pubclean() +
  geom_boxplot(outlier.colour="black", outlier.fill = "black", show.legend = FALSE) +
  expand_limits(x = 0, y = 0) +
  labs(title="", x="Stage", y = "") +
  scale_colour_manual(values=c("#0000FF","#0000FF","#0000FF","#0000FF","#0000FF")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank()) +
  theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="none")

alpha_prokaryote_in
alpha_prokaryote_in$layers
alpha_prokaryote_in$layers <- alpha_prokaryote_in$layers[-1]

# rename facets on the plot 
levels(alpha_prokaryote_in$data$variable)
levels(alpha_prokaryote_in$data$variable)[levels
      (alpha_prokaryote_in$data$variable)=="Observed"] <- "Richness"



alpha_16s_out1_A = ggarrange(alpha_16s_out1,
                              labels = "A",
                              widths = 1,
                              ncol = 1,
                              nrow = 1,
                              legend = "bottom")
alpha_16s_out1_A

alpha_16s_out2_B = ggarrange(alpha_16s_out2,
                              labels = "B",
                              widths = 1,
                              ncol = 1,
                              nrow = 1,
                              legend = "bottom")
alpha_16s_out2_B

alpha_prokaryote_in_C = ggarrange(alpha_prokaryote_in,
                             labels = "C",
                             widths = 1,
                             ncol = 1,
                             nrow = 1,
                             legend = "none")
alpha_prokaryote_in_C


# *** FIGURE 4 ---------------------------------
ggarrange(alpha_16s_out1_A,
          alpha_16s_out2_B,
          alpha_prokaryote_in_C,
          widths = c(0.9, 1.2, 0.89),
          align = "none" ,
          ncol = 3, nrow = 1, 
          legend = c("bottom"))



# create vegan objects ----------------------

otu_fungi_out <- as.data.frame(otu_table(biom_ITS_uparse_out_ev))
taxa_fungi_out <- as.data.frame(as.matrix(tax_table(biom_ITS_uparse_out_ev)))
metadata_fungi_out <- as.data.frame(as.matrix(sample_data(biom_ITS_uparse_out_ev)))
dim(otu_fungi_out)


otu_fungi_in <- as.data.frame(otu_table(biom_ITS_uparse_in_ev))
taxa_fungi_in <- as.data.frame(as.matrix(tax_table(biom_ITS_uparse_in_ev)))
metadata_fungi_in <- as.data.frame(as.matrix(sample_data(biom_ITS_uparse_in_ev)))
dim(otu_fungi_in)


otu_prokaryote_out <- as.data.frame(otu_table(biom_16s_uparse_out_ev))
taxa_prokaryote_out <- as.data.frame(as.matrix(tax_table(biom_16s_uparse_out_ev)))
metadata_prokaryote_out <- as.data.frame(as.matrix(sample_data(biom_16s_uparse_out_ev)))
dim(otu_prokaryote_out)


otu_prokaryote_in <- as.data.frame(otu_table(biom_16s_uparse_in_ev))
taxa_prokaryote_in <- as.data.frame(as.matrix(tax_table(biom_16s_uparse_in_ev)))
metadata_prokaryote_in <- as.data.frame(as.matrix(sample_data(biom_16s_uparse_in_ev)))
dim(otu_prokaryote_in)


# >>> BETA DIVERSITY ----------------------------------------

# MNDS multidimensional scaling -----------------------------

# NMDS prokaryotes ------
sample_data(biom_16s_uparse_in_ev)
nmds_16s_out = ordinate(biom_16s_uparse_out_ev, method ="NMDS", distance="bray", try=200)
nmds_16s_out

p_nmds_16s_out = plot_ordination(biom_16s_uparse_out_ev, nmds_16s_out, color="Growhouse", shape = "Stage") + 
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=1) + # ,aes(shape=Age))
  scale_colour_manual(values=c("cyan","orange")) +
  scale_shape_manual(values=c(16, 17)) +
  #theme(legend.title=element_blank()) +
  theme(legend.position="right")

p_nmds_16s_out

p_nmds_16s_out_Tissue = plot_ordination(biom_16s_uparse_out_ev, nmds_16s_out, color="Origin") + 
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=1) + # ,aes(shape=Age))
  scale_colour_manual(values = c("Cap"="purple", "Stem"="grey", Soil="pink")) +
  theme(legend.position="right") 

p_nmds_16s_out_Tissue


nmds_16s_in = ordinate(biom_16s_uparse_in_ev, method ="NMDS", distance="bray", try=200)
nmds_16s_in

p_nmds_16s_in = plot_ordination(biom_16s_uparse_in_ev, nmds_16s_in, color="Stage", shape = "Stage") + 
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, bg= "#0000FF") + # ,aes(shape=Age))
  scale_colour_manual(values=c("#0000FF","#0000FF","#0000FF","#0000FF","#0000FF")) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  theme(legend.position="right")

p_nmds_16s_in

# *** FIGURE 5 ----------------------------
ggarrange(p_nmds_16s_out,
          p_nmds_16s_out_Tissue,
          p_nmds_16s_in,
          labels = c("A", "B", "C"),
          widths = c(1, 0.95, 1.15),
          align = "none", 
          ncol = 3, nrow = 1, 
          legend = c("right"))

# NMDS fungi -----
sample_data(biom_ITS_uparse_out_ev)

nmds_ITS_out = ordinate(biom_ITS_uparse_out_ev, method ="NMDS", k=2, distance="bray", try=200)
nmds_ITS_out

p_nmds_ITS_out = plot_ordination(biom_ITS_uparse_out_ev, nmds_ITS_out, color="Growhouse", shape = "Stage") + 
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=1) + # ,aes(shape=Age))
  scale_colour_manual(values=c("red","green")) +
  scale_shape_manual(values=c(16, 17)) +
  #theme(legend.title=element_blank())
  theme(legend.position="right") 

p_nmds_ITS_out


sample_data(biom_ITS_uparse_in_ev)

nmds_ITS_in = ordinate(biom_ITS_uparse_in_ev, method ="NMDS", distance="bray", try=200)
nmds_ITS_in

p_nmds_ITS_in = plot_ordination(biom_ITS_uparse_in_ev, nmds_ITS_in, color="Stage", shape = "Stage") + 
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, bg="black") + # ,aes(shape=Age))
  scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000")) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  #theme(legend.title=element_blank()) + 
  theme(legend.position="right")

p_nmds_ITS_in

# *** FIGURE 6 ------------------------------------
ggarrange(p_nmds_ITS_out, 
          p_nmds_ITS_in,
          labels = c("A", "B"),
          widths = c(1,1.1),
          align = "none", 
          ncol = 2, nrow = 1, 
          legend = c("right"))




# Shepard diagrams -----------------

metaMDS_prokaryote_out = metaMDS(t(otu_prokaryote_out), k=2, trymax=200, distance="bray", weakties = TRUE)
metaMDS_prokaryote_out
stressplot(metaMDS_prokaryote_out, main="Shepard diagram\nNMDS Prokaryotes outdoor")

metaMDS_prokaryote_in = metaMDS(t(otu_prokaryote_in), k=2, trymax=200, distance="bray", weakties = TRUE)
metaMDS_prokaryote_in
stressplot(metaMDS_prokaryote_in, main="Shepard diagram\nNMDS Prokaryotes indoor")

metaMDS_fungi_out = metaMDS(t(otu_fungi_out), k=2, trymax=200, distance="bray", weakties = TRUE)
metaMDS_fungi_out
stressplot(metaMDS_fungi_out, main="Shepard diagram\NMDS Fungi outdoor")

metaMDS_fungi_in = metaMDS(t(otu_fungi_in), k=2, trymax=200, distance="bray", weakties = TRUE)
metaMDS_fungi_in
stressplot(metaMDS_fungi_in, main="Shepard diagram\nNMDS Fungi indoor")


# FIGURE S2 ----------------------------------------------------------
par(mfrow=c(2,2)) 
stressplot(metaMDS_prokaryote_out, main="A", adj="0")
stressplot(metaMDS_prokaryote_in, main="B", adj="0")
stressplot(metaMDS_fungi_out, main="C",adj="0")
stressplot(metaMDS_fungi_in, main="D",adj="0")
par(mfrow=c(1,1)) 


# >>> PERMANOVA -----------------------------------------

library("vegan")
library("RVAideMemoire")
metadata_fungi_out

model.matrix(~ Stage * Growhouse, data=metadata_fungi_out)
model.matrix(~ Stage + Growhouse + Stage : Growhouse, data=metadata_fungi_out)

adonis(t(otu_fungi_out) ~ Stage * Growhouse, data=metadata_fungi_out, permutations=9999) # by = "margin"

adonis(t(otu_fungi_out) ~ Stage + Growhouse + Stage : Growhouse, data=metadata_fungi_out, permutations=9999) 
adonis(t(otu_fungi_out) ~ Growhouse + Stage + Stage : Growhouse, data=metadata_fungi_out, permutations=9999)

#adonis.II(t(otu_fungi_out) ~ Stage + Growhouse + Stage : Growhouse, data=metadata_fungi_out, permutations=9999)

vegan::vegdist(t(otu_fungi_out), method="bray") -> dist_otu_fungi_out

permdisp_otu_fungi_out_G <- betadisper(dist_otu_fungi_out, metadata_fungi_out$Growhouse)
permdisp_otu_fungi_out_S <- betadisper(dist_otu_fungi_out, metadata_fungi_out$Stage)

anova(permdisp_otu_fungi_out_G, permutations = 9999)
permutest(permdisp_otu_fungi_out_G, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_out_G)
plot(TukeyHSD(permdisp_otu_fungi_out_G), las=1)
boxplot(permdisp_otu_fungi_out_G)

anova(permdisp_otu_fungi_out_S, permutations = 9999)
permutest(permdisp_otu_fungi_out_S, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_out_S)
plot(TukeyHSD(permdisp_otu_fungi_out_S), las=1)
boxplot(permdisp_otu_fungi_out_S)


metadata_fungi_in
adonis(t(otu_fungi_in) ~ Stage, data=metadata_fungi_in, permutations=9999)

vegan::vegdist(t(otu_fungi_in), method="bray") -> dist_otu_fungi_in
permdisp_otu_fungi_in_S <- betadisper(dist_otu_fungi_in, metadata_fungi_in$Stage)
anova(permdisp_otu_fungi_in_S)
plot(permdisp_otu_fungi_in_S)
boxplot(permdisp_otu_fungi_in_S)


metadata_prokaryote_out

model.matrix(~ Origin * Stage * Growhouse, data=metadata_prokaryote_out)

adonis(t(otu_prokaryote_out) ~ Origin * Stage * Growhouse, data=metadata_prokaryote_out, permutations=9999)
adonis(t(otu_prokaryote_out) ~ Origin + Stage + Growhouse + Origin:Stage + Origin:Growhouse, data=metadata_prokaryote_out, permutations=9999)


vegan::vegdist(t(otu_prokaryote_out), method="bray") -> dist_otu_prokaryote_out
permdisp_otu_prokaryote_out_G <- betadisper(dist_otu_prokaryote_out, metadata_prokaryote_out$Growhouse)
permdisp_otu_prokaryote_out_S <- betadisper(dist_otu_prokaryote_out, metadata_prokaryote_out$Stage)
permdisp_otu_prokaryote_out_T <- betadisper(dist_otu_prokaryote_out, metadata_prokaryote_out$Origin)

anova(permdisp_otu_prokaryote_out_G, permutations = 9999)
permutest(permdisp_otu_prokaryote_out_G, permutations = 9999, pairwise = T)
plot(permdisp_otu_prokaryote_out_G)

anova(permdisp_otu_prokaryote_out_S, permutations = 9999)
permutest(permdisp_otu_prokaryote_out_S, permutations = 9999, pairwise = T)
plot(permdisp_otu_prokaryote_out_S)

anova(permdisp_otu_prokaryote_out_T, permutations = 9999)
permutest(permdisp_otu_prokaryote_out_T, permutations = 9999, pairwise = T)
plot(permdisp_otu_prokaryote_out_T)


metadata_prokaryote_in
adonis(t(otu_prokaryote_in) ~ Stage, data=metadata_prokaryote_in, permutations=9999)

vegan::vegdist(t(otu_prokaryote_in), method="bray") -> dist_otu_prokaryote_in
permdisp_otu_prokaryote_in_S <- betadisper(dist_otu_prokaryote_in, metadata_prokaryote_in$Stage)
anova(permdisp_otu_prokaryote_in_S, permutations=9999)
plot(permdisp_otu_prokaryote_in_S)


# >>> INDICATOR SPECIES ANALYSIS (ISA) -----------------------------

#indicator species analysis (isa) Fungi
library("indicspecies")

# Fungi -----------

otu_fungi_out
otu_fungi_in

isa_ITS_out_stage <- multipatt(as.data.frame(t(otu_fungi_out)),
                              metadata_fungi_out$Stage, control=how(nperm=9999))
summary(isa_ITS_out_stage, indvalcomp=TRUE)

isa_ITS_out_stage -> isa_ITS_out_stage_fdr
isa_ITS_out_stage_fdr$sign$p.value<-p.adjust(isa_ITS_out_stage$sign$p.value, "fdr")
summary(isa_ITS_out_stage_fdr)


isa_ITS_out_grow <- multipatt(as.data.frame(t(otu_fungi_out)),
                               metadata_fungi_out$Growhouse, control=how(nperm=9999))
summary(isa_ITS_out_grow, indvalcomp=TRUE)

isa_ITS_out_grow -> isa_ITS_out_grow_fdr
isa_ITS_out_grow_fdr$sign$p.value<-p.adjust(isa_ITS_out_grow$sign$p.value, "fdr")
summary(isa_ITS_out_grow_fdr)

  
isa_ITS_in_stage <- multipatt(as.data.frame(t(otu_fungi_in)),
                              metadata_fungi_in$Stage, control=how(nperm=9999))
summary(isa_ITS_in_stage, indvalcomp=TRUE)

isa_ITS_in_stage -> isa_ITS_in_stage_fdr
isa_ITS_in_stage_fdr$sign$p.value<-p.adjust(isa_ITS_in_stage$sign$p.value, "fdr")
summary(isa_ITS_in_stage_fdr)


# prokaryote  ----------------

otu_prokaryote_out
otu_prokaryote_in

isa_16s_out_stage <- multipatt(as.data.frame(t(otu_prokaryote_out)),
                               metadata_prokaryote_out$Stage, control=how(nperm=9999))
summary(isa_16s_out_stage, indvalcomp=TRUE)

isa_16s_out_stage -> isa_16s_out_stage_fdr
isa_16s_out_stage_fdr$sign$p.value<-p.adjust(isa_16s_out_stage$sign$p.value, "fdr")
summary(isa_16s_out_stage_fdr)


isa_16s_out_grow <- multipatt(as.data.frame(t(otu_prokaryote_out)),
                              metadata_prokaryote_out$Growhouse, control=how(nperm=9999))
summary(isa_16s_out_grow, indvalcomp=TRUE)

isa_16s_out_grow -> isa_16s_out_grow_fdr
isa_16s_out_grow_fdr$sign$p.value<-p.adjust(isa_16s_out_grow$sign$p.value, "fdr")
summary(isa_16s_out_grow_fdr)

sink(file="isa_16s_out_growhouse_fdr.csv") # A trick to save a super long output!!
summary(isa_16s_out_grow_fdr)
sink()


isa_16s_out_tissue <- multipatt(as.data.frame(t(otu_prokaryote_out)),
                             metadata_prokaryote_out$Origin, control=how(nperm=9999))
summary(isa_16s_out_tissue, indvalcomp=TRUE)

isa_16s_out_tissue -> isa_16s_out_tissue_fdr
isa_16s_out_tissue_fdr$sign$p.value<-p.adjust(isa_16s_out_tissue$sign$p.value, "fdr")
summary(isa_16s_out_tissue_fdr)

sink(file="isa_16s_out_tissue_fdr.csv") 
summary(isa_16s_out_tissue_fdr)
sink()

isa_16s_in_stage <- multipatt(as.data.frame(t(otu_prokaryote_in)),
                              metadata_prokaryote_in$Stage, control=how(nperm=9999))
summary(isa_16s_in_stage, indvalcomp=TRUE)

isa_16s_in_stage -> isa_16s_in_stage_fdr
isa_16s_in_stage_fdr$sign$p.value<-p.adjust(isa_16s_in_stage$sign$p.value, "fdr")
summary(isa_16s_in_stage_fdr)


# extract data from ISA result ---------------------------
library(dplyr)

summary(isa_16s_out_grow_fdr)
head(isa_16s_out_grow_fdr$sign)
str(isa_16s_out_grow_fdr$sign)

result_isa_16s_out_grow_fdr <- subset(isa_16s_out_grow_fdr$sign, p.value <= 0.05, select = c(s.G1, s.G2, p.value))
result_isa_16s_out_grow_fdr <- subset(result_isa_16s_out_grow_fdr, rowSums(result_isa_16s_out_grow_fdr) <= 2, select = c(s.G1, s.G2, p.value))
result_isa_16s_out_grow_fdr
dim(result_isa_16s_out_grow_fdr)

# subset your table accordingly, here a phyloseq object, as an example
biom_16s_uparse_out_ISA_Gh <- biom_16s_uparse_out
otu_table(biom_16s_uparse_out_ISA_Gh) <-otu_table(biom_16s_uparse_out)[row.names(result_isa_16s_out_grow_fdr), ]

#any(colSums(otu_table(biom_16s_uparse_out_ISA_Gh)) == 0)
#otu_table(biom_16s_uparse_out_ISA_Gh) <- otu_table(biom_16s_uparse_out_ISA_Gh)[,which(colSums(otu_table(biom_16s_uparse_out_ISA_Gh)) > 0)] 

biom_16s_uparse_out_ISA_Gh
sample_data(biom_16s_uparse_out_ISA_Gh)
tax_table(biom_16s_uparse_out_ISA_Gh)
otu_table(biom_16s_uparse_out_ISA_Gh)


write.csv(otu_table(biom_16s_uparse_out_ISA_Gh), "otus_16s_uparse_out_ISA_Gh.csv")
write.csv(tax_table(biom_16s_uparse_out_ISA_Gh), "taxa_16s_uparse_out_ISA_Gh.csv")


# plot bars for ISA -------------
source("Rscript_2_plot_ordered_bar.R")

biom_16s_uparse_out_ISA_Gh_class = tax_glom(biom_16s_uparse_out_ISA_Gh, "Class")

biom_16s_uparse_out_ISA_Gh_merged = merge_samples(biom_16s_uparse_out_ISA_Gh_class, "Growhouse")
sample_data(biom_16s_uparse_out_ISA_Gh_merged)$Growhouse <- levels(sample_data(biom_16s_uparse_out_ISA_Gh)$Growhouse)

biom_16s_uparse_out_ISA_Gh_merged_ab = transform_sample_counts(biom_16s_uparse_out_ISA_Gh_merged, function(x) 100 * x/sum(x))

top_30_ITS <- names(sort(taxa_sums(biom_16s_uparse_out_ISA_Gh_merged_ab), TRUE)[1:30]) 
biom_16s_uparse_out_ISA_Gh_merged_ab_30 <- prune_taxa(top_30_ITS, biom_16s_uparse_out_ISA_Gh_merged_ab)
biom_16s_uparse_out_ISA_Gh_merged_ab_30

# setting up color palettes -----------------------------------------------------------

palette_new30 = c("#560d0d","#dba4a4", "#cc1c1c","#a0fffc","#111b77","#283dff","#636bb7","#bfc5ff",
                  "#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121","#ea7f17","#fcb067",
                  "#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000","#a35151", "#5b5b19","#fcfc00",
                  "#ffff9e","#ffb7ef","#fa7efc","#ae09ea","#521899","#1e0047")

pie(rep(1, length(palette_new30)), labels = sprintf("%d (%s)", 
     seq_along(palette_new30),palette_new30), col = palette_new30)


# plot ISA for Growhouse
isa_plot_Gh = plot_ordered_bar(biom_16s_uparse_out_ISA_Gh_merged_ab_30, x = "Growhouse", fill="Class", leg_size = 0.4) +
  theme_classic() +
  labs(x="", y="Relative abundance", title="") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5, colour = "black")) +
  scale_fill_manual(values=palette_new30) +
  theme(legend.position="right")

isa_plot_Gh


summary(isa_16s_out_tissue_fdr)
head(isa_16s_out_tissue_fdr$sign)
str(isa_16s_out_tissue_fdr$sign)

result_isa_16s_out_tissue_fdr <- subset(isa_16s_out_tissue_fdr$sign, p.value <= 0.05, select = c(s.Cap, s.Soil, s.Stem, p.value))
result_isa_16s_out_tissue_fdr <- subset(result_isa_16s_out_tissue_fdr, rowSums(result_isa_16s_out_tissue_fdr) <= 3, select = c(s.Cap, s.Soil, s.Stem, p.value))
result_isa_16s_out_tissue_fdr
dim(result_isa_16s_out_tissue_fdr)

# subset your table accordingly, here a phyloseq object, as an example
biom_16s_uparse_out_ISA_Ts <- biom_16s_uparse_out
otu_table(biom_16s_uparse_out_ISA_Ts) <-otu_table(biom_16s_uparse_out)[row.names(result_isa_16s_out_tissue_fdr), ]

#any(colSums(otu_table(biom_16s_uparse_out_ISA_Ts)) == 0)
#otu_table(biom_16s_uparse_out_ISA_Ts) <- otu_table(biom_16s_uparse_out_ISA_Ts)[,which(colSums(otu_table(biom_16s_uparse_out_ISA_Ts)) > 0)] 

biom_16s_uparse_out_ISA_Ts
sample_data(biom_16s_uparse_out_ISA_Ts)
tax_table(biom_16s_uparse_out_ISA_Ts)
otu_table(biom_16s_uparse_out_ISA_Ts)

write.csv(otu_table(biom_16s_uparse_out_ISA_Ts), "otus_16s_uparse_out_ISA_Ts.csv")
write.csv(tax_table(biom_16s_uparse_out_ISA_Ts), "taxa_16s_uparse_out_ISA_Ts.csv")


# plot bars for ISA -------------

biom_16s_uparse_out_ISA_Ts_class = tax_glom(biom_16s_uparse_out_ISA_Ts, "Class")

biom_16s_uparse_out_ISA_Ts_merged = merge_samples(biom_16s_uparse_out_ISA_Ts_class, "Origin")
sample_data(biom_16s_uparse_out_ISA_Ts_merged)$Origin <- levels(sample_data(biom_16s_uparse_out_ISA_Ts)$Origin)

biom_16s_uparse_out_ISA_Ts_merged_ab = transform_sample_counts(biom_16s_uparse_out_ISA_Ts_merged, function(x) 100 * x/sum(x))

top_30_ITS <- names(sort(taxa_sums(biom_16s_uparse_out_ISA_Ts_merged_ab), TRUE)[1:30]) 
biom_16s_uparse_out_ISA_Ts_merged_ab_30 <- prune_taxa(top_30_ITS, biom_16s_uparse_out_ISA_Ts_merged_ab)
biom_16s_uparse_out_ISA_Ts_merged_ab_30

sample_data(biom_16s_uparse_out_ISA_Ts_merged_ab_30)$Origin <- factor(sample_data(biom_16s_uparse_out_ISA_Ts_merged_ab_30)$Origin, levels=c("Cap","Stem","Soil"))

# Plot ISa for Sample Origin == Tissue
isa_plot_Ts = plot_ordered_bar(biom_16s_uparse_out_ISA_Ts_merged_ab_30, x = "Origin", fill="Class", leg_size = 0.4) +
  theme_classic() +
  labs(x="", y="Relative Abundance", title="") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5, colour = "black")) +
  scale_fill_manual(values=palette_new30) +
  theme(legend.position="right")

isa_plot_Ts

# *** FIGURE S3 ---------------------------------------
ggarrange(isa_plot_Gh, isa_plot_Ts,
          labels = c("A", "B"),
          widths = c(0.88, 1),
          align = "h", ncol = 2, nrow = 1)



# >>> VENN DIAGRAM -------------------------------------
library("limma")

source("Rscript_1_venn_diagram.R")

otu_fungi_out
metadata_fungi_out

biom_ITS_uparse_out_ev_Gh = merge_samples(biom_ITS_uparse_out_ev, "Growhouse")
otu_fungi_out_Gh <- as.data.frame(t(otu_table(biom_ITS_uparse_out_ev_Gh)))
venn_counts_otu_fungi_out_Gh <- vennCounts(otu_fungi_out_Gh, include="both")
venn_counts_otu_fungi_out_Gh

biom_ITS_uparse_out_ev_St = merge_samples(biom_ITS_uparse_out_ev, "Stage")
otu_fungi_out_St <- as.data.frame(t(otu_table(biom_ITS_uparse_out_ev_St)))
venn_counts_otu_fungi_out_St <- vennCounts(otu_fungi_out_St, include="both")
venn_counts_otu_fungi_out_St

biom_ITS_uparse_in_ev_St = merge_samples(biom_ITS_uparse_in_ev, "Stage")
otu_fungi_in_St <- as.data.frame(t(otu_table(biom_ITS_uparse_in_ev_St)))
venn_counts_otu_fungi_in_St <- vennCounts(otu_fungi_in_St, include="both")
venn_counts_otu_fungi_in_St


otu_prokaryote_out
metadata_prokaryote_out

biom_16s_uparse_out_ev_St = merge_samples(biom_16s_uparse_out_ev, "Stage")
otu_prokaryote_out_St <- as.data.frame(t(otu_table(biom_16s_uparse_out_ev_St)))
venn_counts_otu_prokaryote_out_St <- vennCounts(otu_prokaryote_out_St, include="both")
venn_counts_otu_prokaryote_out_St

biom_16s_uparse_out_ev_Ts = merge_samples(biom_16s_uparse_out_ev, "Origin")
otu_prokaryote_out_Ts <- as.data.frame(t(otu_table(biom_16s_uparse_out_ev_Ts)))
venn_counts_otu_prokaryote_out_Ts <- vennCounts(otu_prokaryote_out_Ts, include="both")
venn_counts_otu_prokaryote_out_Ts

biom_16s_uparse_in_ev_St = merge_samples(biom_16s_uparse_in_ev, "Stage")
otu_prokaryote_in_St <- as.data.frame(t(otu_table(biom_16s_uparse_in_ev_St)))
venn_counts_otu_prokaryote_in_St <- vennCounts(otu_prokaryote_in_St, include="both")
venn_counts_otu_prokaryote_in_St


# *** FIGURE 7 --------------------------------

layout(matrix(1:6, ncol=3))
my.venn_Gian(venn_counts_otu_fungi_out_Gh,
            cex=c(0.9),
            circle.col =c("green", "red"),
            mar = c(1,1,1,1),
            lwd = 2, main="A")

my.venn_Gian(venn_counts_otu_prokaryote_out_Ts,
            cex=c(0.9),
            circle.col =c("purple", "grey", "pink"),
            mar = c(1,1,1,1),
            lwd = 2, main="B")


my.venn_Gian(venn_counts_otu_fungi_out_St,
            cex=c(0.9),
            circle.col =c("green", "red"),
            mar = c(1,1,1,1),
            lwd = 2, main="C")

my.venn_Gian(venn_counts_otu_prokaryote_out_St,
            cex=c(0.9),
            circle.col =c("cyan", "orange"),
            mar = c(1,1,1,1),
            lwd = 2, main="D")

my.venn_Gian(venn_counts_otu_fungi_in_St,
            cex=c(0.9), 
            circle.col = "black",
            mar = c(1,1,1,1),
            lwd = 2, main="E")

my.venn_Gian(venn_counts_otu_prokaryote_in_St,
            cex=c(0.9), 
            circle.col = "blue",
            mar = c(1,1,1,1),
            lwd = 2, main="F")

dev.off()








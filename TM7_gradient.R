#clean environment
rm(list=ls())

#import libraries
library(vegan)
# library(ggplot2)
library(EcolUtils)
library(biomformat)
# library(stringr)
library(dplyr)
# library(gridExtra)
# library(tidyverse)
library(ape)

#import data
setwd("C:/Users/Zach/Documents/Microbiome Family Project")

taxonomy <- read.delim("97_otu_taxonomy.txt", header = FALSE, row.names = 1)
taxonomy <- separate(taxonomy, col = V2, sep = ';',
                     into = c('Kingdom', 'Phyla', 'Class',
                              'Order', 'Famiy', 'Genus', 'Species'))
Microbiome_mapping <- read.delim("Microbiome_mapping_14Oct2019.txt")
rownames(Microbiome_mapping) <- Microbiome_mapping[,1]
Microbiome_mapping[,1] <- NULL

my_biom <- read_hdf5_biom("61277_otu_table.biom") #Change here within quotes.
write_biom(my_biom, "formatted_biom.biom")
my_biom <- read_biom("formatted_biom.biom")
OTU_table <- as.data.frame(as.matrix(biom_data(my_biom)))
colnames(OTU_table) <- sub('11959.', '', colnames(OTU_table))

merged_OTU_table <- merge(taxonomy,
                          OTU_table,
                          by.x = "row.names",
                          by.y = "row.names")
rownames(merged_OTU_table) <- merged_OTU_table[,1]
merged_OTU_table[,1] <- NULL


#rarefaction
transOTU <- rowSums(t(OTU_table)) 
Q01 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.0106) #cuts of at ~1350 sample size
rarecurve(t(OTU_table), step = 100, cex = 0.5)
rared_OTU <- as.data.frame((rrarefy.perm(t(OTU_table), sample = Q01, n = 100, round.out = T)))
rared_OTU <- as.data.frame(rared_OTU[rowSums(rared_OTU) >= Q01 - (Q01 * 0.1), colSums(rared_OTU) >= 1])


############### ############### ############### 
###############  Beta diversity ###############
############### ############### ###############


TM7_merged_OTU <- merged_OTU_table[merged_OTU_table$Phyla == 'p__TM7']


TM7_merged_OTU <- merged_OTU_table[grepl('p__TM7',merged_OTU_table$Phyla),]

TM7_scores <- as.data.frame(colSums(TM7_merged_OTU[c(8:103)]))
colnames(TM7_scores) <- c('score')
# chi <- cca(rared_OTU)
# d.chisq <- dist(chi$CA$u[,1:2])
NMDS1 <- metaMDS(rared_OTU, distance = "bray", k = 2)
coordinates <- data.frame(NMDS1$points[,1:2])

plot(x = coordinates$MDS1, y = coordinates$MDS2)

nmds_plus_mapping <- merge(coordinates, Microbiome_mapping, by = 0)
rownames(nmds_plus_mapping) <- nmds_plus_mapping[,1]
nmds_plus_mapping[,1] <- NULL

nmds_plus_mapping_and_scores <- merge(nmds_plus_mapping, TM7_scores, by = 0)
rownames(nmds_plus_mapping_and_scores) <- nmds_plus_mapping_and_scores[,1]
nmds_plus_mapping_and_scores[,1] <- NULL

ggplot(data = nmds_plus_mapping_and_scores) +
  aes(x = MDS1, y = MDS2, shape = Family, color = score)+
  geom_point(size = 3) +
  scale_color_gradient(low = 'blue', high = 'red')
  labs(col = "Family") +  
  theme_bw()

filtersamples <- rownames(rared_OTU)
filtermapping <- subset(Microbiome_mapping, rownames(Microbiome_mapping) %in% filtersamples)

adonis(data = filtermapping, formula = rared_OTU ~ Family 
       / Individual + Day + Day:Family,
       permutations = 999, method = "bray")

adonis(data = filtermapping, formula = rared_OTU ~ Family
       /Individual/Day,
       permutations = 999, method = "bray")

ggplot(data = nmds_plus_Microbiome_mapping) +
  aes(x = MDS1, y = MDS2, color = Family) +
  geom_point() +
  labs(fill = "Family") +
  ggtitle("NMDS of Families on Oral Microbiome") +
  theme_classic(base_size = 14, base_line_size = .5)

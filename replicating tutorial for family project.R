#libraries necessary
library(vegan)
library(ggplot2)
library(EcolUtils)
library(biomformat)
library(stringr)
library(dplyr)
library(gridExtra)
library(tidyverse)


#import data

setwd("C:/Users/Zach/Documents/Microbiome Family Project")

taxonomy <-read.delim("taxonomy.txt")

Microbiome_mapping <- read.delim("Microbiome_mapping.txt")

#OTU-table (needs to have columns in number order)
OTU_table <- read.delim("OTU_table.xlsx")

#did not filter data because it was said to already be clean

#rarefaction

# get quartile ranges for rarefaction
#OTU-table was put in order by columns in excel, then imported to R
transOTU <- rowSums(t(OTU_table[,2:96])) 

Q10 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.10)
Q15 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.15)

# plot 
barplot(sort(transOTU), ylim = c(0, max(transOTU)), 
        xlim = c(0, NROW(transOTU)), col = "Blue", ylab = "Read Depth", xlab = "Sample") 
abline(h = c(Q10, Q15), col = c("red", "pink"))
plot.new()


rarecurve(t(OTU_table[,2:96]), step = 100, cex = 0.5)
abline(v = c(Q10, Q15), col = c("red", "pink"))

# you will need to make a BIG decision here on where to draw the rarefaction cutoff
# samples to the left of the red or pink lines will be thrown out
rared_OTU <- as.data.frame((rrarefy.perm(t(OTU_table[,2:96]), sample = Q10, n = 100, round.out = T)))

# This only keeps the samples that meet the rarefaction cutoff.
rared_OTU <- as.data.frame(rared_OTU[rowSums(rared_OTU) >= Q10 - (Q10 * 0.1), colSums(rared_OTU) >= 1])


############### ############### ############### 
############### Alpha diversity ###############
############### ############### ############### 

# outside of this, we can look at the richness of each sample, or how many taxa are in each sample
richness <- as.data.frame(specnumber(rared_OTU))
colnames(richness) <- c("speciesrich")

# Merge with metadata to create a plot.
nerged_rich <- merge(richness, Microbiome_mapping, by = 0)
rownames(nerged_rich) <- nerged_rich$Row.names
nerged_rich$Row.names <- NULL

# now for alpha diversity
# we will use Shannon diversity as a method of alpha-diversity.
shannon <- as.data.frame(diversity(rared_OTU, index = "shannon"))
colnames(shannon) <- c("alpha_shannon")

# Merge with metadata to create a plot.
nerged_alpha <- merge(nerged_rich, shannon, by = 0)

#filter data frame to only have individuals being compared.
#In this case, fam1 members 
nerged_alpha <- filter(nerged_alpha, Family == 'Fam1')

#Plotting alpha diversity:
#Change 'factor_X' to categorical metadata factor you wish to plot. Press tab after the '$' sign.
factor_X <- as.factor(nerged_alpha$Individual)

#plot richness
p1 <- ggplot(data = nerged_alpha) +
  aes(x = factor_X, y = nerged_alpha$speciesrich, 
      fill = factor_X) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  labs(title = 'Species Richness',
       #Change x-axis label and legend title to metadata factor of interest.
       x = 'Family 1 Individuals', y = 'Species Richness', fill = 'Individual') +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_jitter(width = .2) +
  theme(legend.position = "none")

# preview the figure
p1

# plot alpha diversity
p2 <- ggplot(data = nerged_alpha) +
  aes(x = factor_X, y = nerged_alpha$alpha_shannon, 
      fill = factor_X) +
  geom_boxplot(outlier.shape = NA, lwd = 1) +
  labs(title = 'Alpha Diversity',
       #Change x-axis label and legend title to metadata factor of interest.
       x = 'Family 1 Individuals', y = 'Shannon Diversity Index', fill = 'Individual') +
  theme_classic(base_size = 14, base_line_size = 1) +
  geom_jitter(width = .2) +
  theme(legend.position = "none")

#preview the figure
p2



############### ############### ############### 
###############  Beta diversity ###############
############### ############### ############### 
NMDS1 <- metaMDS(rared_OTU, distance = "bray", k = 2, trymax = 500)

# Extract the two axes of the NMDS to plot the x and y coordinates
coordinates <- data.frame(NMDS1$points[,1:2])

# Quick glance at the NMDS plot
# not too exciting right? well we need to add Microbiome_mapping!!
plot(x = coordinates$MDS1, y = coordinates$MDS2)

# so merge NMDS axes coordinates with Microbiome_mapping
nmds_plus_Microbiome_mapping <- merge(coordinates, Microbiome_mapping, by = 0)

# choose the factor you are interested in (this time let's check the "Family" effects)
Factor_x <- as.factor(nmds_plus_Microbiome_mapping$Family)

# Plot
ggplot(data = nmds_plus_Microbiome_mapping) +
  aes(x = MDS1, y = MDS2, color = Family)+ #Creates and colors legend to match, modify after $ here.
  geom_point(size = 3) +
  labs(col = "Family") + #Renames legend, modifiable. 
  theme_bw()

# for statistical tests, we need to get the data tidied up a bit
# so let's subset the Microbiome_mapping and keeps only the samples that passed filtering and rarefaction
filtersamples <- rownames(rared_OTU)
filtermapping <- subset(Microbiome_mapping, rownames(Microbiome_mapping) %in% filtersamples)

#factors
adonis(data = filtermapping, formula = rared_OTU ~ Family 
       / Individual + Day + Day:Family,
       permutations = 999, method = "bray")

adonis(data = filtermapping, formula = rared_OTU ~ Family
       /Individual/Day,
       permutations = 999, method = "bray")

# again, the adonis test is dependent on what you think is happening in the microbiomes so adjsut accordingly

# if we want we can go back to our plot and add these results in
# To replot: 
# Change the R2 and p-value numbers to match!
ggplot(data = nmds_plus_Microbiome_mapping) +
  aes(x = MDS1, y = MDS2, color = Family) +
  geom_point() +
  labs(fill = "Family") +
  ggtitle("NMDS of Families on Oral Microbiome") +
  theme_classic(base_size = 14, base_line_size = .5)

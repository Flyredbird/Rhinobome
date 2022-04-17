## Install pacakages
install.packages("devtools")
install.packages("rlang")
devtools::install_github("leffj/mctoolsr")
install.packages("multcompView")
install.packages("PMCMRplus")
install.packages("stats")
install.packages("ape")
install.packages("gplots")
install.packages("plotly")
install.packages("tidyr")
install.packages("graphics")
install.packages("phyloseq")
install.packages("vegan")
install.packages("RVAideMemoire")

## Load packages
library(mctoolsr)
library(stats)
library(vegan)
library(ggplot2)
library(multcompView)
library(scico)
library(PMCMRplus)
library(stats)
library(ape)
library(gplots)
library(plotly)
library(tidyr)
library(VennDiagram)
library(graphics)
library(phyloseq)
library(RVAideMemoire)

# Load Data Files
setwd("~/Desktop/bioinformatic_analysis/FullPopulationResults") # set file directory location
input_RhinoAdult <- load_taxa_table("AllDataOTUtable_semicolons.txt","AllRhinoDataMetadataCleaned.txt") #add OTU table and metadata files
#151 samples loaded

#Filter unwanted samples/categories
input_RhinoAdult <- filter_data(input = input_RhinoAdult, filter_cat= "FSJ", filter_vals="N") 
#41 samples remain

#Look at file structure
str(input_RhinoAdult)

# Sort OTUs by DNA reads, select rarefication value
sort <- sort(colSums(input_RhinoAdult$data_loaded))
print(sort)

# Rarefy dataset
input_RhinoAdult = single_rarefy(input_RhinoAdult,29265) 
#40 samples remain, lost 61FFeb

#Check to make sure data is rarefied (i.e. all samples have same number of reads)
sort(colSums(input_RhinoAdult$data_loaded))

### Alpha Diversity ##############################################

# Calculation of Richness

div.rich<-calc_diversity(input_RhinoAdult$data_loaded, metric = 'richness')

# Calculation of Shannon Diversity

div.shan<-calc_diversity(input_RhinoAdult$data_loaded, metric = 'shannon')

## Merge diversity values into original data set

div.rich.map <- merge(div.rich, input_RhinoAdult$map_loaded, by=0, all=TRUE)
div.shan.map <- merge(div.shan, input_RhinoAdult$map_loaded, by=0, all=TRUE)

### rename column "x" with richness values
names(div.rich.map)[2]="richness"

## rename column "x" with shannon values
names(div.shan.map)[2]="shannon"

#Alpha diversity boxplots
boxplot(richness ~ Age, data = div.rich.map, main = "Richness by Age")
boxplot(shannon ~ Age, data = div.shan.map, main = "Shannon Diversity by Age")


##Alpha diversity Kruskal-Wallis##

#Kruskal-Wallis: richness by Age
k.rich <- kruskal.test(`richness`~ Age, data=div.rich.map)
print(k.rich)

#Shannon Div
k.shan <- kruskal.test(`shannon`~ Age, data=div.shan.map)
print(k.shan)

#Pairwise Wilcoxon test for Richness across Age
rich.wcx=pairwise.wilcox.test(div.rich.map$richness, div.rich.map$Age, p.adj = "bonf")
rich.wcx

#Shannon ~ Age
shan.wcx=pairwise.wilcox.test(div.shan.map$shannon, div.shan.map$Age, p.adj = "bonf")
shan.wcx

##You are now entering the Rhinobome##

## Install mctoolsr from GitHub
install.packages("devtools")
install.packages("rlang")
devtools::install_github("leffj/mctoolsr")
#devtools::install_github("r-lib/rlang")
install.packages("multcompView")
install.packages('scico')
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

## Post-intall, load packages
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
library(utils)

#input data
input_RhinoAdult <- load_taxa_table("AllDataOTUtable_semicolons.txt","AllRhinoDataMetadataCleaned.txt") #add OTU table and metadata files
#151 samples loaded

#Filtering for Individual/Full Pop Comparisons#
input_RhinoAdult <- filter_data(input = input_RhinoAdult, filter_cat= "FSJ", filter_vals="N") 
#X samples loaded

#Filtering for Age Comparisons: Only F, S, and J samples from Winter#
#input_RhinoAdult <- filter_data(input = input_RhinoAdult, filter_cat= "FSJNoSummer", filter_vals="N") 
#X samples loaded

#Filtering for Month Comparisons: Only Female and Subadult Female Samples#
#input_RhinoAdult <- filter_data(input = input_RhinoAdult, filter_cat= "FSFrozen", filter_vals="N") 
#51 samples remain

input_RhinoAdult = single_rarefy(input_RhinoAdult,36100) 

#Most abundant phyla 
##2 = Phylum; 3 = Class; 4 = Order; 5 = Family; 6 = Genus
tax_sum_phyla <- summarize_taxonomy(input_RhinoAdult, level = 2, report_higher_tax = FALSE)
#most abundant by metadata characteristic e.g. 'Age'
sum_phy <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_phyla, 'Age', metadata_map = input_RhinoAdult$map_loaded, smry_fun = mean)

taxa.phyla.bars <- plot_taxa_bars(tax_sum_phyla, input_RhinoAdult$map_loaded, type_header = 'Age', num_taxa = 30) + scale_fill_scico_d(
  alpha = 1,
  begin = 0,
  end = 1,
  direction = -1,
  palette = "batlow",
  aesthetics = "fill")

plot(taxa.phyla.bars)

# Calculation of Richness

div.rich<-calc_diversity(input_RhinoAdult$data_loaded, metric = 'richness')

# Calculation of Shannon Diversity

div.shan<-calc_diversity(input_RhinoAdult$data_loaded, metric = 'shannon')

div.rich.map <- merge(div.rich, input_RhinoAdult$map_loaded, by=0, all=TRUE)
div.shan.map <- merge(div.shan, input_RhinoAdult$map_loaded, by=0, all=TRUE)
#div.simp.map <- merge(div.simp, input_RhinoAdult$map_loaded, by=0, all=TRUE)

### rename column "x" with richness values
names(div.rich.map)[2]="richness"
names(div.shan.map)[2]="shannon"

boxplot(richness ~ Month, data = div.rich.map, main = "Richness by Month")
boxplot(shannon ~ Month, data = div.shan.map, main = "Shannon Diversity by Month")

k.rich <- kruskal.test(`richness`~ Individual, data=div.rich.map)
print(k.rich)

#Shannon Div
k.shan <- kruskal.test(`shannon`~ Month, data=div.shan.map)
print(k.shan)

#rich.wcx=pairwise.wilcox.test(div.rich.map$richness, div.rich.map$Month, p.adj = "bonf")
#rich.wcx

#shan.wcx=pairwise.wilcox.test(div.shan.map$shannon, div.shan.map$Month, p.adj = "bonf")
#shan.wcx


dm <- calc_dm(input_RhinoAdult$data_loaded)

ord <- calc_ordination(dm, 'nmds')

plot_ordination(input_RhinoAdult, ord, 'Age', hulls=TRUE)
plot_ordination(input_RhinoAdult, ord, 'Month', hulls=TRUE) + scale_color_manual(values=c("midnightblue","darkcyan","darkolivegreen4","plum1","salmon1","yellow3")) 

#plot_ordination(input_RhinoAdult, ord, 'Age', 'Individual', hulls=TRUE) + scale_shape_manual(values = c(15,16,18,17,1,2,9,13))


otu <- t(input_RhinoAdult$data_loaded)
meta <- input_RhinoAdult$map_loaded
tax <- input_RhinoAdult$taxonomy_loaded

# Relative abundance
rel <- decostand(otu,"total")

# Hellinger transformation
varespec.hel <- decostand(rel, method="hellinger")

# Bray-Curtis dissimilarity
varespec.bray <- vegdist(varespec.hel, method="bray")
varespec.pcoa <- cmdscale(varespec.bray, k=nrow(otu)-1, eig=T)
eig2 <- eigenvals(varespec.pcoa)
eig2 / sum(eig2) 

# Permanova and Permdisp
adonis(varespec.bray ~ meta$Individual, permutations = 1000)
# Significant p = 0.001, R2 = 0.148, F = 5.0115
pwpm = pairwise.perm.manova(varespec.bray, fact = meta$Individual, nperm = 1000) # All sig
print(pwpm)




##***************************
## INTRODUCTION TO SEQUENCE CLUSTERING 
##
## Karl Cottenie
##
## 2025-10-16
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
conflicted::conflicts_prefer(lubridate::setdiff)
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
library(Biostrings)

# Startup ends here

## _ Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

####PART 1 - DATA DESCRIPTION ----

#Info about this specific dataset. Mysis is a genus of freshwater crustaceans. Code I used to download the data, on October 20, 2021
Mysis <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Mysis&format=tsv")
#code used to write data to a file
write_tsv(Mysis, "Mysis_BOLD.tsv")

#Next, we are defining the amount of missing data we will accept for this analysis. This parameter is defined as the proportion of a sequence comprised of Ns, after trimming off terminal Ns and removing all gaps (-). Here, we are using 0.01 as our starting value, which means we will accept up to 1% internal Ns in our DNA sequences.
missing.data <- 0.01

#Next, we are defining the amount of sequence length variability we will accept for our analysis. As a start, I am indicating 50 (nucleotides). We will then use this variable downstream to indicate that we will accept 50 nucleotides above and below the median sequence length.
length.var <- 50

####PART 2 - LOADING AND FILTERING DATA----

#read data. Note here that this line is now flexible for various inputs we might define at the top of the script. Check out the help for the function paste() if this is new.
df_mysis <- read_tsv("../data/Mysis_BOLD.tsv")
df_mysis

#Filtering data: filter out rows missing sequences. Filtering for COI-5P markercode only. End trimming of Ns (any number). Removing gaps (-). Note we aren't removing internal Ns as usually those are a placeholder for an actual nucleotide; so, we would want internal Ns in for downstream alignment. After end trimming, filtering out sequences with more than 1% Ns, as here we want high-quality sequences for clustering. Then, we are keeping data within a certain number of nucleotides (defined above) of the median sequence length. That is a choice that could be reconsidered, particularly if a dataset is small. We are putting our edited (i.e. trimmed) nucleotides in a new column so we can compare. When writing such code, each line should be individually checked before putting this all together!
df_mysis <- df_mysis %>%
  select(processid, genus_name, species_name, country, markercode, nucleotides) %>%
  filter(!is.na(nucleotides)) %>%
  filter(markercode == "COI-5P") %>%
  mutate(nucleotides2 = str_remove_all(nucleotides, "^N+|N+$|-")) %>%
  filter(str_count(nucleotides2, "N") <= (missing.data * str_count(nucleotides))) %>%
  filter(str_count(nucleotides2) >= median(str_count(nucleotides2)) - length.var & str_count(nucleotides2) <= median(str_count(nucleotides2)) + length.var)

#some checks
dim(df_mysis)
unique(df_mysis$markercode)
unique(df_mysis$species_name)
unique(df_mysis$genus_name)
sum(is.na(df_mysis$nucleotides2))
sum(str_count(df_mysis$nucleotides2, "-"))
summary(str_count(df_mysis$nucleotides2))

####PART 3 - FEATURE EXTRACTION----

#Next, we are converting the nucleotides to a DNAStringSet (class) so that we can use functions from the Biostrings package. Note that this time we are editing in place (using same column name before and after assignment operator), as we have already checked out this function before.
df_mysis <- as.data.frame(df_mysis)
df_mysis$nucleotides2 <- DNAStringSet(df_mysis$nucleotides2)
class(df_mysis$nucleotides2)

#Adding dinucleotide frequency (k-mers of length 2), here using proportions as that let's us account for variability in sequence lengths.
df_mysis <- cbind(df_mysis, as.data.frame(dinucleotideFrequency(df_mysis$nucleotides2, as.prob = TRUE)))

#Adding trinucleotide frequency (k-mers of length 3)
df_mysis <- cbind(df_mysis, as.data.frame(trinucleotideFrequency(df_mysis$nucleotides2, as.prob = TRUE)))

df_mysis |> View()

names(df_mysis)

####PART 4 - HIERARCHICAL CLUSTERING AND VISUALIZATION----

# this clustering method starts with a distance matrix 
# there are literally 100's of distance measures
# see for instance the p-distance in the k-mer activity sheet
# here I picked the euclidean distances based on the different k-mer frequencies
---
  ####PART 4 - HIERARCHICAL CLUSTERING AND VISUALIZATION (MODIFIED)----

####PART 4 - HIERARCHICAL CLUSTERING AND VISUALIZATION (MODIFIED)----

# Recalculate the distance matrix based on k-mer frequencies
# Exclude the first 7 columns, which are metadata
dist_mysis = dist(df_mysis[,-(1:7)], method = "euclidean")

# Perform hierarchical clustering using complete linkage
hier_mysis = hclust(dist_mysis, method = "complete")

# Plot the dendrogram labeled by species
plot(hier_mysis, labels = df_mysis$species_name, cex = 0.5)

# Highlight 5 clusters with a red border
rect.hclust(hier_mysis, k = 5, border = "red")

# Optional: highlight clusters using a height cutoff of 0.02 with blue border
rect.hclust(hier_mysis, h = 0.02, border = "blue")

# Save the figure to the figures folder with your UoG username
if(!dir.exists("figures")) dir.create("figures")
ggsave("figures/vlelo_cluster.png", plot = last_plot(), width = 6, height = 4)

####PART 5 - NEXT STEPS----
# code challenges
# use different levels of quality control (e.g., missing data, or length variability)
# use different types of features (e.g., more or fewer k-mer)
# use different types of distance measures
# use different types of clustering algorithms
# use different types of visualizations
# use different types of cluster identifications


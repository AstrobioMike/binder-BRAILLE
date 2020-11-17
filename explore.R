## goes along with our work on OSF here: https://osf.io/uhk48/wiki/home/
  ## write Mike on slack if you don't have and want access to that (not needed here)

## an intro to doing stuff here can be found here: https://hackmd.io/@astrobiomike/BRAILLE-mg-exploring
  ## including some ways to search for KO terms of interest

# this loads the data in R
load("BRAILLE-mg.RData")

# these are all packages of code we need to load also
library(tidyverse)
library(vegan)
library(dendextend)
library(phyloseq)
library(ggrepel)
library(ggpubr)
library(KEGGREST)
library(stringi)

## examples plotting single KO term with coverage normalized to all
plot_single_KO("K02588")

plot_single_KO("K10944")

plot_single_KO("K10944", coverage_round_acc = 10)

## example plotting function normalized to specific metabolism
plot_single_KO("K10944", methane_KO_gene_df_list)

## example plotting most abundant KOs in a metabolism
plot_overview_of_KOs_with_highest_coverages(N_KO_gene_df_list)

## adjusting the total number plotted
plot_overview_of_KOs_with_highest_coverages(N_KO_gene_df_list, number_to_plot = 12)


# N_KO_gene_df_list
# methane_KO_gene_df_list
# sulfur_KO_gene_df_list
# carbon_KO_gene_df_list
# diverse_microbial_metabolism_KO_gene_df_list

# tables with KO information for each group are available in the browser at the bottom right under the "Files" tab
# which can be opened to look across which KOs are associated with each

## e.g. plotting a KO normalized to Carbon metabolism grouping
plot_single_KO("K00196", carbon_KO_gene_df_list)

## e.g. plotting top 12 most abundant in the Diverse microbial metabolism pathway
plot_overview_of_KOs_with_highest_coverages(diverse_microbial_metabolism_KO_gene_df_list, number_to_plot = 12)

## function examples
get_KO_info("K00196")

plot_single_KO("K10944")

plot_overview_of_KOs_with_highest_coverages(all_our_KO_gene_df_list)

plot_overview_of_KOs_with_highest_coverages(all_our_KO_gene_df_list, number_to_plot = 4)


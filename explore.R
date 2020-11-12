## goes along with our work on OSF here: https://osf.io/uhk48/wiki/home/
## write Mike on slack if you don't have access to that

# these all need to be loaded first
library(tidyverse)
library(vegan)
library(dendextend)
library(phyloseq)
library(ggrepel)
library(ggpubr)
library(KEGGREST)
library(stringi)

# this loads the data in R
load("BRAILLE-mg.RData")

# currently can plot any KO term by metabolism group
# ones here now are:
   # nitrogen: N_KO_gene_df_list
   # methane:  methane_KO_gene_df_list
   # sulfur:   sulfur_KO_gene_df_list

# tables with KO information for each group are available in the browser at the bottom right under the "Files" tab

# e.g. to plot K10946	pmoC-amoC	methane/ammonia monooxygenase subunit C from the nitrogem metabolism overview
plot_single_KO("K10946", N_KO_gene_df_list)

# or K11180	dsrA	dissimilatory sulfite reductase alpha subunit from the sulfur metabolism overview
plot_single_KO("K11180", sulfur_KO_gene_df_list)
  # that gives a bunch of warnings and empty plot
  # can check the counts like so:
sulfur_KO_gene_df_list[["counts"]] %>% filter(KO_ID == "K11180")
  # seems there are no hits for that one

# will update this soon with more functionality and examples

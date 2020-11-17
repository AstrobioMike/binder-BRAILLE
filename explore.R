## goes along with our work on OSF here: https://osf.io/uhk48/wiki/home/
## write Mike on slack if you don't have and want access to that (not needed here)
## an intro to doing stuff here can be found here: https://hackmd.io/@astrobiomike/BRAILLE-mg-exploring

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


# currently can plot any KO term by metabolism group
# ones here now are:
   # nitrogen: N_KO_gene_df_list
   # methane:  methane_KO_gene_df_list
   # sulfur:   sulfur_KO_gene_df_list
   # carbon:   carbon_KO_gene_df_list

# tables with KO information for each group are available in the browser at the bottom right under the "Files" tab

# e.g. to plot K10946	pmoC-amoC	methane/ammonia monooxygenase subunit C from the nitrogem metabolism overview
plot_single_KO("K10946", N_KO_gene_df_list)

# if we check for something that wasn't detected, we'll get a message telling us that
# e.g. K02588	nifH	nitrogenase iron protein NifH
plot_single_KO("K02588", N_KO_gene_df_list)


# will update this soon with more functionality and examples

library(tidyverse)
library(vegan)
library(dendextend)
library(phyloseq)
library(ggrepel)
library(ggpubr)
library(KEGGREST)
library(stringi)


KO_tab <- read.table("../combined-outputs/All-KO-function-coverages.tsv", sep="\t", quote="", header=TRUE)
head(KO_tab)
dim(KO_tab)

# some color stuff we'll be re-using
VAL_col <- "#3978e6"
YEL_col <- "#e35b5b"
both_col <- "#b153c2"

VAL_darkest_col <- "#0e5ce6"
VAL_lightest_col <- "#aabfe3"
YEL_darkest_col <- "#e32222"
YEL_lightest_col <- "#e0a6a6"

# vectors of color gradients for each
VAL_gradient <- colorRampPalette(c(VAL_lightest_col, VAL_darkest_col))(5)
YEL_gradient <- colorRampPalette(c(YEL_lightest_col, YEL_darkest_col))(5)


# vectors of color gradients for each
VAL_gradient <- colorRampPalette(c(VAL_lightest_col, VAL_darkest_col))(5)
YEL_gradient <- colorRampPalette(c(YEL_lightest_col, YEL_darkest_col))(5)

## adding some summary stats
KO_norm_sums_across_samples <- rowSums(KO_tab[3:12])
KO_norm_means_across_samples <- rowMeans(KO_tab[3:12])
KO_mean_rank_percentile <- round(rank(KO_norm_means_across_samples) / length(KO_norm_means_across_samples) * 100)
KO_norm_stdev_across_samples <- apply(KO_tab[3:12], MARGIN=1, function(x) sd(x))
KO_norm_coef_of_var_across_samples <- KO_norm_stdev_across_samples / KO_norm_means_across_samples

# sticking to table

KO_tab_w_summary_stats <- data.frame(KO_tab[,c(1,3:12)], "sum"=round(KO_norm_sums_across_samples,2), "mean"=round(KO_norm_means_across_samples,2), "mean_rank_perc"=KO_mean_rank_percentile, "stdev"=round(KO_norm_stdev_across_samples,2), "CoV"=round(KO_norm_coef_of_var_across_samples, 2), KO_tab[,2, drop=F], stringsAsFactors=F)

# starting with hierarchical clustering
all_KO_only_bray_dist <- vegdist(t(KO_tab[, 3:12]))
all_KO_only_bray_hclust <- hclust(all_KO_only_bray_dist, method="ward.D2")
all_KO_only_bray_dendro <- as.dendrogram(all_KO_only_bray_hclust, hang=0.15)

  # coloring labels based on source cave (and creating a vector of source cave)
col_vec <- vector()
source_vec <- vector()
for (label in labels(all_KO_only_bray_dendro)) {
    if (startsWith(label, "VAL")) {
        col_vec <- c(col_vec, "#3978e6")
        source_vec <- c(source_vec, "VAL")
    } else {
        col_vec <- c(col_vec, "#e35b5b")
        source_vec <- c(source_vec, "YEL")
    }
}

labels_colors(all_KO_only_bray_dendro) <- col_vec

plot(all_KO_only_bray_dendro, main="Clustering based on normalized KO-term coverages", ylab="Bray-Curtis dissimilarity")

ggdendrogram(all_KO_only_bray_dendro)

### PCoA view based on these normalized gene-coverages
KO_df <- data.frame(row.names=KO_tab_w_summary_stats$KO_ID, KO_tab_w_summary_stats[,2:11])
KO_norm_phy <- otu_table(KO_df, taxa_are_rows=TRUE)
sample_info_df <- data.frame(row.names=names(KO_tab_w_summary_stats)[2:11], sample_ID=names(KO_tab_w_summary_stats)[2:11], source=c(rep("VAL", 5), rep("YEL", 5)))
sample_info_phy <- sample_data(sample_info_df)
KO_norm_physeq_obj <- phyloseq(KO_norm_phy, sample_info_phy)

KO_bray_pcoa <- ordinate(KO_norm_physeq_obj, method="MDS")
KO_bray_eigen_vals <- KO_bray_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

KO_PCoA_plot <- plot_ordination(KO_norm_physeq_obj, KO_bray_pcoa, color="source") +
    coord_fixed(sqrt(KO_bray_eigen_vals[2]/KO_bray_eigen_vals[1])) +
    scale_color_manual(values=c(VAL_col, YEL_col)) + theme_bw() +
    ggtitle("PCoA of normalized KO terms", subtitle="Based on Bray-Curtis dissimilarities") +
    geom_text_repel(aes(label=sample_ID)) + theme(legend.position="none")

KO_PCoA_plot


## looking for DosSantos et al N-fixation genes
Dos_Santos_N_fixation_KOs <- c("K02585", "K02586", "K02587", "K02588", "K02591", "K02592")

KO_tab %>% filter(KO_ID %in% DosSantos_et_al_2012_N_fixation_KOs)
# none present


## Parsing down to focus on specific KOs based on KEGG metabolism groups

listDatabases()

#### helper functions ####
# helper function to generate KO summary table:
make_KO_summary_tab <- function(target_KO_list) {

      # the kegg API limits individual requests to 10, so breaking the total we want into blocks of 10
    list_of_pathway_KO_blocks <- split(target_KO_list, ceiling(seq_along(target_KO_list)/10))
    num_blocks <- length(list_of_pathway_KO_blocks)

      # initializing some vectors we're going to populate
    entry_vec <- vector()
    name_vec <- vector()
    def_vec <- vector()
    module_ids_vec <- vector()
    module_vec <- vector()
    path_ids_vec <- vector()
    path_vec <- vector()

    # iterating through each block of at most 10 KO terms
    for ( block in seq(1, num_blocks) ) {

        # getting and storing the current block of KO terms' information
        current <- keggGet(list_of_pathway_KO_blocks[[block]])

        # here iterating through that stored block of information, one term at a time to get the info we want for each
        for ( num in seq(1, length(current)) ) {

            # KO ID
            current_entry <- current[[num]]$ENTRY
            if ( length(current_entry) > 0 ) {
                entry_vec <- c(entry_vec, current_entry)
            } else {
                entry_vec <- c(entry_vec, NA)
            }

            # KO name(s)
            current_name <- current[[num]]$NAME
            if ( length(current_name) > 0 ) {
                name_vec <- c(name_vec, current_name)
            } else {
                name_vec <- c(name_vec, NA)
            }

            # KO definition
            current_def <- current[[num]]$DEFINITION
            if ( length(current_def) > 0 ) {
                # getting rid of EC number if present
                current_def <- sub(current_def, pattern=" \\[EC:.*$", replacement="")
                def_vec <- c(def_vec, current_def)
            } else {
                def_vec <- c(def_vec, NA)
            }

            # module ID(s)
            current_module_id <- stri_join(names(current[[num]]$MODULE), collapse="; ")
            if ( length(current_module_id) > 0 ) {
                module_ids_vec <- c(module_ids_vec, current_module_id)
            } else {
                module_ids_vec <- c(module_ids_vec, NA)
            }

            # module definition(s)
            current_module <- stri_join(current[[num]]$MODULE, collapse="; ")
            if ( length(current_module) > 0 ) {
                module_vec <- c(module_vec, current_module)
            } else {
                module_vec <- c(module_vec, NA)
            }

            # pathway ID(s)
            current_path_id <- stri_join(names(current[[num]]$PATHWAY), collapse="; ")
            if ( length(current_path_id) > 0 ) {
                path_ids_vec <- c(path_ids_vec, current_path_id)
            } else {
                path_ids_vec <- c(path_ids_vec, NA)
            }

            # pathway definition(s)
            current_path <- stri_join(current[[num]]$PATHWAY, collapse="; ")
            if ( length(current_path) > 0 ) {
                path_vec <- c(path_vec, current_path)
            } else {
                path_vec <- c(path_vec, NA)
            }

        }
    }

      # now combining into table we are returning
    return(data.frame("KO_ID"=entry_vec, "KO_name"=name_vec, "KO_def"=def_vec, "module_IDs"=module_ids_vec, "module_defs"=module_vec, "pathway_IDs"=path_ids_vec, "pathway_defs"=path_vec, stringsAsFactors=F))

}

# helper function to build table of how many of each target KO term were annotated from each sample and table of normalized coverage
parse_samples_by_target_metabolism <- function( samples, target_metabolism_tab ) {

    target_KOs <- target_metabolism_tab$KO_ID
    # making empty tables to fill
    count_df <- data.frame(matrix(0, ncol=length(all_samples) + 3, nrow=length(target_KOs)), row.names=target_KOs)
    names(count_df) <- c("KO_ID", all_samples, "KO_name", "KO_def")
    count_df$KO_ID <- target_KOs
    count_df$KO_name <- target_metabolism_tab$KO_name
    count_df$KO_def <- target_metabolism_tab$KO_def

    cov_df <- count_df

    for ( sample in samples ) {

        cat("On sample: ", sample, "\n")
        file <- paste0("../annotations-and-taxonomy/", sample, "-gene-coverage-annotation-and-tax.tsv")

        # reading in table
        curr_tab <- read.table(file, header=TRUE, quote="", sep="\t")

        # filtering for target KOs
        curr_sub_tab <- curr_tab %>% filter(KO_ID %in% target_KOs)

        # iterating through table (well as much as we can in R)
        for ( row in 1:nrow(curr_sub_tab) ) {

            curr_KO <- as.vector(curr_sub_tab[row, "KO_ID"])

            # tracking how many genes for each KO
            count_df[curr_KO, sample] <- count_df[curr_KO, sample] + 1
            # building table of coverages summed by KO term
            cov_df[curr_KO, sample] <- cov_df[curr_KO, sample] + curr_sub_tab[row, "coverage"]

        }

    }

    # normalizing coverage table to copies per million
    cov_df[,2:11] <- cov_df[,2:11] %>% apply(2, function(x) ( x / sum(x)) * 1000000 ) %>% as.data.frame()

    # generating some summary stats for each KO across samples
      # for gene copies

    curr_count_sums <- rowSums(count_df[,2:11])
    curr_count_means <- rowMeans(count_df[,2:11])
    curr_count_stdevs <- apply(count_df[,2:11], 1, function(x) sd(x))
    curr_count_coefs_of_var <- curr_count_stdevs / curr_count_means

      # and for normalized coverages
    curr_cov_sums <- rowSums(cov_df[,2:11])
    curr_cov_means <- rowMeans(cov_df[,2:11])
    curr_cov_stdevs <- apply(cov_df[,2:11], 1, function(x) sd(x))
    curr_cov_coefs_of_var <- curr_cov_stdevs / curr_cov_means

    # building final tables
    count_df <- data.frame(count_df[,1:11], "sum"=curr_count_sums, "mean"=curr_count_means, "stdev"=curr_count_stdevs, "CoV"=curr_count_coefs_of_var, count_df[,12:13], row.names=NULL)
    cov_df <- data.frame(cov_df[,1:11], "sum"=curr_cov_sums, "mean"=curr_cov_means, "stdev"=curr_cov_stdevs, "CoV"=curr_cov_coefs_of_var, cov_df[,12:13], row.names=NULL)

    # making one that holds the KEGG annotations for this set
    ko_df <- count_df %>% select(KO_ID, KO_name, KO_def)
    return(list("counts"=count_df, "cov"=cov_df, "annots"=ko_df))

}

# helper function for plotting cov and num unique genes of single KO within a metabolism
plot_single_KO <- function( KO_term, metabolism_df_list, coverage_round_acc = 1000, gene_copy_round_acc = 10) {

    # making long format tables for easier plotting
    cov_long <- metabolism_df_list[["cov"]][,1:11] %>% pivot_longer(-KO_ID, names_to="Sample", values_to="Norm_cov")
    counts_long <- metabolism_df_list[["counts"]][,1:11] %>% pivot_longer(-KO_ID, names_to="Sample", values_to="Num_copies")

    # subsetting to target KO
    target_KO_cov <- cov_long %>% filter(KO_ID == KO_term)
    target_KO_counts <- counts_long %>% filter(KO_ID == KO_term)

    # getting upper bounds for y-axes and transformation to be able to plot on same
       # little hack required here because Hadley doesn't like it: https://stackoverflow.com/questions/3099219/ggplot-with-2-y-axes-on-each-side-and-different-scales/3101876#3101876
       # this site helped me out: https://whatalnk.github.io/r-tips/ggplot2-secondary-y-axis.nb.html
    y_cov_limit <- plyr::round_any(max(target_KO_cov$Norm_cov), coverage_round_acc, ceiling)
    y_gene_count_limit <- plyr::round_any(max(target_KO_counts$Num_copies), gene_copy_round_acc, ceiling)

    # making label
    num_uniq_genes <- metabolism_df_list[["counts"]] %>% filter(KO_ID == KO_term) %>% pull(sum)
    target_KO_def <- metabolism_df_list[["counts"]] %>% filter(KO_ID == KO_term) %>% pull(KO_def)
    label <- paste0(KO_term, " - ", target_KO_def, " (", num_uniq_genes, " copies recovered)")
    names(label) <- KO_term

    # plotting
    my_plot <- ggplot() +
        geom_bar(data=target_KO_cov, aes(y=Norm_cov, x=Sample, fill=Sample), stat="identity", width=0.75) +
        scale_fill_manual(values = c(VAL_gradient, YEL_gradient)) + facet_wrap(. ~ KO_ID, scales="free_y", labeller = labeller(KO_ID=label)) + theme_bw() +
        geom_point(data=target_KO_counts, aes(x=Sample, y = Num_copies * y_cov_limit / y_gene_count_limit ), shape = 3, size=3, stroke = 1) +
        scale_y_continuous(sec.axis = sec_axis(~ . * y_gene_count_limit / y_cov_limit , name = "Num. gene-copies recovered (+)"), limits = c(0,y_cov_limit)) +
        theme(axis.title.y = element_text(face = "bold"), axis.text.y = element_text(size=12)) +
        theme(strip.text = element_text(size = 12, face="bold"), axis.title.x=element_blank()) +
        theme(axis.text.x = element_text(size = 10, face="bold", angle = 45, hjust = 1)) +
        labs(fill = NULL, y="Norm. Cov.") + theme(legend.position = "none")


    return(my_plot)
}

# helper functions for plotting overview of many KO terms
plot_single_for_plotting_multiple_KOs <- function( KO_term, metabolism_df_list, coverage_round_acc = 1000, gene_copy_round_acc = 10) {

    # making long format tables for easier plotting
    cov_long <- metabolism_df_list[["cov"]][,1:11] %>% pivot_longer(-KO_ID, names_to="Sample", values_to="Norm_cov")
    counts_long <- metabolism_df_list[["counts"]][,1:11] %>% pivot_longer(-KO_ID, names_to="Sample", values_to="Num_copies")

    # subsetting to target KO
    target_KO_cov <- cov_long %>% filter(KO_ID == KO_term)
    target_KO_counts <- counts_long %>% filter(KO_ID == KO_term)

    # getting upper bounds for y-axes and transformation to be able to plot on same
       # little hack required here because Hadley doesn't like it: https://stackoverflow.com/questions/3099219/ggplot-with-2-y-axes-on-each-side-and-different-scales/3101876#3101876
       # this site helped me out: https://whatalnk.github.io/r-tips/ggplot2-secondary-y-axis.nb.html
    y_cov_limit <- plyr::round_any(max(target_KO_cov$Norm_cov), coverage_round_acc, ceiling)
    y_gene_count_limit <- plyr::round_any(max(target_KO_counts$Num_copies), gene_copy_round_acc, ceiling)

    # making label
    num_uniq_genes <- metabolism_df_list[["counts"]] %>% filter(KO_ID == KO_term) %>% pull(sum)
    label <- paste0(KO_term, " (", num_uniq_genes, ")")
    names(label) <- KO_term

    # plotting
    my_plot <- ggplot() +
        geom_bar(data=target_KO_cov, aes(y=Norm_cov, x=Sample, fill=Sample), stat="identity", width=0.75) +
        scale_fill_manual(values = c(VAL_gradient, YEL_gradient)) + facet_wrap(. ~ KO_ID, scales="free_y", labeller = labeller(KO_ID=label)) + theme_bw() +
        geom_point(data=target_KO_counts, aes(x=Sample, y = Num_copies * y_cov_limit / y_gene_count_limit ), shape = 3, size=2) +
        theme(axis.title.y = element_blank(), axis.text.y = element_text(size=12)) +
        theme(strip.text = element_text(size = 10, face="bold"), axis.title.x=element_blank()) +
        theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position="bottom") +
        labs(fill = NULL, y="Norm. Cov.") + guides(fill = guide_legend(nrow = 1))


    return(my_plot)
}
plot_overview <- function( metabolism_df_list, mean_coverage, title="Gene coverages and relative copy numbers") {

    # getting target KOs that have >= the specified mean coverage
    target_KOs <- metabolism_df_list[["cov"]] %>% filter(mean >= mean_coverage) %>% arrange(desc(CoV)) %>% pull(KO_ID)
    target_KOs <- factor(target_KOs, levels = target_KOs)

    # initializing plotting list
    curr_plots <- vector("list", length(target_KOs))

    # iterating through generating each
    for ( num in 1:length(target_KOs) ) {
        my_plots[[num]] <- plot_single_for_plotting_multiple_KOs( target_KOs[num], metabolism_df_list)
    }

    # making combined plot
    my_plot <- ggarrange(plotlist=my_plots, common.legend=TRUE, legend = "bottom")

    my_plot <- annotate_figure(my_plot, top = text_grob(title, face="bold", size=14), left = text_grob("Norm. Cov. (CPM)", face="bold", size=11, rot=90))

    return(my_plot)
}

##################

## nitrogen metabolism ##
# searching KEGG pathways for Nitrogen
N_pathway_info <- keggFind(database = "pathway", query = "nitrogen")

# getting "Nitrogen metabolism" pathway ID
N_pathway_ID <- sub(names(N_pathway_info), pattern = "path:", replacement = "")

# getting info for all KOs associated with "Nitrogen metabolism" pathway
N_pathway_KOs_info <- keggLink("ko", N_pathway_ID)

# getting KO IDs alone
N_pathway_KOs <- as.character(sub(N_pathway_KOs_info, pattern="ko:", replacement=""))

length(N_pathway_KOs) # 65

Nitrogen_metabolism_KEGG_tab <- make_KO_summary_tab(N_pathway_KOs)

all_samples <- names(KO_tab[3:12])

N_KO_gene_df_list <- parse_samples_by_target_metabolism(all_samples, Nitrogen_metabolism_KEGG_tab)

# making directory and writing out
dir.create("Nitrogen-focus")

write.table(N_KO_gene_df_list[["counts"]], "Nitrogen-focus/N-metabolism-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(N_KO_gene_df_list[["cov"]], "Nitrogen-focus/N-metabolism-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(N_KO_gene_df_list[["annots"]], "Nitrogen-focus/N-metabolism-KO-info.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")

head(N_KO_gene_df_list[["counts"]])

tab <- N_KO_gene_df_list[["counts"]] %>% select(KO_ID, KO_name, KO_def)

## which vary the least
N_KO_gene_df_list[["cov"]] %>% arrange(CoV) %>% head()

## plotting this KO term
plot_single_KO("K00266", N_KO_gene_df_list)

# histogram of all CoVs
ggplot(N_KO_gene_df_list[["cov"]] %>% filter(! is.na(CoV)), aes(x=CoV)) +
    geom_histogram(fill=both_col, color="black", boundary = 1, binwidth=0.1) + theme_bw() +
    ggtitle("N-KO-term Coefficients of Variation") + ylab(label="Num. KO terms") +
    scale_y_continuous(limit = c(0,7), oob = function(x, limits) x) +
    scale_x_continuous(limit = c(0,3.5), oob = function(x, limits) x)

## highest CoVs
N_KO_gene_df_list[["cov"]] %>% arrange(desc(CoV)) %>% head()

## let's filter by mean coverage
summary(N_KO_gene_df_list[["cov"]]$mean)

N_KO_gene_df_list[["cov"]] %>% arrange(desc(CoV)) %>% filter(mean >= 4000) %>% head()

## looking at top 3 all methane/ammonia monooxygenase subunits A/B/C
plot_single_KO("K10946", N_KO_gene_df_list)
plot_single_KO("K10945", N_KO_gene_df_list)
plot_single_KO("K10944", N_KO_gene_df_list)


### Overview
# making overview of those with greater than 10000 mean coverage per million ordered by decreasing CoV
plot_overview(N_KO_gene_df_list, 10000, "N-related gene coverages and relative copy numbers")


## methane metabolism ##
# searching KEGG pathways for methane
methane_pathway_info <- keggFind(database = "pathway", query = "methane")

# getting "methane metabolism" pathway ID
methane_pathway_ID <- sub(names(methane_pathway_info), pattern = "path:", replacement = "")

# getting info for all KOs associated with "methane metabolism" pathway
methane_pathway_KOs_info <- keggLink("ko", methane_pathway_ID)

# getting KO IDs alone
methane_pathway_KOs <- as.character(sub(methane_pathway_KOs_info, pattern="ko:", replacement=""))

length(methane_pathway_KOs) # 190

methane_metabolism_KEGG_tab <- make_KO_summary_tab(methane_pathway_KOs)

all_samples <- names(KO_tab[3:12])

methane_KO_gene_df_list <- parse_samples_by_target_metabolism(all_samples, methane_metabolism_KEGG_tab)

# making directory and writing out
dir.create("Methane-focus")

write.table(methane_KO_gene_df_list[["counts"]], "Methane-focus/Methane-metabolism-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(methane_KO_gene_df_list[["cov"]], "Methane-focus/Methane-metabolism-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(methane_KO_gene_df_list[["annots"]], "Methane-focus/Methane-metabolism-KO-info.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")


## which vary the least
methane_KO_gene_df_list[["cov"]] %>% arrange(CoV) %>% head()

## plotting this KO term
plot_single_KO("K00600", methane_KO_gene_df_list)

# histogram of all CoVs
ggplot(methane_KO_gene_df_list[["cov"]] %>% filter(! is.na(CoV)), aes(x=CoV)) +
    geom_histogram(fill=both_col, color="black", boundary = 1, binwidth=0.1) + theme_bw() +
    ggtitle("Methane-KO-term Coefficients of Variation") + ylab(label="Num. KO terms") +
    scale_y_continuous(limit = c(0,18), oob = function(x, limits) x) +
    scale_x_continuous(limit = c(0,3.5), oob = function(x, limits) x)

## highest CoVs
methane_KO_gene_df_list[["cov"]] %>% arrange(desc(CoV)) %>% head()

## let's filter by mean coverage
summary(methane_KO_gene_df_list[["cov"]]$mean)

methane_KO_gene_df_list[["cov"]] %>% arrange(desc(CoV)) %>% filter(mean >= 5000) %>% head()

## looking at top 3
plot_single_KO("K17068", methane_KO_gene_df_list)
plot_single_KO("K00148", methane_KO_gene_df_list)
plot_single_KO("K05979", methane_KO_gene_df_list)


### Overview
# making overview of those with greater than 14200 mean coverage per million ordered by decreasing CoV (just so we don't over over 20 in the summary fig)
plot_overview(methane_KO_gene_df_list, 14200, "Methane-related gene coverages and relative copy numbers")


## sulfur metabolism ##
# searching KEGG pathways for sulfur
sulfur_pathway_info <- keggFind(database = "pathway", query = "sulfur")

# getting "sulfur metabolism" pathway ID
sulfur_pathway_ID <- sub(names(sulfur_pathway_info), pattern = "path:", replacement = "")

# getting info for all KOs associated with "sulfur metabolism" pathway
sulfur_pathway_KOs_info <- keggLink("ko", sulfur_pathway_ID)

# getting KO IDs alone
sulfur_pathway_KOs <- unique(as.character(sub(sulfur_pathway_KOs_info, pattern="ko:", replacement="")))

length(sulfur_pathway_KOs) # 137

sulfur_metabolism_KEGG_tab <- make_KO_summary_tab(sulfur_pathway_KOs)

all_samples <- names(KO_tab[3:12])

sulfur_KO_gene_df_list <- parse_samples_by_target_metabolism(all_samples, sulfur_metabolism_KEGG_tab)

# making directory and writing out
dir.create("Sulfur-focus")

write.table(sulfur_KO_gene_df_list[["counts"]], "Sulfur-focus/Sulfur-metabolism-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(sulfur_KO_gene_df_list[["cov"]], "Sulfur-focus/Sulfur-metabolism-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(sulfur_KO_gene_df_list[["annots"]], "Sulfur-focus/Sulfur-metabolism-KO-info.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")

## which vary the least
sulfur_KO_gene_df_list[["cov"]] %>% arrange(CoV) %>% head()

## plotting this KO term
plot_single_KO("K03635", sulfur_KO_gene_df_list)

# histogram of all CoVs
ggplot(sulfur_KO_gene_df_list[["cov"]] %>% filter(! is.na(CoV)), aes(x=CoV)) +
    geom_histogram(fill=both_col, color="black", boundary = 1, binwidth=0.1) + theme_bw() +
    ggtitle("Sulfur-KO-term Coefficients of Variation") + ylab(label="Num. KO terms") +
    scale_y_continuous(limit = c(0,15), oob = function(x, limits) x) +
    scale_x_continuous(limit = c(0,3.5), oob = function(x, limits) x)

## highest CoVs
sulfur_KO_gene_df_list[["cov"]] %>% arrange(desc(CoV)) %>% head()

## let's filter by mean coverage
summary(sulfur_KO_gene_df_list[["cov"]]$mean)

sulfur_KO_gene_df_list[["cov"]] %>% arrange(desc(CoV)) %>% filter(mean >= 7000) %>% head()

## looking at top 3
plot_single_KO("K22622", sulfur_KO_gene_df_list)
plot_single_KO("K21148", sulfur_KO_gene_df_list)
plot_single_KO("K23163", sulfur_KO_gene_df_list)


### Overview
# making overview of those with greater than 10000 mean coverage per million ordered by decreasing CoV
plot_overview(sulfur_KO_gene_df_list, 15000, "Sulfur-related gene coverages and relative copy numbers")


#########

## carbon metabolism ##
# searching KEGG pathways for carbon
carbon_pathway_info <- keggFind(database = "pathway", query = "carbon")

# getting "carbon metabolism" pathway ID
carbon_pathway_ID <- sub(names(carbon_pathway_info), pattern = "path:", replacement = "")

# getting info for all KOs associated with "carbon metabolism" pathway
carbon_pathway_KOs_info <- keggLink("ko", carbon_pathway_ID)

# getting KO IDs alone
carbon_pathway_KOs <- unique(as.character(sub(carbon_pathway_KOs_info, pattern="ko:", replacement="")))

length(carbon_pathway_KOs) # 494

carbon_metabolism_KEGG_tab <- make_KO_summary_tab(carbon_pathway_KOs)

all_samples <- names(KO_tab[3:12])

carbon_KO_gene_df_list <- parse_samples_by_target_metabolism(all_samples, carbon_metabolism_KEGG_tab)

# making directory and writing out
dir.create("Carbon-focus")

write.table(carbon_KO_gene_df_list[["counts"]], "Carbon-focus/Carbon-metabolism-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(carbon_KO_gene_df_list[["cov"]], "Carbon-focus/Carbon-metabolism-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")

## which vary the least
carbon_KO_gene_df_list[["cov"]] %>% arrange(CoV) %>% head()

## plotting this KO term
plot_single_KO("K03635", carbon_KO_gene_df_list)

# histogram of all CoVs
ggplot(carbon_KO_gene_df_list[["cov"]] %>% filter(! is.na(CoV)), aes(x=CoV)) +
    geom_histogram(fill=both_col, color="black", boundary = 1, binwidth=0.1) + theme_bw() +
    ggtitle("Carbon-KO-term Coefficients of Variation") + ylab(label="Num. KO terms") +
    scale_y_continuous(limit = c(0,15), oob = function(x, limits) x) +
    scale_x_continuous(limit = c(0,3.5), oob = function(x, limits) x)

## highest CoVs
carbon_KO_gene_df_list[["cov"]] %>% arrange(desc(CoV)) %>% head()

## let's filter by mean coverage
summary(carbon_KO_gene_df_list[["cov"]]$mean)

carbon_KO_gene_df_list[["cov"]] %>% arrange(desc(CoV)) %>% filter(mean >= 7000) %>% head()

## looking at top 3
plot_single_KO("K22622", carbon_KO_gene_df_list)
plot_single_KO("K21148", carbon_KO_gene_df_list)
plot_single_KO("K23163", carbon_KO_gene_df_list)



### Overview
# making overview of those with greater than 10000 mean coverage per million ordered by decreasing CoV
plot_overview(carbon_KO_gene_df_list, 15000, "Carbon-related gene coverages and relative copy numbers")
















save.image("BRAILLE-mg.RData")


## this is not meant to be run in the binder environment, large input files are not present
## it is here as reference as this is the code that produced the data and functions that are used in the binder
## see our OSF for more details: https://osf.io/uhk48/wiki/home/

# load("BRAILLE-mg.RData")

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

all_samples <- names(KO_tab)[3:12]

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

#### helper functions ####

# helper function to generate KO summary table:
make_KO_summary_tab <- function(target_KOs) {

    # removing "Not annotated" if provided, will be added as needed below
    target_KOs <- target_KOs[target_KOs != "Not annotated"]

    # the kegg API limits individual requests to 10, so breaking the total we want into blocks of 10
    list_of_pathway_KO_blocks <- split(target_KOs, ceiling(seq_along(target_KOs)/10))
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

    # now combining into table
    out_tab <- data.frame("KO_ID"=entry_vec, "KO_name"=name_vec, "KO_def"=def_vec, "module_IDs"=module_ids_vec, "module_defs"=module_vec, "pathway_IDs"=path_ids_vec, "pathway_defs"=path_vec, stringsAsFactors=F)

    # some KO terms may have been removed from KEGG since annotation was performed, and therefore wouldn't be found
    # putting a check and reporting if that's the case for any, adding to table to keep track of them, but making all values set to 'Not found at KEGG'
    got_KOs <- out_tab %>% pull(KO_ID)

    missed_KOs <- setdiff(target_KOs, got_KOs)

    if ( length(missed_KOs) != 0 ) {
        num_missed <- length(missed_KOs)

        cat(" ", num_missed, "KOs were not found at KEGG, they were likely removed after the utilized annotation db was created:\n\n")

        for ( term in missed_KOs ) {
            cat("\t", term, "\n")
            out_tab <- rbind(out_tab, c(term, rep("Not found at KEGG", 6)))
        }

        cat("\n")
        cat("  These have been added to the output table with 'Not found at KEGG' in all other fields besides the 'KO_ID' one.\n\n")

    }

    # sorting out table
    out_tab <- out_tab %>% arrange(KO_ID)

    return(out_tab)

}

# helper function to build table of how many of each target KO term were annotated from each sample and table of normalized coverage
parse_samples_by_target_KO_tab <- function( samples, target_KO_tab ) {

    target_KOs <- target_KO_tab$KO_ID
    # making empty tables to fill
    count_df <- data.frame(matrix(0, ncol=length(samples) + 3, nrow=length(target_KOs)), row.names=target_KOs)
    names(count_df) <- c("KO_ID", samples, "KO_name", "KO_def")
    count_df$KO_ID <- target_KOs
    count_df$KO_name <- target_KO_tab$KO_name
    count_df$KO_def <- target_KO_tab$KO_def

    cov_df <- count_df

    for ( sample in samples ) {

        cat("On sample: ", sample, "\n")
        file <- paste0("../annotations-and-taxonomy/", sample, "-gene-coverage-annotation-and-tax.tsv")

        # reading in table
        curr_tab <- read.table(file, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)

        # setting all KO_IDs that are NA to "Not annotated" so they are captured and counted in our building table if we are running this on everything
        curr_tab$KO_ID[is.na(curr_tab$KO_ID)] <- "Not annotated"

        # filtering for target KOs
        curr_sub_tab <- curr_tab %>% filter(KO_ID %in% target_KOs)

        # iterating through table (well as much as we can in R)
        for ( row in 1:nrow(curr_sub_tab) ) {

            curr_KO <- curr_sub_tab[row, "KO_ID"]

            # tracking how many genes for each KO (or not annotated)
            count_df[curr_KO, sample] <- count_df[curr_KO, sample] + 1
            # building table of coverages summed by KO term (or not annotated)
            cov_df[curr_KO, sample] <- cov_df[curr_KO, sample] + curr_sub_tab[row, "coverage"]

        }

    }

    # normalizing coverage table to copies per million
    cov_df[,2:(1 + length(samples))] <- cov_df[,2:(1 + length(samples))] %>% apply(2, function(x) ( x / sum(x)) * 1000000 ) %>% as.data.frame()

    # generating some summary stats for each KOs gene-counts across samples
    curr_count_sums <- rowSums(count_df[,2:(1 + length(samples))])
    curr_count_means <- rowMeans(count_df[,2:(1 + length(samples))])
    curr_count_stdevs <- apply(count_df[,2:(1 + length(samples))], 1, function(x) sd(x))
    curr_count_coefs_of_var <- curr_count_stdevs / curr_count_means

    # and for normalized coverages
    curr_cov_sums <- rowSums(cov_df[,2:(1 + length(samples))])
    curr_cov_means <- rowMeans(cov_df[,2:(1 + length(samples))])
    curr_cov_stdevs <- apply(cov_df[,2:(1 + length(samples))], 1, function(x) sd(x))
    curr_cov_coefs_of_var <- curr_cov_stdevs / curr_cov_means

    # building final tables
    count_df <- data.frame(count_df[,1:(1 + length(samples))], "sum"=curr_count_sums, "mean"=curr_count_means, "stdev"=curr_count_stdevs, "CoV"=curr_count_coefs_of_var, count_df[,c("KO_name", "KO_def")], row.names=NULL)
    cov_df <- data.frame(cov_df[,1:(1 + length(samples))], "sum"=curr_cov_sums, "mean"=curr_cov_means, "stdev"=curr_cov_stdevs, "CoV"=curr_cov_coefs_of_var, cov_df[,c("KO_name", "KO_def")], row.names=NULL)

    # making one that holds the KEGG annotations for this set
    ko_df <- count_df %>% select(KO_ID, KO_name, KO_def)
    return(list("counts"=count_df, "cov"=cov_df, "annots"=ko_df))

}

# helper function for plotting cov and num unique genes of single KO
plot_single_KO <- function( KO_term, target_KO_df_list = all_our_KO_gene_df_list, coverage_round_acc = 100, gene_copy_round_acc = 10) {

    # making sure it present/was assembled/annotated
    if ( ! KO_term %in% ( target_KO_df_list[['cov']] %>% pull(KO_ID) ) ) {

        cat("\n   ", KO_term, "was not detected.\n\n")

        # couldn't figure out a better way to just report this and not give an error while exiting
        stop_no_error <- function() {
            opt <- options(show.error.messages = FALSE)
            on.exit(options(opt))
            stop()
        }
        stop_no_error()
    }

    N_KO_gene_df_list[["cov"]] %>% filter(KO_ID == "K02588") %>% select(2:11) %>% sum()
    # making long format tables for easier plotting
    cov_long <- target_KO_df_list[["cov"]][,1:11] %>% pivot_longer(-KO_ID, names_to="Sample", values_to="Norm_cov")
    counts_long <- target_KO_df_list[["counts"]][,1:11] %>% pivot_longer(-KO_ID, names_to="Sample", values_to="Num_copies")

    # subsetting to target KO
    target_KO_cov <- cov_long %>% filter(KO_ID == KO_term)
    target_KO_counts <- counts_long %>% filter(KO_ID == KO_term)

    # getting upper bounds for y-axes and transformation to be able to plot on same
       # little hack required here because Hadley doesn't like it: https://stackoverflow.com/questions/3099219/ggplot-with-2-y-axes-on-each-side-and-different-scales/3101876#3101876
       # this site helped me out: https://whatalnk.github.io/r-tips/ggplot2-secondary-y-axis.nb.html
    y_cov_limit <- plyr::round_any(max(target_KO_cov$Norm_cov), coverage_round_acc, ceiling)
    y_gene_count_limit <- plyr::round_any(max(target_KO_counts$Num_copies), gene_copy_round_acc, ceiling)

    # making label
    num_uniq_genes <- target_KO_df_list[["counts"]] %>% filter(KO_ID == KO_term) %>% pull(sum)
    target_KO_def <- target_KO_df_list[["counts"]] %>% filter(KO_ID == KO_term) %>% pull(KO_def)
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
plot_single_for_plotting_multiple_KOs <- function( KO_term, target_KO_df_list, coverage_round_acc = 1000, gene_copy_round_acc = 10) {

    # making long format tables for easier plotting
    cov_long <- target_KO_df_list[["cov"]][,1:11] %>% pivot_longer(-KO_ID, names_to="Sample", values_to="Norm_cov")
    counts_long <- target_KO_df_list[["counts"]][,1:11] %>% pivot_longer(-KO_ID, names_to="Sample", values_to="Num_copies")

    # subsetting to target KO
    target_KO_cov <- cov_long %>% filter(KO_ID == KO_term)
    target_KO_counts <- counts_long %>% filter(KO_ID == KO_term)

    # getting upper bounds for y-axes and transformation to be able to plot on same
       # little hack required here because Hadley doesn't like it: https://stackoverflow.com/questions/3099219/ggplot-with-2-y-axes-on-each-side-and-different-scales/3101876#3101876
       # this site helped me out: https://whatalnk.github.io/r-tips/ggplot2-secondary-y-axis.nb.html
    y_cov_limit <- plyr::round_any(max(target_KO_cov$Norm_cov), coverage_round_acc, ceiling)
    y_gene_count_limit <- plyr::round_any(max(target_KO_counts$Num_copies), gene_copy_round_acc, ceiling)

    # making label
    num_uniq_genes <- target_KO_df_list[["counts"]] %>% filter(KO_ID == KO_term) %>% pull(sum)
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

plot_overview_of_KOs_with_highest_coverages <- function( target_KO_df_list, number_to_plot=25, title="Gene coverages and relative copy numbers ('+' signs)") {

    # keeping the top number_to_plot KOs based on coverage, sorting by CoV and pulling target KOs
    target_KOs <- target_KO_df_list[["cov"]] %>% slice_max(mean, n = number_to_plot) %>% arrange(desc(CoV)) %>% pull(KO_ID)

    # printing out KO info for these:
    print(target_KO_df_list[['annots']] %>% filter(KO_ID %in% target_KOs) %>% arrange(match(KO_ID, target_KOs)))

    # initializing plotting list
    curr_plots <- vector("list", length(target_KOs))


    # iterating through generating each
    for ( num in 1:length(target_KOs) ) {
        curr_plots[[num]] <- plot_single_for_plotting_multiple_KOs( target_KOs[num], target_KO_df_list)
    }

    # making combined plot
    my_plot <- ggarrange(plotlist=curr_plots, common.legend=TRUE, legend = "bottom")

    my_plot <- annotate_figure(my_plot, top = text_grob(title, face="bold", size=14), left = text_grob("Norm. Cov. (CPM)", face="bold", size=11, rot=90))

    return(my_plot)
}

plot_overview_filtering_by_coverage <- function( target_KO_df_list, mean_coverage, title="Gene coverages and relative copy numbers") {

    # getting target KOs that have >= the specified mean coverage
    target_KOs <- target_KO_df_list[["cov"]] %>% filter(mean >= mean_coverage) %>% arrange(desc(CoV)) %>% pull(KO_ID)

    # printing out KO info for these:
    print(target_KO_df_list[['annots']] %>% filter(KO_ID %in% target_KOs) %>% arrange(match(KO_ID, target_KOs)))

    # initializing plotting list
    curr_plots <- vector("list", length(target_KOs))

    # iterating through generating each
    for ( num in 1:length(target_KOs) ) {
        my_plots[[num]] <- plot_single_for_plotting_multiple_KOs( target_KOs[num], target_KO_df_list)
    }

    # making combined plot
    my_plot <- ggarrange(plotlist=my_plots, common.legend=TRUE, legend = "bottom")

    my_plot <- annotate_figure(my_plot, top = text_grob(title, face="bold", size=14), left = text_grob("Norm. Cov. (CPM)", face="bold", size=11, rot=90))

    return(my_plot)
}

plot_overview_filtering_by_CoV <- function( target_KO_df_list, min_CoV, title="Gene coverages and relative copy numbers") {

    # getting target KOs that have >= the specified mean coverage
    target_KOs <- target_KO_df_list[["cov"]] %>% filter(CoV >= min_CoV) %>% arrange(desc(CoV)) %>% pull(KO_ID)

    # printing out KO info for these:
    print(target_KO_df_list[['annots']] %>% filter(KO_ID %in% target_KOs) %>% arrange(match(KO_ID, target_KOs)))

    # initializing plotting list
    curr_plots <- vector("list", length(target_KOs))

    # iterating through generating each
    for ( num in 1:length(target_KOs) ) {
        my_plots[[num]] <- plot_single_for_plotting_multiple_KOs( target_KOs[num], target_KO_df_list)
    }

    # making combined plot
    my_plot <- ggarrange(plotlist=my_plots, common.legend=TRUE, legend = "bottom")

    my_plot <- annotate_figure(my_plot, top = text_grob(title, face="bold", size=14), left = text_grob("Norm. Cov. (CPM)", face="bold", size=11, rot=90))

    return(my_plot)
}

# this function prints out a histogram and info of the coefficient of variation of the provided KO terms
plot_and_summarize_coverage_CoVs <- function( target_KO_df_list, plot_title="", min_mean_coverage=0) {

    # filtering by specified minimum mean coverage
    curr_cov_tab <- target_KO_df_list[["cov"]] %>% filter(mean >= min_mean_coverage)

    # getting mean coverages with and without un-annotated included
    curr_mean_coverages <- curr_cov_tab %>% filter(! is.na(mean)) %>% pull(mean)
    curr_mean_coverages_no_unannotated <- curr_cov_tab %>% filter(! is.na(mean)) %>% filter(! KO_ID == "Not annotated") %>% pull(mean)
    # getting number of entries
    n = length(curr_mean_coverages)
    n_no_unannotated = length(curr_mean_coverages_no_unannotated)

    cat("\n--------------------------------------------------------------------------------\n")
    cat("\n Summary of mean coverages:\n")
    cat("\n   Including those not annotated:\n")
    cat("    N =", n, "\n")
    print(summary(curr_mean_coverages))
    cat("\n   Excluding those not annotated:\n")
    cat("    N =", n_no_unannotated, "\n")
    print(summary(curr_mean_coverages_no_unannotated))


    # getting CoVs
    curr_CoVs <- curr_cov_tab %>% filter(! is.na(CoV)) %>% pull(CoV)
    # getting number of entries
    n = length(curr_CoVs)
    # printing out summary

    cat("\n--------------------------------------------------------------------------------\n")
    cat("\n Summary of CoVs:\n")
    cat("   N =", n, "\n")
    print(summary(curr_CoVs))

    cat("\n--------------------------------------------------------------------------------\n")
    cat("\n Top least varied:\n")
    print(curr_cov_tab %>% arrange(CoV) %>% head(n=10))

    cat("\n--------------------------------------------------------------------------------\n")
    cat("\n Top most varied:\n")
    print(curr_cov_tab %>% arrange(desc(CoV)) %>% head(n=10))

    # plotting histogram
    plot <- ggplot(curr_cov_tab %>% filter(! is.na(CoV)), aes(x=CoV)) +
        geom_histogram(fill=both_col, color="black", boundary = 1, binwidth=0.1) + theme_bw() +
        ggtitle(paste0(plot_title, " KO-term Coefficients of Variation")) + ylab(label="Num. KO terms") +
        scale_y_continuous(oob = function(x, limits) x) +
        scale_x_continuous(oob = function(x, limits) x)

    return(plot)
}

# This function will pull out the info for any given KO once we've made our overview of all genes tables, and add the link to the KO page
get_KO_info <- function( target_KO, KO_KEGG_tab = all_our_KOs_KEGG_tab ) {

    # bulding link
    curr_KO_link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ko:", target_KO)

    # getting info on KEGG term
    try(curr_info <- keggGet(target_KO), silent = TRUE)

    # checking if it was found at KEGG
    if ( ! exists("curr_info") ) {

        cat(paste0("\n  It seems '", target_KO, "' wasn't found at KEGG. It may have been removed, or may have never existed.\n\n  This should be the link if it were there if you wanna take a look:\n\n      ", curr_KO_link, "\n\n"))
        # couldn't figure out a better way to just report this and not give an error while exiting
        stop_no_error <- function() {
            opt <- options(show.error.messages = FALSE)
            on.exit(options(opt))
            stop()
        }

        stop_no_error()
    }

    # parsing some info
    if ( length(curr_info[[1]]$ENTRY) > 0 ) { curr_KO_ID <- curr_info[[1]]$ENTRY %>% as.vector } else { curr_KO_ID <- NA }
    if ( length(curr_info[[1]]$NAME) > 0 ) { curr_KO_name <- curr_info[[1]]$NAME %>% as.vector } else { curr_KO_name <- NA }
    if ( length(curr_info[[1]]$DEFINITION) > 0 ) { curr_KO_def <- curr_info[[1]]$DEFINITION %>% as.vector } else { curr_KO_def <- NA }

    # seeing if it's in our table
    if ( target_KO %in% (KO_KEGG_tab %>% pull(KO_ID)) ) {

      curr_in_data <- "Yes"

    } else {

      curr_in_data <- "No"

    }

    # reporting info
    cat("\n  KO ID          :  ", curr_KO_ID, "\n")
    cat("  KO name        :  ", curr_KO_name, "\n")
    cat("  KO definition  :  ", curr_KO_def, "\n")
    cat("  KO link        :  ", curr_KO_link, "\n")
    cat("  In our data?   :  ", curr_in_data, "\n\n")

}

#### beta-stuff ####
# starting with hierarchical clustering
all_KO_only_bray_dist <- vegdist(t(KO_tab[, 3:12]))
all_KO_only_bray_hclust <- hclust(all_KO_only_bray_dist, method="ward.D2")
all_KO_only_bray_dendro <- as.dendrogram(all_KO_only_bray_hclust, hang=0.15)

  # coloring labels based on source cave (and creating a vector of source cave)
col_vec <- vector()
source_vec <- vector()
for (label in labels(all_KO_only_bray_dendro)) {
    if (startsWith(label, "VAL")) {
        col_vec <- c(col_vec, VAL_col)
        source_vec <- c(source_vec, "VAL")
    } else {
        col_vec <- c(col_vec, YEL_col)
        source_vec <- c(source_vec, "YEL")
    }
}

labels_colors(all_KO_only_bray_dendro) <- col_vec

plot(all_KO_only_bray_dendro, main="Clustering based on normalized KO-term coverages", ylab="Bray-Curtis dissimilarity")

### PCoA view based on these normalized gene-coverages
KO_df <- data.frame(row.names=KO_tab$KO_ID, KO_tab[,3:12])
KO_norm_phy <- otu_table(KO_df, taxa_are_rows=TRUE)
sample_info_df <- data.frame(row.names=names(KO_tab)[3:12], sample_ID=names(KO_tab)[3:12], source=c(rep("VAL", 5), rep("YEL", 5)))
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


#### searching for specific things ####

## looking for DosSantos et al N-fixation genes
Dos_Santos_et_al_2012_N_fixation_KOs <- c("K02585", "K02586", "K02587", "K02588", "K02591", "K02592")

KO_tab %>% filter(KO_ID %in% Dos_Santos_et_al_2012_N_fixation_KOs)
# none present

#### overview of how we can use KEGGREST ####
#########
# With the KEGGREST package we loaded above, we can see the different KEGG databases accessible to us with:
listDatabases()

# And we can see which pathways are available with:
keggList("pathway") # which prints out a lot

# And we can search for those with "metabolism" in their name like so:
keggFind(database = "pathway", query = "metabolism")
  # Giving us an overview of those we might want to look into.


#### all KOs overview ####

all_our_KOs <- KO_tab %>% pull(KO_ID) %>% as.vector()

# took 35 minutes
all_our_KOs_KEGG_tab <- make_KO_summary_tab(all_our_KOs)

# adding "Not annotated" to the KO tab so our next function pulls that info too
all_our_KOs_KEGG_tab <- rbind(all_our_KOs_KEGG_tab, rep("Not annotated", 7))

# took 2 hours
all_our_KO_gene_df_list <- parse_samples_by_target_KO_tab(all_samples, all_our_KOs_KEGG_tab)

# making directory and writing out
dir.create("All-our-KO-terms")

write.table(all_our_KO_gene_df_list[["counts"]], "All-our-KO-terms/All-our-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(all_our_KO_gene_df_list[["cov"]], "All-our-KO-terms/All-our-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(all_our_KO_gene_df_list[["annots"]], "All-our-KO-terms/All-our-KO-info.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")

# looking at those that vary the least and most across all our recovered KO terms
plot_and_summarize_coverage_CoVs(all_our_KO_gene_df_list)

plot_and_summarize_coverage_CoVs(all_our_KO_gene_df_list, plot_title = "Mean coverage >= 200", min_mean_coverage=200)

plot_single_KO("K12546", all_our_KO_gene_df_list, coverage_round_acc=100)

plot_overview_of_KOs_with_highest_coverages(all_our_KO_gene_df_list)

#### nitrogen metabolism ####

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

N_KO_gene_df_list <- parse_samples_by_target_KO_tab(all_samples, Nitrogen_metabolism_KEGG_tab)

# making directory and writing out
dir.create("Nitrogen-focus")

write.table(N_KO_gene_df_list[["counts"]], "Nitrogen-focus/N-metabolism-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(N_KO_gene_df_list[["cov"]], "Nitrogen-focus/N-metabolism-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(N_KO_gene_df_list[["annots"]], "Nitrogen-focus/N-metabolism-KO-info.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")

## looking at how much coverage varies across the samples
plot_and_summarize_coverage_CoVs(N_KO_gene_df_list, "Nitrogen-related")

plot_and_summarize_coverage_CoVs(N_KO_gene_df_list, "Mean coverage >= 4000 Nitrogen-related", min_mean_coverage = 4000)

## plotting the most varied KO terms
plot_single_KO("K10946", N_KO_gene_df_list)
plot_single_KO("K10945", N_KO_gene_df_list)
plot_single_KO("K10944", N_KO_gene_df_list)

plot_overview_of_KOs_with_highest_coverages(N_KO_gene_df_list)


#### methane metabolism ####
# searching KEGG pathways for Methane
methane_pathway_info <- keggFind(database = "pathway", query = "methane")

# getting "Methane metabolism" pathway ID
methane_pathway_ID <- sub(names(methane_pathway_info), pattern = "path:", replacement = "")

# getting info for all KOs associated with "Methane metabolism" pathway
methane_pathway_KOs_info <- keggLink("ko", methane_pathway_ID)

# getting KO IDs alone
methane_pathway_KOs <- as.character(sub(methane_pathway_KOs_info, pattern="ko:", replacement=""))

length(methane_pathway_KOs) # 190

methane_metabolism_KEGG_tab <- make_KO_summary_tab(methane_pathway_KOs)

methane_KO_gene_df_list <- parse_samples_by_target_KO_tab(all_samples, methane_metabolism_KEGG_tab)

# making directory and writing out
dir.create("Methane-focus")

write.table(methane_KO_gene_df_list[["counts"]], "Methane-focus/Methane-metabolism-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(methane_KO_gene_df_list[["cov"]], "Methane-focus/Methane-metabolism-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(methane_KO_gene_df_list[["annots"]], "Methane-focus/Methane-metabolism-KO-info.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")

## looking at how much coverage varies across the samples
plot_and_summarize_coverage_CoVs(methane_KO_gene_df_list, "Methane-related")

plot_and_summarize_coverage_CoVs(methane_KO_gene_df_list, "Mean coverage >= 5000 Methane-related", min_mean_coverage = 5000)

## plotting the most varied KO terms
plot_single_KO("K17068", methane_KO_gene_df_list)
plot_single_KO("K00148", methane_KO_gene_df_list)
plot_single_KO("K05979", methane_KO_gene_df_list)

plot_overview_of_KOs_with_highest_coverages(methane_KO_gene_df_list)


#### sulfur metabolism ####
# searching KEGG pathways for Sulfur
sulfur_pathway_info <- keggFind(database = "pathway", query = "sulfur")

# setting "Sulfur metabolism" pathway ID
sulfur_pathway_ID <- "map04122"

# getting info for all KOs associated with "Sulfur metabolism" pathway
sulfur_pathway_KOs_info <- keggLink("ko", sulfur_pathway_ID)

# getting KO IDs alone
sulfur_pathway_KOs <- as.character(sub(sulfur_pathway_KOs_info, pattern="ko:", replacement=""))

length(sulfur_pathway_KOs) # 29

sulfur_metabolism_KEGG_tab <- make_KO_summary_tab(sulfur_pathway_KOs)

sulfur_KO_gene_df_list <- parse_samples_by_target_KO_tab(all_samples, sulfur_metabolism_KEGG_tab)

# making directory and writing out
dir.create("Sulfur-focus")

write.table(sulfur_KO_gene_df_list[["counts"]], "Sulfur-focus/Sulfur-metabolism-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(sulfur_KO_gene_df_list[["cov"]], "Sulfur-focus/Sulfur-metabolism-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(sulfur_KO_gene_df_list[["annots"]], "Sulfur-focus/Sulfur-metabolism-KO-info.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")

## looking at how much coverage varies across the samples
plot_and_summarize_coverage_CoVs(sulfur_KO_gene_df_list, "Sulfur-related")

plot_and_summarize_coverage_CoVs(sulfur_KO_gene_df_list, "Mean coverage >= 30000 Sulfur-related", min_mean_coverage = 30000)

## plotting the most varied KO terms
plot_single_KO("K21029", sulfur_KO_gene_df_list)
plot_single_KO("K04085", sulfur_KO_gene_df_list)
plot_single_KO("K21147", sulfur_KO_gene_df_list)

plot_overview_of_KOs_with_highest_coverages(sulfur_KO_gene_df_list)



#### carbon metabolism ####
# searching KEGG pathways for Carbon
carbon_pathway_info <- keggFind(database = "pathway", query = "carbon")

# setting "Carbon metabolism" pathway ID
carbon_pathway_ID <- "map01200"

# getting info for all KOs associated with "Carbon metabolism" pathway
carbon_pathway_KOs_info <- keggLink("ko", carbon_pathway_ID)

# getting KO IDs alone
carbon_pathway_KOs <- as.character(sub(carbon_pathway_KOs_info, pattern="ko:", replacement=""))

length(carbon_pathway_KOs) # 356

carbon_metabolism_KEGG_tab <- make_KO_summary_tab(carbon_pathway_KOs)

carbon_KO_gene_df_list <- parse_samples_by_target_KO_tab(all_samples, carbon_metabolism_KEGG_tab)

# making directory and writing out
dir.create("Carbon-focus")

write.table(carbon_KO_gene_df_list[["counts"]], "Carbon-focus/Carbon-metabolism-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(carbon_KO_gene_df_list[["cov"]], "Carbon-focus/Carbon-metabolism-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(carbon_KO_gene_df_list[["annots"]], "Carbon-focus/Carbon-metabolism-KO-info.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")

## looking at how much coverage varies across the samples
plot_and_summarize_coverage_CoVs(carbon_KO_gene_df_list, "Carbon-related")

plot_and_summarize_coverage_CoVs(carbon_KO_gene_df_list, "Mean coverage >= 1000 Carbon-related", min_mean_coverage = 1000)

## plotting the most varied KO terms
plot_single_KO("K00196", carbon_KO_gene_df_list)
plot_single_KO("K00245", carbon_KO_gene_df_list)
plot_single_KO("K17067", carbon_KO_gene_df_list)

plot_overview_of_KOs_with_highest_coverages(carbon_KO_gene_df_list)




#### diverse microbial metabolism metabolism ####
# searching KEGG pathways for Carbon
diverse_microbial_metabolism_pathway_info <- keggFind(database = "pathway", query = "diverse microbial metabolism")

# setting "Carbon metabolism" pathway ID
diverse_microbial_metabolism_pathway_ID <- "map01120"

# getting info for all KOs associated with "Carbon metabolism" pathway
diverse_microbial_metabolism_pathway_KOs_info <- keggLink("ko", diverse_microbial_metabolism_pathway_ID)

# getting KO IDs alone
diverse_microbial_metabolism_pathway_KOs <- as.character(sub(diverse_microbial_metabolism_pathway_KOs_info, pattern="ko:", replacement=""))

length(diverse_microbial_metabolism_pathway_KOs) # 1201

diverse_microbial_metabolism_metabolism_KEGG_tab <- make_KO_summary_tab(diverse_microbial_metabolism_pathway_KOs)

diverse_microbial_metabolism_KO_gene_df_list <- parse_samples_by_target_KO_tab(all_samples, diverse_microbial_metabolism_metabolism_KEGG_tab)

# making directory and writing out
dir.create("Diverse-microbial-metabolism-focus")

write.table(diverse_microbial_metabolism_KO_gene_df_list[["counts"]], "Diverse-microbial-metabolism-focus/Diverse-microbial-metabolism-KO-gene-copy-counts.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(diverse_microbial_metabolism_KO_gene_df_list[["cov"]], "Diverse-microbial-metabolism-focus/Diverse-microbial-metabolism-KO-coverages-CPM.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")
write.table(diverse_microbial_metabolism_KO_gene_df_list[["annots"]], "Diverse-microbial-metabolism-focus/Diverse-microbial-metabolism-KO-info.tsv", row.names=FALSE, sep="\t", quote=FALSE, na="NA")

## looking at how much coverage varies across the samples
plot_and_summarize_coverage_CoVs(diverse_microbial_metabolism_KO_gene_df_list, "Diverse metabolism-related")

plot_and_summarize_coverage_CoVs(diverse_microbial_metabolism_KO_gene_df_list, "Mean coverage >= 2000 Diverse-metabolism-related", min_mean_coverage = 2000)

## plotting the most varied KO terms
plot_single_KO("K01053", diverse_microbial_metabolism_KO_gene_df_list)
plot_single_KO("K17068", diverse_microbial_metabolism_KO_gene_df_list)
plot_single_KO("K11178", diverse_microbial_metabolism_KO_gene_df_list)

plot_overview_of_KOs_with_highest_coverages(diverse_microbial_metabolism_KO_gene_df_list)


save.image("BRAILLE-mg.RData")

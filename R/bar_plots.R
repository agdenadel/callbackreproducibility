library(ggplot2)
library(ggpattern)
library(reshape2)


# hacky way to get rid of notes in R CMD CHECK
# todo: switch to aes_string?
# no visible binding for global variable
base_algorithm = NULL
variable = NULL
value = NULL
is_knockoff = NULL
method = NULL
num_clusters = NULL
base_algorithm = NULL
method_clustering_algo = NULL




#' @title Computes standard clustering metrics for clustering assignments vs true group labels.
#'
#' @description Simulates single-cell RNAseq data for three groups with splat, clusters with Seurat and with KCC, and saves UMAP scatterplots of the clusters determined.
#'
#' @param true The true group labels.
#' @param clusters The clustering assignments.
#'
#' @returns A dataframe with the metrics.
#'
#' @name get_cluster_metrics
#' @export
get_cluster_metrics <- function(true, clusters) {
  metrics_name <- c("Completeness",
                    "Homogeneity",
                    "V-measure",
                    "ARI",
                    "Jaccard",
                    "NMI",
                    "FMI")
  
  metrics <- c(clevr::completeness,
               clevr::homogeneity,
               clevr::v_measure,
               clevr::adj_rand_index,
               clusteval::jaccard,
               aricode::NMI,
               clevr::fowlkes_mallows)
  
  results <- data.frame()
  

  
  for (i in 1:length(metrics)) {
    metric <- metrics[[i]]
    metric_name <- metrics_name[i]
    
    res <- metric(true, clusters)

    results[1, metric_name] <- res
    
  }
  
  return(results)
  
}



get_metrics_seurat <- function(seurat_obj, ground_truth_name) {
  method_names <- c("leiden_idents", # todo convert this to louvain idents
                    "knockoff_louvain_idents",
                    "leiden_idents",               
                    "knockoff_leiden_idents",
                    "scSHC_clusters",
                    "knockoff_hierarchical_clustering_idents",
                    "gap_k_means_idents",
                    "knockoff_k_means_idents"
                    )
  
  base_algorithm <- c("Louvain",
                      "Louvain",
                      "Leiden",               
                      "Leiden",
                      "Hierarchical",
                      "Hierarchical",
                      "K-means",
                      "K-means"
  )
  
  is_knockoff <- rep(c("Default", "KCC"), 4)
  
  true <- seurat_obj@meta.data[[ground_truth_name]]
  
  metrics_df <- data.frame()

  for (method_name in method_names) {
    clusters <- seurat_obj@meta.data[[method_name]]
    metrics <- get_cluster_metrics(true, clusters)
    

    metrics_df <- rbind(metrics_df, metrics)
  }
  
  metrics_df$method <- method_names
  metrics_df$base_algorithm <- base_algorithm
  metrics_df$is_knockoff <- is_knockoff
  
  return(metrics_df)
}


if (FALSE) {

setwd("~/Code/repos/PCKnockoffs_paper/tabula_muris/knockoff_clustering_output/")

tissue_seurat_files <- list.files(pattern = "rds")

metrics_df <- data.frame()

for (tissue_seurat_file in tissue_seurat_files) {
  tissue <- readRDS(tissue_seurat_file)
  
  tissue_name <- sub("_seurat_obj.rds", "", tissue_seurat_file)
  
  print(tissue_name)

  tissue_metrics_df <- get_metrics_seurat(tissue, "cell_ontology_class")
  tissue_metrics_df$tissue <- tissue_name

  metrics_df <- rbind(metrics_df, tissue_metrics_df)
  
}
}





#' @title Creates a facet plot (by tissue) displaying a clustering evaluation metric for KCC and default clustering for various clustering algorithlms.
#'
#' @description Creates a facet plot (by tissue) displaying a clustering evaluation metric for KCC and default clustering for Louvain, Leiden, Hierarchical, and K-means clustering algorithms.
#'
#' @param metrics_df The dataframe containing the metrics.
#' @param statistic The clustering evaluation metric to plot.
#' @param y_label A string with the name of the metric to be the y-axis label.
#'
#' @returns A ggplot2 object with the facet bar plot.
#'
#' @name tissue_facet_bar_plot
#' @export

library(ggpattern)
tissue_facet_bar_plot <- function(cluster_metrics_df, statistic, y_label) {

  # re-order factor levels
  cluster_metrics_df$method_clustering_algo <- factor(cluster_metrics_df$method_clustering_algo, levels=c("Louvain", "Leiden", "Hierarchical", "K-means", "CHOIR"))

  
  
  pattern_bar_plot <- ggplot2::ggplot(cluster_metrics_df) + 
    ggpattern::geom_bar_pattern(ggplot2::aes(x=method_clustering_algo, y=!! rlang::sym(statistic), fill=method_clustering_algo, pattern=is_knockoff, alpha=is_knockoff), stat='identity', position="dodge",
                     color = "black", 
                     pattern_fill = "black",
                     pattern_angle = 45,
                     pattern_density = 0.05,
                     pattern_spacing = 0.05,
                     pattern_key_scale_factor = 0.6) +
    ggpattern::scale_pattern_manual(values = c(KCC = "stripe", Default = "none")) +
    ggplot2::facet_wrap(~ tissue_name, ncol=5) +
    ggplot2::xlab("Tabula Muris Tissues") + 
    ggplot2::ylab(y_label) + 
    #labs(fill='Clustering Algorithm\n(KCC percent improved)', pattern='Selection Method') + 
    ggplot2::labs(fill='Clustering Algorithm', pattern='Selection Method') + 
    ggplot2::scale_fill_brewer(palette = "Set1", 
                      guide = ggplot2::guide_legend(override.aes = list(pattern = "none"))) + #,
                      #labels = c(paste0("Louvain (", louvain_percent_improved, ")"),
                      #           paste0("Leiden (", leiden_percent_improved, ")"),
                      #           paste0("Hierarchical (", hierarchical_percent_improved, ")"),
                      #           paste0("K-means (", k_means_percent_improved, ")"))) +
    ggplot2::scale_alpha_manual(values = c(KCC = 1.0, Default = 0.5), guide=FALSE) +
    ggplot2::guides(pattern = ggplot2::guide_legend(override.aes = list(fill = "white", pattern_spacing = 0.015))) + 
    ggplot2::theme_bw() +
    ggplot2::theme(#axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      #axis.text.x = ggplot2::element_text(angle = 90, margin = margin(t = 25)),
      axis.text.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size=13),
      axis.title = ggplot2::element_text(size=18),
      strip.text = ggplot2::element_text(size=13),
      legend.text = ggplot2::element_text(size=13),
      legend.title = ggplot2::element_text(size=13)
    )
  
  return(pattern_bar_plot)
}



metrics_bar_plot <- function(metrics_df, algorithm) {
  subsetted_df <- subset(metrics_df, base_algorithm == algorithm)
  
  text_size <- 24
  
  pattern_bar_plot <- ggplot2::ggplot(subsetted_df) + 
    #geom_bar(aes(x=variable, y=value, fill=base_algorithm), stat='identity', position="dodge", color="black") + 
    ggpattern::geom_bar_pattern(ggplot2::aes(x=variable, y=value, fill=base_algorithm, pattern=is_knockoff, alpha=is_knockoff), stat='identity', position="dodge",
                     color = "black", 
                     pattern_fill = "black",
                     pattern_angle = 45,
                     pattern_density = 0.05,
                     pattern_spacing = 0.025,
                     pattern_key_scale_factor = 0.6) +
    ggpattern::scale_pattern_manual(values = c(KCC = "stripe", Default = "none")) +
    ggplot2::scale_alpha_manual(values = c(KCC = 1.0, Default = 0.5), guide=FALSE) +
    #xlab("Clustering Metrics") + 
    #ylab("Value") + 
    ggplot2::labs(fill='Method') + #, pattern='Selection Method') + 
    ggplot2::theme_bw() +
    ggplot2::theme(#axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, margin = ggplot2::margin(t = 50)),
      axis.text = ggplot2::element_text(size=text_size),
      axis.title=ggplot2::element_blank(),
      #axis.title = ggplot2::element_text(size=15),
      strip.text = ggplot2::element_text(size=text_size),
      legend.text = ggplot2::element_text(size=text_size),
      #legend.background = ggplot2::element_rect(fill="grey", colour = "black"),
      legend.position = c(0.7, 0.85),
      legend.title=ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_brewer(palette = "Set1", limits = c("Louvain", "Leiden", "Hierarchical", "K-means", "CHOIR")) + #, 
                      #guide = guide_legend(override.aes = list(pattern = "none"))) +
    ggplot2::guides(pattern = ggplot2::guide_legend(override.aes = list(fill = "white"))) + 
    ggplot2::guides(fill="none")
  return(pattern_bar_plot)
}






num_clusters_bar_plot <- function(num_clusters_df) {
  text_size <- 24
  
  p <- ggplot2::ggplot(num_clusters_df) +
    #geom_bar(aes(x=method, y=num_clusters, fill=is_knockoff), stat='identity', position="dodge", color = "black") +
    ggpattern::geom_bar_pattern(ggplot2::aes(x=method, y=num_clusters, fill=method, pattern=is_knockoff, alpha=is_knockoff), stat='identity', position="dodge",
                     color = "black", 
                     pattern_fill = "black",
                     pattern_angle = 45,
                     pattern_density = 0.05,
                     pattern_spacing = 0.025,
                     pattern_key_scale_factor = 0.6) +
    ggpattern::scale_pattern_manual(values = c(Knockoff = "stripe", Default = "none", FACS = "none")) + 
    ggplot2::scale_alpha_manual(values = c(KCC = 1.0, Default = 0.5), guide=FALSE) +
  
    ggplot2::xlab("Metric") +
    ggplot2::ylab("Number of Clusters") + 
    ggplot2::theme_bw() + 
    #geom_abline(slope=0, intercept=10,  col = "black",lty=2) +
    ggplot2::theme(#axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, margin = ggplot2::margin(t = 45)),
      axis.text = ggplot2::element_text(size=text_size),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size=text_size),
      #axis.title = ggplot2::element_text(size=15),
      strip.text = ggplot2::element_text(size=text_size),
      legend.text = ggplot2::element_text(size=text_size + 8),
      legend.position = c(0.85, 0.9),
      #legend.background = ggplot2::element_rect(fill="grey", colour = "black"),
      legend.title=ggplot2::element_blank()
    ) +
    ggplot2::guides(fill=ggplot2::guide_legend(title="Method", override.aes = list(pattern = "none"))) + 
    ggplot2::scale_fill_brewer(palette = "Set1", limits = c("Louvain", "Leiden", "Hierarchical", "K-means", "FACS")) +
    ggplot2::guides(pattern = ggplot2::guide_legend(override.aes = list(fill = "white", size=24)), fill = "none")
  return(p)
}









get_num_clusters_df <- function(seurat_obj) {
  
  louvain_kcc_num_clusters <- length(levels(as.factor(seurat_obj@meta.data$knockoff_louvain_idents)))
  louvain_num_clusters <- length(levels(as.factor(seurat_obj@meta.data$seurat_clusters)))
  
  leiden_kcc_num_clusters <- length(levels(as.factor(seurat_obj@meta.data$knockoff_leiden_idents)))
  leiden_num_clusters <- length(levels(as.factor(seurat_obj@meta.data$leiden_idents)))
  
  hier_kcc_num_clusters <- length(levels(as.factor(seurat_obj@meta.data$knockoff_hierarchical_clustering_idents)))
  hier_num_clusters <- length(levels(as.factor(seurat_obj@meta.data$scSHC_clusters)))
  
  kmeans_kcc_num_clusters <- length(levels(as.factor(seurat_obj@meta.data$knockoff_k_means_idents)))
  kmeans_num_clusters <- length(levels(as.factor(seurat_obj@meta.data$gap_k_means_idents)))
  
  
  num_clusters <- c(10,
                    louvain_num_clusters, louvain_kcc_num_clusters,
                    leiden_num_clusters, leiden_kcc_num_clusters,
                    hier_num_clusters, hier_kcc_num_clusters,
                    kmeans_num_clusters, kmeans_kcc_num_clusters
  )
  
  method <- c("FACS",
               "Louvain", "Louvain",
               "Leiden", "Leiden",
               "Hierarchical", "Hierarchical",
               "K-means", "K-means"
  )
  
  is_knockoff <- c("Default",
                   "Default", "Knockoff",
                   "Default", "Knockoff",
                   "Default", "Knockoff",
                   "Default", "Knockoff"
  )
  
  
  
  num_clusters_df <- data.frame(method, num_clusters, is_knockoff)
  
  num_clusters_df$method <- factor(num_clusters_df$method, levels=c("FACS", "Louvain", "Leiden", "Hierarchical", "K-means"))
  num_clusters_df$is_knockoff <- factor(num_clusters_df$is_knockoff, levels=c("Default", "Knockoff", "FACS"))
  
  return(num_clusters_df)

}


if (FALSE) {

ari_bar_plot <- tissue_facet_bar_plot(metrics_df, "adj_rand_index", "Adjusted Rand Index (ARI)")
ari_bar_plot

v_measure_bar_plot <- tissue_facet_bar_plot(metrics_df, "v_measure", "v-measure")
v_measure_bar_plot

homogeneity_bar_plot <- tissue_facet_bar_plot(metrics_df, "homogeneity", "Homogeneity")
homogeneity_bar_plot

completeness_bar_plot <- tissue_facet_bar_plot(metrics_df, "completeness", "Completeness")
completeness_bar_plot

v_measure_bar_plot <- tissue_facet_bar_plot(metrics_df, "v_measure", "v-measure")
v_measure_bar_plot

#nmi_bar_plot <- tissue_facet_bar_plot(metrics_df, "nmi_results", "Normalized Mutual Information (NMI)")
#nmi_bar_plot

fm_bar_plot <- tissue_facet_bar_plot(metrics_df, "Fowlkes_Mallows", "Fowlkes Mallows Index")
fm_bar_plot

jaccard_bar_plot <- tissue_facet_bar_plot(metrics_df, "jaccard", "Jaccard Index")
jaccard_bar_plot

}








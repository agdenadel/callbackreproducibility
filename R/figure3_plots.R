#todo library(tidyverse) check if i need to add this to DESCRIPTION


# hacky way to get rid of notes in R CMD CHECK
# todo: switch to aes_string?
# no visible binding for global variable
avg_log2FC = NULL
Name = NULL
color = NULL
log10pval = NULL
cluster = NULL
p_val_adj = NULL

#' @title Loads a Tabula Muris tissue seurat object from an RDS file.
#'
#' @description Loads a Tabula Muris tissue seurat object from an RDS file.
#'
#' @param tissue_name The name of the Tabula Muris tissue of interest.
#' @param directory The directory that the Tabula Muris tissue rds file is in.
#'
#' @returns A seurat object for the Tabula Muris tissue of interest.
#'
#' @name get_tissue
#' @export
get_tissue <- function(tissue_name, directory) {
  tissue <- readRDS(paste0(directory, tissue_name, "_seurat_obj.rds"))
  
  tissue_ch_index_seurat_obj <- readRDS(paste0(directory, "ch_index_seurat_objs/", tissue_name, "_seurat_obj_k_means_by_ch_index_.rds"))
  tissue@meta.data$ch_k_means_idents <- tissue_ch_index_seurat_obj@meta.data$ch_k_means_idents
  
  tissue <- Seurat::ScaleData(tissue, features = rownames(tissue))
  return(tissue)
}



#' @title Saves umap scatterplots of the cell ontology labels, Seurat default clusters, and callback clusters for the Tabula Muris tissue of interest.
#'
#' @description Saves umap scatterplots of the cell ontology labels, Seurat default clusters, and callback clusters for the Tabula Muris tissue of interest.
#'
#' @param tissue The Seurat object for the Tabula Muris tissue of interest.
#' @param tissue_name The name of the Tabula Muris tissue of interest.
#' @param legend_pos The location to place the legend in the cell ontology plot.
#'
#' @name save_scatter_plots
#' @export
fig3_scatter_plots <- function(tissue, tissue_name, legend_pos=c(0.6, 0.2)) {  
  louvain_default <- custom_scatter(tissue, "umap", group_by = "seurat_clusters", x_title = "UMAP 1", y_title = "UMAP 2", pt.size = 2, label=FALSE) + Seurat::NoLegend()
  louvain_callback <- custom_scatter(tissue, "umap", group_by = "callback_idents", x_title = "UMAP 1", y_title = "UMAP 2", pt.size = 2) + Seurat::NoLegend() 
  cell_ontology <- custom_scatter(tissue, "umap", group_by = "cell_ontology_class", x_title = "UMAP 1", y_title = "UMAP 2", pt.size = 2) + 
    ggplot2::theme(legend.position = legend_pos,
                   legend.text = ggplot2::element_text(size=20)) + 
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=6), ncol = 1)) + 
    ggplot2::scale_colour_discrete(na.translate = F)

  column_label_1 <- patchwork::wrap_elements(panel = grid::textGrob('Cell Ontology', gp = grid::gpar(fontsize = 64)))
  column_label_2 <- patchwork::wrap_elements(panel = grid::textGrob('Seurat Default', gp = grid::gpar(fontsize = 64)))
  column_label_3 <- patchwork::wrap_elements(panel = grid::textGrob('callback', gp = grid::gpar(fontsize = 64, fontfamily = "Courier")))


  umap_grid <- column_label_1 + column_label_2 + column_label_3 +
    cell_ontology + louvain_default + louvain_callback +
    patchwork::plot_layout(widths = c(5, 5, 5),
                           heights = c(1,3))

  return(umap_grid)
}



#' @title Saves RNA expression heatmaps with group labels for the Seurat default clusters and callback clusters for the Tabula Muris tissue of interest.
#'
#' @description Saves RNA expression heatmaps with group labels for the Seurat default clusters and callback clusters for the Tabula Muris tissue of interest.
#'
#' @param tissue The Seurat object for the Tabula Muris tissue of interest.
#' @param tissue_name The name of the Tabula Muris tissue of interest.
#'
#' @name save_heatmaps
#' @export
#' @importFrom magrittr %>%
fig3_heatmaps <- function(tissue, tissue_name){
  # re-index clusters to be one-based
  tissue@meta.data$callback_idents <- as.factor(as.numeric(tissue@meta.data$callback_idents))
  
  Seurat::Idents(tissue) <- tissue@meta.data$callback_idents
  
  markers <- Seurat::FindAllMarkers(tissue, only.pos = TRUE, logfc.threshold = 0.25)
  
  markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup() -> top10
  
  clusters_picked <- levels(Seurat::Idents(tissue))
  levels(tissue) <- clusters_picked
  callback_heatmap <- Seurat::DoHeatmap(tissue, features = top10$gene,
                                        label = FALSE,
                                        size = 8,
                                        angle = 0,
                                        hjust = 0.5,
                                        disp.max = 2.5,
                                        disp.min = -2.5,
                                        group.bar.height = 0.06) +
    ggplot2::scale_fill_gradientn(colors = c("grey", "white", "blue"), na.value = "white", limits = c(-2.5,2.5)) +
    #NoLegend() + 
    ggplot2::xlab("Cells") + 
    ggplot2::ylab("Marker Genes") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(size = 48),
                   legend.text= ggplot2::element_text(size=32),
                   legend.title = ggplot2::element_text(size=36),
                   legend.key.width = ggplot2::unit(1,"cm"),
                   legend.key.height = ggplot2::unit(2, "cm"),
                   plot.title = ggplot2::element_text(size = 48, hjust = 0.5)) +
                   ggplot2::guides(colour="none")
  
  
  
  # re-index clusters to be one-based
  tissue@meta.data$seurat_clusters <- as.factor(as.numeric(tissue@meta.data$seurat_clusters))
  Seurat::Idents(tissue) <- tissue@meta.data$seurat_clusters
  
  markers <- Seurat::FindAllMarkers(tissue, only.pos = TRUE)
  #markers <- FindAllMarkers(krasnow_umi_seurat, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1.0, features = VariableFeatures(krasnow_umi_seurat))
  
  markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup() -> top10  
  
  clusters_picked <- levels(Seurat::Idents(tissue))
  levels(tissue) <- clusters_picked
  default_heatmap <- Seurat::DoHeatmap(tissue, features = top10$gene, 
                                       label = FALSE,
                                       size = 8,
                                       angle = 0,
                                       hjust = 0.5,
                                       disp.max = 2.5, 
                                       disp.min = -2.5,
                                       group.bar.height = 0.06) +
    ggplot2::scale_fill_gradientn(colors = c("grey", "white", "blue"), na.value = "white", limits = c(-2.5,2.5)) +
    #NoLegend() +
    ggplot2::xlab("Cells") + 
    ggplot2::ylab("Marker Genes") +
    #ggtitle("Default Clusters") + 
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.title = ggplot2::element_text(size = 48),
                   legend.text= ggplot2::element_text(size=32),
                   legend.title = ggplot2::element_text(size=36),
                   legend.key.width = ggplot2::unit(1,"cm"),
                   legend.key.height = ggplot2::unit(2, "cm"),
                   plot.title = ggplot2::element_text(size = 48, hjust = 0.5)) +
    ggplot2::guides(colour="none")
  
  
  return(list("callback_heatmap" = callback_heatmap, "default_heatmap" = default_heatmap))
}




#' @title Saves volcano plots for particular Seurat default clusters and particular callback clusters for the Tabula Muris tissue of interest.
#'
#' @description Saves volcano plots for particular Seurat default clusters and particular callback clusters for the Tabula Muris tissue of interest.
#'
#' @param tissue The Seurat object for the Tabula Muris tissue of interest.
#' @param tissue_name The name of the Tabula Muris tissue of interest.
#' @param callback_cluster1 The name of the first callback cluster of interest.
#' @param callback_cluster2 The name of the second callback cluster of interest.
#' @param default_cluster1 The name of the first default Seurat cluster of interest.
#' @param default_cluster2 The name of the second default Seurat cluster of interest.
#' @param ymax The maximum of the y-axis.
#' @param y_increment The increment of the y-axis tick marks.
#' @param genes_to_label_left The names of the genes to label on the left hand side of plot.
#' @param genes_to_label_right The names of the genes to label on the right hand side of plot.
#'
#' @name save_heatmaps
#' @export
fig3_volcano_plots <- function(tissue, tissue_name,
                               callback_cluster1, callback_cluster2,
                               default_cluster1, default_cluster2,
                               ymax=150,
                               y_increment = 10,
                               genes_to_label_left=c(),
                               genes_to_label_right=c()) {
  
  volcano_plot <- function(seurat_obj, markers, cluster1, cluster2, title, logfc_thresh=1.0,
                           genes_to_label_left=c(),
                           genes_to_label_right=c()) {
    markers$log10pval <- -log10(markers$p_val)
    
    markers$Name <- rownames(markers)
    
    p_val_thresh <- -log10(0.05 / dim(seurat_obj)[1])
    
    
    markers$log10pval[markers$log10pval > ymax] <- ymax 
    
    markers$color <- "grey"
    
    markers[(markers$avg_log2FC > 1) & (markers$log10pval > p_val_thresh),]$color <- "red"
    markers[(markers$avg_log2FC < -1) & (markers$log10pval > p_val_thresh),]$color <- "blue"
    
    
    p <- ggplot2::ggplot(markers, ggplot2::aes(x=avg_log2FC, y=log10pval, color = color)) + 
      ggplot2::geom_point() +
      ggplot2::scale_colour_identity() +
      ggrepel::geom_label_repel(data = markers %>% dplyr::filter(Name %in% genes_to_label_right), ggplot2::aes(label = Name),
                                min.segment.length = 0,
                                box.padding = 1.5,
                                point.size = 2,
                                size = 10,
                                force = 12,
                                xlim  = c(14,19),
                                ylim  = c(50,160),
                                hjust=0,
                                direction = "y",
                                max.overlaps = Inf) + # right side isn't showing half of the labels
      ggrepel::geom_label_repel(data = markers %>% dplyr::filter(Name %in% genes_to_label_left), ggplot2::aes(label = Name),
                                min.segment.length = 0,
                                box.padding = 1.5,
                                point.size = 2,
                                size = 10,
                                force = 12,
                                seed = 123,
                                xlim  = c(-14,-19),
                                ylim  = c(50,160),
                                hjust=1,
                                direction = "y",
                                max.overlaps = Inf) + # left side isn't showing half of the labels
      #gghighlight::gghighlight(log10pval > p_val_thresh) + 
      #gghighlight::gghighlight(abs(avg_log2FC) > logfc_thresh) + 
      ggplot2::xlim(-20,20) + 
      ggplot2::xlab("Average Log2-Fold Change") +
      ggplot2::ggtitle(title) + 
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     axis.title = ggplot2::element_text(size=40),
                     axis.text = ggplot2::element_text(size=32),
                     title = ggplot2::element_text(size=44),
                     legend.title = ggplot2::element_blank(),
                     legend.background = ggplot2::element_blank(),
                     legend.box.background = ggplot2::element_rect(colour = "black"), 
                     legend.text = ggplot2::element_text(size=32),
                     legend.position = c(0.8, 0.6),
                     legend.key.size = ggplot2::unit(3, "cm")) +
      ggplot2::geom_hline(ggplot2::aes(yintercept=p_val_thresh), size = 1, linetype = 'dashed', show.legend = TRUE) +
      ggplot2::geom_vline(ggplot2::aes(xintercept=logfc_thresh), linetype = 'dashed', size = 1) +
      ggplot2::geom_vline(ggplot2::aes(xintercept=-logfc_thresh), linetype = 'dashed', size = 1) +
      ggplot2::scale_linetype_manual(values=c("dashed")) +
      #scale_color_manual(values=c("blue", "red")) +
      # annotate p-value threshold
      ggplot2::annotate("segment", x = 13, xend = 15, y = 18, yend = p_val_thresh, colour = "black", linetype = "dashed", size = 2) +
      ggplot2::annotate("label", x = 13, y = 25, label = "Adj. P-value = 0.05", size = 12) +
      # annotate logFC thresholds
      ggplot2::annotate("segment", x = 11.6, xend = 1.0, y = 40, yend = 60, colour = "black", linetype = "dashed", size = 2) +
      ggplot2::annotate("label", x = 13, y = 40, label = "Avg. Log2-FC = 1.0", size = 12) +
      ggplot2::annotate("segment", x = -11.6, xend = -1.0, y = 40, yend = 60, colour = "black", linetype = "dashed", size = 2) +
      ggplot2::annotate("label", x = -13, y = 40, label = "Avg. Log2-FC = -1.0", size = 12) +
      #geom_text_repel(data=subset(markers, abs(avg_log2FC) > logfc_thresh & p_val_adj < 0.05), aes(label = Name), size = 8)
      ggplot2::scale_y_continuous(name = "-log10 P-value", limits = c(0, ymax))#,
                       #breaks = seq(from = 0, to = ymax, by = y_increment)
                       #labels = c(seq(from = 0, to = ymax - 100, by = y_increment), paste0('\u2265', ymax)))
    return(p)
  }
  
  
  
  
  
  # re-index clusters to be one-based
  tissue@meta.data$callback_idents <- as.factor(as.numeric(tissue@meta.data$callback_idents))
  
  Seurat::Idents(tissue) <- tissue@meta.data$callback_idents
  #DimPlot(tissue)
  
  callback_markers <- Seurat::FindMarkers(tissue,
    ident.1 = callback_cluster1,
    ident.2 = callback_cluster2,
    logfc.threshold = 0.0,
    )
  # highlighted markers
  subset(callback_markers, abs(avg_log2FC) > 1.0 & p_val_adj < 0.05)
  
  callback_title <- paste0("callback: Cluster ", callback_cluster1, " vs Cluster ", callback_cluster2)
  volcano_callback <- volcano_plot(tissue, callback_markers, callback_cluster1, callback_cluster2, callback_title,
                              genes_to_label_left = genes_to_label_left, genes_to_label_right = genes_to_label_right)
  
  # re-index clusters to be one-based
  tissue@meta.data$seurat_clusters <- as.factor(as.numeric(tissue@meta.data$seurat_clusters))
  Seurat::Idents(tissue) <- tissue@meta.data$seurat_clusters
  
  default_markers <- Seurat::FindMarkers(tissue,
    ident.1 = default_cluster1,
    ident.2 = default_cluster2,
    logfc.threshold = 0.0,
    )
  # highlighted markers
  subset(default_markers, abs(avg_log2FC) > 1.0 & p_val_adj < 0.05)
  
  default_title <- paste0("Default: Cluster ", default_cluster1, " vs Cluster ", default_cluster2)
  volcano_default <- volcano_plot(tissue, default_markers, default_cluster1, default_cluster2, default_title)
  

  return(list("volcano_callback" = volcano_callback, "volcano_default" = volcano_default))

}








#' @title Plots P-values and LFC for over-clustered clusters.
#'
#' @description Plots test statistics for clusters that are (putatively) of the same type when each cluster is compared to a third cluster.
#'
#' @param tissue The Seurat object for the Tabula Muris tissue of interest.
#' @param default_cluster1a The first over-clustered cluster of interest.
#' @param default_cluster1b The second over-clustered cluster of interest.
#' @param default_cluster2 The cluster to compare the other two clusters against.

#'
#' @name fig3_p_value_scatterplot
#' @export
fig3_p_value_scatterplot <- function(tissue,
                                     default_cluster1a,
                                     default_cluster1b,
                                     default_cluster2) {
  
  # re-index clusters to be one-based
  tissue@meta.data$seurat_clusters <- as.factor(as.numeric(tissue@meta.data$seurat_clusters))
  Seurat::Idents(tissue) <- tissue@meta.data$seurat_clusters
  
  default_markers_a <- Seurat::FindMarkers(tissue,
                                         ident.1 = default_cluster1a,
                                         ident.2 = default_cluster2,
                                         logfc.threshold = 0.0,
                                         min.pct = 0)
  
  default_markers_b <- Seurat::FindMarkers(tissue,
                                           ident.1 = default_cluster1b,
                                           ident.2 = default_cluster2,
                                           logfc.threshold = 0.0,
                                           min.pct = 0)

  # sort genes alphabetically
  default_markers_a <- default_markers_a[ order(row.names(default_markers_a)), ]
  default_markers_b <- default_markers_b[ order(row.names(default_markers_b)), ]
  
  log_p_val_a <- -log(default_markers_a$p_val)
  log_p_val_b <- -log(default_markers_b$p_val)
  
  
  df <- data.frame(log_p_val_a, log_p_val_b)

  
  r_p_val <- cor(log_p_val_a, log_p_val_b)

  large_text_size <- 32
  small_text_size <- 24
  
  correlation_text_size <- 16
  
  p_value_scatterplot <- ggplot2::ggplot(df, ggplot2::aes(x = log_p_val_a, y = log_p_val_b)) + 
    ggplot2::geom_point() + 
    ggplot2::theme_bw() +
    xlim(0, 175) +
    ylim(0, 175) +
    ggplot2::xlab(paste0("Default Clusters ", default_cluster1a, " vs ", default_cluster2, "\n-log P-values")) +
    ggplot2::ylab(paste0("Default Clusters ", default_cluster1b, " vs ", default_cluster2, "\n-log P-values")) +
    ggplot2::annotate("text", x = 25, y = 125, label = paste0("r = ", round(r_p_val, 3)), size = correlation_text_size) + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = small_text_size),
                   axis.title = ggplot2::element_text(size = large_text_size))
  
  
  
  lfc_a <- default_markers_a$avg_log2FC
  lfc_b <- default_markers_b$avg_log2FC
  
  df <- data.frame(lfc_a, lfc_b)
  
  
  r_lfc <- cor(lfc_a, lfc_b)
  
  lfc_scatterplot <- ggplot2::ggplot(df, ggplot2::aes(x = lfc_a, y = lfc_b)) + 
    ggplot2::geom_point() + 
    ggplot2::theme_bw() +
    ggplot2::xlab(paste0("Default Clusters ", default_cluster1a, " vs ", default_cluster2, "\nLog-Fold Change")) +
    ggplot2::ylab(paste0("Default Clusters ", default_cluster1b, " vs ", default_cluster2,"\nLog-Fold Change")) +
    ggplot2::annotate("text", x = -10, y = 5, label = paste0("r = ", round(r_lfc, 3)), size = correlation_text_size) + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = small_text_size),
                   axis.title = ggplot2::element_text(size = large_text_size))
 
 return(list("p_value_scatterplot" = p_value_scatterplot, "lfc_scatterplot" = lfc_scatterplot))
}


#' @title Scatterplot with a custom theme.
#'
#' @description Given a Seurat object and plotting parameters, returns a scatterplot with a custom theme based on those parameters.
#'
#' @param seurat_obj The Seurat object that data will be taken from.
#' @param group_name The name of the metadata column to group by.
#' @param title The plot title.
#' @param add_legend Whether or not to include a legend.
#' @returns Returns a ggplot2 object containing the scatterplot.

#' @name one_umap
#' @export
one_umap <- function(seurat_obj, group_name, title, add_legend=FALSE) {
  p <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = group_name, pt.size = 1) + 
    ggplot2::theme(axis.line = ggplot2::element_line(size=3),
                   axis.title.x = ggplot2::element_text(size = 40, vjust=0), 
                   axis.title.y = ggplot2::element_text(size = 40), 
                   axis.ticks = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 24),
                   plot.title = ggplot2::element_text(size = 48, family = "Courier", face = "plain")) +
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2") +
    ggplot2::ggtitle(title) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size=8), nrow = 15))


  if (!add_legend) {
    p <- p + Seurat::NoLegend()
  }
  
  return(p)
}



#' @title UMAP grid for the Cell Ontology Classes, callback clusters, sc-SHC
#' clusters, and CHOIR clusters
#'
#' @description Given a Seurat object and plotting parameters, returns a grid
#' of UMAP scatterplots for the Cell Ontology Classes, callback clusters, 
#' sc-SHC clusters, and CHOIR clusters
#'
#' @param tissue The Seurat object that data will be taken from
#' @param tissue_name The name of the tissue being plotted.
#' @returns Returns a ggplot2 object containing the grid of UMAP plots.

#' @name get_umap
#' @export
get_umap <- function(tissue, tissue_name) {
    
  umap_fig <- one_umap(tissue, "cell_ontology_class", "Cell Ontology Class", add_legend=TRUE) + 
    one_umap(tissue, "callback_idents", "callback", add_legend=TRUE) + 
    one_umap(tissue, "scSHC_clusters", "sc-SHC", add_legend=TRUE) + 
    one_umap(tissue, "CHOIR_clusters_0.05", "CHOIR", add_legend=TRUE) +
      patchwork::plot_layout(widths = c(1, 1),
                             heights = c(1,1)) + 
    plot_annotation(title = tissue_name,
                    theme = theme(plot.title = element_text(size = 64, hjust = 0.5, vjust = 1.0)))

  return(umap_fig)
}

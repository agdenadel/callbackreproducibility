
#' @title Scatterplot with a custom theme.
#'
#' @description Given a Seurat object and plotting parameters, returns a scatterplot with a custom theme based on those parameters.
#'
#' @param seurat_obj The Seurat object that data will be taken from
#' @param reduction The dimensionality reduction type to plot (i.e. 'umap' or 'tsne')
#' @param group_by The groupings to color by.
#' @param x_title The title of the x-axis.
#' @param y_title The title of the y-axis.
#' @param pt.size The size of points in the scatterplot.
#' @param label Whether or not the groups should be labelled.
#' @returns Returns a ggplot2 object containing the scatterplot.

#' @name custom_scatter
#' @export
custom_scatter <- function(seurat_obj, reduction, group_by, x_title, y_title, pt.size, label=FALSE) {
  p <- Seurat::DimPlot(seurat_obj, reduction = reduction, group.by = group_by, pt.size = pt.size, label=label) + 
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black", linewidth=3),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = 48, vjust=0),
                   axis.title.y = ggplot2::element_text(size = 48), 
                   axis.ticks = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 32),
                   plot.title = ggplot2::element_blank())  +
    #theme(legend.position = "none") + 
    ggplot2::xlab(x_title) + 
    ggplot2::ylab(y_title)
  
  return(p)
}


#' @title Barplot faceted by tissue with a custom theme.
#'
#' @description Given a Seurat object and plotting parameters, returns a scatterplot with a custom theme based on those parameters.
#'
#' @param cluster_metrics_df The dataframe containing the clustering metrics.
#' @param statistic The clustering metric to plot.
#' @param y_label The name of the clustering metric
#' @returns Returns a ggplot2 object containing the facet plot.

#' @name clustering_metrics_facet_plot
#' @export
clustering_metrics_facet_plot <- function (cluster_metrics_df, statistic, y_label) 
{
  # remove underscores from tissue names
  cluster_metrics_df$tissue_name <- gsub("_", "\n", cluster_metrics_df$tissue_name)
  
  # re-order factor levels
  cluster_metrics_df$method <- factor(cluster_metrics_df$method, levels = c("callback", "sc-SHC", "CHOIR"))
  
  
  small_text_size <- 16
  large_text_size <- 22
  
  bar_plot <- ggplot2::ggplot(cluster_metrics_df, ggplot2::aes(x = method, 
                                                             y = !!rlang::sym(statistic),
                                                             fill = method,
                                                             #label=sprintf("%0.3f", round(!!rlang::sym(statistic), digits = 3)))) + 
                                                             label=sprintf("%0.2f", round(!!rlang::sym(statistic), digits = 2)))) + 
  ggplot2::geom_bar(
                      stat = "identity", 
                      position = "dodge",
                      color = "black",
                      alpha = 0.7) + 
    ggplot2::facet_wrap(~tissue_name, ncol = 4) + 
    ggplot2::scale_y_continuous(breaks=c(0, 0.5, 1.00), limits = c(-0.05, 1.1)) +
    ggplot2::xlab("Tabula Muris Tissues") +
    ggplot2::ylab(y_label) + 
    ggplot2::ggtitle(y_label) + 
    ggplot2::labs(fill = "Method") + 
    #ggplot2::scale_fill_brewer(palette = "Set1", labels = c("callback", "sc-SHC", "CHOIR")) + 
    ggplot2::scale_fill_manual(values = c("red", "grey", "black"), labels = c("callback", "sc-SHC", "CHOIR")) + 
    ggplot2::theme_bw() +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = small_text_size),
                   axis.title = ggplot2::element_text(size = large_text_size),
                   strip.text = ggplot2::element_text(size = small_text_size), 
                   legend.text = ggplot2::element_text(size = small_text_size, family = "Courier"),
                   legend.title = ggplot2::element_text(size = small_text_size),
                   plot.title = ggplot2::element_text(size = large_text_size, hjust = 0.5)) + 
    ggplot2::geom_text(vjust = -0.2)
    
  return(bar_plot)
}
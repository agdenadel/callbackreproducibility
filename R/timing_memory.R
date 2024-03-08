
#' @title Timing barplot for benchmarking a dataset.
#'
#' @description Creates a barplot showing the timing data for callback, sc-SHC,
#' and CHOIR benchmarked on a dataset.
#'
#' @param timing_memory_df The dataframe containing the timing data
#' @param small_text_size The text size for all but the title and axis titles.
#' @param large_text_size The title and axis titles text size.
#' @param title The title for the plot.
#' @returns Returns a ggplot2 object containing the timing bar plot.

#' @name timing_barplot
#' @export
timing_barplot <- function(timing_memory_df, small_text_size=16, large_text_size=22, title) {
  timing_summary_df <- timing_memory_df %>% dplyr::group_by(method) %>% dplyr::summarize(
    mean = mean(time), 
    sd = stats::sd(time))
  
  p <- ggplot2::ggplot(timing_summary_df, ggplot2::aes(x=method, y=mean, fill=method)) + #, label=sprintf("%0.1f", round(mean, digits = 1)))) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.7, colour = 'black') +
    ggplot2::geom_errorbar(ggplot2::aes(x=method, ymin=mean-sd, ymax=mean+sd), width=.5) + 
    #geom_text(vjust = -0.2) +
    ggplot2::ggtitle(title) + 
    ggplot2::xlab("Method") + 
    ggplot2::ylab("Time Taken (minutes)") +
    ggplot2::ylim(0, 120) +
    ggplot2::labs(fill = "Method") + 
    #ggplot2::ylim(c(0,18)) +
    ggplot2::scale_fill_manual(values = c("red", "grey", "black")) + 
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = large_text_size + 6, hjust = 0.5),
                   axis.ticks.x = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = small_text_size),
                   axis.title = ggplot2::element_text(size = large_text_size),
                   strip.text = ggplot2::element_text(size = small_text_size), 
                   legend.text = ggplot2::element_text(size = small_text_size, family = "Courier"),
                   legend.title = ggplot2::element_text(size = small_text_size))
  
  return(p)
}



#' @title Memory barplot for benchmarking a dataset.
#'
#' @description Creates a barplot showing the peak memory usage data for
#' callback, sc-SHC, and CHOIR benchmarked on a dataset.
#'
#' @param timing_memory_df The dataframe containing the memory usage data
#' @param small_text_size The text size for all but the title and axis titles.
#' @param large_text_size The title and axis titles text size.
#' @param title The title for the plot.
#' @returns Returns a ggplot2 object containing the memory usage bar plot.

#' @name memory_barplot
#' @export
memory_barplot <- function(timing_memory_df, small_text_size=16, large_text_size=22, title) {
  timing_memory_df$memory <- timing_memory_df$memory * 0.00104858 # convert from mebibytes to gigabytes
  
  memory_summary_df<- timing_memory_df %>% dplyr::group_by(method) %>% dplyr::summarize(
    mean = mean(memory), 
    sd = stats::sd(memory))
  
  p <- ggplot2::ggplot(memory_summary_df, ggplot2::aes(x=method, y=mean, fill=method)) + #, label=sprintf("%0.1f", round(mean, digits = 1)))) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.7, colour = 'black') + 
    ggplot2::geom_errorbar(ggplot2::aes(x=method, ymin=mean-sd, ymax=mean+sd), width=.5) + 
    #geom_text(vjust = -0.2) +
    ggplot2::ggtitle(title) + 
    ggplot2::xlab("Method") + 
    ggplot2::ylim(0, 13) +
    ggplot2::ylab("Peak Memory (GB)") +
    ggplot2::labs(fill = "Method") + 
    #ggplot2::ylim(c(0,18)) +
    ggplot2::scale_fill_manual(values = c("red", "grey", "black")) + 
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = large_text_size + 6, hjust = 0.5),
                   axis.ticks.x = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = small_text_size),
                   axis.title = ggplot2::element_text(size = large_text_size),
                   strip.text = ggplot2::element_text(size = small_text_size), 
                   legend.text = ggplot2::element_text(size = small_text_size, family = "Courier"),
                   legend.title = ggplot2::element_text(size = small_text_size))
  return(p)
}




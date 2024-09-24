`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' Save Figures
#'
#' @param name Name of figure
#' @param fig ggplot or similar figure object
#' @param width Width of plot in inches. Default = 6
#' @param height Height of plot in inches. Default = 6
#' @param plot_dir Plotting directory. Defaults to "analysis/plots"
#' @param pdf_plot Logical for plotting pdf too. Default = TRUE
#' @param font_family If specified, sets all font family. Default = NULL
#' @param res Image resolution in dpi. Default = 300
#' @param ... Other parameters to pass to ragg::agg_png
#' @importFrom grDevices dev.off pdf
save_figs <- function(name,
                      fig,
                      width = 6,
                      height = 6,
                      plot_dir = file.path(here::here(), "analysis/plots"),
                      pdf_plot = TRUE,
                      font_family = "Helvetica",
                      res = 300,
                      ...) {

  if(!is.null(font_family)) {
    fig <- fig + ggplot2::theme(text = ggplot2::element_text(family = font_family))
  }

  dir.create(plot_dir, showWarnings = FALSE)
  fig_path <- function(name) {paste0(plot_dir, "/", name)}

  ragg::agg_png(fig_path(paste0(name,".png")),
                width = width,
                height = height,
                units = "in",
                res = res,
                ...)
  print(fig)
  dev.off()

  if(pdf_plot) {
    pdf(file = fig_path(paste0(name,".pdf")), width = width, height = height)
    print(fig)
    dev.off()
  }

}

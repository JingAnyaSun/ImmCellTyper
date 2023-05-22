#' Plot median marker expression heatmap of binary classification results
#'
#' Visaulation of the median marker expression of each auto-annotated cluster
#' @param sce a SingleCellExperiment object
#' @param binary.results list, results from binaryClust function
#' @param remove.unclass logical vector, option to remove unclassified cells from plot, default is F
#' @return a ggplot object
#' @export
#'
#' @examples
#' plotbcHeatmap(sce,binary.results,remove.unclass=F)
plotbcHeatmap<-function (sce,binary.results, remove.unclass = FALSE)
{
  #define the function of binary_summary
  binary_summary<-function (data = NULL, binary.results = NULL)
  {
    if (is.null(data)) {
      return("Please provide data. Use load_data")
    }
    if (is.null(binary.results)) {
      return("Please provide binary classification results. Use binary_class")
    }
    frequencies <- data.frame(table(binary.results[, "Cell.Type"]))
    colnames(frequencies) <- c("Cell.Type", "Frequency")
    frequencies[, "Percentage"] <- frequencies[, "Frequency"] *
      100/sum(frequencies[, "Frequency"])
    medians <- data.frame(matrix(nrow = dim(frequencies)[1],
                                 ncol = dim(data)[2]))
    rownames(medians) <- frequencies[, 1]
    colnames(medians) <- colnames(data)
    for (i in frequencies[, 1]) {
      medians[i, ] <- as.vector(apply(data[binary.results[,
                                                          "Cell.Type"] == i, ], 2, median))
    }
    medians[, "Cell.Type"] <- rownames(medians)
    to.return <- merge(medians, frequencies, by = "Cell.Type")
    return(to.return)
  }
  exprs_n <- t(assay(sce, 'exprs'))

  binary.summary <- binary_summary(exprs_n, binary.results)

  if (is.null(binary.summary)) {
    return("Error: Please provide binary classification summary (use binary_summary).")
  }
  if (remove.unclass == FALSE) {
    to.plot <- binary.summary
  }
  else if (remove.unclass == TRUE) {
    to.plot <- binary.summary[binary.summary[, "Cell.Type"] !=
                                "Unclassified", ]
  }
  rownames(to.plot) <- to.plot[, "Cell.Type"]
  to.plot <- subset(to.plot, select = -c(Cell.Type, Frequency,
                                         Percentage))
  hm_pal = rev(brewer.pal(11, "RdYlBu"))
  medians.heatmap<-pheatmap(to.plot, color = colorRampPalette(hm_pal)(100),
                            cluster_cols = F,
                            cluster_rows = T,
                            border_color="white",
                            cellwidth = 10,
                            cellheight = 17)

  grid.newpage()
  plot <- grid.draw(medians.heatmap$gtable)
  return(plot)
}

#' Plot total cell population frequencies after binary classification
#'
#'Visualisation of cell population frequencies/composition after binary classification, in the format of bar plot
#' @param sce a SingleCellExperiment object
#' @param binary_results list, results of binary classification
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' plotbcFreq(sce,binary.results)
#'
plotbcFreq<-function(sce,binary.results){

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
  p<-ggplot(binary.summary,  aes("", Percentage, fill = Cell.Type)) +
    geom_bar(stat="identity", position = position_dodge())  +
    theme_minimal()+
    scale_y_continuous()+
    geom_text(aes(label=round(Percentage,2)), position=position_dodge(width=0.9), vjust=0.1,size=3.5)+
    labs(title = "Cell Abundance", x = "Cell types", y = "Percentage(%)")+
    coord_flip()
  return(p)
}

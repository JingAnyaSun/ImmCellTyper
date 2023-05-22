#' Performing binary classification
#'
#' Binary classification and annotation will be performed based on user-defined classification matrix
#' @param data data frame of marker expression, which can be exported from sce object
#' @param class.file file path to user pre-defined classification matrix(a csv file)
#' @param scale logical vector, standarisation of data, default is F
#'
#' @return table of the classification/defination results of each cell
#' @export
#'
#' @examples
#' binary_results <- binaryClass(exprs_n, class_dir)
binaryClass <- function ( data = NULL,
                           class.file = NULL,
                           scale = FALSE ) {

  if ( is.null(data) ) {
    return ('Please provide data. Use load_data')
  }

  if ( is.null(class.file) ) {
    return ('Please provide a classification file.')
  } else {

    data <- data[, colSums(data != 0, na.rm = TRUE) > 0]

    if ( scale == TRUE ) {
      data <- data.frame(scale(data))
    } else if ( scale == FALSE ) {
      data <- data
    } else {
      return ('Please choose either TRUE or FALSE for scale.')
    }

    kmeans.results <- data.frame(apply(data, 2, run_kmeans))
    colnames(kmeans.results) <- colnames(data)
    kmeans.limits <- data.frame(matrix(nrow = 0, ncol = 3))
    colnames(kmeans.limits) <- c('pos_min', 'Q1', 'Q3')

    for ( marker in colnames(kmeans.results) ) {

      if ( length(unique(kmeans.results[,marker])) == 1 & unique(kmeans.results[,marker])[1] == '-') {
        pos_min <- 999
        Q1 <- 0
        Q3 <- 0
      } else {
        pos_min <- as.numeric(min(data[kmeans.results[,marker] == '+',marker]))
        Q1 <- as.numeric(quantile(data[kmeans.results[,marker] == '+',marker], 0.25))
        Q3 <- as.numeric(quantile(data[kmeans.results[,marker] == '+',marker], 0.75))
      }
      kmeans.results[,marker] <- as.character(kmeans.results[,marker])
      kmeans.results[data[,marker] <= Q1 & kmeans.results[,marker] == '+',marker] <- 'L'
      kmeans.results[data[,marker] > Q1 & data[,marker] <= Q3,marker] <- 'I'
      kmeans.results[data[,marker] > Q3,marker] <- 'H'
      kmeans.limits[marker,] <- c(pos_min, Q1, Q3)

    }


    kmeans.results[,'Number'] <- rownames(kmeans.results)
    class.results <- data.frame(matrix(nrow = dim(data)[1], ncol = 1))
    colnames(class.results) <- c('Cell.Type')
    class.cat <- read.table(class.file, header = TRUE, sep = ',', row.names = 1, check.names = FALSE)
    class.cat <- apply(class.cat, c(1,2), function (x) { gsub(' ', '', x) })
    class.cat[class.cat == ''] <- 'A'

    for ( cell.type in rownames(class.cat) ) {

      temp.row <- class.cat[cell.type,, drop = FALSE]
      temp.row <- temp.row[,temp.row != 'A', drop = FALSE]

      temp.row.fine <- temp.row[,temp.row != '+', drop = FALSE]
      temp.row.all <- temp.row[,temp.row == '+', drop = FALSE]

      rows.required <- merge(kmeans.results, temp.row.fine)
      if ( length(colnames(temp.row.all)) == 1 ) {
        rows.required <- rows.required[rows.required[,colnames(temp.row.all)] %in% c('L', 'I', 'H'),][,'Number']
      } else {
        rows.required <- rows.required[apply(rows.required[,colnames(temp.row.all)], 1, function (x) { all(x %in% c('L', 'I', 'H')) }),][,'Number']
      }
      class.results[rows.required, 'Cell.Type'] <- cell.type

    }

    class.results[is.na(class.results[,'Cell.Type']),'Cell.Type'] <- 'Unclassified'

    return ( class.results )

  }

}

run_kmeans <- function ( data,
                         cutoff.limit = 2 ) {

  k.means <- kmeans(data, 2, iter.max = 500)
  pos <- which(k.means$centers == max(k.means$centers))
  neg <- which(k.means$centers == min(k.means$centers))
  results <- k.means$cluster

  if ( min(k.means$centers) > cutoff.limit | max(k.means$centers) < cutoff.limit | abs(max(k.means$centers) - min(k.means$centers)) < 1 | abs(max(data) - max(k.means$centers)) > cutoff.limit ) {

    results[data > cutoff.limit] <- '+'
    results[data <= cutoff.limit] <- '-'

  } else {

    results[results == pos] <- '+'
    results[results == neg] <- '-'

  }

  return ( results )

}



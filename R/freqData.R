#' Extract raw data(table) of cell frequencies from a sce object
#'
#'Directly extract raw data of cell frequencies from a sce object in order to perform more customized analysis dependent on
#'the research questions
#' @param sce a SingleCellExperiment objcet
#' @param metadata logical vector, default is T which means meatadata will be added into the table
#'
#' @return a table with sample_id, cluster_id, frequency, and metadata(e.g.condition, outcome) if 'metadata=T'.
#' @export
#'
#' @examples
#'
#' fq<-freqData(sce,metadata=T)
freqData<-function(sce,metadata=T){
  if(is.null(cluster_ids(sce))){
    stop('Error:Input cluster_id value is null')
  }
  n <- table(
    sample_id = sce$sample_id,
    cluster_id = sce$cluster_id)
  p <- prop.table(n, 1)
  # turn into tidy table
  df <- as.data.frame(p)
  # add relevant metadata
  if(metadata==T){
    md <- colnames(colData(sce))
    md <- md[!md %in% c('sample_id','cluster_id','batch')]
    message('The following metadata will be added:')
    print(md)
    Sys.sleep(2)
    i <- match(df$sample_id, sce$sample_id)
    cd <- colData(sce)[i, md]
    df <- cbind(df, cd)
    return(df)
  }else{
    message('No metadata is added')
    Sys.sleep(2)
    return(df)
  }
}

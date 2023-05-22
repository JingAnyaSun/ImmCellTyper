#'
#' Extract a certain cell cluster of interest from the total population
#' @param sce a SingleCellExperiment object
#' @param cluster character string specifying the cell cluster to extract
#' @param null_cluster_id logical vector to determine whether to clear out the cluster_id in the cluster to be extracted,
#' default is 'T', which is also recommended, because further interrogation/clustering might be done in the new cluster
#'
#' @return a SingleCellExperiment object for the extracted cluster
#' @export
#'
#' @examples
#' # to extract B cells from total cells, and define it to a new sce object named b_sce
#' b_sce<-extractCluster(sce,cluster='B Cells',null_cluster_id = F)

extractCluster<-function(sce,cluster,null_cluster_id=T){
  clu<-unique(sce$cluster_id)
  if(!is(sce,'SingleCellExperiment')){
    message('Error: please input an sce object')
  }
  if(is.null(cluster_ids(sce))){
    stop('Error: Input cluster id value is null, please do binaryClust first!')
  }
  if(!cluster %in% clu){
    message('Error: Input cluster name must be consistent with types file! Please select one of the following:')
    print(unique(cluster_ids(sce)))
  }else{
    clu_sce <- sce[,sce$cluster_id == cluster]
    if(null_cluster_id==T){
      clu_sce$cluster_id <- NULL
      return(clu_sce)}else{
        return(clu_sce)
      }
  }
}

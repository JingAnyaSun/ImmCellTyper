#' Downsample sce object
#'
#' Decrease the number of cells per sample by random selection to help reduce the computational burden
#' which may incur during data processing or clustering
#' @param sce a SingleCellExperiment object
#' @param n integer, number of cells to be downsampled
#' @param seed numeric vector; the seed set for each run of clustering to keep the results consistent between runs
#' @return a sce object
#' @export
#'
#' @examples
#' sce_d <- downsample(sce, n=2000)
#'
downsample<-function(sce,n,seed=1234){

  if(is.numeric(n)==T & n>0){
    message('Perform downsampling ',n, ' cells per sample')
    #take the same number of cells from each patient sample
    inds <- split(1:length(colData(sce)[,'sample_id']), colData(sce)[,'sample_id'])
    down_cells <- pmin(table(colData(sce)[,'sample_id']), n)
    #set the seed for downsampling
    set.seed(seed)
    down_inds <- lapply(names(inds), function(i){
      s <- sample(inds[[i]], down_cells[i], replace = FALSE)
    })
    down_inds <- unlist(down_inds)
    sce<-sce[,down_inds]
    return(sce)
  }else{
    stop('Error: downsample has to be a positive interger!')
  }
}

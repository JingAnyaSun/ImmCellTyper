#' Plot expression intensity of the state markers
#'
#' Plot state marker expression intensity of each cluster
#' @param sce a SingleCellExperiment object
#' @param group character string specifying the parameter for comparison
#'
#' @return a ggplot object and table for statistical analysis
#' @export
#'
#' @examples
#' plotStateMarkerExprs(sce,group='condition')
plotStateMarkerExprs <- function(sce,group){
  if(is.null(sce)){
    message('sce is Null, please input an sce object')
  }else if(!is(sce,'SingleCellExperiment')){
    stop('Error:Input data has to be a sce object')
  }
  s_marker<-state_markers(sce)
  exprs <- as.data.frame(t(assay(sce)))
  stat_exprs <-select(exprs,s_marker)
  stat_exprs <- cbind(stat_exprs, sample_id = sce$sample_id,cluster_id=cluster_ids(sce))
  suppressMessages(medexprs<-stat_exprs%>%group_by(cluster_id,sample_id)%>% summarize(across(where(is.numeric), median)))
  meta<-select(as.data.frame(metadata(sce)$experiment_info),-n_cells)
  medexprs<-merge(medexprs,meta,by='sample_id')
  suppressMessages(medexprs_n_l <- melt(medexprs))
  names <- colnames(meta)[!colnames(meta)%in%c('sample_id','patient_id')]
  if(!group %in% names){
    stop('Please choose the right group names for comparison as follows:',paste0(names,', '))
  }

  ####plotting#####
  p <-ggplot(medexprs_n_l,aes_string('variable','value',fill=group))+
    geom_boxplot(position=position_dodge(width =0.6),width=0.4,outlier.size=0.1,size=0.2)+
    theme_bw()+labs(y='Expression')+
    facet_wrap(~cluster_id,scales = 'free')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8))
  print(p)
}

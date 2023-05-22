#' Plot of cell population frequencies across samples, conditions or clusters
#'
#' Visualization of cell population frequencies across samples or conditions after clustering
#' in the format of boxplot or stacked histogram
#' @param sce a SingleCellExperiment object
#' @param type character string specifying the type of plot, options are 'stacked' indicating
#' stacked histogram, and 'box' indicating boxplot;default is 'box'
#' @param title character string, title of the plot; default is 'Cell Abundance'
#' @param group_by character string specifying the metadata used to group by; valid values are for instance 'condition'etc.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' #plot a boxplot grouped by condition
#' plotFreq(sce,title='Cell Abundance',group_by='condition')
#' #plot a stacked histogram grouped by batch
#' plotFreq(sce,type='stacked',title='Cell Abundance',group_by='batch')
#'
plotFreq<-function(sce,type='box',title='Cell abundances',group_by){
  if(!is(sce,'SingleCellExperiment')){
    stop('Error: Input data is not an sce object')
  }
  if(is.null(cluster_ids(sce))){
    stop('Error:Input cluster_id value is null')
  }
  #stacked bar plot
  n <- table(
    sample_id = sce$sample_id,
    cluster_id = sce$cluster_id)
  p <- prop.table(n, 1)
  # turn into tidy table
  df <- as.data.frame(p)
  md <- colnames(colData(sce))
  md <- md[!md %in% c('sample_id','cluster_id')]
  i <- match(df$sample_id, sce$sample_id)
  cd <- colData(sce)[i, md]
  df <- cbind(df, cd)
  fp <- colnames(colData(sce))
  fp<-fp[!fp%in%c('sample_id','cluster_id')]
  if(!group_by%in%fp){
    stop('Error: group_by has to be one of the below:',paste0(fp,', '))
  }
  if(type=='stacked'){
    p <-ggplot(df,aes(x=patient_id,y=Freq,fill=cluster_id))+
      geom_bar(stat ="identity",alpha=0.9)+
      scale_y_continuous(expand = c(0,0),label=scales::percent)+
      labs(y="Cell Proportion(%)")+
      theme_bw()+
      theme(axis.title.x =element_blank(),
            legend.key.size = unit(0.5,'cm'),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
      theme(strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.spacing.y = unit(2, "lines"),
            panel.border = element_blank())+
      ggtitle(title)
    p+facet_wrap(group_by, scales = "free_x")
  } else if(type == 'box'){
    ggplot(df,aes_string(x=group_by,y='Freq', color=group_by))+
      geom_boxplot() + facet_wrap(.~cluster_id,scales="free")+
      labs(y='Cell proportion(%)')+
      theme_bw()+
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA,size = 0.5),
            strip.background = element_rect(size = 0.5),
            panel.background = element_blank()
      )+ggtitle(title)
  } else{
    stop('Please choose the right plot type:stacked or box')
  }
}

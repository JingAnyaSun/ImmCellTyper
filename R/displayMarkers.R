#' Display marker intensity distribution selected for Binary clustering
#'
#' @param sce a SingleCellExpriment object
#' @param type a data frame of user pre-defined classification matrix
#' @param color_by character string specifying the color coding, default is 'sample_id'
#' @param ncol a numeric vector specifying the number of subplot per row, default is 3
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' #input/read types file as a data frame
#' class_dir <- 'types.csv'
#' types <- read_csv('types.csv')
#'
#' #construct a sce object
#' sce <- prepData(fs, panel, md, features = panel$fcs_colname)
#'
#' #show marker intensity distribution of the type file(classification matrix) with 4 subplots per row, colored by study conditions
#' displayMarkers(sce,types,color_by='condition',n=4)
#'
displayMarkers <- function(sce,type,color_by='sample_id',ncol=3){
  if(is(sce,'SingleCellExperiment')){
    m <- colnames(type)[colnames(type)!='Cell Type']
    if(!all(m %in% rowData(sce)$marker_name)){
      message('Error: Please check types file to make sure all markers listed are in the panel')
    } else{
      df.in<- data.frame(t(assay(sce,'exprs'))[,m], colData(sce), check.names = FALSE)
      gg_df <- melt(df.in,
                    variable.name = "antigen",
                    id.vars = names(colData(sce)))

      p<-ggplot(gg_df, fill = NULL,
                aes_string(
                  x = 'value', y = "after_stat(ndensity)",
                  col = color_by, group = "sample_id"))+
        facet_wrap(~ antigen, scales = "free_x") +
        geom_density() +
        ylab("normalized density")+
        theme_classic() + theme(
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"))}
    p$facet$params$ncol <- ncol

    return(p)
  } else{
    message('Error: Input data is ', class(sce),', please input sce data')
  }
}

#' Prepare data for batch effect examination
#'
#' By loading sample metadata, panel metadata and fcs files to construct a sce object for batch effect examination and normalisation
#' @param fcs_dir file path for fcs files
#' @param md_dir  file path for sample metadata
#' @param panel_dir file path for panel metadata
#' @return a sce object
#' @export
#'
#' @examples
#' fcs_dir <- 'fcs'
#' metadata_dir <- 'fcs/metadata3.xlsx'
#' panel_dir <- 'fcs/panel_metadata2.xlsx'
#' prepsce <- batchPrep(fcs_dir,metadata_dir,panel_dir)
batchPrep<-function(fcs_dir,md_dir,panel_dir){
  #check input data
  if(is.null(md_dir)){
    stop('Error:md_dir is missing, with no default')
  }
  if(is.null(panel_dir)){
    stop('Error:panel_dir is missing, with no default')
  }
  if(is.null(fcs_dir)){
    stop('Error:fcs_dir is missing, with no default')
  }
  message('Loading metadata and fcs files...')
  md <- read_excel(md_dir)
  md$batch <- factor(md$batch)
  md$condition  <- factor(md$condition)
  md$sample_id  <- factor(md$sample_id, levels = md$sample_id[order(md$condition)])
  md$patient_id <- factor(md$patient_id)
  panel <- read_excel(panel_dir)
  #remove possible problematic strings
  panel$antigen <- gsub("-", "_", panel$antigen)
  #read fcs files
  fcs_files <- list.files(fcs_dir, pattern = "fcs$")
  fs <- read.flowSet(file.path(fcs_dir, fcs_files))
  #construct a sce object
  message('Constructing sce object...')
  sce <- prepData(fs, panel, md, md_cols = list(file='file_name',id='sample_id',factors=c('batch','condition','patient_id')), features = panel$fcs_colname)
  return(sce)
}

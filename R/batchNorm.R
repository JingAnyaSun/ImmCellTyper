#' Batch correction to remove unwanted non-biological variation among dataset
#'
#' Batch normalization using CytoNorm or CytofRUV to remove unwanted variation arising from
#' experimental or instrumental artifacts. Please refer to Mtrussart et al.(CytofRUV) and Gassen et al.(CytoNorm)
#' to set relevant parameters.
#' @param fcs_dir file path of the fcs files
#' @param md_dir file path of the metadata file
#' @param panel_dir file path of the panel data file
#' @param method character string specifying the correction method to be used, default is 'CytoNorm', other option is 'CytofRUV'
#' @param nCells integer, CytoNorm parameter to specify the number of cells for training and model,default is 6000
#' @param xdim integer, CytoNorm parameter to specify the xdim of flowSOM clustering, default is 5
#' @param ydim integer,CytoNorm parameter to specify the ydim of flowSOM clustering, default is 5
#' @param nClus integer, CytoNorm and CytofRUV parameter to specify the k number of flowSOM clustering, default is 10
#' @param nQ integer, CytoNorm parameter to specify the number of quantiles, default is 101
#' @param goal character string, CytoNorm parameter, default is 'mean'
#' @param k integer,CytofRUV parameters,default is 5
#' @param output_dir file path of the output files(after normalisation),default is 'Norm_output'folder under working directory
#'
#' @return fcs files after normalisation
#' @export
#'
#' @examples
#' #set the path for fcs files, sample metadata and panel metadata
#' fcs_dir <- 'fcs'
#' metadata_dir <- 'fcs/metadata2.xlsx'
#' panel_dir <- 'fcs/panel_metadata2.xlsx'
#' #Use the default CytoNorm method and parameters for batch normalisation
#' batchNorm(fcs_dir,metadata_dir,panel_dir)
#' #Still use CytoNorm method, but change k number(nClus) as 20 for flowSOM clustering
#' batchNorm(fcs_dir,metadata_dir,panel_dir,nClus=20)
#' #Use CytofRUV method, set nClus as 20 and k value as 6
#' batchNorm(fcs_dir,metadata_dir,panel_dir,method='CytofRUV',nClus=20,k=6)
#'
batchNorm<-function(fcs_dir=NULL,md_dir,panel_dir,method='CytoNorm',nCells=6000,xdim=5,ydim=5,nClus=10,nQ=101,
                    goal='mean',k=5,output_dir="Norm_output"){

  if(method=='CytoNorm'){
    #Identifying the data
    md <- read_excel(md_dir)
    panel <- read_excel(panel_dir)
    fcs_files <- list.files(fcs_dir, pattern = "fcs$")
    if(length(fcs_files)==0){
      stop('Error: Please point the right path for your fcs files')
    }
    data  <- data.frame(File = fcs_files,
                        Path = file.path(fcs_dir, fcs_files),
                        Type = md[match(fcs_files, md$file_name), ]$condition,
                        Batch = md[match(fcs_files, md$file_name), ]$batch,
                        stringsAsFactors = FALSE)
    data$Type <- ifelse(data$Type == 'Ref', 'Train', 'Validation')
    train_data <- dplyr::filter(data, Type == "Train")
    validation_data <- dplyr::filter(data, Type == "Validation")
    channels <- dplyr:: filter(panel, marker_class == 'type')$fcs_colname

    transformList <- flowCore::transformList(channels,cytofTransform)
    transformList.reverse <- flowCore::transformList(channels,
                                                     cytofTransform.reverse)
    #Testing whether clustering is appropriate
    fsom <- prepareFlowSOM(train_data$Path,
                           channels,
                           nCells = nCells,
                           FlowSOM.params = list(xdim = xdim,
                                                 ydim = ydim,
                                                 nClus = nClus,
                                                 scale = FALSE),
                           transformList = transformList,
                           seed = 1)
    cvs <- testCV(fsom,
                  cluster_values = c(5, 10, 15))
    message('CV of each cluster is below:')
    Sys.sleep(1)
    print(cvs$pctgs$`10`)
    message('Building the model...')
    #Training the model
    model <- CytoNorm.train(files = train_data$Path,
                            labels = train_data$Batch,
                            channels = channels,
                            transformList = transformList,
                            FlowSOM.params = list(nCells = nCells,
                                                  xdim = xdim,
                                                  ydim = ydim,
                                                  nClus = nClus,
                                                  scale = FALSE),
                            normMethod.train = QuantileNorm.train,
                            normParams = list(nQ = nQ,
                                              goal = goal),
                            seed = 1,
                            verbose = TRUE)
    #Normalizing data
    message('Normalizing the fcs files...')
    CytoNorm.normalize(model = model,
                       files = validation_data$Path,
                       labels = validation_data$Batch,
                       transformList = transformList,
                       transformList.reverse = transformList.reverse,
                       normMethod.normalize = QuantileNorm.normalize,
                       outputDir = output_dir,
                       prefix = "Norm_",
                       clean = TRUE,
                       verbose = TRUE)
    #copy Ref/control file into the output_dir
    file.copy(train_data$Path,output_dir)
    #save panel and sample metadata after normalisation
    message('Writing sample and panel metadata after normalisation...')
    write_xlsx(panel, path = file.path(output_dir,'Norm_panel.xlsx'))
    md$file_name[md$condition!='Ref'] <- paste("Norm_", md$file_name[md$condition!='Ref'], sep = "")
    write_xlsx(md, path = file.path(output_dir,'Norm_md.xlsx'))


    ##############Please note the markers for normalisation here are set as type markers
  } else if(method=='CytofRUV'){
    wd_data <-getwd()
    seed=1234
    clusters_nb=nClus
    md <- read_excel(md_dir)
    ## Loading the data
    data=load_data(wd_data,md_dir,panel_dir)
    ## Cluster the data
    data$daf=cluster_data(data$daf,seed,markers_to_use=data$lineage_markers,clusters_nb)
    ## prep for normalisation
    dir_name_norm_data=output_dir
    raw_data <- data.frame(sample = data$daf$sample_id, cluster=cluster_ids(data$daf,paste0('meta',nClus)), t(SummarizedExperiment::assay(data$daf,"exprs")))
    colnames(raw_data) <- gsub("^X", "",  colnames(raw_data))

    rep_samples=list(c(filter(md,condition=='Ref')$sample_id))
    cluster_list_rep_samples <- list(seq(1,20))
    k_value <- k
    seed=1234
    normalise_data(data=data,raw_data=raw_data,rep_samples=rep_samples, norm_clusters=cluster_list_rep_samples, k=k_value, num_clusters=clusters_nb,wd_data=wd_data,dir_norm_data=dir_name_norm_data)
    message('Done!')
  } else{
    message('Error: please choose the right method for batch normalisation: CytoNorm, CytofRUV')
  }
}

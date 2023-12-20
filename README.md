# ImmCellTyper
 Cytof data analysis toolkit with cell population auto-annotation function  
 
## Introduction
 
ImmCellTyper was developed based on the framework and structure of CyTOF data analysis workflow described by Nowicka et al (CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets), and has a novel in house semi-supervised clustering tool named BinaryClust, which enables automatic classification and annotation of main cell types as well as in-depth interrogation of certain subpopulation of interest. By using k-means, BinaryClust takes advantage of the behaviour of CyTOF markers, which is log normal distributed with zero inflation in most of cases, to separate positive and negative cell populations of each protein marker. Then align the results with user-defined classification matrix, cells can be classified and annotated as different populations.The accuracy and performance of the automatic annotation is comparable to manual gating. Users can further extract a certain population of interest, to perform further clustering and investigation. Moreover, a variety of tools are integrated to support data quality check, batch effect correction, dimension reduction, unsupervised clustering(flowSOM and Rphenograph), and sophisticated statistical testing for multiple group comparison etc. 

<img width="1180" alt="BinaryClust diagram2" src="https://github.com/JingAnyaSun/BinaryClust2/assets/106811059/2fb5e5c2-210f-45b4-a0ac-451ebda201df">


## Installation
The package can be installed via running the below commands:
```
if(!require(devtools, quietly = TRUE))
  install.packages("devtools")
  devtools::install_github("JingAnyaSun/ImmCellTyper")
```
## Required files 
ImmCellTyper pipeline requires fcs files(.fcs), sample metadata(.xlsx), panel metadata(.xlsx), and classification matrix(.csv). It is recommended to put all the files needed under the same folder. 
### Sample metadata
1. For running the batch evaluation/normalisation session, an example of sample metadata is listed below, which should contain the information of file names('file_name'), sample ids('sample_id'), grouping('condition'), patient ids('patient_id') and batches('batch'). And for the reference/control/anchor sample in each run, it should be labelled as 'Ref' in 'condition' column:
```
# A tibble: 4 × 5                                                                                 
  file_name   sample_id condition patient_id batch
  <chr>       <chr>     <chr>     <chr>      <dbl>
1 A1.fcs      HC1_B1    Ref       HC1            1
2 A2.fcs      HC2_B1    HC        HC2            1
3 Run3_A1.fcs HC1_B2    Ref       HC1            2
4 Run3_A2.fcs HC2_B2    HC        HC2            2
```
2. For running the differential analysis session, an example of sample metadata is as below. Apart from 'file_name', 'sample_id', 'condition', and'patient_id', users can add more grouping information like 'outcome' etc. But remember to specify it correspondingly in 'prepData' function. 
Batch information is not needed in this metadata file. 
```
# A tibble: 11 × 5                                                                                
   file_name        sample_id condition patient_id outcome
   <chr>            <chr>     <chr>     <chr>        <dbl>
 1 AM_Non_JAK2.fcs  pt1       Non_JAK2  pt1              1
 2 CJ_JAK2.fcs      pt2       JAK2      pt2              1
 3 FD_Non_JAK2.fcs  pt3       Non_JAK2  pt3              2
 4 JB_Non_JAK2.fcs  pt4       Non_JAK2  pt4              3
 5 JY_JAK2.fcs      pt5       JAK2      pt5              1
 6 KA_JAK2.fcs      pt6       JAK2      pt6              2
 7 LC_Non_JAK2.fcs  pt7       Non_JAK2  pt7              3
 8 LDD_Non_JAK2.fcs pt8       Non_JAK2  pt8              1
 9 PC_JAK2.fcs      pt9       JAK2      pt9              3
10 SBB_JAK2.fcs     pt10      JAK2      pt10             3
11 YE_Non_JAK2.fcs  pt11      Non_JAK2  pt11             2
```
### Panel metadata
The requirement for panel metadata is that it shall comprise the column names of the fcs files('fcs_colname'), protein names('antigen') and the status of markers(marker_class):'type' means phenotyping markers whereas 'state' indicates functional proteins. An example is shown below:  
```
# A tibble: 36 × 3                                                                                
   fcs_colname antigen      marker_class
   <chr>       <chr>        <chr>       
 1 Y89Di       CD45         type       
 2 Pr141Di     CD196_CCR6   type        
 3 Nd142Di     CD19         type        
 4 Nd143Di     CD127_IL-7Ra state        
 5 Nd144Di     CD38         type        
 6 Nd145Di     CD33         type        
 7 Nd146Di     IgD          type        
 8 Sm147Di     CD11c        type        
 9 Nd148Di     CD16         type        
10 Sm149Di     CD194_CCR4   state       
# ℹ 26 more rows
```
### Classification matrix 
The classification matrix is defined by users to specify the marker expression of each cell type, based on their biological knowledge and gating strategy. This file should be in the format of csv. As exemplified below: `+` means positively expressed, `-` means negatively expressed and `A` means "any". The marker name should be consistent with the panel metadata 'antigen' column. 
```
# A tibble: 7 × 11
  `Cell Type`          CD14  CD16  CD161 CD19  CD20  CD3   CD4   CD56_NCAM CD8   TCRgd
  <chr>                <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr>     <chr> <chr>
1 NK Cells             -     A     A     -     A     -     A     +         A     A    
2 Dendritic Cells      -     A     A     -     A     -     A     -         A     A    
3 Monocytes            +     A     A     -     -     -     A     A         A     A    
4 B Cells              -     -     -     +     A     -     A     A         A     A    
5 T Cells, Gamma Delta -     A     A     -     A     +     A     A         A     +    
6 T Cells, CD4         -     A     A     -     A     +     +     A         -     -    
7 T Cells, CD8         -     A     A     -     A     +     -     A         +     -   
```
## Procedure
If you have prepared all the required files mentioned above, then you are ready to start the workflow. It is suggested to first do the data quality and batch effects check(please refer to vignettes/Intro_to_batch_exam_corrct), then start the differential analysis(please see vignettes/Intro_to_data_analysis for detailed instructions). 

Example dataset used can be downloaded from Zenodo [10.5281/zenodo.7982165](https://doi.org/10.5281/zenodo.7982165).


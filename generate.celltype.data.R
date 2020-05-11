#' generate.celltype.data
#'
#' \code{generate.celltype.data} Takes expression & cell type annotations and creates celltype_data files which contain the mean and specificity matrices
#'
#' @param exp Numerical matrix with row for each gene and column for each cell. Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be cross referenced against the annot data frame.
#' @param annot dataframe with columns containing the cell type names associated with each column in exp. Requires "patient_id" and "pathology" columns.
#' @param groupName A human readable name for refering to the dataset being loaded
#' @param numCores Number of cores that should be used to speedup the computation
#' @return Filenames for the saved celltype_data files
#' @examples
#' # Load the single cell data
#' data("cortex_mrna")
#' expData = cortex_mrna$exp
#' fNames_ALLCELLS = generate.celltype.data(exp=expData,annotLevels,"allKImouse")
#' @export
#' @import parallel
#' @import future
#' @import ggdendro
#' @import gridExtra

generate.celltype.data <- function(exp,annot, GroupName, run_sct_transform = TRUE, numCores=1){
  library(Matrix)
  library(future)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  options(datatable.optimize=1)
  exp = Matrix::as.matrix(exp, sparse = TRUE)
  # Convert to data.table and inner join so it's clear what expression is from which cell type (this can be very slow with base R so use data.table)
  # - guide to data.table merge functions is here: https://gist.github.com/nacnudus/ef3b22b79164bbf9c0ebafbf558f22a0
  if(class(exp)=="matrix"){
    exp_dt = data.table(exp,keep.rownames = TRUE)
    exp_sparse = Matrix(exp)
  }else{
    stop("Expected exp to be a normal matrix... need to write functions for other instances")
  }
  # Set parameters to enable scTransform to use multiple cores
  #these steps are for all the cells at all levels
  future::plan(strategy = 'multicore', workers = numCores)
  options(future.globals.maxSize = 10 * 1024 ^ 3)
  
  # Variable transform the data using scTransform
  #devtools::install_github(repo = 'ChristophH/sctransform')
 
  if(run_sct_transform == FALSE){
	print("run_sct_transform set to FALSE. Skipping sct transform.")
	normalized_exp_umi = exp_sparse
  }else{
    print("run_sct_transform set to TRUE. Performing sct transform can take some time.")
	normalized_exp <- sctransform::vst(exp_sparse,return_corrected_umi=TRUE) #$umi_corrected
	normalized_exp_umi = normalized_exp$umi_corrected

  }
  

  normalized_exp_dt = data.table(as.matrix(normalized_exp_umi),keep.rownames = TRUE) # Need to convert into a normal matrix here.... this will slow the function down  when using large datasets
  rm(normalized_exp_umi)
  
  annot = data.table(annot)
  
  exp_long   = data.table::melt(normalized_exp_dt) %>% dplyr::rename(GeneSymbol=rn,cell_id=variable,Exp=value)
  rm(normalized_exp_dt)
  
  #if(all.equal(exp_long$cell_id, annot$cell_id, ignore.row.order = TRUE)){
  #  print('super duper')
    exp_merged = merge.data.table(exp_long,annot, on="cell_id", nomatch=0)
  #}else{
  #  print('Cell_ids in annot and exp_long do not match (replace - with . in annot).')
  #}
  
  rm(exp_long)
  
  ct_mean_level1_long = exp_merged[,list(mean=mean(Exp)),by = list(GeneSymbol,level1class)]
  ct_mean_level1 = acast(ct_mean_level1_long, GeneSymbol~level1class, value.var="mean")
  
  ct_totalExp_level1_long = ct_mean_level1_long[,list(totalExp=sum(mean)),by = list(GeneSymbol)]
  
  ct_meansSummed_level1 = ct_mean_level1_long[ct_totalExp_level1_long, on="GeneSymbol", nomatch=0]
  ct_totalExp_level1 = acast(ct_meansSummed_level1, GeneSymbol~level1class, value.var="totalExp")
  
  ct_specificity_level1_long = ct_meansSummed_level1[,list(specificity=mean/totalExp),by = list(GeneSymbol,level1class)]
  ct_specificity_level1 = acast(ct_specificity_level1_long, GeneSymbol~level1class, value.var="specificity")
  
  #generate list of matrices
  ct_level1 = list(annot=annot$level1class, mean_exp = ct_mean_level1,  
                   specificity=ct_specificity_level1)
  
  #level2
  ct_mean_level2_long = exp_merged[,list(mean=mean(Exp)),by = list(GeneSymbol,level2class)]
  ct_mean_level2 = acast(ct_mean_level2_long, GeneSymbol~level2class, value.var="mean")
  
  ct_totalExp_level2_long = ct_mean_level2_long[,list(totalExp=sum(mean)),by = list(GeneSymbol)]
  
  ct_meansSummed_level2 = ct_mean_level2_long[ct_totalExp_level2_long, on="GeneSymbol", nomatch=0]
  ct_totalExp_level2 = acast(ct_meansSummed_level2, GeneSymbol~level2class, value.var="totalExp")
  
  ct_specificity_level2_long = ct_meansSummed_level2[,list(specificity=mean/totalExp),by = list(GeneSymbol,level2class)]
  ct_specificity_level2 = acast(ct_specificity_level2_long, GeneSymbol~level2class, value.var="specificity")
  
  #generate list of matrices
  ct_level2 = list(annot=annot$level2class, mean_exp = ct_mean_level2,  
                   specificity=ct_specificity_level2)
  
  ctd = list(ct_level1,ct_level2)
 
  #add dendrogram data
  ctd = lapply(ctd,bin.specificity.into.quantiles,numberOfBins=40)
  library(ggdendro)
  ctd = lapply(ctd,prep.dendro)
    
  #save the data
  fNames=sprintf("CellTypeData_%s.rda", GroupName)
  save(ctd,file=fNames)
  
  options(datatable.optimize=0)
  
  return(fNames)
}





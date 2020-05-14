#' generate.individual.ctds
#'
#' \code{generate.individual.ctds} Splits expression & cell type annotations per individual and creates celltype_data files which contain the mean and specificity matrices
#'
#' @param exp Numerical matrix with row for each gene and column for each cell. Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be cross referenced against the annot data frame.
#' @param annot dataframe with columns containing the cell type names associated with each column in exp. Annot requires "patient_id" and "pathology" columns.
#' @param numCores Number of cores that should be used to speedup the computation
#' @return Filenames for the saved celltype_data files
#' @examples
#' # Load the single cell data
#' data("TM2_rawcounts")
#' data("TM2_metadata") #annot!!!
#' exp = TM2_rawcounts
#' annot = annot
#' ctd_individuals = generate.individual.ctds(exp,annot, run_sctransform = TRUE, numCores = 4)
#' @export
#' @import parallel
#' @import future
#' @import ggdendro
#' @import gridExtra
#' @import EWCE

generate.individual.ctds <- function(exp,annot,run_sctransform = TRUE, numCores = 4){
  library(EWCE)
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(reshape2)
  library(Matrix)
  library(stringr)

  annot$patient_info <- paste(annot$patient_id,annot$pathology, sep = "_")
  annot_split <- split(annot, annot$patient_info)
  
  #convert exp to matrix
  mat = as.matrix(exp)
  
  ## vector of donors
  #patient_id <- unique(annot$patient_id)
  
  #create list of split matrices based on annotation
  mat_list <- list()
  for(i in annot_split){
    #for item in annot_split, select cell id
    case = list(i$cell_id)
    #then split the data from the expression matrix based on the cell id
    mat_split = mat[,colnames(mat) %in% unlist(case)]
    #add each matrix to a list
    mat_list[[length(mat_list)+1]] <- mat_split
  }
  
  ctd_loc <- list()
  library(foreach)
  #remove_short!!!
  foreach(i = annot_split, name = names(annot_split), j = mat_list) %do% {
    print(i[1:3,1:3])
    print(j[1:3,1:3])
    GroupName <- name
    print(GroupName)
    j = drop.uninformative.genes(exp=j,level2annot = annot[["level2class"]])
  	annotLevels = list(level1class=i[["level1class"]],level2class=i[["level2class"]])
  	celltype_data_file = EWCE::generate.celltype.data(exp= j, annotLevels = annotLevels, groupName = GroupName)
    ctd_loc[[length(ctd_loc)+1]] <- celltype_data_file
  }
  save(ctd_loc, file = "ctd_loc.rda")
  return(ctd_loc) 
}

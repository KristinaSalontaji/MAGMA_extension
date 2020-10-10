################################################################################
#' generate.individual.ctds
#'
#' \code{generate.individual.ctds} Splits expression & cell type annotations per
#' individual and creates celltype_data files which contain the mean and
#'  specificity matrices
#'
#' @param exp Numerical matrix with row for each gene and column for each cell.
#' Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be
#' cross referenced against the annot data frame.
#' @param annot dataframe with columns containing the cell type names
#' associated with each column in exp. Annot requires "patient_id" and
#' "pathology" columns.
#' @param numCores Number of cores that should be used to speedup the computation.
#' Use 1 for using this package in windows system.
#' @param savePath Path where to save the celltype_data files.
#' @return Filenames for the saved celltype_data files
#' @examples
#' # Load the single cell data
#' data("TM2_rawcounts")
#' data("TM2_metadata") #annot!!!
#' exp = TM2_rawcounts
#' annot = annot
#' ctd_individuals = generate.individual.ctds(exp,annot, numCores = 1, savePath=getwd())
#' @import foreach
#' @importFrom Matrix as.matrix
#' @importFrom EWCE drop.uninformative.genes generate.celltype.data
#' @export

generate.individual.ctds <- function(exp,annot,numCores = 1,savePath=getwd()){

  annot$patient_info <- paste(annot$patient_id,annot$pathology, sep = "_")
  annot_split <- split(annot, annot$patient_info)

  #convert exp to matrix
  mat = Matrix::as.matrix(exp)

  ## vector of donors
  #patient_id <- unique(annot$patient_id)

  #create list of split matrices based on annotation
  mat_list <- list()
  for(i in names(annot_split)){
    #for item in annot_split, select cell id
    case = annot_split[[i]][["cell_id"]] %>% as.character()
    #then split the data from the expression matrix based on the cell id and
    #add each matrix to a list
    mat_list[[i]] = mat[,colnames(mat) %in% case]
  }

  outdir <- file.path(savePath, "ctd_out")
  dir.create(path = outdir)

  ctd_loc <- list()
  foreach::foreach(i = annot_split, name = names(annot_split), j = mat_list) %do% {
    GroupName <- name
    print(GroupName)
    print(i[1:3,1:3])
    print(j[1:3,1:3])
    j = EWCE::drop.uninformative.genes(exp=j,level2annot = i[["level2class"]])
    print("Uninformative genes are dropped")
  	annotLevels = list(level1class=i[["level1class"]],level2class=i[["level2class"]])
  	print("Generating celltype data file...")
  	celltype_data_file = EWCE::generate.celltype.data(exp= j,
  	                                                  annotLevels = annotLevels,
  	                                                  groupName = GroupName,
  	                                                  no_cores = numCores,
  	                                                  savePath = outdir)
    print("Saving ctd files...")
  	ctd_loc[[length(ctd_loc)+1]] <- celltype_data_file
  }
  save(ctd_loc, file = file.path(outdir, "ctd_loc.rda"))
  return(ctd_loc)
}

#' magma.run.individual.analysis
#'
#' \code{magma.run.individual.analysis} Takes cell type data files from individual cases, merges all MAGMA results and returns graphs based on diagnosis.
#'
#' @param List of cell type data file locations, ctd_loc.
#' @param Specificity species as the species of the cell type data provided
#' @param GWAS species as the species of the GWAS
#' @param numCores Number of cores that should be used to speedup the computation
#' @return Filenames for the saved celltype_data files
#' @examples
#' # Load the single cell data
#' data("ctd_loc")
#' ctd_loc = ctd_loc
#' ctd_individuals = magma.run.individual.analysis(ctd_loc,specificity_species, gwas_species, gwas_sumstats_path, genome_ref_path, EnrichmentMode, numberOfBins = 41)
#' @export
#' @import parallel
#' @import future
#' @import ggdendro
#' @import gridExtra
#' @import EWCE
#' @import tidyr
#' @import dplyr
#' @import tibble
#' @import stringr


magma.run.individual.analysis <- function(ctd_loc,specificity_species, gwas_species, gwas_sumstats_path, genome_ref_path, EnrichmentMode, numberOfBins = 41){
    
  # load the generated rda files
  result_files <- list()
  for(i in ctd_loc){
    fNames= i
    load(fNames)
    celltypedat = ctd
    celltypedat = prepare.quantile.groups(celltypedat,specificity_species=deparse(substitute(specificity_species)), gwas_species = deparse(substitute(gwas_species)), numberOfBins=41)
    print("prepared plots")
    #run magma analysis (input option about linear or Top10)
    #save result as an rda file based on name
    if(EnrichmentMode == "Top 10%") {
      result = calculate_celltype_associations(celltypedat,gwas_sumstats_path,genome_ref_path=genome_ref_path,specificity_species = specificity_species, EnrichmentMode = EnrichmentMode)
      fResult=sprintf("MAGMA_Top10_%s", i)
      save(result, file = fResult)
      result_files[[length(result_files)+1]] <- fResult
	  save(result_files, file= "result_files_Top10.rda")
    } else if (EnrichmentMode == "Linear") {
      result = calculate_celltype_associations(celltypedat,gwas_sumstats_path,genome_ref_path=genome_ref_path,specificity_species = specificity_species, EnrichmentMode = EnrichmentMode)
      fResult=sprintf("MAGMA_Linear_%s", i)
      save(result, file = fResult)
      result_files[[length(result_files)+1]] <- fResult
	  save(result_files, file= "result_files_linear.rda")
    }
  }
  
  
  #process the MAGMA output files
  level_1_list <- list()
  level_2_list <- list()
  
  library(stringr)
  
  for(i in result_files){
    fResult= i
    load(fResult)
    
    #create df per levell: select each results and pvalue by celltype level
    level_1 = select(result[[1]][["results"]], Celltype, log10p)
    level_2 = select(result[[2]][["results"]], Celltype, log10p)
    
    #assign diagnosis and ID
    level_1$ID = str_extract(fResult, "[A-Z][0-9]")
    level_2$ID = str_extract(fResult, "[A-Z][0-9]")
    
    #res <- strapplyc(fResult, "((?:[^_]*_){2})(.*).rda", simplify = TRUE)[2,]
    library(gsubfn)
    level_1$diagnosis = strapplyc(fResult, "((?:[^_]*_){4})(.*).rda", simplify = TRUE)[2,]
    level_2$diagnosis = strapplyc(fResult, "((?:[^_]*_){4})(.*).rda", simplify = TRUE)[2,]
  
    
    level_1_list[[length(level_1_list)+1]] <- level_1
    level_2_list[[length(level_2_list)+1]] <- level_2
  }
  
  #create list 
  
  level_1_results = do.call("rbind", level_1_list)
  level_2_results = do.call("rbind", level_2_list)
  full_results = list(level_1_results = level_1_results, level_2_results = level_2_results)

  log10p_adj_l1 = log10(0.5/(n_distinct(level_1_results$Celltype)*2))
  log10p_adj_l2 = log10(0.5/(n_distinct(level_2_results$Celltype)*2))

  library(ggplot2)
  library(cowplot)
  
  png(file = "results_level1.png", width = 662, height = 515)
  ggplot(level_1_results) +
    geom_boxplot(aes(x=log10p,y=Celltype,fill=diagnosis)) +
    geom_vline(xintercept = log10p_adj_l1) +
    theme_cowplot()
  dev.off()
  
  png(file = "results_level2.png", width = 662, height = 1022)
  ggplot(level_2_results) +
    geom_boxplot(aes(x=log10p,y=Celltype,fill=diagnosis)) +
    geom_vline(xintercept = log10p_adj_l2) +
    theme_cowplot()
  dev.off()
 
 
  save(full_results, file = "full_results_per_individual.rda")
  return(full_results)

}

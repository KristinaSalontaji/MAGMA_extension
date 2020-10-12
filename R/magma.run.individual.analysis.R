################################################################################
#' magma.run.individual.analysis
#'
#' \code{magma.run.individual.analysis} Takes cell type data files from
#' individual cases, merges all MAGMA results and returns graphs based on diagnosis.
#'
#' @param List of cell type data file locations, ctd_loc.
#' @param Specificity species as the species of the cell type data provided.
#' @param GWAS species as the species of the GWAS.
#' @param numCores Number of cores that should be used to speedup the computation.
#' @return Filenames for the saved celltype_data files
#' @examples
#' # Load the single cell data
#' data("ctd_loc")
#' ctd_loc = ctd_loc
#' ctd_individuals = magma.run.individual.analysis(ctd_loc,
#' specificity_species,
#' gwas_species,
#' gwas_sumstats_path,
#' genome_ref_path,
#' EnrichmentMode,
#' numberOfBins = 41,
#' savePath)
#' @importFrom stringr str_extract
#' @importFrom magrittr %>%
#' @importFrom dplyr select n_distinct
#' @importFrom gsubfn strapplyc
#' @importFrom ggplot2 ggplot geom_boxplot geom_hline aes ggsave ylab coord_flip
#' @importFrom cowplot theme_cowplot
#' @importFrom ggpubr stat_compare_means
#' @importFrom MAGMA.Celltyping prepare.quantile.groups calculate_celltype_associations
#' @export


magma.run.individual.analysis <- function(ctd_loc,
                                          specificity_species,
                                          gwas_species,
                                          gwas_sumstats_path,
                                          genome_ref_path,
                                          EnrichmentMode,
                                          numberOfBins = 41,
                                          savePath = getwd()){

  # load the generated rda files
  result_files <- list()
  for(i in ctd_loc){
    fNames= i
    print(fNames)
    load(fNames)
    celltypedat = ctd
    celltypedat = MAGMA.Celltyping::prepare.quantile.groups(celltypedat,
                                                            specificity_species=deparse(substitute(specificity_species)),
                                                            gwas_species = deparse(substitute(gwas_species)),
                                                            numberOfBins=41)
    #print("prepared plots")

    magma_out <- file.path(savePath, "magma_out")
    dir.create(magma_out)

    #run magma analysis (input option about linear or Top10)
    #save result as an rda file based on name
    if(EnrichmentMode == "Top 10%") {
      result = MAGMA.Celltyping::calculate_celltype_associations(celltypedat,
                                                                 gwas_sumstats_path,
                                                                 genome_ref_path=genome_ref_path,
                                                                 specificity_species = specificity_species,
                                                                 EnrichmentMode = EnrichmentMode)

      #fResult=sprintf("MAGMA_Top10_%s", i)
      fResult = file.path(magma_out, sprintf("MAGMA_Top10_%s", basename(i)))
      save(result, file = fResult)
      result_files[[length(result_files)+1]] <- fResult
	  save(result_files, file= file.path(magma_out, "result_files_Top10.rda"))
    } else if (EnrichmentMode == "Linear") {
      result = MAGMA.Celltyping::calculate_celltype_associations(celltypedat,
                                                                 gwas_sumstats_path,
                                                                 genome_ref_path=genome_ref_path,
                                                                 specificity_species = specificity_species,
                                                                 EnrichmentMode = EnrichmentMode)
      #fResult=sprintf("MAGMA_Linear_%s", i)
      fResult = file.path(magma_out, sprintf("MAGMA_Linear_%s", basename(i)))
      save(result, file = fResult)
      result_files[[length(result_files)+1]] <- fResult
	  save(result_files, file= file.path(magma_out, "result_files_linear.rda"))
    }
  }


  #process the MAGMA output files
  level_1_list <- list()
  level_2_list <- list()

  for(i in result_files){
    fResult= i
    load(fResult)

    #create df per level: select each results and pvalue by celltype level
    level_1 = dplyr::select(result[[1]][["results"]], Celltype, log10p)
    level_2 = dplyr::select(result[[2]][["results"]], Celltype, log10p)

    #assign diagnosis and ID
    level_1$ID = stringr::str_extract(basename(fResult), "[A-Z][0-9]")
    level_2$ID = stringr::str_extract(basename(fResult), "[A-Z][0-9]")

    #res <- strapplyc(fResult, "((?:[^_]*_){2})(.*).rda", simplify = TRUE)[2,]
    level_1$diagnosis = gsubfn::strapplyc(basename(fResult), "((?:[^_]*_){4})(.*).rda", simplify = TRUE)[2,]
    level_2$diagnosis = gsubfn::strapplyc(basename(fResult), "((?:[^_]*_){4})(.*).rda", simplify = TRUE)[2,]


    level_1_list[[length(level_1_list)+1]] <- level_1
    level_2_list[[length(level_2_list)+1]] <- level_2
  }

  #create list

  level_1_results = do.call("rbind", level_1_list)
  level_2_results = do.call("rbind", level_2_list)
  full_results = list(level_1_results = level_1_results, level_2_results = level_2_results)

  log10p_adj_l1 = log10(0.05/(dplyr::n_distinct(level_1_results$Celltype)*2))
  log10p_adj_l2 = log10(0.05/(dplyr::n_distinct(level_2_results$Celltype)*2))

  ggplot2::ggplot(level_1_results, ggplot2::aes(x= Celltype, y = abs(log10p))) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=diagnosis)) +
    ggplot2::geom_hline(yintercept = abs(log10p_adj_l1)) +
    ggplot2::ylab( "-log10(p)") +
    cowplot::theme_cowplot() +
    ggpubr::stat_compare_means(ggplot2::aes(group = diagnosis),
                               label.y = max(abs(level_1_results$log10p)),
                               label = "p.format") +
    ggplot2::coord_flip(ylim = c(min(abs(level_1_results$log10p)),
                                 max(abs(level_1_results$log10p))+1))
  ggplot2::ggsave(file = file.path(magma_out, "results_level1.png"), width = 6, height = 5, units = "in", dpi = 300)


  ggplot2::ggplot(level_2_results, ggplot2::aes(x= Celltype, y = abs(log10p))) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=diagnosis)) +
    ggplot2::geom_hline(yintercept = abs(log10p_adj_l2)) +
    ggplot2::ylab( "-log10(p)") +
    cowplot::theme_cowplot() +
    ggpubr::stat_compare_means(ggplot2::aes(group = diagnosis),
                               label.y = max(abs(level_2_results$log10p)),
                               label = "p.format") +
    ggplot2::coord_flip(ylim = c(min(abs(level_2_results$log10p)),
                                 max(abs(level_2_results$log10p))+1))
  ggplot2::ggsave(file = file.path(magma_out, "results_level2.png"), width = 6, height = 5, units = "in", dpi = 300)


  save(full_results, file = file.path(magma_out, "full_results_per_individual.rda"))
  return(full_results)

}

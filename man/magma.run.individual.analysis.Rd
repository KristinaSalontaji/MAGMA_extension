% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/magma.run.individual.analysis.R
\name{magma.run.individual.analysis}
\alias{magma.run.individual.analysis}
\title{magma.run.individual.analysis}
\usage{
magma.run.individual.analysis(
  ctd_loc,
  specificity_species,
  gwas_species,
  gwas_sumstats_path,
  genome_ref_path,
  EnrichmentMode,
  numberOfBins = 41,
  savePath = getwd()
)
}
\arguments{
\item{List}{of cell type data file locations, ctd_loc.}

\item{Specificity}{species as the species of the cell type data provided.}

\item{GWAS}{species as the species of the GWAS.}

\item{numCores}{Number of cores that should be used to speedup the computation.}
}
\value{
Filenames for the saved celltype_data files
}
\description{
\code{magma.run.individual.analysis} Takes cell type data files from
individual cases, merges all MAGMA results and returns graphs based on diagnosis.
}
\examples{
# Load the single cell data
data("ctd_loc")
ctd_loc = ctd_loc
ctd_individuals = magma.run.individual.analysis(ctd_loc,
specificity_species,
gwas_species,
gwas_sumstats_path,
genome_ref_path,
EnrichmentMode,
numberOfBins = 41,
savePath)
}

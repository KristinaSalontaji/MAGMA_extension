% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.individual.ctds.R
\name{generate.individual.ctds}
\alias{generate.individual.ctds}
\title{generate.individual.ctds}
\usage{
generate.individual.ctds(exp, annot, numCores = 1, savePath = getwd())
}
\arguments{
\item{exp}{Numerical matrix with row for each gene and column for each cell.
Row names are MGI/HGNC gene symbols. Column names are cell IDs which can be
cross referenced against the annot data frame.}

\item{annot}{dataframe with columns containing the cell type names
associated with each column in exp. Annot requires "patient_id" and
"pathology" columns.}

\item{numCores}{Number of cores that should be used to speedup the computation.
Use 1 for using this package in windows system.}

\item{savePath}{Path where to save the celltype_data files.}
}
\value{
Filenames for the saved celltype_data files
}
\description{
\code{generate.individual.ctds} Splits expression & cell type annotations per
individual and creates celltype_data files which contain the mean and
 specificity matrices
}
\examples{
# Load the single cell data
data("TM2_rawcounts")
data("TM2_metadata") #annot!!!
exp = TM2_rawcounts
annot = annot
ctd_individuals = generate.individual.ctds(exp,annot, numCores = 1, savePath=getwd())
}

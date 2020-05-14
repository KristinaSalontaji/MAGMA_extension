# MAGMA_extension
MAGMA celltyping analysis based on disease state of subjects

## Introduction

This R package can be used to generate cell type data files based on patient_id of single cell RNA sequencing data and to run MAGMA celltyping on each individual, allowing a comparison of cell type association results between disease and control subjects.

The package is an extension of the R packages "EWCE" (https://github.com/NathanSkene/EWCE) and MAGMA_Celltyping (https://github.com/NathanSkene/MAGMA_Celltyping)

For the installation and further details, please check these packages.


## Installation

Before installing this package it is neccesary to install the magma software package. Please download it from https://ctg.cncr.nl/software/magma. Please do note, the magma software which forms the backend of this package was developed by Christian de Leeuw from Daniella Posthuma's lab. If you use this package to generate publishable results then you must cite their publication (listed below).

The executable should be copied to /usr/local/bin so that R can find it. Then install this package as follows:

For the installation and further details, please check the above mentioned packages.

```
install.packages("devtools")
library(devtools)
install_github("nathanskene/ewce")
install_github("ChristophH/sctransform")
install_github("nathanskene/MAGMA_Celltyping")
install_github("KristinaSalontaji/MAGMA_extension")
```

## Using the package (basic usage)

### Set parameters to be used for the analysis

```{r, eval=FALSE}
# Set path the 1000 genomes reference data.
genome_ref_dir = "~/Downloads/g1000_eur"
if(!file.exists(sprintf("%s/g1000_eur.bed",genome_ref_dir))){
    download.file("https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip",destfile=sprintf("%s.zip",genome_ref_dir))
    unzip(sprintf("%s.zip",genome_ref_dir),exdir=genome_ref_dir)
}
genome_ref_path = sprintf("%s/g1000_eur",genome_ref_dir)
```

### Install and load all the required R packages and data

The MAGMA_extension package comes with a single cell RNA-seq dataset with cell IDs as columns and genes as rows 
The package also comes with annotation data which contains "cell_id", "level1class", "level2_class", "patient_id" and "pathology" columns. 
For the tutorial, we first need to load these datasets. 

```{r, eval=FALSE }
library(EWCE)
library(MAGMA.Celltyping)
library(MAGMA_extension)

# Load the raw data
data(TM2_rawcounts)

# Load the annotation data
data(TM2_annot)
```

### Generate the cell type data files based on individual subjects

After loading the raw data and the annotation data with the previously mentioned column specifications, we can run the generate.individual.ctds function.

```{r, eval=FALSE }
ctd_loc = generate.individual.ctds(exp,annot, numCores = 4)
```
The function returns a list containing the names of each celltype data file and will be used in the next analyses steps.


### Download and format the summary statistics file

The GWAS summary statistics need to meet the following requirements in order to be processed. You may have to do some pre-processing to make sure that these requirements are met.

- SNP, CHR, BP as first three columns, in this exact order.
- It has at least one of these columns: ("Z","OR","BETA","LOG_ODDS","SIGNED_SUMSTAT")
- It has all of these columns: ("SNP","CHR","BP","P","A1","A2")
- If no "N" column exists in the summary statistics, we need to input N into the call.

```{r, eval=FALSE }
# specify the directory of the GWAS summary statistics:
gwas_sumstats_path = "~/AD_sumstats_Jansenetal_MOD.txt"


# Format it (i.e. column headers, provide N argument if not contained withing the summary statistics fle)
col_headers = format_sumstats_for_magma(gwas_sumstats_path, N=369385)

```

### Map SNPs to Genes

Next to the summary statistics file, we also need to map SNPs to genes.

```{r, eval=FALSE}
genesOutPath = map.snps.to.genes(gwas_sumstats_path, genome_ref_path=genome_ref_path, N=369385)
```

### Run the main cell type association analysis per individual and generate plots per diagnosis

The to run the cell type association analysis per individual, a list of file locations with the individual cell type data files has to be provided ("ctd_loc"). This list was returned during the generate.individual.ctds function.
The analyses can be run in either linear or top10% enrichment modes. Which analysis is run is determined by the input "EnrichmentMode = "Linear" or "EnrichmentMode = "Top 10%".

The magma.run.individual.analysis function will automatically prpare quantile groups for each celltype and run the desired celltype association analysis for each individual cell type data file.
In our example, we are using an RNA-seq dataset from human tissue, so we will specify the "specificity_species" argument to "human". The "gwas_species" is also set to "human" in this example.
If we have another specificity species, we need to specify this with the "specificity_species" argument. To check species that can are compatible with MAGMA_Celltyping and the MAGMA_extension, we can use the One2One package (https://github.com/NathanSkene/One2One).

For the linear association analysis we will run


```{r, eval=FALSE }
ctd_individuals = magma.run.individual.analysis(ctd_loc = ctd_loc,specificity_species = "human", gwas_species = "human", gwas_sumstats_path = gwas_sumstats_path = gwas_sumstats_path, genome_ref_path = genome_ref_path, EnrichmentMode = "Linear")

```

For the Top 10% mode we will run

```{r, eval=FALSE }
ctd_individuals = magma.run.individual.analysis(ctd_loc,specificity_species = "human", gwas_species = "human", gwas_sumstats_path = gwas_sumstats_path_formatted, genome_ref_path = genome_ref_path, EnrichmentMode = "Top 10%")

```

The generated boxplots based on diagnosis will be saved in these files:

     "results_level1_boxplot.png"          "results_level2_boxplot.png"

## References

Please cite the following if you use the MAGMA_celltyping:

[Skene, et al. Genetic identification of brain cell types underlying schizophrenia.
Nature Genetics, 2018.](https://www.nature.com/articles/s41588-018-0129-5)

The package utilises the MAGMA package developed in the Complex Trait Genetics lab at VU university (not us!) so please also cite their work:

[de Leeuw, et al. MAGMA: Generalized gene-set analysis of GWAS data.
PLoS Comput Biol, 2015.](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004219)

The MAGMA_extension package also depends on the EWCE package, so if you use MAGMA_extension, please cite:

[Skene, et al. Identification of Vulnerable Cell Types in Major Brain Disorders Using Single Cell Transcriptomes and Expression Weighted Cell Type Enrichment.
Front. Neurosci, 2016.](https://www.frontiersin.org/articles/10.3389/fnins.2016.00016/full) 

If you use the cell type dataset from University College London, please cite:


If you use the GWAS on Alzheimer's disease, please cite:

[Jansen, et al. Genome-wide meta-analysis identifies new loci and functional pathways influencing Alzheimer’s disease risk.
Nature Genetics, 2019.] (https://www.nature.com/articles/s41588-018-0311-9)



# Composite mutations

![alt text](https://github.com/taylor-lab/composite-mutations/blob/master/data/fig1a.png "Hi there!")


### Instructions
Code in this repository can be used to recreate the essential parts of the main figures in the manuscript. All code is in R (>=3.5.0).

Install required R packages from CRAN and Bioconductor:
```r
## check for missing required packages, install them
required.packages <- c('data.table','ggplot2','cowplot','RColorBrewer',
                       'parallel','ggsignif','binom','scales',
                       'MASS','ggrepel','Hmisc','Rcpp','pheatmap',
                       'car','here','magrittr','knitr','rmarkdown')

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

## install DESeq2 from Bioconductor if not installed
if(!'DESeq2' %in% installed.packages()) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("DESeq2"))
}
```

(Optional) Several datasets used in this manuscript can take multiple hours to run. For convenience, precalculated results from these tests are available in the `data/` directory. These can be regenerated with the following commands:
```shell

## generate data/gene_enrichment.txt
Rscript r/run_gene_enrichment_test.R

## generate data/residue_enrichment.txt
Rscript r/run_residue_enrichment_test.R

## generate data/observed_vs_expected_compounds_impact.rds; data/observed_vs_expected_compounds_tcga.rds
## warning: this can take hours to finish
Rscript r/run_cohort_permutation_tests.R

## generate data/observed_vs_expected_composites_per_tmb_impact.rds
## warning: this can take hours to finish
Rscript r/run_tmb_permutation_tests.R

## generate data/cancertype_enrichment_permutation_test.txt
## warning: this can only be run on a linux-based compute cluster 
## with LSF installed for job scheduling 
## (default resources requested per LSF job are 32GB RAM and 24 available cores)
Rscript r/run_cancertype_permutation_tests.R
```

Generate HTML files containing graphics for each main figure from the command line, which can be viewed in your web browser:
```shell
R -e "rmarkdown::render('html/figure-1.Rmd', output_file = 'figure-1.html')"
R -e "rmarkdown::render('html/figure-2.Rmd', output_file = 'figure-2.html')"
R -e "rmarkdown::render('html/figure-3.Rmd', output_file = 'figure-3.html')"
R -e "rmarkdown::render('html/figure-4.Rmd', output_file = 'figure-4.html')"
```

External Data Figures (EDF) can similarly be generated:
```shell
R -e "rmarkdown::render('html/figure-extended.Rmd', output_file = 'figure-extended.html')"
```

To generate an HTML file containing values references in the main text:
```shell
R -e "rmarkdown::render('html/text_values.Rmd', output_file = 'text_values.html')"
```

### Citation
- URL: https://www.nature.com/articles/s41586-020-2315-8 
- DOI: https://doi.org/10.1038/s41586-020-2315-8

### Contact
E-mail any questions to [gorelica@stanford.edu](mailto:gorelica@stanford.edu?subject=[GitHub]%20Composite-Mutations%20paper).

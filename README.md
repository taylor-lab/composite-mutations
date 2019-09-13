# compound-mutations

Companion code to the paper: \
**__Pervasive compound oncogenic mutations in cancer__**

Code in this repository can be used to recreate the essential parts of the main figures in the manuscript. All code is in R (>=3.5.0).

### Instructions
Install required R packages from CRAN:
```r
install.packages(c('data.table','ggplot2','cowplot','RColorBrewer',
                   'ggsignif','binom','scales','MASS',
                   'ggrepel','here','Hmisc','pheatmap','car'))
```
Install additional required R packages from Bioconductor:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2"))
```
Download RNA-Seq count file GSE136295_FSR_RNAseq.featureCounts.cnt.csv.gz from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136295) (requires secure token until Aug 23, 2020).

Extract compressed datasets in the `data/` directory:
```shell
gunzip data/data_mutations.txt.gz
gunzip data/data_mutations_phased.txt.gz
gunzip data/data_mutations_tcga.txt.gz
gunzip data/GSE136295_FSR_RNAseq.featureCounts.cnt.csv.gz
```

Several datasets used in this manuscript can take multiple hours to run. For convenience, precalculated results from these tests are available in the `data/` directory. These can be regenerated with the following commands:
```shell
Rscript r/run_gene_enrichment_test.R                    ## generates data/gene_enrichment.txt
Rscript r/run_residue_enrichment_test.R                 ## generates data/residue_enrichment.txt
Rscript r/run_compound_rate_permutation_test.R          ## generates data/observed_vs_expected_compounds_impact.rds
Rscript r/run_compound_rate_permutation_test_per_tmb.R  ## generates data/observed_vs_expected_compounds_impact_per_tmb.rds
```

HTML files containing the figures can be generated from the command-line for each main figure as such:
```shell
R -e "rmarkdown::render('do/figure-1.Rmd', output_file = 'figure-1.html')"
R -e "rmarkdown::render('do/figure-2.Rmd', output_file = 'figure-2.html')"
R -e "rmarkdown::render('do/figure-3.Rmd', output_file = 'figure-3.html')"
R -e "rmarkdown::render('do/figure-4.Rmd', output_file = 'figure-4.html')"
```

An HTML file containing values references in the main text can also be generated:
```shell
R -e "rmarkdown::render('do/text_values.Rmd', output_file = 'text_values.html')"
```

### Citation
- URL: **pending** 
- DOI: **pending**

### Contact
E-mail any questions to [gorelica@mskcc.org](mailto:gorelica@mskcc.org?subject=[GitHub]%20Compound-Mutations%20paper).

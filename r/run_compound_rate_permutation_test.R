## % of tumors with any compound: expected vs observed
## fig 1b/c of paper; also extended with tcga permutation
## also extended figure-1 (expected vs observed at each TMB)

source(here::here('r/prerequisites.R'))
cpus <- 12 ## number of CPUs available on this machine for parallelization 

## C++ function to execute permutations
cppFunction('IntegerVector shuffle( IntegerVector samples, IntegerVector genes, int iterations) {
            IntegerVector possible_samples = unique(samples);
            IntegerVector possible_genes = unique(genes);
            int possible_samples_length = possible_samples.length();    
            int possible_genes_length = possible_genes.length();    
            int n = samples.length();
            bool replace = false;
            IntegerMatrix sample_compound_matrix(possible_samples_length, iterations);
            for (int i=0; i<iterations; ++i) {
                // repeat this per iteration 
                IntegerMatrix sample_gene_matrix(possible_samples_length,possible_genes_length);
                IntegerVector shuffled_genes = sample(genes, n, replace); // shuffle all mutations gene labels
                // initialize the sample/gene matrix with 0 mutation counts
                for (int s=0; s<possible_samples_length; ++s) { 
                    sample_compound_matrix(s,i) = 0;
                    for (int g=0; g<possible_genes_length; ++g) { 
                        sample_gene_matrix(s,g) = 0; 
                    }
                }
                // loop over the data and populate matrix of samples with the count of shuffled mutations in each gene
                for(int r=0; r<n; ++r) {
                    int current_sample = samples(r);
                    int current_gene = shuffled_genes(r);
                    sample_gene_matrix(current_sample, current_gene) = sample_gene_matrix(current_sample,current_gene) + 1;
                }
                // populate the matrix of sample|iteration, where cells are 0=no compound, 1=any compound mutations.
                for (int s=0; s<possible_samples_length; ++s) { 
                    for (int g=0; g<possible_genes_length; ++g) { 
                        if(sample_gene_matrix(s,g) > 1) {
                            sample_compound_matrix(s,i) = 1;
                            g = possible_genes_length;
                        }
                    }
                }
            }
            return sample_compound_matrix;
}')


## function to parallelize permutations across batches of 100 permutations each
run_for_batch <- function(batch,m,reps_per_batch=100) {
    message(batch)
    N <- length(unique(m$Tumor_Sample_Barcode))
    qc <- shuffle(m$Tumor_Sample_Barcode,m$Hugo_Symbol,reps_per_batch)
    out <- data.table(x=colSums(qc))
    out$batch <- batch
    out
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run permutation test for all samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## load sampleid and gene-symbols from mutation data
d <- fread(here('data/data_mutations.txt'),select=c('Tumor_Sample_Barcode','Hugo_Symbol'))
d <- d[order(Tumor_Sample_Barcode, Hugo_Symbol),]

## get the total number of samples (including any with 0 observed mutations)
samples <- fread(here('data/data_clinical.txt'),select=c('Tumor_Sample_Barcode'))[[1]]
N <- length(unique(samples))

## convert the data to a 0-indexed integer matrix
## Tumor_Sample_Barcode = 0,...,N-1 for N unique samples
## Hugo_Symbol = 0,...,467 (468 genes in panel)
m <- d
m$Tumor_Sample_Barcode <- as.integer(factor(m$Tumor_Sample_Barcode)) - 1
m$Hugo_Symbol <- as.integer(factor(m$Hugo_Symbol)) - 1

## run 1e5 iterations across 1000 parallelized batches of 100 
set.seed(42)
batches <- 1:1000
reps_per_batch <- 100
l <- mclapply(batches, run_for_batch, m, reps_per_batch, mc.cores=cpus)
ll <- rbindlist(l)
ll$prop <- ll$x / N 

## get the observed proportion of compound-mutant samples
d$id <- paste(d$Tumor_Sample_Barcode,d$Hugo_Symbol)
tbl <- table.freq(d$id)
d <- merge(d, tbl, by.x='id', by.y='value')
d$compound <- d$N > 1
summarize_sample <- function(d) {
    any.compound <- any(d$compound)
    list(any.compound=any.compound)
}
res <- d[,summarize_sample(.SD),by=Tumor_Sample_Barcode]
obs <- sum(res$any.compound) / N

## save the results to data/
results_all <- list(dat_all=ll,obs_all=obs)
saveRDS(results_all,file=here('data/observed_vs_expected_compounds_impact.rds'))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run permutation tests for tmb
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_permutation_per_tmb <- function(burden, d) {
    message(burden)
    m <- d[tmb==burden,c('Tumor_Sample_Barcode','Hugo_Symbol'),with=F]
    m$Tumor_Sample_Barcode <- as.integer(factor(m$Tumor_Sample_Barcode)) - 1
    m$Hugo_Symbol <- as.integer(factor(m$Hugo_Symbol)) - 1

    ## run 1e5 iterations across 1000 parallelized batches of 100 
    set.seed(42)
    batches <- 1:1000
    reps_per_batch <- 100
    l <- mclapply(batches, run_for_batch, m, reps_per_batch, mc.cores=cpus)
    ll <- rbindlist(l)
    Ntmb <- length(unique(m$Tumor_Sample_Barcode))
    ll$prop <- ll$x / Ntmb 

    ## get the observed proportion of compound-mutant samples
    m <- d[tmb==burden,c('Tumor_Sample_Barcode','Hugo_Symbol'),with=F]
    m$id <- paste(m$Tumor_Sample_Barcode,m$Hugo_Symbol)
    tbl <- table.freq(m$id)
    m <- merge(m, tbl, by.x='id', by.y='value')
    m$compound <- m$N > 1
    summarize_sample <- function(d) {
        any.compound <- any(d$compound)
            list(any.compound=any.compound)
    }
    res <- m[,summarize_sample(.SD),by=Tumor_Sample_Barcode]
    obs <- sum(res$any.compound) / Ntmb
    list(dat_all=ll, obs_all=obs, burden=burden, samples=Ntmb)
}

## load sampleid and gene-symbols from mutation data
d <- fread(here('data/data_mutations.txt'),select=c('Tumor_Sample_Barcode','Hugo_Symbol','tmb'))
d <- d[order(Tumor_Sample_Barcode, Hugo_Symbol),]
d$tmb <- round(d$tmb)

## run permutations for each TMB from 1:50
prev_run <- -1
burdens <- 1:50
burdens <- burdens[burdens %nin% prev_run]
l <- lapply(burdens, run_permutation_per_tmb, d)

## annotate results from each permutation with the mean, CIs and P-value
get_results <- function(qc) { 
    mu <- mean(qc$dat_all$prop)
    mid <- median(qc$dat_all$prop)
    lwr <- as.numeric(quantile(qc$dat_all$prop,0.025))
    upr <- as.numeric(quantile(qc$dat_all$prop,0.975))
    obs <- qc$obs_all
    R <- nrow(qc$dat_all)
    p <- sum(qc$dat_all$prop >= obs)/R
    samples <- qc$samples
    out <- list(burden=qc$burden,samples=samples,obs=obs,mu=mu,mid=mid,lwr=lwr,upr=upr,data=qc$dat_all,p=p)
}

l2 <- lapply(l, get_results)
saveRDS(l2,file=here('data/observed_vs_expected_compounds_per_tmb_impact.rds'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run permutation test for all samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## prep TCGA mutation data
clin <- fread(here('data/data_clinical_tcga.txt'))

## load and prep TCGA mutation data
dat <- fread(here('data/data_mutations_tcga.txt'))

## load sampleid and gene-symbols from mutation data
d <- dat[impact468_gene==T,]
d <- d[order(Tumor_Sample_Barcode, Hugo_Symbol),]

## get the total number of samples (including any with 0 observed mutations)
samples <- clin$Tumor_Sample_Barcode
N <- length(unique(samples))

## convert the data to a 0-indexed integer matrix
## Tumor_Sample_Barcode = 0,...,N-1 for N unique samples
## Hugo_Symbol = 0,...,467 (468 genes in panel)
m <- d
m$Tumor_Sample_Barcode <- as.integer(factor(m$Tumor_Sample_Barcode)) - 1
m$Hugo_Symbol <- as.integer(factor(m$Hugo_Symbol)) - 1

## run 1e5 iterations across 1000 parallelized batches of 100 
set.seed(42)
batches <- 1:1000
reps_per_batch <- 100
l <- mclapply(batches, run_for_batch, m, reps_per_batch, mc.cores=cpus)
ll <- rbindlist(l)
ll$prop <- ll$x / N 

## get the observed proportion of compound-mutant samples
d$id <- paste(d$Tumor_Sample_Barcode,d$Hugo_Symbol)
tbl <- table.freq(d$id)
d <- merge(d, tbl, by.x='id', by.y='value')
d$compound <- d$N > 1
summarize_sample <- function(d) {
    any.compound <- any(d$compound)
    list(any.compound=any.compound)
}
res <- d[,summarize_sample(.SD),by=Tumor_Sample_Barcode]
obs <- sum(res$any.compound) / N

## save the results to data/
results_all <- list(dat_all=ll,obs_all=obs)
saveRDS(results_all,file=here('data/observed_vs_expected_compounds_tcga.rds'))

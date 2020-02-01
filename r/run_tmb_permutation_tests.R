## % of tumors with any composite: expected vs observed
## fig 1b/c of paper; also extended with tcga permutation
## also extended figure-1 (expected vs observed at each TMB)

source(here::here('r/prerequisites.R'))
cpus <- cpu_prompt()

## C++ function to execute permutations
cppFunction('IntegerVector shuffle( IntegerVector samples, IntegerVector genes, int iterations) {
            IntegerVector possible_samples = unique(samples);
            IntegerVector possible_genes = unique(genes);
            int possible_samples_length = possible_samples.length();    
            int possible_genes_length = possible_genes.length();    
            int n = samples.length();
            bool replace = false;
            IntegerMatrix sample_composite_matrix(possible_samples_length, iterations);
            for (int i=0; i<iterations; ++i) {
                // repeat this per iteration 
                IntegerMatrix sample_gene_matrix(possible_samples_length,possible_genes_length);
                IntegerVector shuffled_genes = sample(genes, n, replace); // shuffle all mutations gene labels
                // initialize the sample/gene matrix with 0 mutation counts
                for (int s=0; s<possible_samples_length; ++s) { 
                    sample_composite_matrix(s,i) = 0;
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
                // populate the matrix of sample|iteration, where cells are 0=no composite, 1=any composite mutations.
                for (int s=0; s<possible_samples_length; ++s) { 
                    for (int g=0; g<possible_genes_length; ++g) { 
                        if(sample_gene_matrix(s,g) > 1) {
                            sample_composite_matrix(s,i) = 1;
                            g = possible_genes_length;
                        }
                    }
                }
            }
            return sample_composite_matrix;
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
# run permutation tests per tmb for IMPACT samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_permutation_per_tmb <- function(burden, d) {
    message('Running for TMB=',burden)
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

    ## get the observed proportion of composite-mutant samples
    m <- d[tmb==burden,c('Tumor_Sample_Barcode','Hugo_Symbol'),with=F]
    m$id <- paste(m$Tumor_Sample_Barcode,m$Hugo_Symbol)
    tbl <- table.freq(m$id)
    m <- merge(m, tbl, by.x='id', by.y='value')
    m$composite <- m$N > 1
    summarize_sample <- function(d) {
        any.composite <- any(d$composite)
            list(any.composite=any.composite)
    }
    res <- m[,summarize_sample(.SD),by=Tumor_Sample_Barcode]
    obs <- sum(res$any.composite) / Ntmb
    list(dat_all=ll, obs_all=obs, burden=burden, samples=Ntmb)
}

## load sampleid and gene-symbols from mutation data
d <- fread(here('data/data_mutations.txt.gz'),select=c('Tumor_Sample_Barcode','Hugo_Symbol','tmb'))
d <- d[order(Tumor_Sample_Barcode, Hugo_Symbol),]
d$tmb <- round(d$tmb)

## run permutations for each TMB from 1:50
prev_run <- -1
burdens <- 1:50
burdens <- burdens[burdens %nin% prev_run]

message('Will now run 100,000 permutations for samples with TMB from 1-50 (individually) to the determine proportion of samples with composite mutations expected per TMB in IMPACT dataset.\n')
message('This will run 1000 batches (100 permutations each) for samples with TMB from 1-50. Batches and TMB will be counted below as they run. Higher TMBs will run faster because these have fewer samples. This will take a while ...\n')
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
saveRDS(l2,file=here('data/observed_vs_expected_composites_per_tmb_impact.rds'))




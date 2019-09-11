## helper script to run permutations via Rscript for LSF parallelization 
## for figs 1c/1d of paper

setwd('/home/ang46/lab/projects/compound_mutations')
source('scripts/func.R')

library(Rcpp)
library(parallel)
args = commandArgs(trailingOnly=TRUE)

## function to execute permutations
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
    N <- length(unique(m$pid))
    qc <- shuffle(m$pid,m$Hugo_Symbol,reps_per_batch)
    out <- data.table(x=colSums(qc))
    out$prop <- out$x / N
    out$batch <- batch
    out
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run permutations on input data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

input_txt <- args[1]
output_rds <- args[2]
cores <- as.integer(args[3])
burden <- args[4]

d0 <- fread(input_txt)
samples <- length(unique(d0$pid))

d <- d0[,c('pid','Hugo_Symbol'),with=F]
d <- d[order(pid, Hugo_Symbol),]

## the data for shuffling is an integer matrix, pid|Hugo_Symbol, 
## where each is a 0-indexed integer for its corresponding factor
## i.e. 0 = genes[1]; 340=genes[341] (same for samples)
m <- d
m$pid <- as.integer(factor(m$pid)) - 1
m$Hugo_Symbol <- as.integer(factor(m$Hugo_Symbol)) - 1


## run 1e5 iterations across parallelized batches of 100 
set.seed(42)
batches <- 1:1000
reps_per_batch <- 100
cpus <- cores
l <- mclapply(batches, run_for_batch, m, reps_per_batch, mc.cores=cpus)
ll <- rbindlist(l)

d$id <- paste(d$pid,d$Hugo_Symbol)
tbl <- table.freq(d$id)
d <- merge(d, tbl, by.x='id', by.y='value')
d$compound <- d$N > 1
summarize_sample <- function(d) {
    any.compound <- any(d$compound)
    list(any.compound=any.compound)
}
res <- d[,summarize_sample(.SD),by=pid]
obs <- sum(res$any.compound) / nrow(res)
results_all <- list(dat_all=ll,obs_all=obs,samples=samples,burden=burden)
saveRDS(results_all,file=output_rds)




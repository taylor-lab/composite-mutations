
require(vegan)
require(data.table)
require(rpkg)
require(here)

args = commandArgs(trailingOnly=TRUE)
cancertype <- noquote(args[1])
dir_out <- here('data/cancertype_permutation_tests')

## load sampleid and gene-symbols from mutation data
message('loading data ...')
d <- fread(here('data/data_mutations.txt.gz'),
           select=c('Tumor_Sample_Barcode','Hugo_Symbol','exclude','putative_resistance_mutation','metamaintype'))
d <- d[order(Tumor_Sample_Barcode, Hugo_Symbol),]

message('Starting permutations for [',cancertype,'] at ',as.character(date()))
clin <- fread(here('data/data_clinical.txt.gz'))
clin <- clin[exclude==F & metamaintype %in% cancertype,]
N_samples <- nrow(clin)
d <- d[exclude==F & putative_resistance_mutation==F,]
d[,c('exclude','putative_resistance_mutation'):=NULL]
dd <- d[metamaintype==cancertype]
label <- gsub(' ','_',tolower(cancertype))


run_permutations <- function(batch,label,d,N_samples,R) {
    message('Batch ',batch)

    ## save result from this batch
    batchid <- as.character(1e5 + batch)
    batchid <- substr(batchid,2,nchar(batchid))
    fout <- file.path(dir_out,paste0(label,'_',batchid,'_batch.txt'))

    if(!file.exists(fout)) {
        x <- xtabs(~Tumor_Sample_Barcode + Hugo_Symbol,data=dd)
        x <- as.data.frame.matrix(x)
        ps <- permatswap(x, times = R, burnin = 100, thin = 10, mtype = "count", fixedmar='both', shuffle='both')
        f = function(m) {
            composites <- sum(rowSums(m > 1) > 0)
            list(composites=composites)
        }
        l <- rbindlist(lapply(ps$perm, f))
        l$prop <- l$composites / N_samples
        l$batch <- batch
        write.tsv(l, fout)
        rm(ps)
        trash <- gc()
    } else {
        trash <- sample(1:100)
        message('File exists: ',fout)
    }
}

batches <- 1:100
R <- 100
set.seed(42)
l <- lapply(batches, run_permutations, label, d, N_samples, R)
message('Finished permutations at ',as.character(date()))





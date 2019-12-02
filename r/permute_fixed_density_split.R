
require(vegan)
require(data.table)
require(rpkg)

args = commandArgs(trailingOnly=TRUE)
cancertype <- noquote(args[1])
N_samples <- as.integer(args[2])


## load sampleid and gene-symbols from mutation data
message('loading data ...')
d <- fread('/home/ang46/lab/repos/compound-mutations/data/data_mutations.txt',select=c('Tumor_Sample_Barcode','Hugo_Symbol','exclude','putative_resistance_mutation','metamaintype'))
d <- d[exclude==F & putative_resistance_mutation==F,]
d[,c('exclude','putative_resistance_mutation'):=NULL]
d <- d[order(Tumor_Sample_Barcode, Hugo_Symbol),]
message('Starting permutations for [',cancertype,'] at ',as.character(date()))

run_cancertype <- function(batch,cancertype,d,N_samples,R) {
    message('Batch ',batch)

    ## save result from this batch
    batchid <- as.character(1e5 + batch)
    batchid <- substr(batchid,2,nchar(batchid))
    ctype <- gsub(' ','_',tolower(cancertype))
    fout <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches/',ctype,'_',batchid,'_batch.txt')

    if(!file.exists(fout)) {
        dd <- d[metamaintype==cancertype]
        x <- xtabs(~Tumor_Sample_Barcode + Hugo_Symbol,data=dd)
        x <- as.data.frame.matrix(x)
        ps <- permatswap(x, times = R, burnin = 100, thin = 10, mtype = "count", fixedmar='both', shuffle='both')
        f = function(m) {
            compounds <- sum(rowSums(m > 1))
            list(compounds=compounds)
        }
        l <- rbindlist(lapply(ps$perm, f))
        l$prop <- l$compounds / N_samples
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
#batches <- 1:1000
#R <- 10


set.seed(42)
l <- lapply(batches, run_cancertype, cancertype, d, N_samples, R)
message('Finished permutations at ',as.character(date()))





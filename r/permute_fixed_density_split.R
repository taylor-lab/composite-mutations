
require(vegan)
require(data.table)
require(rpkg)
require(here)

args = commandArgs(trailingOnly=TRUE)
cancertype <- noquote(args[1])
tmb <- as.integer(noquote(args[2]))


## load sampleid and gene-symbols from mutation data
message('loading data ...')
d <- fread('/home/ang46/lab/repos/compound-mutations/data/data_mutations.txt',select=c('Tumor_Sample_Barcode','Hugo_Symbol','exclude','putative_resistance_mutation','metamaintype'))
d <- d[order(Tumor_Sample_Barcode, Hugo_Symbol),]

if(!is.na(cancertype) & cancertype %in% d$metamaintype) {
    message('Starting permutations for [',cancertype,'] at ',as.character(date()))

    clin <- fread(here('data/data_clinical.txt'))
    clin <- clin[exclude==F & metamaintype %in% cancertype,]
    N_samples <- nrow(clin)

    d <- d[exclude==F & putative_resistance_mutation==F,]
    d[,c('exclude','putative_resistance_mutation'):=NULL]
    dd <- d[metamaintype==cancertype]
    mode <- 'cancertype'
    label <- gsub(' ','_',tolower(cancertype))
} else if(!is.na(tmb) & tmb > 0) {
    ## permutation by TMB
    message('Starting permutations for [',tmb,'] at ',as.character(date()))
    clin <- fread(here('data/data_clinical.txt'))
    clin$tmb <- round(clin$tmb)
    retained_samples <- clin$Tumor_Sample_Barcode[clin$tmb %in% tmb]
    N_samples <- length(retained_samples)
    dd <- d[Tumor_Sample_Barcode %in% retained_samples,]
    mode <- 'tmb'
    label <- tmb
} else {
    message('Using complete MAF')
    clin <- fread(here('data/data_clinical.txt'))
    N_samples <- nrow(clin)
    dd <- d
    mode <- 'fullmaf'
    label <- 'fullmaf'
}


run_permutations <- function(batch,label,mode,d,N_samples,R) {
    message('Batch ',batch)

    ## save result from this batch
    batchid <- as.character(1e5 + batch)
    batchid <- substr(batchid,2,nchar(batchid))

    if(mode=='cancertype') {
        fout <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches_v2/',label,'_',batchid,'_batch.txt')
    } else if(mode=='tmb') {
        fout <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/tmb_permutations_batches/tmb-',label,'_',batchid,'_batch.txt')   
    } else if(mode=='fullmaf') {
        fout <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/fullmaf_permutations_batches/',batchid,'_batch.txt')      
    }

    if(!file.exists(fout)) {
        x <- xtabs(~Tumor_Sample_Barcode + Hugo_Symbol,data=dd)
        x <- as.data.frame.matrix(x)
        ps <- permatswap(x, times = R, burnin = 100, thin = 10, mtype = "count", fixedmar='both', shuffle='both')
        f = function(m) {
            compounds <- sum(rowSums(m > 1) > 0)
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
set.seed(42)
l <- lapply(batches, run_permutations, label, mode, d, N_samples, R)
message('Finished permutations at ',as.character(date()))





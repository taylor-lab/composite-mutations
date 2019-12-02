
require(vegan)
require(data.table)
require(rpkg)
require(here)

clin <- fread(here('data/data_clinical.txt'))
clin <- clin[exclude==F,]
tumortype_tbl <- table.freq(clin$metamaintype)
tumortypes <- tumortype_tbl$value
tumortypes <- tumortypes[tumortypes!='Other']

launch_job <- function(type) {
    N_samples <- tumortype_tbl$N[tumortype_tbl$value==type]
    id <- gsub(' ','_',tolower(type))
    rscript <- '/home/ang46/lab/repos/compound-mutations/r/permute_fixed_density_split.R'
    cmd <- paste0('Rscript ',rscript," '",type,"' ",N_samples)
    mem <- 16
    cores <- 16

    job_id <- paste0('perm_',id)
    out_file <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches/',id,'.out')
    err_file <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches/',id,'.err')

    bsub_cmd <- paste0('bsub -R rusage[mem=',mem,'] -J ',job_id,' -oo ',
                       out_file,' -eo ',err_file,' -We 24:00 -n ',cores,' ',cmd)

    message(bsub_cmd)
    system(bsub_cmd, intern=F, wait=T)
}
trash <- lapply(tumortypes[1:length(tumortypes)], launch_job)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge/save results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clin <- fread(here('data/data_clinical.txt'))
clin <- clin[exclude==F,]
tumortype_tbl <- table.freq(clin$metamaintype)
tumortypes <- tumortype_tbl$value
tumortypes <- tumortypes[tumortypes!='Other']


prev_runs <- dir('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches',full.names=T)
prev_runs <- grep('[.]txt',prev_runs,value=T)
f <- function(file) {
    x <- fread(file)
    x$type <- file
    x$type <- gsub('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches/','',x$type)
    x
}
prev <- lapply(prev_runs, f)
prev <- rbindlist(prev)
f <- function(s) strsplit(s,'_00')[[1]][1]
prev$type <- sapply(prev$type, f)


ctypes <- gsub(' ','_',tolower(tumortypes))
map <- data.table(tumortype=tumortypes,type=ctypes)
prev <- merge(prev, map, by='type', all.x=T)

tbl <- table.freq(prev$tumortype)
donetypes <- tbl$value[tbl$N==10000]
prev <- prev[tumortype %in% donetypes,]

summarize_type <- function(query.tumortype,prev,tumortype_tbl,d) {
    message(query.tumortype)
    N_samples <- tumortype_tbl$N[tumortype_tbl$value==query.tumortype]
    m <- as.data.frame.matrix(xtabs(~Tumor_Sample_Barcode + Hugo_Symbol,data=d[metamaintype==query.tumortype]))
    compounds <- sum(rowSums(m > 1))
    obs <- compounds / N_samples
    res <- prev[tumortype==query.tumortype,] 
    mu <- mean(res$prop)
    qs <- as.numeric(quantile(res$prop,c(0.025,0.975)))
    list(tumortype=query.tumortype,x=compounds,N=N_samples,obs=obs,mu=mu,lwr=qs[1],upr=qs[2])
}

l <- lapply(donetypes, summarize_type, prev, tumortype_tbl, d)
ll <- rbindlist(l)
ll <- ll[order(obs)]
ll$tumortype <- factor(ll$tumortype, levels=rev(ll$tumortype))
write.tsv(ll,'~/lab/repos/compound-mutations/data/cancertype_enrichment_permutation_test.txt')





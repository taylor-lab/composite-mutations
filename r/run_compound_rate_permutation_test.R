
require(vegan)
require(data.table)
require(rpkg)
require(here)


## run cancertype permutation tests with fixed density, but here we actually get the prop compound among samples overall
clin <- fread(here('data/data_clinical.txt'))
clin <- clin[exclude==F,]
tumortype_tbl <- table.freq(clin$metamaintype)
tumortypes <- tumortype_tbl$value
tumortypes <- tumortypes[tumortypes!='Other']

launch_job <- function(type) {
    id <- gsub(' ','_',tolower(type))
    rscript <- '/home/ang46/lab/repos/compound-mutations/r/permute_fixed_density_split.R'
    cmd <- paste0('Rscript ',rscript," '",type,"' NA ")
    mem <- 32
    cores <- 24

    job_id <- paste0('perm_',id)
    out_file <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches_v2/',id,'.out')
    err_file <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches_v2/',id,'.err')

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

prev_runs <- dir('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches_v2',full.names=T)
prev_runs <- grep('[.]txt',prev_runs,value=T)
f <- function(file) {
    x <- fread(file)
    x$type <- file
    x$type <- gsub('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/cancertype_permutations_batches_v2/','',x$type)
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
    compounds <- sum(rowSums(m > 1) > 0)
    obs <- compounds / N_samples
    x <- prev$prop[prev$tumortype==query.tumortype] 
    p.value = sum(x >= obs)/length(x)
    mu <- mean(x)
    qs <- as.numeric(quantile(x,c(0.025,0.975)))
    lor <- log2(obs/mu)
    list(tumortype=query.tumortype,x=compounds,N=N_samples,obs=obs,mu=mu,lwr=qs[1],upr=qs[2],p.value=p.value,logOR=lor)
}

d <- fread(here('data/data_mutations.txt'))
d <- d[exclude==F & putative_resistance_mutation==F,]

l <- lapply(donetypes, summarize_type, prev, tumortype_tbl, d)
ll <- rbindlist(l)
ll <- ll[order(obs)]
ll$tumortype <- factor(ll$tumortype, levels=rev(ll$tumortype))
ll[p.value==0,p.value:=1e-4]
ll$q.value <- p.adjust(ll$p.value,method='BH')
ll[q.value < 0.01]
write.tsv(ll,'~/lab/repos/compound-mutations/data/cancertype_enrichment_permutation_test.txt')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run permutations per TMB
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## run cancertype permutation tests with fixed density, but here we actually get the prop compound among samples overall
tmbs <- 1:40
type <- NA
launch_job <- function(tmb, type) {
    id <- tmb
    rscript <- '/home/ang46/lab/repos/compound-mutations/r/permute_fixed_density_split.R'
    cmd <- paste0('Rscript ',rscript," '",type,"' ",tmb)
    mem <- 16
    cores <- 16

    job_id <- paste0('perm_',id)
    out_file <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/tmb_permutations_batches/',id,'.out')
    err_file <- paste0('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/tmb_permutations_batches/',id,'.err')
    bsub_cmd <- paste0('bsub -R rusage[mem=',mem,'] -J ',job_id,' -oo ',
                       out_file,' -eo ',err_file,' -We 24:00 -n ',cores,' ',cmd)

    message(bsub_cmd)
    system(bsub_cmd, intern=F, wait=T)
}
trash <- lapply(tmbs, launch_job, type)
#trash <- lapply(1, launch_job, type)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge the TMB permutations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prev_runs <- dir('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/tmb_permutations_batches',full.names=T)
prev_runs <- grep('[.]txt',prev_runs,value=T)
f <- function(file) {
    x <- fread(file)
    x$tmb <- file
    x$tmb <- gsub('/home/ang46/lab/projects/compound_mutations/data/processed_data/ifs/reviewer_response/tmb_permutations_batches/tmb-','',x$tmb)
    f <- function(s) strsplit(s,'_')[[1]][1]
    x$tmb <- sapply(x$tmb, f)
    x
}
prev <- lapply(prev_runs, f)
prev <- rbindlist(prev)
prev$tmb <- as.integer(prev$tmb)
tbl <- table.freq(prev$tmb)
donetmbs <- as.integer(tbl$value[tbl$N==10000])
prev <- prev[tmb %in% donetmbs,]

d <- fread(here('data/data_mutations.txt'),select=c('Tumor_Sample_Barcode','Hugo_Symbol'))
clin <- fread(here('data/data_clinical.txt'))

summarize_type <- function(query.tmb,prev,clin,d) { 
    message(query.tmb)
    tmp <- clin[round(tmb) %in% query.tmb,]
    valid_samples <- tmp$Tumor_Sample_Barcode
    N_samples <- length(unique(valid_samples))
    m <- as.data.frame.matrix(xtabs(~Tumor_Sample_Barcode + Hugo_Symbol,data=d[Tumor_Sample_Barcode %in% valid_samples]))
    compounds <- sum(rowSums(m > 1) > 0)
    obs <- compounds / N_samples
    x <- prev$prop[prev$tmb==query.tmb] 
    p.value = sum(x >= obs)/length(x)
    mu <- mean(x)
    qs <- as.numeric(quantile(x,c(0.025,0.975)))
    lor <- log2(obs/mu)
    list(tmb=query.tmb,x=compounds,N=N_samples,obs=obs,mu=mu,lwr=qs[1],upr=qs[2],p.value=p.value,logOR=lor)
}

l <- lapply(donetmbs, summarize_type, prev, clin, d)
ll <- rbindlist(l)
ll[p.value==0,p.value:=1e-4]
ll <- ll[order(obs)]
ll$q.value <- p.adjust(ll$p.value,method='BH')
ll[q.value < 0.01]
write.tsv(ll,'~/lab/repos/compound-mutations/data/tmb_enrichment_permutation_test.txt')




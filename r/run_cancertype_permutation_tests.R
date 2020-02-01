source(here::here('r/prerequisites.R'))

## run cancertype permutation tests with fixed density
clin <- fread(here('data/data_clinical.txt.gz'))
clin <- clin[exclude==F,]
tumortype_tbl <- table.freq(clin$metamaintype)
tumortypes <- tumortype_tbl$value
tumortypes <- tumortypes[tumortypes!='Other']

dir.create(here('data/lsf'),showWarnings=F)
dir.create(here('data/cancertype_permutation_tests'),showWarnings=F)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# launch LSF jobs to run permutation test
# per tumortype. This can take days!
# The while() loop below will merge/save results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

launch_job <- function(type) {
    id <- gsub(' ','_',tolower(type))
    rscript <- here('r/permute_cancertype_data.R')
    cmd <- paste0('Rscript ',rscript," '",type,"'")
    mem <- 32
    cores <- 24
    job_id <- paste0('perm_',id)
    out_file <- paste0(here('data/lsf/',id,'.out'))
    err_file <- paste0(here('data/lsf/',id,'.err'))
    bsub_cmd <- paste0('bsub -R rusage[mem=',mem,'] -J ',job_id,' -oo ',
                       out_file,' -eo ',err_file,' -We 24:00 -n ',cores,' ',cmd)
    message(bsub_cmd)
    system(bsub_cmd, intern=F, wait=T)
}

## load and merge the results, check if any cancer types did not finish
clin <- fread(here('data/data_clinical.txt.gz'))
clin <- clin[exclude==F,]
tumortype_tbl <- table.freq(clin$metamaintype)
tumortypes <- tumortype_tbl$value
tumortypes <- tumortypes[tumortypes!='Other']

prev_runs <- dir(here('data/cancertype_permutation_tests'),full.names=T)
prev_runs <- grep('[.]txt',prev_runs,value=T)
f <- function(file) {
    x <- fread(file)
    x$type <- file
    x$type <- gsub(here('data/cancertype_permutation_tests/'),'',x$type)
    x
}
prev <- lapply(prev_runs, f)
prev <- rbindlist(prev)
f <- function(s) strsplit(s,'_00')[[1]][1]
prev$type <- sapply(prev$type, f)
ctypes <- gsub(' ','_',tolower(tumortypes))
map <- data.table(tumortype=tumortypes,type=ctypes)
prev <- merge(prev, map, by='type', all.x=T)
if(nrow(prev)==0) {
    todo <- tumortypes
} else { 
    tbl <- table.freq(prev$tumortype)
    donetypes <- tbl$value[tbl$N==10000]
    prev <- prev[tumortype %in% donetypes,]
    types_not_finished <- tbl$value[tbl$N < 10000]
    types_not_yet_run <- tumortypes[!tumortypes %in% tbl$value]
    todo <- c(types_not_finished, types_not_yet_run)
}


## If there were any tumortypes which either were not started, or did not finish, run them here.
if(length(todo)>0) {
    msg <- paste('These cancer types did not finish:',paste(sort(todo),collapse=', '),'. Re-launching LSF jobs to finish them.')
    trash <- lapply(todo, launch_job)
    stop('Quitting. Resume when all LSF jobs complete.')
}


## once permutations tests for each cancer type are finished, merge and annotate the results
summarize_type <- function(query.tumortype,prev,tumortype_tbl,d) {
    message(query.tumortype)
    N_samples <- tumortype_tbl$N[tumortype_tbl$value==query.tumortype]
    m <- as.data.frame.matrix(xtabs(~Tumor_Sample_Barcode + Hugo_Symbol,data=d[metamaintype==query.tumortype]))
    composites <- sum(rowSums(m > 1) > 0)
    obs <- composites / N_samples
    x <- prev$prop[prev$tumortype==query.tumortype] 
    p.value = sum(x >= obs)/length(x)
    mu <- mean(x)
    qs <- as.numeric(quantile(x,c(0.025,0.975)))
    lor <- log2(obs/mu)
    list(tumortype=query.tumortype,x=composites,N=N_samples,obs=obs,mu=mu,lwr=qs[1],upr=qs[2],p.value=p.value,logOR=lor)
}

d <- fread(here('data/data_mutations.txt.gz'))
d <- d[exclude==F & putative_resistance_mutation==F,]
l <- lapply(donetypes, summarize_type, prev, tumortype_tbl, d)
ll <- rbindlist(l)
ll <- ll[order(obs)]
ll$tumortype <- factor(ll$tumortype, levels=rev(ll$tumortype))
ll[p.value==0,p.value:=1e-4]
ll$q.value <- p.adjust(ll$p.value,method='BH')
ll[q.value < 0.01]
write.tsv(ll,here('data/cancertype_enrichment_permutation_test.txt'))





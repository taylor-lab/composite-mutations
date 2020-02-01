## This script will test individually recurrent variants in the mutation data
## for enriched incidence as part of a composite mutation.
## output: data/residue_enrichment.txt


source(here::here('r/prerequisites.R'))
cpus <- cpu_prompt()

test_mutation <- function(query_mutation,dd) {
    ## function to test a mutation for enrichment in composite mutations
    
    message(query_mutation)
    query_gene <- strsplit(query_mutation,' ')[[1]][1]
    prop_resistance <- mean(dd$putative_resistance_mutation[dd$tm %in% query_mutation],na.rm=T)
    prop_truncating <- mean(dd$truncating[dd$tm %in% query_mutation],na.rm=T)
    hotspot_class <- paste(sortunique(dd$hotspot_class[dd$tm %in% query_mutation]),collapse=',')
    oncogenic <- paste(sortunique(dd$oncogenic[dd$tm %in% query_mutation]),collapse=',')

    samples_with_composite_mutation_in_gene <- sortunique(dd$Tumor_Sample_Barcode[dd$Hugo_Symbol %in% query_gene & dd$composite==T])
    samples_with_composite_mutation_including_residue <- sortunique(dd$Tumor_Sample_Barcode[dd$tm %in% query_mutation & dd$composite==T])
    samples_with_composite_mutation_not_including_residue <- samples_with_composite_mutation_in_gene[samples_with_composite_mutation_in_gene %nin% samples_with_composite_mutation_including_residue] 

    samples_with_singleton_mutation_in_gene <- sortunique(dd$Tumor_Sample_Barcode[dd$Hugo_Symbol %in% query_gene & dd$composite==F])
    samples_with_singleton_mutation_including_residue <- sortunique(dd$Tumor_Sample_Barcode[dd$tm %in% query_mutation & dd$composite==F])
    samples_with_singleton_mutation_not_including_residue <- samples_with_singleton_mutation_in_gene[samples_with_singleton_mutation_in_gene %nin% samples_with_singleton_mutation_including_residue] 

    n.samples_with_singleton_mutation_including_residue = length(unique(strsplit(samples_with_singleton_mutation_including_residue,9)))
    n.samples_with_singleton_mutation_not_including_residue = length(unique(strsplit(samples_with_singleton_mutation_not_including_residue,9)))
    n.samples_with_composite_mutation_including_residue = length(unique(strsplit(samples_with_composite_mutation_including_residue,9)))
    n.samples_with_composite_mutation_not_including_residue = length(unique(strsplit(samples_with_composite_mutation_not_including_residue,9)))

    m <- rbind(c(n.samples_with_singleton_mutation_not_including_residue,n.samples_with_singleton_mutation_including_residue),
               c(n.samples_with_composite_mutation_not_including_residue,n.samples_with_composite_mutation_including_residue))

    ## two-sided p-values for text
    tst <- fisher.test(m,alternative='two.sided')
    p.value.two.sided <- tst$p.value
    OR <- tst$estimate
    
    ## enrichment p-value for plot
    tst <- fisher.test(m,alternative='greater')
    p.value.enriched <- tst$p.value

    cis_samples <- dd$Tumor_Sample_Barcode[dd$cis==T & dd$tm %in% query_mutation]
    trans_samples <- dd$Tumor_Sample_Barcode[dd$trans==T & dd$tm %in% query_mutation]
    trans_or_separate_cells_samples <- dd$Tumor_Sample_Barcode[dd$trans_or_separate_cell==T & dd$tm %in% query_mutation]
    ambiguous_samples <- dd$Tumor_Sample_Barcode[dd$ambiguous==T & dd$tm %in% query_mutation]
    unknown_samples <- dd$Tumor_Sample_Barcode[dd$unknown==T & dd$tm %in% query_mutation]
    not_phased_samples <- dd$Tumor_Sample_Barcode[dd$not_phased==T & !is.na(dd$not_phased) & dd$tm %in% query_mutation]

    n_cis <- 0; n_trans <- 0; n_ambiguous <- 0; n_trans_or_separate_cells <- 0; n_unknown <- 0; n_not_phased <- 0
    if(length(cis_samples)>0) n_cis <- length(unique(strtrim(cis_samples,9)))
    if(length(trans_samples)>0) n_trans <- length(unique(strtrim(trans_samples,9)))
    if(length(ambiguous_samples)>0) n_ambiguous <- length(unique(strtrim(ambiguous_samples,9)))
    if(length(trans_or_separate_cells_samples)>0) n_trans_or_separate_cells <- length(unique(strtrim(trans_or_separate_cells_samples,9)))
    if(length(unknown_samples)>0) n_unknown <- length(unique(strtrim(unknown_samples,9)))
    if(length(not_phased_samples)>0) n_not_phased <- length(unique(strtrim(not_phased_samples,9)))

    list(query_gene=query_gene,
         query_mutation=query_mutation,
         n.samples_with_singleton_mutation_including_residue=n.samples_with_singleton_mutation_including_residue,
         n.samples_with_singleton_mutation_not_including_residue=n.samples_with_singleton_mutation_not_including_residue,
         n.samples_with_composite_mutation_including_residue=n.samples_with_composite_mutation_including_residue,
         n.samples_with_composite_mutation_not_including_residue=n.samples_with_composite_mutation_not_including_residue,
         OR=OR,
         p.value.two.sided=p.value.two.sided,
         p.value.enriched=p.value.enriched,
         prop_resistance=prop_resistance,
         prop_truncating=prop_truncating,
         hotspot_class=hotspot_class,
         oncogenic=oncogenic,
         n_cis=n_cis, n_trans=n_trans, n_ambiguous=n_ambiguous, n_trans_or_separate_cells=n_trans_or_separate_cells, 
         n_unknown=n_unknown, n_not_phased=n_not_phased)    
}


## load/format mutation data
d_mutations <- fread(here('data/data_mutations.txt.gz'),
                     select=c('Tumor_Sample_Barcode','exclude','Hugo_Symbol','putative_resistance_mutation',
                              'Variant_Classification','Variant_Type','tcn', 'mutationid', 'high_tmb','hotspot_class',
                              'tm','hotspotid','HGVSp_Short','truncating','Start_Position'))
d <- d_mutations[exclude==F & putative_resistance_mutation==F]
d[hotspotid!='',tm:=hotspotid] ## in cases of INDEL hotspots, group mutation IDs together into the hotspot cluster
d[Variant_Classification=='TERT promoter', tm:=paste('TERT',gsub('p[.]','',HGVSp_Short))]
d <- d[,c('Tumor_Sample_Barcode','Hugo_Symbol','tm','putative_resistance_mutation','mutationid','exclude','truncating','hotspotid','Variant_Classification','Start_Position'),with=F]


## load/format phasing data for annotation the results
d_phased <- fread(here('data/data_mutations_phased.txt.gz'))
d_phased <- d_phased[exclude==F & putative_resistance_mutation.1==F & putative_resistance_mutation.2==F]
d_phased$id <- paste(d_phased$mutationid.1,d_phased$mutationid.2,sep=' + ')
d_phased <- d_phased[!duplicated(id),]
phased <- d_phased
phased[hotspotid.1!='',tm.1:=hotspotid.1]
phased[hotspotid.2!='',tm.2:=hotspotid.2]
phased[Variant_Classification.1=='TERT promoter', tm.1:=paste('TERT',gsub('p[.]','',HGVSp_Short.1))]
phased[Variant_Classification.2=='TERT promoter', tm.2:=paste('TERT',gsub('p[.]','',HGVSp_Short.2))]
phased <- phased[,c('mutationid.1','mutationid.2','phase'),with=F]
phased <- melt(phased,id.var='phase')
phased[,variable:=NULL]
setnames(phased,'value','mutationid')



## here we annotate the rate of cis/trans composites with each mutation in the table
## annotate phase of all mutations (NB: a single mutation in a sample with 
## 3+ mutations in composite can have multiple phases)
annotate_phases_of_variant <- function(query.mutationid,phased) {
    phased <- phased[mutationid==query.mutationid,]
    if(nrow(phased) > 0) {
        cis <- any(phased$phase=='cis',na.rm=T)
        trans <- any(phased$phase=='trans',na.rm=T)
        trans_or_separate_cells <- any(phased$phase=='trans_or_separate_cells',na.rm=T)
        ambiguous <- any(phased$phase=='ambiguous: 3+ cis, 3+ trans',na.rm=T)
        unknown <- any(phased$phase=='unknown',na.rm=T)
        not_phased <- F
    } else {
        cis <- F
        trans <- F
        trans_or_separate_cells <- F
        ambiguous <- F
        unknown <- F
        not_phased <- T
    }
    list(mutationid=query.mutationid, cis=cis, trans=trans, trans_or_separate_cells=trans_or_separate_cells, ambiguous=ambiguous, unknown=unknown, not_phased=not_phased)
}

message('Annotating phase info for mutations ...')
l <- mclapply(d$mutationid, annotate_phases_of_variant, phased, mc.cores=cpus)
ll <- rbindlist(l)
d <- merge(d, ll, by='mutationid', all.x=T)
d$tsbgene <- paste(d$Tumor_Sample_Barcode,d$Hugo_Symbol)

## annotate composites (including the phased and not-phased)
tbl <- table.freq(d$tsbgene)
d <- merge(d, tbl, by.x='tsbgene', by.y='value', all.x=T)
d$composite <- d$N > 1
d[,N:=NULL]

## annotate non-composites as NA for not phased
d[composite==F,not_phased:=NA]
d <- d[order(Hugo_Symbol, Start_Position),]




## test SNP residues and in-frame INDELS mutated in 5+ patients
message('Testing residues, in-frame indels and TERT promoter alleles mutated in 5+ patients for enrichment ...')
tmp <- d[,c('tm','Tumor_Sample_Barcode'),with=F]
tmp$pid <- strtrim(tmp$Tumor_Sample_Barcode,9)
summarize_tm <- function(d) {
    n_pts <- length(unique(d$pid))
    list(n_pts=n_pts)
}
info <- tmp[,summarize_tm(.SD),by=tm]
test_mutations <- info$tm[info$n_pts >= 5]
l <- mclapply(test_mutations, test_mutation, d, mc.cores=cpus)
ll <- rbindlist(l)
ll <- ll[order(p.value.enriched,decreasing=F),]
ll$q.value.enriched <- p.adjust(ll$p.value.enriched,method='BH')
ll$q.value.two.sided <- p.adjust(ll$p.value.two.sided,method='BH')
ll <- merge(ll, d[!duplicated(tm),c('tm','hotspotid'),with=F], by.x='query_mutation', by.y='tm', all.x=T)

## save result
write.tsv(ll,here('data/residue_enrichment.txt'))
message('Done!')





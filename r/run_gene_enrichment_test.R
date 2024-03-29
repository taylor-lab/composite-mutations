## this script will test for genes significantly enriched for composite mutations
## output: data/gene_enrichment.txt

source(here::here('r/prerequisites.R'))

## load mutation data
d_mutations <- fread(here('data/data_mutations.txt.gz'),
                     select=c('Tumor_Sample_Barcode','exclude','Hugo_Symbol','Role','putative_resistance_mutation',
                              'Variant_Classification','Variant_Type','tcn', 'mutationid', 'high_tmb','hotspot_class',
                              'tm','hotspotid','HGVSp_Short','truncating','Reference_Allele','Tumor_Seq_Allele2',
                              'Start_Position','End_Position'))


## load gene annotations including covariates
gene_info <- fread(here('data/gene_info.txt'),
                   select=c('Hugo_Symbol','Role','panel_introduced','reptime',
                            'percentage_gene_gc_content','cds_length','hic'))


## merge gene annotations to mutation data, drop excluded samples
d <- d_mutations
d <- d[exclude==F & putative_resistance_mutation==F,]
d <- merge(d, gene_info, by='Hugo_Symbol', all.x=T)


## consider the TERT promoter as a unique gene, where the cds_length is now the total length of the promoter region covered by IMPACT
tert_promoter <- which(d$Hugo_Symbol=='TERT' & d$Variant_Classification=='TERT promoter')
d[tert_promoter, Hugo_Symbol:='TERT promoter']
d[tert_promoter, cds_length:=499] ## this is the length of the promoter region of TERT covered by IMPACT468
d[tert_promoter, percentage_gene_gc_content:=0.7896] ## this is the %GC content of TERT promoter region


## get observed number of composite mutated samples and mutated samples per gene
get_composite_rate_per_gene <- function(d) {

    ## get sample barcodes for singletons/composites/overall
    tbl <- table.freq(d$Tumor_Sample_Barcode)
    samples <- tbl$value
    composite_samples <- tbl$value[tbl$N > 1]
    singleton_samples <- tbl$value[tbl$N == 1]

    ## don't overcount multiple samples per patient
    samples <- length(unique(strtrim(samples, 9)))
    composites <- length(unique(strtrim(composite_samples, 9)))
    singletons <- length(unique(strtrim(singleton_samples, 9)))

    ## calculate the average TCN to include as covariate for this gene
    tcn <- mean(d$tcn,na.rm=T) 
    list(
         singletons=singletons,composites=composites,samples=samples,
         tcn=tcn,reptime=unique(d$reptime),percentage_gene_gc_content=unique(d$percentage_gene_gc_content),cds_length=unique(d$cds_length),
         hic=unique(d$hic),
         panel_introduced=unique(d$panel_introduced))
}
res <- d[,get_composite_rate_per_gene(.SD),by=Hugo_Symbol]


## impute missing values
res[is.nan(tcn),tcn:=2]
res$hic[is.na(res$hic)] <- round(mean(res$hic,na.rm=T))


## prep for test
res$panel_introduced <- factor(res$panel_introduced)
res$Hugo_Symbol <- gsub('-','_',res$Hugo_Symbol)
res$Hugo_Symbol <- gsub(' ','_',res$Hugo_Symbol)


## use negative-binomial regression to model the predicted number of composite-mutant samples per gene by chance
m <- glm.nb(composites ~ offset(log(samples)) + reptime + cds_length + percentage_gene_gc_content + panel_introduced + tcn + hic, data=res)
predictions <- predict(m,newdata=res,type='response')
res$predicted_proportion <- predictions / res$samples
res$observed_proportion <- res$composites / res$samples
dat <- res[,c('Hugo_Symbol','composites','samples','predicted_proportion','observed_proportion'),with=F]


## test each gene for enriched or depleted composite mutant samples
test_gene <-function(dat) {
    x <- tryCatch({
        tst <- binom.test(dat$composites,dat$samples,dat$predicted_proportion,alternative='greater')
        p_enriched <- tst$p.value
        tst2 <- binom.test(dat$composites,dat$samples,dat$predicted_proportion,alternative='two.sided')
        p_twosided <- tst2$p.value
        list(p_enriched=p_enriched, p_twosided=p_twosided)
    },error=function(e) {
        list(p_enriched=as.numeric(NA), p_twosided=as.numeric(NA))
    })
    dat$p_enriched <- x$p_enriched
    dat$p_twosided <- x$p_twosided
    dat
}
dat <- dat[,test_gene(.SD),by=Hugo_Symbol]
dat$logOR <- log2(dat$observed_proportion/dat$predicted_proportion)
x <- gene_info[,c('Hugo_Symbol','Role'),with=F]
x$Hugo_Symbol <- gsub('-','_',x$Hugo_Symbol)
dat <- merge(dat, x,by='Hugo_Symbol', all.x=T)
dat <- dat[order(p_enriched,decreasing=F),]
dat$q_enriched <- p.adjust(dat$p_enriched,method='BH')
dat$q_twosided <- p.adjust(dat$p_twosided,method='BH')

## save results
write.tsv(dat, here('data/gene_enrichment.txt'))




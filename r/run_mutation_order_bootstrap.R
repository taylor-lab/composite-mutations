source(here::here('r/prerequisites.R'))
require(boot)


## load starting datasets
d_mutations <- fread(here('data/data_mutations.txt'),select=c('Tumor_Sample_Barcode','PATIENT_ID','exclude','Hugo_Symbol','Role','putative_resistance_mutation',
                                                              'Variant_Classification','Variant_Type','tcn', 'mutationid', 'high_tmb','hotspot_class',
                                                              'tm','hotspotid','HGVSp_Short','truncating','Reference_Allele','Tumor_Seq_Allele2',
                                                              'Start_Position','End_Position'))
d_clinical <- fread(here('data/data_clinical.txt'))
d_phased <- fread(here('data/data_mutations_phased.txt'))

## exclude duplicates from phasing data
d_phased$id <- paste(d_phased$mutationid.1,d_phased$mutationid.2,sep=' + ')
d_phased <- d_phased[!duplicated(id),]


get_ordered_variants <- function(d_phased, d_mutations) {
    d_phased <- d_phased[putative_resistance_mutation.1==F & putative_resistance_mutation.2==F & exclude==F,]
    d_phased[hotspotid.1!='',tm.1:=hotspotid.1] ## this will group IF-INDEL hotspots
    d_phased[hotspotid.2!='',tm.2:=hotspotid.2]
    d_phased$functional.1 <- d_phased$hotspotid.1!='' | d_phased$oncogenic.1 %in% c('Oncogenic','Likely Oncogenic','Predicted Oncogenic') |
    (d_phased$Role %in% c('TSG','Oncogene/TSG') & d_phased$truncating.1==T) | d_phased$Variant_Classification.1=='TERT promoter'
    d_phased$functional.2 <- d_phased$hotspotid.2!='' | d_phased$oncogenic.2 %in% c('Oncogenic','Likely Oncogenic','Predicted Oncogenic') |
    (d_phased$Role %in% c('TSG','Oncogene/TSG') & d_phased$truncating.2==T) | d_phased$Variant_Classification.2=='TERT promoter'

    d_phased[Variant_Classification.1 %in% c('TERT promoter'),tm.1:=paste0('TERT ',gsub('p[.]','',HGVSp_Short.1))]
    d_phased[Variant_Classification.1 %in% c('Splice_Site','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins'), tm.1:=paste(Hugo_Symbol,'Truncating')]
    d_phased[Variant_Classification.2 %in% c('TERT promoter'),tm.2:=paste0('TERT ',gsub('p[.]','',HGVSp_Short.2))]
    d_phased[Variant_Classification.2 %in% c('Splice_Site','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins'), tm.2:=paste(Hugo_Symbol,'Truncating')]

    d_phased <- d_phased[,c('Tumor_Sample_Barcode','Role','Hugo_Symbol','hotspot.1','hotspot.2','functional.1','functional.2',
              'tm.1','tm.2','phase','common_reads_alt1_alt2',
              'common_reads_alt1_ref2','common_reads_ref1_alt2','ccf_Mcopies_lower.1','ccf_Mcopies_upper.1',
              'ccf_Mcopies_lower.2','ccf_Mcopies_upper.2','Variant_Classification.1','Variant_Classification.2','phase','same_cell_known'),with=F]

    d_phased$firstallele <- ''
    d_phased$secondallele <- ''
    d_phased[ccf_Mcopies_lower.1 > ccf_Mcopies_upper.2,firstallele:=tm.1]
    d_phased[ccf_Mcopies_lower.1 > ccf_Mcopies_upper.2,secondallele:=tm.2]
    d_phased[ccf_Mcopies_upper.1 < ccf_Mcopies_lower.2,firstallele:=tm.2]
    d_phased[ccf_Mcopies_upper.1 < ccf_Mcopies_lower.2,secondallele:=tm.1]
    dd <- d_phased[,c('Tumor_Sample_Barcode','Role','Hugo_Symbol','firstallele','secondallele','phase','same_cell_known','hotspot.1','hotspot.2',
               'functional.1','functional.2','Variant_Classification.1','Variant_Classification.2'),with=F]
    dd <- dd[firstallele!='' & secondallele!='',]

    ## get number of unique patients with each mutation
    d_mutations <- d_mutations[exclude==F & putative_resistance_mutation==F,]
    d_mutations$functional <- d_mutations$hotspot==T | d_mutations$oncogenic %in% c('Oncogenic','Likely Oncogenic','Predicted Oncogenic') |
    (d_mutations$Role %in% c('TSG','Oncogene/TSG') & d_mutations$truncating==T) | d_mutations$Variant_Classification=='TERT promoter'
    d_mutations[Variant_Classification %in% c('TERT promoter'),tm:=paste0('TERT ',gsub('p[.]','',HGVSp_Short))]
    d_mutations[Variant_Classification %in% c('Splice_Site','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins'), tm:=paste(Hugo_Symbol,'Truncating')]
    d_mutations[hotspotid!='',tm:=hotspotid] ## this will group IF-INDEL hotspots

    getN <- function(d_mutations) {
        N <- length(unique(strtrim(d_mutations$Tumor_Sample_Barcode,9)))
        list(N=N)
    }
    tbl <- d_mutations[,getN(.SD),by=tm]
    setnames(tbl,'tm','value')

    ## annotate first and second hits with their prevalence across full dataset
    dd2 <- merge(dd, tbl, by.x='firstallele', by.y='value', all.x=T)
    setnames(dd2,'N','firstallele_N')
    dd2 <- merge(dd2, tbl, by.x='secondallele', by.y='value', all.x=T)
    setnames(dd2,'N','secondallele_N')
    tbl_gene <- table.freq(d_mutations$Hugo_Symbol)
    dd2 <- merge(dd2, tbl_gene, by.x='Hugo_Symbol', by.y='value', all.x=T)
    setnames(dd2,'N','totalN')

    ## compare prevalence of each hit for signficantly more prevalent:
    first_ci <- binom::binom.confint(dd2$firstallele_N, dd2$totalN, method='exact')
    second_ci <- binom::binom.confint(dd2$secondallele_N, dd2$totalN, method='exact')
    dd2$prop_samples_with_first_allele <- first_ci$mean
    dd2$prop_samples_with_first_allele_lwr <- first_ci$lower
    dd2$prop_samples_with_first_allele_upr <- first_ci$upper
    dd2$prop_samples_with_second_allele <- second_ci$mean
    dd2$prop_samples_with_second_allele_lwr <- second_ci$lower
    dd2$prop_samples_with_second_allele_upr <- second_ci$upper
    dd2$more_common_allele <- 'no difference'
    dd2[prop_samples_with_first_allele_lwr > prop_samples_with_second_allele_upr, more_common_allele:='first']
    dd2[prop_samples_with_first_allele_upr < prop_samples_with_second_allele_lwr, more_common_allele:='second']
    dd2$Hugo_Symbol <- as.character(dd2$Hugo_Symbol)
    dd2$firstallele <- as.character(dd2$firstallele)
    dd2$secondallele <- as.character(dd2$secondallele)
    dd2
}

## get contingency table of 1st/2nd more prevalent vs Oncogene/TSG
dat <- get_ordered_variants(d_phased, d_mutations)
dat <- dat[same_cell_known==T | phase=='cis',] ## only consider compounds definitely in same cell
dat$n_functional <- dat$functional.1 + dat$functional.2
dat <- dat[n_functional == 2,]
dat <- dat[more_common_allele!='no difference',]
m <- as.data.frame.matrix(xtabs(~Role + more_common_allele, data=dat[Role %in% c('Oncogene','TSG')]))
m$overall <- m$first + m$second

## bootstap the ratio to get CIs
getratio <- function(dat,index) {
    dat <- dat[index,]
    m <- as.data.frame.matrix(xtabs(~Role + more_common_allele, data=dat[Role %in% c('Oncogene','TSG')]))
    m$overall <- m$first + m$second
    l <- m$first/m$overall
    names(l) <- c('first_onc','first_tsg')
    l
}
set.seed(42)
R <- 1e5
b <- boot(data=dat, statistic=getratio, R=R, parallel='multicore', ncpus=12)

## format the resulting data
onc_mu <- m$first[1]/m$overall[1]
onc_p <- sum(0.5 >= b$t[,1])/R;
onc_ci <- as.numeric(quantile(b$t[,1],c(0.025,0.975)))
tsg_mu <- m$first[2]/m$overall[2]
tsg_p <- sum(0.5 >= b$t[,2])/R;
tsg_ci <- as.numeric(quantile(b$t[,2],c(0.025,0.975)))
dat <- data.table(role=c('Oncogene','TSG'))
dat$mu <- c(onc_mu, tsg_mu)
dat$lwr_boot <- c(onc_ci[1], tsg_ci[1])
dat$upr_boot <- c(onc_ci[2], tsg_ci[2])
setnames(dat,'mu','prop_first')
dat$prop_second <- 1-dat$prop_first

## add in binomial CIs for comparison
binomialCIs <- binom::binom.confint(x=m[,1],n=m[,3],method='exact')
dat$lwr_binom <- binomialCIs$lower
dat$upr_binom <- binomialCIs$upper
dat <- dat[,c('role','prop_first','prop_second','lwr_boot','upr_boot','lwr_binom','upr_binom'),with=F]

## save resulting data
write.tsv(dat,here('data/rebuttal_only/mutation_order_bootstrapped_precalculated.txt'))




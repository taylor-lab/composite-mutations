require(data.table)
require(ggplot2)
require(cowplot)
require(RColorBrewer)
require(parallel)
require(ggsignif)
require(binom)
require(scales)
require(MASS)
require(ggrepel)
require(Hmisc)
require(Rcpp)
require(here)

important_classes <- c('Missense_Mutation','Nonsense_Mutation','Splice_Site','In_Frame_Del','In_Frame_Ins','Frame_Shift_Del','Frame_Shift_Ins','TERT promoter','Translation_Start_Site')
truncating_classes <- c('Nonsense_Mutation','Splice_Site','Frame_Shift_Del','Frame_Shift_Ins')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define commonly used colors for figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

colors <- list()

## 2-level enrichment (significant, NS)
sig2 <- c('#e50000','#a6a6a6')
names(sig2) <- c('Enriched','Not significant')
colors$sig2 <- sig2

## 3-level enrichment (enriched/NS/depleted)
sig3 <- c('#B2182B','#BFBFBF','#2166AC')
names(sig3) <- c('Enriched','Not significant','Depleted')
colors$sig3 <- sig3

## compounds
compound <- c('#cccccc','#005b96') ## gray / blue
names(compound) <- c('Singleton','Compound')
colors$compound <- compound

## expected vs observed
observed <- c('#003366','#bea893')
names(observed) <- c('Observed','Expected')
colors$observed <- observed

## 2-level phase
phase2 <- c('#CE1256','#0057e7')
names(phase2) <- c('Cis','Trans')
colors$phase2 <- phase2

## oncogene/TSG/na 
role <- c('#5a7f48','#91b4ff','#D9D9D9','black')
names(role) <- c('Oncogene','TSG','Oncogene/TSG','n/a')
colors$role <- role

## 2-color bar plot
barplot2col <- c('#000000','#D9D9D9')
colors$barplot2col <- barplot2col


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define helper functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ggplot theme
theme_std <- function(base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    #theme_classic(base_size = base_size, base_family = 'ArialMT', base_line_size = base_line_size, base_rect_size = base_rect_size)  %+replace%
    theme_classic(base_size = base_size, base_family = 'ArialMT')  %+replace%
    theme(
          line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "round"),
          text = element_text(family = 'ArialMT', face = "plain",
                              colour = "black", size = base_size, lineheight = 0.9,
                              hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug=F),
          axis.text = element_text(colour = "black", family='ArialMT', size=rel(0.8)),
          axis.ticks = element_line(colour = "black", size=rel(1)),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size = rel(1)),
          legend.key = element_blank(),
          strip.background = element_blank())
}


## function to force axis break in ggplot
break_axis <- function(y, maxlower, minupper=NA, lowerticksize, upperticksize, ratio_lower_to_upper) {
    if(is.na(minupper)) {
        breakpos <- maxlower
        lowerticklabels <- seq(0,breakpos,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(breakpos+upperticksize,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- breakpos + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > breakpos
        newy[ind] <- breakpos + uppertickspacing*((newy[ind]-breakpos) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    } else {
        lowerticklabels <- seq(0,maxlower,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(minupper,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- maxlower + 0.5*lowerticksize + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > maxlower
        newy[ind] <- maxlower + 0.5*lowerticksize + 1*uppertickspacing + uppertickspacing*((newy[ind]-minupper) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    }
}


## function to convert ggplot object into separate objects for plot vs legend
extract_gglegend <- function(p){
    require(ggplot2)
    require(cowplot)

    ## extract the legend from a ggplot object
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if(length(leg) > 0) leg <- tmp$grobs[[leg]]
    else leg <- NULL
    leg

    ## return the legend as a ggplot object
    legend <- cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(leg))
    plot <- p + theme(legend.position='none')
    list(plot=plot,legend=legend)
}


## function to split HGVSp_Short field in MAF into Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, Amino_acids
HGVSp_Short_parse <- function(x,cpus=1) {
    require(parallel)

    x1 <- gsub('p.','',x)
    m <- gregexpr('[0-9]+',x1)
    nums <- regmatches(x1,m)
    get.aa <- function(num) num[1]
    aas <- sapply(nums,get.aa,USE.NAMES=F)
    tmp <- data.table(HGVSp_Short=x1,aa=aas)
    tmp$i <- 1:nrow(tmp)

    split <- function(d) {
        s <- strsplit(d$HGVSp_Short,d$aa)[[1]]
        s[2] <- gsub('_sice','splice',s[2])
        list(Reference_Amino_Acid=s[1],Variant_Amino_Acid=s[2])
    }
    info <- tmp[,split(.SD),by=i]
    info$HGVSp_Short <- x
    info$Amino_Acid_Position=as.integer(aas)
    info <- info[,c(4,2,5,3),with=F]
    info$Amino_acids <- paste(info$Reference_Amino_Acid,info$Variant_Amino_Acid,sep='/')
    info$Amino_acids[is.na(info$Amino_Acid_Position)] <- NA
    info
}


## creates a data.table with histogram of values in a vector 
table.freq <- function(value) {
    if(is.null(value) | length(value)==0) {
        tbl <- data.table(value=NA,N=NA)
    } else {
        tbl <- adt(table(value))
        tbl <- tbl[order(tbl$N,decreasing=T),]
    }
    tbl
}


## shortcut for sort(unique(...))
sortunique <- function(x,...) {
        sort(unique(x),na.last=T,...)
}


## generate a qq-plot in ggplot2
gg_qqplot <- function(xs, labels=NULL, ci=0.95, show.lambda=T) {

    # inspired by gg_qqplot.R
    # Kamil Slowikowski
    # February 16, 2014

    ## remove NAs from pvalues
    xs <- xs[!is.na(xs)]

    ## generate qqplot
    N = length(xs)
    df = data.frame(observed=-log10(sort(xs)),
                    expected=-log10(1:N / N),
                    p=sort(xs),
                    cupper=-log10(qbeta(ci,     1:N, N - 1:N + 1)),
                    clower=-log10(qbeta(1 - ci, 1:N, N - 1:N + 1)))
    if(!is.null(labels)) df$label <- labels

    ## make plot
    log10Pe = expression(paste("Expected -log"[10], plain(P)))
    log10Po = expression(paste("Observed -log"[10], plain(P)))
    p <- ggplot(df,aes(x=expected,y=observed)) +
    geom_point(pch=21, size=2, color='black') +
    geom_abline(intercept=0, slope=1, alpha=0.5) +
    geom_line(aes(expected, cupper), linetype=2) +
    geom_line(aes(expected, clower), linetype=2) +
    xlab(log10Pe) +
    ylab(log10Po)

    if(!is.null(labels)) p <- p + geom_text_repel(aes(label=label), color='blue',size=4, na.rm=T)
    if(show.lambda) {
        ## calculate lambda
        set.seed(1234)
        chisq <- qchisq(1 - xs, 1)
        lambda <- paste('lambda =',round(median(chisq) / qchisq(0.5, 1),4))
        message(lambda)
        lambda_dat <- data.table(x=0.20*max(df$expected),y=0.93*max(df$observed),label=lambda)
        p <- p + geom_text(data=lambda_dat,aes(x=x,y=y,label=label),color='red',size=6)
    }
    p
}

## shortcut to write tab-delimited data with consistent format 
write.tsv <- function(d, file, sep = "\t", quote = F, row.names = F, ...) {
    write.table(d, file = file, sep = sep, quote = quote, row.names = row.names, ...)
}

## shortcut for as.data.table
adt <- function(d) as.data.table(d)

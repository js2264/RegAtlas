#!/usr/bin/env Rscript

# Load full data
load('./shared/data/classification_tissue-spe-genes-REs.RData')
# Load required libraries
library(GenomicRanges)
library(magrittr)
# Load custom functions
source("./shared/bin/R/custom_R_functions.R")
#
all$regulatory_class <- readr::read_tsv('./shared/data/reg_elements_dev_ageing_tissues_20190824.tsv')$annot_dev_age_tissues
all.deconv <- cbind(
    all, 
    ATAC, 
    max.tissue.df[, grepl('max.tissue|ratio', colnames(max.tissue.df))]
) %>% 
    dplyr::mutate(uniqueWormBaseID = strsplit(as.character(WormBaseID),',')) %>% 
    tidyr::unnest(uniqueWormBaseID)
atac.dt <- cbind(
    all[, c('chr','start','stop','gene_name','regulatory_class','which.tissues')], 
    round(ATAC, 3), 
    max.tissue.df[, grepl('max.tissue$', colnames(max.tissue.df))]
)
row.names(atac.dt) <- all$coords
colnames(atac.dt)[colnames(atac.dt) == 'gene_name'] <- "geneID"
colnames(atac.dt)[grep(paste(order.tissues, collapse = '|'), colnames(atac.dt))] <- paste0(
    order.tissues[1:length(grep(paste(order.tissues, collapse = '|'), colnames(atac.dt)))], 
    '_TPM'
)
atac.dt$regulatory_class <- factor(atac.dt$regulatory_class)
lcap.dt <- cbind(
    as.data.frame(genes.gtf)[,c(1:3,5,10,11,18)], 
    round(LCAP, 3),
    max.tissue.df.LCAP[, grepl('max.tissue$', colnames(max.tissue.df.LCAP))]
)
colnames(lcap.dt)[1:3] <- c('chr', 'start', 'stop')
colnames(lcap.dt)[colnames(lcap.dt) == 'gene_id'] <- "WormBaseID"
colnames(lcap.dt)[colnames(lcap.dt) == 'gene_name'] <- "geneID"
genes.gtf2 <- as.data.frame(genes.gtf)
# Save data
save(file = "./shared/data/minimal-data.RData", list = c(
    # color/names vars
    'color.tissues',
    'order.tissues',
    # data vars
    'mat.proms.gene',
    'cao03',
    'ATAC',
    'all',
    'max.tissue.df',
    'LCAP',
    'LCAP_normalized',
    'max.tissue.df.LCAP',
    'genes.gtf',
    'genes.gtf2',
    'convID',
    'LCAPdev',
    'all.deconv', 
    'atac.dt', 
    'lcap.dt',
    # Functions
    'getGeneInfos',
    'name2WB',
    'WB2name',
    'fetchWBinfos'
    )
)
# Quit
q('no')

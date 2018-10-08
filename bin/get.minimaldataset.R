#!/usr/bin/env Rscript

#=============================================================================#
#                                                                             #
#     USAGE: Run from the app folder (./dashboard.Ahringer/releases/dev/)     #
#     This script takes data from ../../_tissuesTSNE and minimize it          #
#                                                                             #
#=============================================================================#

message('\n\n-- Getting minimal datasets object --\n\n')

load('~/Documents/PhD/__Bioinfo/_tissuesTSNE/.tSNE.RData')
library(GenomicRanges)
source("~/Documents/PhD/__Bioinfo/_scripts/bin/useful_R_functions.R")

valid.groups = c(1:33)
CLASSES.proms <- colortransp[factor(order.tissues[all[which(all$is.prom & all$clustersATAC.tissues %in% valid.groups),]$clustersATAC.tissues], levels = order.tissues)]
all.proms.valid <- all[all$is.prom & all$clustersATAC.tissues %in% valid.groups,]
valid.groups = order.tissues[c(1:5, 33)]
CLASSES.genes = colortransp[c(1:5, 33)][factor(order.tissues[c(1:5, 33)][max.tissue.df.LCAP[genes.gtf$is.prot.cod & max.tissue.df.LCAP$which.tissues %in% valid.groups,]$which.tissues], levels = order.tissues[c(1:5, 33)])]
genes.valid <- max.tissue.df.LCAP$which.tissues %in% order.tissues[1:33]

save(file = "data/tSNE.minimal.RData", list = c(
    'colortransp',
    'color.tissues',
    'order.tissues',

    'mat.proms.gene',
    'cao03',
    'ATAC',
    'all',
    'max.tissue.df',
    'LCAP',
    'max.tissue.df.LCAP',
    'genes.gtf',
    'SUMM',
    'convID',

    'tSNE.df.proms',
    'tSNE.df.LCAP',
    'LCAPdev',
    'ATACdev',
    'CLASSES.proms',
    'all.proms.valid',
    'CLASSES.genes',
    'genes.valid',

    'getGeneInfos',
    'name2WB',
    'WB2name',
    'fetchWBinfos'
    )
)

q('no')

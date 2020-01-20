#!/usr/bin/env Rscript

library(optparse)
option_list = list(
  make_option("--input", type="character", default = '',
              help="Input bed file"),
  make_option("--length", type="integer", default = '',
              help="Final size of genomic loci"),
  make_option("--sizes", type="character", default = '',
              help="chrom.sizes file", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#
#
#

input_file <- opt$input
length <- opt$length
chromsize_file <- opt$sizes

#
#
#

require(tidyverse)
require(rtracklayer)

g <- import(input_file)
seqinfo <- read.table(chromsize_file, stringsAsFactor = FALSE) %>% filter(V1 %in% seqlevels(g)) %>% '['(match(seqlevels(g), .$V1), )
seqinfo <- Seqinfo(seqnames = seqinfo$V1, seqlengths = seqinfo$V2, isCircular = rep(F, nrow(seqinfo)), genome = rep(NA, nrow(seqinfo)))
seqinfo(g) <- seqinfo

#
#
#

g2 <- resize(g, fix = 'center', width = length) %>% trim()
export(g2, gsub('.bed', paste0('.l', length, '.bed'), input_file))

q('no')
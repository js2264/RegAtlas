#!/usr/bin/env Rscript

# This script reads metadata from the 2nd tab of a metadata excel datasheet containing the info on which sequencing runs to pool together for each replicate of each ATAC seq sample.
# It generates a 'combined' fastq file of all the re-sequenced replicates.

require(readxl)
dir.create('_raw_data', showWarnings = F)

#########################################################
#                                                       #
# Required functions                                    #
#                                                       #
#########################################################

getMetadata <- function(metadata.path = '_metadata.ATAC-seq.xlsx') {
    if (file.exists(metadata.path)) {
        metadata <- read_excel('_metadata.ATAC-seq.xlsx', sheet = 2) # Read data from the excel metadata file
    } else { stop('\nMetadata file not found. Aborting now...\n') }
}

checkSeqIDs <- function(metadata) {
    seqIDs <- na.omit(unlist(strsplit(c(metadata$rep1, metadata$rep2), ';')))
    if (length(seqIDs) != length(unique(seqIDs))) { stop('\nSome seqIDs are in duplicate. Check metadata file. Aborting now...\n') } else { message(paste0(length(seqIDs), ' seqIDs found: ', paste(seqIDs, collapse = ', '))) }
    for (seqID in seqIDs) {
        #message(seqID)
        if (!length(list.files(list.files('../__ATAC_mapping/_raw_data', full.names = T, pattern = seqID), full.names = T)) %in% c(2,4)) {
            stop('\nSome seqIDs raw files are not found. Check metadata file. Aborting now...\n')
        }
    }
    #message('\nAll associated seqIDs files are found. Proceeding to combining data...\n')
    return(T)
}

combineReplicates <- function(metadata) {
    replicates <- unlist(lapply(metadata$tissue.reporter, function(x) { c(paste(x, 'YA_ATAC_rep1', sep = '_'), paste(x, 'YA_ATAC_rep2', sep = '_')) }))
    message(paste0(length(replicates), ' replicates found: ', paste(replicates, collapse = ', ')))

    for (replicate in replicates) {
        message("Combining files for: ", replicate)
        raw.files <- c()
        tissue <- unlist(strsplit(replicate, '_'))[1]
        stage <- unlist(strsplit(replicate, '_'))[2]
        assay <- unlist(strsplit(replicate, '_'))[3]
        rep <- unlist(strsplit(replicate, '_'))[4]
        seqIDs <- unlist(strsplit(as.character(metadata[metadata$tissue.reporter == tissue, rep]), ';'))
        raw.files <- paste(unlist( lapply(seqIDs, function(x) {list.files(list.files('../__ATAC_mapping/_raw_data', full.names = T, pattern = x), full.names = T)} ) ), collapse = ' ')

        system(sprintf('cat %s > %s', raw.files, paste0('_raw_data/', replicate, '.combined.fq.gz')), wait = T)

    }
}

#########################################################
#                                                       #
# EXECUTE COMMANDS TO GROUP SEQ. RESULTS BY REPLICATES  #
#                                                       #
#########################################################

metadata <- getMetadata(metadata.path = '_metadata.ATAC-seq.xlsx')

if (checkSeqIDs(metadata) == T) { message("All files found. proceeding to combining...\n") ; combineReplicates(metadata) }

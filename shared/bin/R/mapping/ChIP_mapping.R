#!/usr/bin/env Rscript

# Raw single-end reads must be located in individual subfolders in the _raw_data subfolder.
# You must be in a folder directly created in ~/, e.g. ~/_JSbiXXX-HSQXXX-xxxxxx/
# which.genome files located in ~/_data.

#########################################################
#                                                       #
# Read arguments (smaples and steps to perform)         #
#                                                       #
#########################################################

library(optparse)
option_list = list(
  make_option("--samples", type="character", default='all',
              help="Which samples (JSsXXX) to map (separated by commas and between single quotes)", metavar="character"),
  make_option("--steps", type="character", default='combine,aln,tracks,beads,norm',
              help="Which steps to perform (separated by commas and between single quotes)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

SAMPLES <- if (opt$samples == 'all') { gsub('-.*', '', list.files(path = '_raw_data/', pattern = '-')) } else { unlist(strsplit(opt$samples, ',')) }
which.steps <- if (opt$steps == 'all') { unlist(strsplit('combine,aln,tracks,beads,norm', ',')) } else { unlist(strsplit(opt$steps, ',')) }

#########################################################
#                                                       #
# Master function                                       #
#                                                       #
#########################################################

chip.se_raw2bw <- function(sample, list.samples = list.samples, steps = c('combine', 'aln', 'tracks', 'beads', 'norm', 'rename', 'stats'), ncore = 8, which.genome = '../_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa') {

    suppressWarnings(suppressMessages(library(rtracklayer)))
    suppressWarnings(suppressMessages(library(GenomicAlignments)))
    suppressWarnings(suppressMessages(library(Rsamtools)))
    suppressWarnings(suppressMessages(library(BSgenome)))
    suppressWarnings(suppressMessages(library(rbeads)))

    write(paste0('\n\n\n--------- Script ran @ ', date(), ' ---------'), file = paste0(sample, '.log'), append = TRUE)
    write(paste0('--------- Now starting mapping of sample ', sample, ' ---------'), file = paste0(sample, '.log'), append = TRUE)
    write(paste0('The sample that is going to be mapped is: \n', sample, '\nThe steps that are going to be exectuted are: \n', paste(which.steps, collapse = ', '), '\n'), file = paste0(sample, '.log'), append = TRUE)

    if ('combine' %in% steps) {
        write('\nCombining reads...', file = paste0(sample, '.log'), append = TRUE)
            combine.reads(sample, list.samples)
    } else { write('\nSkipping reads combining...', file = paste0(sample, '.log'), append = TRUE) }

    if ('aln' %in% steps) {
        write('\nAligning reads...', file = paste0(sample, '.log'), append = TRUE)
            se.bwa_aln(sample, list.samples, ncore = ncore, genome = which.genome)
    } else { write('\nSkipping reads alignment...', file = paste0(sample, '.log'), append = TRUE) }

    if ('tracks' %in% steps) {
        write('\nGenerating alignments tracks...', file = paste0(sample, '.log'), append = TRUE)
            addAlnTracks(sample, list.samples, genome = 'ce11')
    } else { write('\nSkipping alignments tracks...', file = paste0(sample, '.log'), append = TRUE) }

    if ('beads' %in% steps) {
        write('\nNormalizing alignments tracks by BEADS...', file = paste0(sample, '.log'), append = TRUE)
            addBeadsTracks(sample, list.samples, genome = which.genome)
    } else { write('\nSkipping normalizing alignments tracks by BEADS...', file = paste0(sample, '.log'), append = TRUE) }

    if ('norm' %in% steps) {
        write('\nNormalizing alignments tracks by zscore/log2...', file = paste0(sample, '.log'), append = TRUE)
            addScaledTrack(sample, list.samples)
    } else { write('\nSkipping normalizing alignments tracks by zscore/log2...', file = paste0(sample, '.log'), append = TRUE) }

    if (file.exists('_metadata.ChIP-seq.csv')) {
        if ('rename' %in% steps & length(unique(read.csv('_metadata.ChIP-seq.csv', header=T)$SampleName[(SAMPLES %in% read.csv('_metadata.ChIP-seq.csv', header=T)$libraries)])) == length(SAMPLES)) {
            write('\nRenaming files...', file = paste0(sample, '.log'), append = TRUE)
                renameFiles(sample, list.samples)
        }
    } else { write('\nSkipping files renaming...', file = paste0(sample, '.log'), append = TRUE) }

    if ('stats' %in% steps) {
        write('\nGetting bam stats for the generated bam files', file = paste0(sample, '.log'), append = TRUE)
            bamstats <- bamStats(list.samples[[sample]]$aligned)
            write(bamstats, file = paste0(sample, '.log'), append = TRUE)
    } else { write('\nSkipping tracks stats...', file = paste0(sample, '.log'), append = TRUE) }

}

#########################################################
#                                                       #
# Other functions                                       #
#                                                       #
#########################################################

# Combine reads for each sample
combine.reads <- function(sample, list.samples) {

    sample.files <- grep('combined', grep('.gz', list.files(pattern = sample, recursive = T), value = T), invert = T, value = T)
    cmd <- sprintf('cat %s %s > %s', sample.files[1], sample.files[2], list.samples[[sample]]$combined)
    system(cmd, wait = TRUE)

}

# Align reads to ce11
se.bwa_aln <- function(sample, list.samples = list.samples, genome = '../_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa', ncore = 8) {

    cmd1 <- sprintf("bwa samse %s <(bwa aln -t %i %s %s) %s | samtools view -@ %i -bSu - | samtools sort -@ %i -m 44294967296 - -o %s", genome, ncore, genome, list.samples[[sample]]$combined, list.samples[[sample]]$combined, ncore, ncore, list.samples[[sample]]$aligned)
    cmd2 <- sprintf('echo "%s" | bash', cmd1)
    system(cmd2, wait = TRUE)

}

# Add alnNQNU/alnQ10NU bigwig tracks
addAlnTracks <- function(sample, list.samples = list.samples, genome = 'ce11') {

    what <- c("rname", "strand", "pos", "qwidth", "mapq")
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(flag = flag, simpleCigar = FALSE, what = what)

    a <- scanBam(list.samples[[sample]]$aligned, param = param)[[1]]
    grng <- GRanges(seqnames = a$rname, ranges = IRanges(a$pos, width = a$qwidth), strand = a$strand, seqinfo=SeqinfoForUCSCGenome(genome), mapq=a$mapq)
    grng <- trim(resize(grng, 200L))
    export.bw(coverage(grng), list.samples[[sample]]$aligned.alnNQNU.bw)

    grngQ10 <- grng[grng$mapq >= 10]
    export.bw(coverage(grngQ10), list.samples[[sample]]$aligned.alnQ10NU.bw)

}

# Add BEADS tracks
addBeadsTracks <- function(sample, list.samples, genome = '../_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa', input = '../_data/ce11/Input^NA_fb23^F^N2^L3_aligned^NA^ce11_AA617^A9a44476.bam', map = '../_data/ce11/ce11_gem-mappability_36bp.bw') {

    # Generate Q0-non-unique
    beadsobj <- beads(
        list.samples[[sample]]$aligned,
        input,
        NULL,
        genome,
        uniq = FALSE, insert = 200L, mapq_cutoff = 0, export = "BEADS",
        rdata = FALSE, export_er = TRUE, quickMap = TRUE
    )
    file.rename(path(beadsobj), list.samples[[sample]]$aligned.BEADSNQNU.linear.bw)

    # Generate Q10-non-unique
    beadsobj <- beads(
        list.samples[[sample]]$aligned,
        input,
        map,
        genome,
        uniq = FALSE, insert = 200L, mapq_cutoff = 10, export = "BEADS",
        rdata = FALSE, export_er = TRUE, quickMap = TRUE
    )
    file.rename(path(beadsobj), list.samples[[sample]]$aligned.BEADSQ10NU.linear.bw)

    # Generate Q10-unique
    beadsobj <- beads(
        list.samples[[sample]]$aligned,
        input,
        map,
        genome,
        uniq = TRUE, insert = 200L, mapq_cutoff = 10, export = "BEADS",
        rdata = FALSE, export_er = TRUE, quickMap = TRUE
    )
    file.rename(path(beadsobj), list.samples[[sample]]$aligned.BEADSQ10U.linear.bw)

    # Remove junk
    cmd <- sprintf("rm %s %s", '*pdf', '*bed')
    system(cmd, wait = T)

}

# Add scaled tracks
addScaledTrack <- function(sample, list.samples) {

    require(rtracklayer)

    # Generate Q0-non-unique norm. tracks
    cov <- import.bw(list.samples[[sample]]$aligned.BEADSNQNU.linear.bw, as='RleList')
        ucov <- unlist(cov)
        mi <- mean(ucov)
        mu <- sd(ucov)
        zsc <- (cov-mi)/mu
        export.bw(zsc, list.samples[[sample]]$aligned.BEADSNQNU.zscore.bw)

        log2sc <- log2(cov)
        log2sc_gr <- as(log2sc, 'GRanges')
        log2sc_gr <- log2sc_gr[is.finite(log2sc_gr$score)]
        export.bw(log2sc_gr, list.samples[[sample]]$aligned.BEADSNQNU.log2.bw)

    # Generate Q10-non-unique norm. tracks
    cov <- import.bw(list.samples[[sample]]$aligned.BEADSQ10NU.linear.bw, as='RleList')
        ucov <- unlist(cov)
        mi <- mean(ucov)
        mu <- sd(ucov)
        zsc <- (cov-mi)/mu
        export.bw(zsc, list.samples[[sample]]$aligned.BEADSQ10NU.zscore.bw)

        log2sc <- log2(cov)
        log2sc_gr <- as(log2sc, 'GRanges')
        log2sc_gr <- log2sc_gr[is.finite(log2sc_gr$score)]
        export.bw(log2sc_gr, list.samples[[sample]]$aligned.BEADSQ10NU.log2.bw)

    # Generate Q0-unique norm. tracks
    cov <- import.bw(list.samples[[sample]]$aligned.BEADSQ10U.linear.bw, as='RleList')
        ucov <- unlist(cov)
        mi <- mean(ucov)
        mu <- sd(ucov)
        zsc <- (cov-mi)/mu
        export.bw(zsc, list.samples[[sample]]$aligned.BEADSQ10U.zscore.bw)

        log2sc <- log2(cov)
        log2sc_gr <- as(log2sc, 'GRanges')
        log2sc_gr <- log2sc_gr[is.finite(log2sc_gr$score)]
        export.bw(log2sc_gr, list.samples[[sample]]$aligned.BEADSQ10U.log2.bw)

}

# Rename files of the sample
renameFiles <- function(sample, list.samples) {
    metadata <- read.csv('_metadata.ChIP-seq.csv', header=T, stringsAsFactors = F)
    new.id <- unlist(lapply(SAMPLES, function(x) (metadata$SampleName[metadata$libraries == x])))
    names(new.id) <- SAMPLES
    old.names <- unlist(list.samples)
    new.names <- unlist(lapply(SAMPLES, function(x) (gsub(x, new.id[names(new.id) == x], list.samples[[x]]))))
    if ( length(unique(new.id)) == length(SAMPLES) & length(unique(new.names)) == length(unique(old.names)) ) {
        if (sample == 'all') {
            for (k in 1:length(old.names)) { file.rename(old.names[k], new.names[k]) }
        } else if (length(sample) == 1) {
            for (k in grep(sample, old.names)) { file.rename(old.names[k], new.names[k]) }
        }
    } else { 'No renaming was performed'}
}
reverse.renameFiles <- function(sample, list.samples) {
    metadata <- read.csv('_metadata.ChIP-seq.csv', header=T, stringsAsFactors = F)
    new.id <- unlist(lapply(SAMPLES, function(x) (metadata$SampleName[metadata$libraries == x])))
    names(new.id) <- SAMPLES
    old.names <- unlist(list.samples)
    new.names <- unlist(lapply(SAMPLES, function(x) (gsub(x, new.id[names(new.id) == x], list.samples[[x]]))))
    if ( length(unique(new.id)) == length(SAMPLES) & length(unique(new.names)) == length(unique(old.names)) ) {
        if (sample == 'all') {
            for (k in 1:length(old.names)) { file.rename(old.names[k], new.names[k]) }
        } else if (length(sample) == 1) {
            for (k in grep(sample, old.names)) { file.rename(new.names[k], old.names[k]) }
        }
    } else { 'No renaming was performed'}
}

# Get some bam stats from bam file f (from JADBtools)
bamStats <- function(files) {

    load('../_data/rrnamodel.rda')

    out <- list()

    for (f in files) {
        what <- c("rname", "strand", "pos", "mapq", "qwidth")
        a <- scanBam(f, param = ScanBamParam(what = what))[[1]]
        r <- length(a[[1]])
        align <- !is.na(a$pos)
        lg <- a$mapq[align] >= 10L
        grng <- GRanges(
            seqnames = a$rname[align],
            ranges = IRanges(a$pos[align], width = a$qwidth[align]),
            strand = a$strand[align]
        )
        aa <- length(grng)
        q <- sum(lg)
        u <- length(unique(grng[lg]))
        rr <- sum(overlapsAny(grng, rrnamodel, ignore.strand=TRUE))

        out[f] <- sprintf(
            'all=%.2fM, aligned=%.2fM[%.0f%%], mapq10=%.2fM[%.0f%%], unique10=%.2fM[%.0f%%], rRNA=%.2fM[%.0f%%]',
            r/10^6, aa/10^6, (aa/r)*100, q/10^6, (q/r)*100, u/10^6, (u/r)*100, rr/10^6, (rr/r)*100
        )
    }

    out <- apply(matrix(c(names(out), unname(unlist(lapply(out, '[[', 1)))), ncol = 2), 1, function(x) {gsub('_bam-files/', '', paste(x, collapse = ":   "))} )

    return(out)

}


#########################################################
#                                                       #
# LAUNCH MASTER FUNCTION IN PARALLEL                    #
#                                                       #
#########################################################

cat(paste0('The samples that are going to be mapped are: \n', paste(SAMPLES, collapse = ', '), '\nThe steps that are going to be exectuted are: \n', paste(which.steps, collapse = ', '), '\n\nIs that correct? (y/n)'))
answer <- readLines("stdin",n=1)
if (answer == 'y') { cat('\nGood! Processing now!\n') } else { stop('\nAborting mapping...\n') }

dir.create('_bam-files/')
dir.create('_bw-files/')
dir.create('_bw-files/alnQ0nu')
dir.create('_bw-files/alnQ10nu')
dir.create('_bw-files/BEADSnuQ0')
dir.create('_bw-files/BEADSnuQ10')
dir.create('_bw-files/BEADSuQ10')

list.samples <- list()
for (sample in SAMPLES) {
    list.samples[[sample]] <- list(
        combined = paste0("_raw_data/", sample, '.combined.fq.gz'),
        aligned = paste0('_bam-files/', sample, '.map_se.bam'),
        aligned.alnNQNU.bw = paste0('_bw-files/alnQ0nu/', sample, '.map_se^alnQ0nu.bw'),
        aligned.alnQ10NU.bw = paste0('_bw-files/alnQ10nu/', sample, '.map_se^alnQ10nu.bw'),
        aligned.BEADSNQNU.linear.bw = paste0('_bw-files/BEADSnuQ0/', sample, '.map_se^BEADSnuQ0.linear.bw'),
        aligned.BEADSNQNU.log2.bw = paste0('_bw-files/BEADSnuQ0/', sample, '.map_se^BEADSnuQ0.log2.bw'),
        aligned.BEADSNQNU.zscore.bw = paste0('_bw-files/BEADSnuQ0/', sample, '.map_se^BEADSnuQ0.zscore.bw'),
        aligned.BEADSQ10NU.linear.bw = paste0('_bw-files/BEADSnuQ10/', sample, '.map_se^BEADSnuQ10.linear.bw'),
        aligned.BEADSQ10NU.log2.bw = paste0('_bw-files/BEADSnuQ10/', sample, '.map_se^BEADSnuQ10.log2.bw'),
        aligned.BEADSQ10NU.zscore.bw = paste0('_bw-files/BEADSnuQ10/', sample, '.map_se^BEADSnuQ10.zscore.bw'),
        aligned.BEADSQ10U.linear.bw = paste0('_bw-files/BEADSuQ10/', sample, '.map_se^BEADSuQ10.linear.bw'),
        aligned.BEADSQ10U.log2.bw = paste0('_bw-files/BEADSuQ10/', sample, '.map_se^BEADSuQ10.log2.bw'),
        aligned.BEADSQ10U.zscore.bw = paste0('_bw-files/BEADSuQ10/', sample, '.map_se^BEADSuQ10.zscore.bw')
    )
}

library(parallel)
mclapply(SAMPLES, chip.se_raw2bw, list.samples = list.samples, steps = which.steps, ncore = 4, which.genome = '../_data/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa', mc.cores=6)

# Generate report
#library(knitr)
#library(rmarkdown)
#params <- list(SAMPLES = SAMPLES, steps = which.steps)
#render('../_scripts/bin/ChIP_mapping_report.Rmd', output_file = paste0(getwd(), '/_log-files/report.pdf'), output_dir = paste0(getwd(), '/_log-files'), output_format = 'pdf_document', params = params)

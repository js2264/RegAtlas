#!/usr/bin/env Rscript

# Raw single-end reads must be located in individual subfolders in the _raw_data subfolder.
# You must be in a folder directly created in ~/, e.g. ~/_JSbiXXX-HSQXXX-xxxxxx/
# which.genome files located in ~/_data.

#########################################################
#                                                       #
# Read arguments (samples and steps to perform)         #
#                                                       #
#########################################################

library(optparse)
option_list = list(
  make_option("--samples", type="character", default='all',
              help="Which samples (JSsXXX) to map (separated by commas and between single quotes)", metavar="character"),
  make_option("--steps", type="character", default='combine,trim,aln,chrM,blacklist,q10,tracks,rename,cleaning,stats',
              help="Which steps to perform (separated by commas and between single quotes)", metavar="character"),
  make_option("--cores", type="integer", default=4,
              help="How many cores to use?", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

SAMPLES <- if (opt$samples == 'all') { gsub('-.*', '', list.files(path = '_raw_data/', pattern = '-')) } else { unlist(strsplit(opt$samples, ',')) }
which.steps <- if (opt$steps == 'all') { unlist(strsplit('combine,trim,aln,chrM,blacklist,q10,tracks,rename,cleaning,stats', ',')) } else { unlist(strsplit(opt$steps, ',')) }
cores <- opt$cores

#########################################################
#                                                       #
# Master function                                       #
#                                                       #
#########################################################

atac.se_raw2bw <- function(sample, list.samples = list.samples, steps, ncore = 8, which.genome = '../_data/ce10.fa') {

    write(paste0('\n\n\n--------- Script ran @ ', date(), ' ---------'), file = list.samples[[sample]]$log, append = TRUE)
    write(paste0('--------- Now starting mapping of sample ', sample, ' ---------'), file = list.samples[[sample]]$log, append = TRUE)
    write(paste0('The sample that is going to be mapped is: \n', sample, '\nThe steps that are going to be exectuted are: \n', paste(steps, collapse = ', '), '\n'), file = list.samples[[sample]]$log, append = TRUE)

    if ('combine' %in% steps) {
        write('\nCombining reads...', file = list.samples[[sample]]$log, append = TRUE)
            combine.reads(sample, list.samples)
    } else { write('\nSkipping reads combining...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('trim' %in% steps) {
        write('\nTrimming reads...', file = list.samples[[sample]]$log, append = TRUE)
            trim.reads(sample, list.samples, l = 50)
    } else { write('\nSkipping reads trimming...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('aln' %in% steps) {
        write('\nAligning reads...', file = list.samples[[sample]]$log, append = TRUE)
            se.bwa_aln(sample, list.samples, genome = which.genome, ncore = ncore)
    } else { write('\nSkipping reads alignment...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('chrM' %in% steps) {
        write('\nFiltering of chrM reads...', file = list.samples[[sample]]$log, append = TRUE)
            filter.chrM(sample, list.samples, ncore = ncore)
    } else { write('\nSkipping filtering of chrM reads...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('blacklist' %in% steps) {
        write('\nFiltering of blacklist regions...', file = list.samples[[sample]]$log, append = TRUE)
            filter.blacklist(sample, list.samples, ncore = ncore, blacklist = '../_data/ce10_blacklist.bed')
    } else { write('\nSkipping filtering of blacklist regions...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('q10' %in% steps) {
        write('\nFiltering of <q10 reads...', file = list.samples[[sample]]$log, append = TRUE)
            filter.q10(sample, list.samples, ncore = ncore)
    } else { write('\nSkipping filtering of <q10 reads...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('tracks' %in% steps) {
        write('\nGenerating alignments tracks...', file = list.samples[[sample]]$log, append = TRUE)
            getBwTrack(sample, list.samples, genome.chrom.sizes = '../_data/ce10.chrom.sizes')
    } else { write('\nSkipping alignments tracks...', file = list.samples[[sample]]$log, append = TRUE) }

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

# Trim reads
trim.reads <- function(sample, list.samples, l = 20) {

    cmd <- sprintf('zcat -d -c %s | fastx_trimmer -z -l %i -Q 33 -o %s', list.samples[[sample]]$combined, l, list.samples[[sample]]$trimmed)
    system(cmd, wait = TRUE)

}

# Align reads to ce10
se.bwa_aln <- function(sample, list.samples, genome = '../_data/ce10.fa', ncore = 8) {

    cmd1 <- sprintf("bwa samse %s <(bwa aln -t %i %s %s) %s | samtools view -@ %i -bSu - | samtools sort -@ %i -m 44294967296 - -o %s", genome, ncore, genome, list.samples[[sample]]$trimmed, list.samples[[sample]]$trimmed, ncore, ncore, list.samples[[sample]]$aligned)
    cmd2 <- sprintf('echo "%s" | bash', cmd1)
    system(cmd2, wait = TRUE)

}

# Filter ChrM
filter.chrM <- function(sample, list.samples, ncore = 8) {

    cmd <- sprintf('samtools index -@ %i %s ; samtools view -@ %i -b %s chrI chrII chrIII chrIV chrV chrX > %s ; rm %s', ncore, list.samples[[sample]]$aligned, ncore, list.samples[[sample]]$aligned, list.samples[[sample]]$rm_chrM, paste0(list.samples[[sample]]$aligned, '.bai'))
    system(cmd, wait = TRUE)

}

# Filter blacklisted regions
filter.blacklist <- function(sample, list.samples, ncore = 8, blacklist = '../_data/ce10_blacklist.bed') {

    cmd <- sprintf('samtools view -@ %i -b -L %s -U %s %s > /dev/null', ncore, blacklist, list.samples[[sample]]$rm_blacklist, list.samples[[sample]]$rm_chrM)
    system(cmd, wait = TRUE)

}

# Filter Q10 reads
filter.q10 <- function(sample, list.samples, ncore = 8) {

    cmd <- sprintf('samtools view -@ %i -b -q 10 %s > %s', ncore, list.samples[[sample]]$rm_blacklist, list.samples[[sample]]$q10)
    system(cmd, wait = TRUE)

}

# Generate BigWig track
getBwTrack <- function(sample, list.samples, genome.chrom.sizes = '../_data/ce10.chrom.sizes') {


    cmd1 <- sprintf('macs2 callpeak -t %s --format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 150 --shift -75 --keep-dup all --name %s', list.samples[[sample]]$q10, paste0(sample, '.tmp.bg'))
    cmd2 <- sprintf('sort -k1,1 -k2,2n %s > %s', paste0(sample, '.tmp.bg_treat_pileup.bdg'), paste0(sample, '.tmp.2.bdg'))
    cmd3 <- sprintf('bedGraphToBigWig %s %s %s', paste0(sample, '.tmp.2.bdg'), genome.chrom.sizes, list.samples[[sample]]$track)
    cmd4 <- sprintf('mv %s %s ; mv %s %s', paste0(sample, '.tmp.bg_peaks.narrowPeak'), list.samples[[sample]]$peaks, paste0(sample, '.tmp.bg_peaks.xls'), list.samples[[sample]]$peaks_xls)

    system(cmd1, wait = T)
    system(cmd2, wait = T)
    system(cmd3, wait = T)
    system(cmd4, wait = T)

}

# Clean up temporary files (must be done outside of master function)
cleanTemp <- function(sample) {


    cmd <- sprintf("rm %s %s", paste0('*', sample, '*tmp*'), paste0('_raw_data/*', sample, '*combined*'))
    system(cmd, wait = T)
    cmd <- sprintf("mkdir _log-files ; mv *log _log-files")
    system(cmd, wait = T)

}

# Rename all files (must be done outside of the master function)
renameFiles <- function(sample, list.samples, metadata.path = '../__ATAC_mapping/_metadata.ATAC-seq.csv', direction = 'JSIDtoNAME') {
    metadata <- read.csv(metadata.path, header=T, stringsAsFactors = F)
    old.names <- unlist(list.samples[[sample]])
    if (length(metadata$SampleName[metadata$libraries == sample]) > 0) { new.names <- gsub(sample, metadata$SampleName[metadata$libraries == sample], old.names) } else { new.names = '' }
    if (any(file.exists(old.names)) & all(nchar(new.names) > 0) & direction == 'JSIDtoNAME') {
        for (k in 1:length(old.names)) {
            file.rename(old.names[k], new.names[k])
            message(paste0(old.names[k], '  --> ', new.names[k]))
        }
    }
    if (any(file.exists(new.names)) & all(nchar(old.names) > 0) & direction == 'NAMEtoJSID') {
        for (k in 1:length(old.names)) {
            file.rename(new.names[k], old.names[k])
            message(paste0(new.names[k], '  --> ', old.names[k]))
        }
    }
}

# Get some bam stats for all bam files (from JADBtools ; must be done outside of master function)
bamStats <- function(sample, list.samples, metadata.path = '../__ATAC_mapping/_metadata.ATAC-seq.csv', plot.path = 'bamstats.pdf') {

    metadata <- read.csv(metadata.path, header=T, stringsAsFactors = F)
    files <- if (file.exists(list.samples[[sample]]$aligned)) {
        list.samples[[sample]]$aligned
    } else if (file.exists(metadata.path) & file.exists(gsub(sample, metadata$SampleName[metadata$libraries == sample], list.samples[[sample]]$aligned))) {
        gsub(sample, metadata$SampleName[metadata$libraries == sample], list.samples[[sample]]$aligned)
    }

    peak.files <- if (file.exists(list.samples[[sample]]$peaks)) {
        list.samples[[sample]]$peaks
    } else if (file.exists(metadata.path) & file.exists(gsub(sample, metadata$SampleName[metadata$libraries == sample], list.samples[[sample]]$peaks))) {
        gsub(sample, metadata$SampleName[metadata$libraries == sample], list.samples[[sample]]$peaks)
    }

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

        if (file.exists(peak.files)) { peaknb <- length(readLines(peak.files)) } else { peaknb <- 0 }

        out[f] <- sprintf(
            'all=%.2fM, aligned=%.2fM[%.0f%%], mapq10=%.2fM[%.0f%%], unique10=%.2fM[%.0f%%], rRNA=%.2fM[%.0f%%], peaks=%i',
            r/10^6, aa/10^6, (aa/r)*100, q/10^6, (q/r)*100, u/10^6, (u/r)*100, rr/10^6, (rr/r)*100, peaknb
        )
    }

    df <- as.data.frame(matrix(NA, ncol = 2), byrow=T)
    df[,1] <- names(out)
    df[,2] <- out

    return(df)

}

# Plot some QC from deeptools
plotDeepToolsQC <- function() {

    cmd1 <- "multiBigwigSummary bins -b $(ls _bw-files/*bw | tr '\n' ' ') -out deeptools.1000.summary -bs 1000"
    cmd2 <- "plotCorrelation -in deeptools.1000.summary -o deeptools.1000.corr.pearson.heatmap.pdf -c pearson -p heatmap --plotNumbers --removeOutliers --labels $(ls _bw-files/* | sed 's/.map.*//' | sed 's,_bw-files/,,' | tr '\n' ' ')"
    cmd3 <- "plotPCA -in deeptools.1000.summary -o deeptools.1000.pca.pdf --labels $(ls _bw-files/* | sed 's/.map.*//' | sed 's,_bw-files/,,' | tr '\n' ' ')"

    system(cmd1, wait = TRUE)
    system(cmd2, wait = TRUE)
    system(cmd3, wait = TRUE)

}

#########################################################
#                                                       #
# LAUNCH MASTER FUNCTION IN PARALLEL                    #
#                                                       #
#########################################################

cat(paste0('The samples that are going to be mapped are: (', length(SAMPLES), ') \n', paste(SAMPLES, collapse = ', '), '\nThe steps that are going to be exectuted are: \n', paste(which.steps, collapse = ', '), '\n\nIs that correct? (y/n)'))
answer <- readLines("stdin",n=1)
if (answer == 'y') { cat('\nGood! Processing now!\n') } else { stop('\nAborting mapping...\n') }

suppressWarnings(suppressMessages(library(rtracklayer)))
suppressWarnings(suppressMessages(library(GenomicAlignments)))
suppressWarnings(suppressMessages(library(Rsamtools)))
suppressWarnings(suppressMessages(library(BSgenome)))
if (!dir.exists('_bam-files/')) {dir.create('_bam-files/')}
if (!dir.exists('_log-files/')) {dir.create('_log-files/')}
if (!dir.exists('_bw-files/')) {dir.create('_bw-files/')}
if (!dir.exists('_peaks-files/')) {dir.create('_peaks-files/')}

list.samples <- list()
for (sample in SAMPLES) {
    list.samples[[sample]] <- list(
        log = paste0('_log-files/', sample, '.log'),
        combined = paste0('_raw_data/', sample, '.combined.fq.gz'),
        trimmed = paste0('_raw_data/', sample, '.combined.trimmed.fq.gz'),
        aligned = paste0('_bam-files/', sample, '.map_se.bam'),
        rm_chrM = paste0('_bam-files/', sample, '.map_se^rm_chrM.bam'),
        rm_blacklist = paste0('_bam-files/', sample, '.map_se^rm_chrM^rm_blacklist.bam'),
        q10 = paste0('_bam-files/', sample, '.map_se^rm_chrM^rm_blacklist^q10.bam'),
        track = paste0('_bw-files/', sample, '.map_se^rm_chrM^rm_blacklist^q10.bw'),
        peaks = paste0('_peaks-files/', sample, '.peaks.bed'),
        peaks_xls = paste0('_peaks-files/', sample, '.peaks.xls')
    )
}

library(parallel)
mclapply(SAMPLES, atac.se_raw2bw, list.samples = list.samples, steps = which.steps, ncore = 4, which.genome = '../_data/ce10.fa', mc.cores = cores)

# Then finish up the script (rename, stats, cleanup)
if ('cleaning' %in% which.steps) {
    write('\nCleaning-up temp. files...', file = list.samples[[sample]]$log, append = TRUE)
    for (sample in SAMPLES) { cleanTemp(sample) }
} else { write('\nSkipping cleaning up...', file = list.samples[[sample]]$log, append = TRUE) }

if (file.exists('../__ATAC_mapping/_metadata.ATAC-seq.csv') & 'rename' %in% which.steps) {
    write('\nRenaming files...', file = list.samples[[sample]]$log, append = TRUE)
    for (sample in SAMPLES) { renameFiles(sample, list.samples, metadata.path = '../__ATAC_mapping/_metadata.ATAC-seq.csv', direction = 'JSIDtoNAME') }
} else if (file.exists('../__ATAC_mapping/_metadata.ATAC-seq.csv') & 'revrename' %in% which.steps) {
    for (sample in SAMPLES) { renameFiles(sample, list.samples, metadata.path = '../__ATAC_mapping/_metadata.ATAC-seq.csv', direction = 'NAMEtoJSID') }
} else { write('\nSkipping files renaming...', file = list.samples[[sample]]$log, append = TRUE) }

if ('stats' %in% which.steps) {
    write('\nGetting bam stats for the generated bam files...', file = list.samples[[sample]]$log, append = TRUE)
    write(paste0('\n\n--------- ', date(), ' ---------\n'), file = 'bamstats.log', append = TRUE)
    for (sample in SAMPLES) {
        bamstats <- bamStats(sample, list.samples, metadata.path = '../__ATAC_mapping/_metadata.ATAC-seq.csv')
        write.table(bamstats, file = 'bamstats.log', quote = F, col.names = F, row.names = F, sep = ': ', append = TRUE)
    }
} else { write('\nSkipping tracks stats...', file = list.samples[[sample]]$log, append = TRUE) }

if ('deeptools' %in% which.steps) {
    write('\nQCing the generated files...', file = list.samples[[sample]]$log, append = TRUE)
    plotDeepToolsQC()
} else { write('\nSkipping QC...', file = list.samples[[sample]]$log, append = TRUE) }



# Generate report
#library(knitr)
#library(rmarkdown)
#params <- list(SAMPLES = SAMPLES, steps = which.steps)
#render('../_scripts/bin/ATAC_mapping_report.Rmd', output_file = paste0(getwd(), '/_log-files/report.pdf'), output_dir = paste0(getwd(), '/_log-files'), output_format = 'pdf_document', params = params)

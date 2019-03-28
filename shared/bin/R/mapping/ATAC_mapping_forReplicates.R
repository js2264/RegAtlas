#!/usr/bin/env Rscript

# This script is fully automated (tested on 13/07/2018)

# which.genome files located in ~/_data.
# Raw pooled replicates single-end sequencing results should be in _raw_data, labelled myo3_repX.combined.fq.gz

#########################################################
#                                                       #
# Read arguments (samples and steps to perform)         #
#                                                       #
#########################################################

library(optparse)
option_list = list(
  make_option("--samples", type="character", default = 'all',
              help="Which samples (xxx_repX) to map (separated by commas and between single quotes)", metavar="character"),
  make_option("--steps", type="character", default = 'all',
              help="Which steps to perform (separated by commas and between single quotes)", metavar="character"),
  make_option("--cores", type="integer", default = 5,
              help="How many cores to use?", metavar="character"),
  make_option("--doCombineReps", type="logical", default = F,
              help="Should the raw files be combined first?", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

doCombineReps <- opt$doCombineReps
if (doCombineReps == T) { system("Rscript ../_scripts/bin/master-scripts/ATAC_pooling.R", wait = T) }

SAMPLES <- if (opt$samples == 'all') { gsub('.combined.fq.gz', '', list.files(path = '_raw_data/', pattern = 'rep..combined.fq.gz')) } else { unlist(strsplit(opt$samples, ',')) }
which.steps <- if (opt$steps == 'all') { unlist(strsplit('trim,aln,chrM,blacklist,q10,tracks,cleaning,stats,combinetracks', ',')) } else { unlist(strsplit(opt$steps, ',')) }
cores <- opt$cores

#########################################################
#                                                       #
# Master function                                       #
#                                                       #
#########################################################

atac.se_raw2bw <- function(sample, list.samples = list.samples, steps, ncore = 8, which.genome = '../_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa') {

    write(paste0('\n\n\n--------- Script ran @ ', date(), ' ---------'), file = list.samples[[sample]]$log, append = TRUE)
    write(paste0('--------- Now starting mapping of sample ', sample, ' ---------'), file = list.samples[[sample]]$log, append = TRUE)
    write(paste0('The sample that is going to be mapped is: \n', sample, '\nThe steps that are going to be exectuted are: \n', paste(steps, collapse = ', '), '\n'), file = list.samples[[sample]]$log, append = TRUE)

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
            filter.blacklist(sample, list.samples, ncore = ncore, blacklist = '../_data/ce11/ce10_blacklist.bed')
    } else { write('\nSkipping filtering of blacklist regions...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('q10' %in% steps) {
        write('\nFiltering of <q10 reads...', file = list.samples[[sample]]$log, append = TRUE)
            filter.q10(sample, list.samples, ncore = ncore)
    } else { write('\nSkipping filtering of <q10 reads...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('tracks' %in% steps) {
        write('\nGenerating alignments tracks...', file = list.samples[[sample]]$log, append = TRUE)
            getBwTrack(sample, list.samples, genome.chrom.sizes = '../_data/ce11/ce11.chrom.sizes')
    } else { write('\nSkipping alignments tracks...', file = list.samples[[sample]]$log, append = TRUE) }

}

#########################################################
#                                                       #
# Other functions                                       #
#                                                       #
#########################################################

# Trim reads
trim.reads <- function(sample, list.samples, l = 20) {

    cmd <- sprintf('zcat -d -c %s | fastx_trimmer -z -l %i -Q 33 -o %s', list.samples[[sample]]$combined, l, list.samples[[sample]]$trimmed)
    system(cmd, wait = TRUE)

}

# Align reads to ce11
se.bwa_aln <- function(sample, list.samples, genome = '../_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa', ncore = 8) {

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
filter.blacklist <- function(sample, list.samples, ncore = 8, blacklist = '../_data/ce11/ce10_blacklist.bed') {

    cmd <- sprintf('samtools view -@ %i -b -L %s -U %s %s > /dev/null', ncore, blacklist, list.samples[[sample]]$rm_blacklist, list.samples[[sample]]$rm_chrM)
    system(cmd, wait = TRUE)

}

# Filter Q10 reads
filter.q10 <- function(sample, list.samples, ncore = 8) {

    cmd <- sprintf('samtools view -@ %i -b -q 10 %s > %s', ncore, list.samples[[sample]]$rm_blacklist, list.samples[[sample]]$q10)
    system(cmd, wait = TRUE)

}

# Generate BigWig track
getBwTrack <- function(sample, list.samples, genome.chrom.sizes = '../_data/ce11/ce11.chrom.sizes') {


    cmd1 <- sprintf('macs2 callpeak -t %s --format BAM --bdg --SPMR --gsize ce --nolambda --nomodel --extsize 150 --shift -75 --keep-dup all --name %s', list.samples[[sample]]$q10, paste0(sample, '.tmp.bg'))
    cmd2 <- sprintf('sort -k1,1 -k2,2n %s > %s', paste0(sample, '.tmp.bg_treat_pileup.bdg'), paste0(sample, '.tmp.2.bdg'))
    cmd3 <- sprintf('bedGraphToBigWig %s %s %s', paste0(sample, '.tmp.2.bdg'), genome.chrom.sizes, list.samples[[sample]]$track)
    cmd4 <- sprintf('mv %s %s ; mv %s %s', paste0(sample, '.tmp.bg_peaks.narrowPeak'), list.samples[[sample]]$peaks, paste0(sample, '.tmp.bg_peaks.xls'), list.samples[[sample]]$peaks_xls)

    system(cmd1, wait = T)
    system(cmd2, wait = T)
    system(cmd3, wait = T)
    system(cmd4, wait = T)

}

# Get some bam stats for all bam files (from JADBtools ; must be done outside of master function)
bamStats <- function(list.samples) {

    files <- unlist(lapply(list.samples, function(x) {'[['(x,'aligned')}))
    peaks <- unlist(lapply(list.samples, function(x) {'[['(x,'peaks')}))

    load('../_data/rrnamodel.rda')

    out <- list()

    for (f in files) {
        
        message('\t>> Bamstats for: ', f)
        what <- c("rname", "strand", "pos", "mapq", "qwidth")
        a <- scanBam(f, param = ScanBamParam(what = what))[[1]]
        align <- !is.na(a$pos)
        grng <- GRanges(
            seqnames = a$rname[align],
            ranges = IRanges(a$pos[align], width = a$qwidth[align]), 
            qual = a$mapq[align]
        )
        map_pe <- a[[1]] # All individual reads
        rm_unmapped <- grng # All individual reads mapping the genome
        rm_chrM <- rm_unmapped[seqnames(rm_unmapped) != 'chrM'] # All individual reads mapping the genome but not on chrM 
        rm_rRNA <- rm_chrM[!(rm_chrM %over% rtracklayer::import.bed('../_data/ce11/WS253_ce10.rRNA_broad.bed'))] # All individual reads mapping the genome but not on chrM nor on rRNA clusters nor on blacklisted regions
        rm_blacklist <- rm_rRNA[!(rm_rRNA %over% rtracklayer::import.bed('../_data/ce11/ce11_blacklist.bed'))] # All individual reads mapping the genome but not on chrM nor on rRNA clusters nor on blacklisted regions
        q10 <- rm_blacklist[rm_blacklist$qual >= 10L] # All individual reads mapping the genome but not on chrM nor on rRNA clusters nor on blacklisted regions, with a auqlity score > 10
        peaknb <- length(readLines(peaks[names(files[files %in% f])]))

        out[f] <- sprintf(
            'all=%.2fM, -unmapped=%.2fM [%.0f%% of total], -chrM=%.2fM [%.0f%% of total], -rRNA=%.2fM [%.0f%% of total], -blacklist=%.2fM [%.0f%% of total], q10=%.2fM [%.0f%% of total], peaks=%i', 
            length(map_pe)/10^6, length(rm_unmapped)/10^6, (length(rm_unmapped)/length(map_pe))*100, length(rm_chrM)/10^6, (length(rm_chrM)/length(map_pe))*100, length(rm_rRNA)/10^6, (length(rm_rRNA)/length(map_pe))*100, length(rm_blacklist)/10^6, (length(rm_blacklist)/length(map_pe))*100, length(q10)/10^6, (length(q10)/length(map_pe))*100, peaknb
        )

        df <- as.data.frame(matrix(NA, ncol = 2, nrow = length(out)), byrow=T)
        df[,1] <- names(out)
        df[,2] <- unlist(out)
        
    }

    df <- as.data.frame(matrix(NA, ncol = 2, nrow = length(out)), byrow=T)
    df[,1] <- names(out)
    df[,2] <- unlist(out)

    return(df)

}

# Clean up temporary files (must be done outside of master function)
cleanTemp <- function(sample) {


    cmd <- sprintf("rm %s %s", paste0('*', sample, '*tmp*'), paste0('_raw_data/*', sample, '*trimmed*'))
    system(cmd, wait = T)

}

# Plot some QC from deeptools
plotDeepToolsQC <- function() {

    cmd1 <- "multiBigwigSummary bins -b $(ls _bw-files/*_rep*bw | tr '\n' ' ') -out deeptools.1000.summary -bs 1000"
    cmd2 <- "plotCorrelation -in deeptools.1000.summary -o deeptools.1000.corr.pearson.heatmap.pdf -c pearson -p heatmap --plotNumbers --removeOutliers --labels $(ls _bw-files/*_rep*bw | sed 's/.map.*//' | sed 's,_bw-files/,,' | tr '\n' ' ')"
    cmd3 <- "plotCorrelation -in deeptools.1000.summary -o deeptools.1000.corr.spearman.heatmap.pdf -c spearman -p heatmap --plotNumbers --removeOutliers --labels $(ls _bw-files/*_rep*bw | sed 's/.map.*//' | sed 's,_bw-files/,,' | tr '\n' ' ')"
    system(cmd1, wait = TRUE)
    system(cmd2, wait = TRUE)
    system(cmd3, wait = TRUE)
    cmd1 <- "multiBigwigSummary BED-file --BED ../_20180427_devATAC-JJ-paper/20180427.annot_Apr27/Fig2D2_regulatory_annotation_Apr27.bed -b $(ls _bw-files/*_rep*bw | tr '\n' ' ') -out deeptools.1000.bed.summary -bs 1000"
    cmd2 <- "plotCorrelation -in deeptools.1000.bed.summary -o deeptools.1000.bed.corr.pearson.heatmap.pdf -c pearson -p heatmap --plotNumbers --removeOutliers --labels $(ls _bw-files/*_rep*bw | sed 's/.map.*//' | sed 's,_bw-files/,,' | tr '\n' ' ')"
    cmd3 <- "plotCorrelation -in deeptools.1000.bed.summary -o deeptools.1000.bed.corr.spearman.heatmap.pdf -c spearman -p heatmap --plotNumbers --removeOutliers --labels $(ls _bw-files/*_rep*bw | sed 's/.map.*//' | sed 's,_bw-files/,,' | tr '\n' ' ')"
    system(cmd1, wait = TRUE)
    system(cmd2, wait = TRUE)
    system(cmd3, wait = TRUE)

    system('rm *deeptools*summary*')
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
if (!dir.exists('_log-files/')) {dir.create('_log-files/')}
if (!dir.exists('_bam-files/')) {dir.create('_bam-files/')}
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
        q10 = paste0('_bam-files/',  sample, '.map_se^rm_chrM^rm_blacklist^q10.bam'),
        track = paste0('_bw-files/', sample, '.map_se^rm_chrM^rm_blacklist^q10.bw'),
        peaks = paste0('_peaks-files/', sample, '.peaks.bed'),
        peaks_xls = paste0('_peaks-files/', sample, '.peaks.xls')
    )
}

library(parallel)
mclapply(SAMPLES, atac.se_raw2bw, list.samples = list.samples, steps = which.steps, ncore = 4, which.genome = '../_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa', mc.cores = cores)

# Then finish up the script (rename, stats, cleanup)
if ('cleaning' %in% which.steps) {
    write('\nCleaning-up temp. files...', file = list.samples[[sample]]$log, append = TRUE)
    for (sample in SAMPLES) { cleanTemp(sample) }
} else { write('\nSkipping cleaning up...', file = list.samples[[sample]]$log, append = TRUE) }

if ('stats' %in% which.steps) {
    write('\n------------------\n', file = "bamstats.log", append = TRUE)
    bamstats <- bamStats(list.samples)
    write.table(bamstats, file = 'bamstats.log', quote = F, col.names = F, row.names = F, sep = ': ', append = TRUE)
}

if ('combinetracks' %in% which.steps) {
    message('\nCombining the bigwig tracks for replicates...')
    for (sample in unique(gsub('_rep.', '', SAMPLES))) {
        message('Merging ', sample, '...')
        files <- grep(paste0(sample, '_rep'), list.files('_bw-files', full.names = T), value = T)
        cmd1 <- sprintf('bigWigMerge %s %s %s', files[1], files[2], paste0(sample, '.tmp.bg'))
        cmd2 <- sprintf('bedGraphToBigWig %s %s %s', paste0(sample, '.tmp.bg'), '../_data/ce11/ce11.chrom.sizes', paste0('_bw-files/', sample, '_combined.bw'))
        system(cmd1, wait = T)
        system(cmd2, wait = T)
    }
    system('rm *.tmp.bg', wait = T)
}

if ('deeptools' %in% which.steps) {
    write('\nQCing the generated files...', file = list.samples[[sample]]$log, append = TRUE)
    plotDeepToolsQC()
} else { write('\nSkipping QC...', file = list.samples[[sample]]$log, append = TRUE) }

#!/usr/bin/env Rscript

# Raw paired-end reads must be located immediatly in the _raw_data subfolder. They should be named XX.combined.rX.fq.gz.
# You must be in a folder directly created in ~/.
# which.genome, which.rrna, which.blacklist and which.genome.chrom.sizes are all files located in ~/_data.

#########################################################
#                                                       #
# Read arguments (smaples and steps to perform)         #
#                                                       #
#########################################################

library(optparse)
option_list = list(
  make_option("--samples", type="character", default='all',
              help="Which samples (JSsXXX) to map (separated by commas and between single quotes)", metavar="character"),
  make_option("--steps", type="character", default='combine,trim,aln,tracks,stats,cleaning',
              help="Which steps to perform (separated by commas and between single quotes)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

SAMPLES <- if (opt$samples == 'all') { gsub('-.*', '', list.files(path = '_raw_data/', pattern = '-')) } else { unlist(strsplit(opt$samples, ',')) }
which.steps <- if (opt$steps == 'all') { unlist(strsplit('combine,trim,aln,tracks,stats,cleaning', ',')) } else { unlist(strsplit(opt$steps, ',')) }

#########################################################
#                                                       #
# Master function                                       #
#                                                       #
#########################################################

lcap_raw2bw <- function(sample, list.samples, steps, l = 20, ncore = 8, which.genome = '../_data/ce10.fa', which.rrna = '../_data/WS253_ce10.rRNA_broad.bed', which.blacklist = '../_data/ce10_blacklist.bed', which.genome.chrom.sizes = '../_data/ce10.chrom.sizes') {

    suppressWarnings(suppressMessages(library(rtracklayer)))
    suppressWarnings(suppressMessages(library(GenomicAlignments)))
    suppressWarnings(suppressMessages(library(Rsamtools)))
    suppressWarnings(suppressMessages(library(BSgenome)))

    write(paste0('\n\n\n--------- Script ran @ ', date(), ' ---------'), file = paste0(sample, '.log'), append = TRUE)
    write(paste0('--------- Now starting mapping of sample ', sample, ' ---------'), file = paste0(sample, '.log'), append = TRUE)
    write(paste0('The sample that is going to be mapped is: \n', sample, '\nThe steps that are going to be exectuted are: \n', paste(steps, collapse = ', '), '\n'), file = paste0(sample, '.log'), append = TRUE)

    if ('combine' %in% steps) {
        write('\nCombining reads...', file = paste0(sample, '.log'), append = TRUE)
            combine.reads(sample, list.samples)
    } else { write('\nSkipping reads combining...', file = paste0(sample, '.log'), append = TRUE) }

    if ('trim' %in% steps) {
        write('\nTrimming reads...', file = paste0(sample, '.log'), append = TRUE)
            trim.reads(sample, list.samples, l = 50)
    } else { write('\nSkipping reads trimming...', file = paste0(sample, '.log'), append = TRUE) }

    if ('aln' %in% steps) {
        write('\nAligning reads...', file = paste0(sample, '.log'), append = TRUE)
            lcap.pe.bwa_aln(sample, list.samples, ncore = ncore, genome = which.genome, rrna = which.rrna, blacklist = which.blacklist)
    } else { write('\nSkipping reads alignment...', file = paste0(sample, '.log'), append = TRUE) }

    if ('tracks' %in% steps) {
        write('\nGenerating alignments tracks...', file = paste0(sample, '.log'), append = TRUE)
            pe.bam2filled.bw(sample, list.samples, genome.chrom.sizes = '../_data/ce10.chrom.sizes')
    } else { write('\nSkipping alignments tracks...', file = paste0(sample, '.log'), append = TRUE) }

    if ('stats' %in% steps) {
        write('\nGetting bam stats for the generated bam files', file = paste0(sample, '.log'), append = TRUE)
            bamstats <- bamStats(list.samples[[sample]]$aligned)
            write.table(bamstats, file = paste0(sample, '.log'), append = TRUE, quote = F, col.names = F, row.names = F, sep = ': ')
    } else { write('\nSkipping tracks stats...', file = paste0(sample, '.log'), append = TRUE) }

    if ('cleaning' %in% steps) {
        write('\nCleaning up temp. files', file = paste0(sample, '.log'), append = TRUE)
            cleanTemp()
    } else { write('\nSkipping cleaning up...', file = paste0(sample, '.log'), append = TRUE) }


}

#########################################################
#                                                       #
# Other functions                                       #
#                                                       #
#########################################################

# Combine reads for each sample
combine.reads <- function(sample, list.samples) {

    sample.files <- grep('combined', grep('R1_001.fastq.gz', list.files(pattern = sample, recursive = T), value = T), invert = T, value = T)
    cmd <- sprintf('cat %s %s > %s', sample.files[1], sample.files[2], list.samples[[sample]]$combined.r1)
    system(cmd, wait = TRUE)

    sample.files <- grep('combined', grep('R2_001.fastq.gz', list.files(pattern = sample, recursive = T), value = T), invert = T, value = T)
    cmd <- sprintf('cat %s %s > %s', sample.files[1], sample.files[2], list.samples[[sample]]$combined.r2)
    system(cmd, wait = TRUE)

}

# Trim paired reads to length n
trim.reads <- function(sample, list.samples, l = 20) {

    cmd <- sprintf('zcat -d -c %s | fastx_trimmer -z -l %i -Q 33 -o %s', list.samples[[sample]]$combined.r1, l, list.samples[[sample]]$trimmed.r1)
    system(cmd, wait = TRUE)

    cmd <- sprintf('zcat -d -c %s | fastx_trimmer -z -l %i -Q 33 -o %s', list.samples[[sample]]$combined.r2, l, list.samples[[sample]]$trimmed.r2)
    system(cmd, wait = TRUE)
}

# Align reads to ce10 and clean bam file (remove chrM, rRNA and blacklist regions)
lcap.pe.bwa_aln <- function(sample, list.samples = list.samples, genome = '../_data/ce10.fa', rrna = '../_data/WS253_ce10.rRNA_broad.bed', blacklist = '../_data/ce10_blacklist.bed', ncore = 8) {

    cmd1 <- sprintf("bwa aln -t %i %s %s > %s", ncore, genome, list.samples[[sample]]$trimmed.r1, paste0(sample, '.r1.sai'))
    cmd2 <- sprintf("bwa aln -t %i %s %s > %s", ncore, genome, list.samples[[sample]]$trimmed.r2, paste0(sample, '.r2.sai'))
    cmd3 <- sprintf("bwa sampe %s %s %s %s %s | samtools view -b -@ %i - | samtools sort -@ %i -O sam -T %s -n - | samtools fixmate -O sam - - | samtools sort -@ %i -O bam -T %s - > %s", genome, paste0(sample, '.r1.sai'), paste0(sample, '.r2.sai'), list.samples[[sample]]$trimmed.r1, list.samples[[sample]]$trimmed.r2, ncore, ncore, paste0(list.samples[[sample]]$trimmed.r1, '.tmp'), ncore, paste0(list.samples[[sample]]$trimmed.r2, '.tmp'), list.samples[[sample]]$aligned)
    cmd4 <- sprintf("samtools view -@ %i -b -f 3 -F 12 %s > %s", ncore, list.samples[[sample]]$aligned, list.samples[[sample]]$no_unmapped)
    cmd5 <- sprintf("samtools index -@ %i %s ; samtools view -@ %i -b %s chrI chrII chrIII chrIV chrV chrX > %s ; rm %s", ncore, list.samples[[sample]]$no_unmapped, ncore, list.samples[[sample]]$no_unmapped, list.samples[[sample]]$noChrM, paste0(list.samples[[sample]]$no_unmapped, '.bai'))
    cmd6 <- sprintf("samtools view -@ %i -u %s | samtools view -b -L %s -U %s - > /dev/null", ncore, list.samples[[sample]]$noChrM, rrna, list.samples[[sample]]$norRNA)
    cmd7 <- sprintf("samtools view -@ %i -b -L %s -U %s %s > /dev/null", ncore, blacklist, list.samples[[sample]]$noBlacklist, list.samples[[sample]]$norRNA)
    cmd8 <- sprintf("samtools view -@ %i -b -q 10 %s > %s", ncore, list.samples[[sample]]$noBlacklist, list.samples[[sample]]$q10)

    system(cmd1, wait = T)
    system(cmd2, wait = T)
    system(cmd3, wait = T)
    system(cmd4, wait = T)
    system(cmd5, wait = T)
    system(cmd6, wait = T)
    system(cmd7, wait = T)
    system(cmd8, wait = T)

}

# Generate spmr-normalized filled fwd and rev tracks
pe.bam2filled.bw <- function(sample, list.samples = list.samples, genome.chrom.sizes = '../_data/ce10.chrom.sizes') {

        cmd1 <- sprintf("samtools view -b -u %s | bedtools genomecov -ibam stdin -bg -pc -strand - > %s", list.samples[[sample]]$q10, paste0(sample, '.fwd.bg'))
        cmd2 <- sprintf("sort -k1,1 -k2,2n %s -o %s", paste0(sample, '.fwd.bg'), paste0(sample, '.fwd.bg'))
        cmd3 <- sprintf("bedGraphToBigWig %s %s %s", paste0(sample, '.fwd.bg'), genome.chrom.sizes, list.samples[[sample]]$fwd.bw)
        system(cmd1, wait = T)
        system(cmd2, wait = T)
        system(cmd3, wait = T)

        cmd1 <- sprintf("samtools view -b -u %s | bedtools genomecov -ibam stdin -bg -pc -strand + > %s", list.samples[[sample]]$q10, paste0(sample, '.rev.bg'))
        cmd2 <- sprintf("sort -k1,1 -k2,2n %s -o %s", paste0(sample, '.rev.bg'), paste0(sample, '.rev.bg'))
        cmd3 <- sprintf("bedGraphToBigWig %s %s %s", paste0(sample, '.rev.bg'), genome.chrom.sizes, list.samples[[sample]]$rev.bw)
        system(cmd1, wait = T)
        system(cmd2, wait = T)
        system(cmd3, wait = T)

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

    df <- as.data.frame(matrix(NA, ncol = 2), byrow=T)
    df[,1] <- names(out)
    df[,2] <- out

    return(out)

}

# Clean up temporary files
cleanTemp <- function() {

    cmd <- sprintf("rm %s %s %s %s %s %s", paste0('*', sample, '*bg*'), paste0('*', sample, '*sai*'), paste0('*', sample, '*tmp*'), paste0('_raw_data/*', sample, '*tmp*'), paste0('_raw_data/*', sample, '*trimmed*'), paste0('_raw_data/*', sample, '*fq.gz*'), paste0('_raw_data/*', sample, '*combined*'))
    system(cmd, wait = T)

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
dir.create('_bw-files/rev_tracks/')
dir.create('_bw-files/fwd_tracks/')

list.samples <- list()
for (sample in SAMPLES) {
    dir.create(paste0('_bam-files/', sample))
    list.samples[[sample]] <- list(
        combined.r1 = paste0("_raw_data/", sample, '.combined.r1.fq.gz'),
        combined.r2 = paste0("_raw_data/", sample, '.combined.r2.fq.gz'),
        trimmed.r1 = paste0("_raw_data/", sample, '.combined^trimmed.r1.fq'),
        trimmed.r2 = paste0("_raw_data/", sample, '.combined^trimmed.r2.fq'),
        aligned = paste0("_bam-files/", sample, '/', sample, '.combined^trimmed^map_pe.bam'),
        no_unmapped = paste0("_bam-files/", sample, '/', sample, '.combined^trimmed^map_pe^rm_unmapped.bam'),
        noChrM = paste0("_bam-files/", sample, '/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM.bam'),
        norRNA = paste0("_bam-files/", sample, '/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA.bam'),
        noBlacklist = paste0("_bam-files/", sample, '/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist.bam'),
        q10 = paste0("_bam-files/", sample, '/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist^q10.bam'),
        fwd.bw = paste0('_bw-files/fwd_tracks/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist^q10.fwd.bw'),
        rev.bw = paste0('_bw-files/rev_tracks/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist^q10.rev.bw')
    )
}

library(parallel)
mclapply(SAMPLES, lcap_raw2bw, list.samples = list.samples, steps = which.steps, l = 20, ncore = 4, which.genome = '../_data/ce10.fa', which.rrna = '../_data/WS253_ce10.rRNA_broad.bed', which.blacklist = '../_data/ce10_blacklist.bed', which.genome.chrom.sizes = '../_data/ce10.chrom.sizes', mc.cores=6)

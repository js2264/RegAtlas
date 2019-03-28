#!/usr/bin/env Rscript

# which.genome files located in ~/_data.
# Raw pooled replicates paired-end sequencing results should be in _raw_data, labelled myo3_repX_read1.combined.fq.gz / myo3_repX_read2.combined.fq.gz

#########################################################
#                                                       #
# Read arguments (samples and steps to perform)         #
#                                                       #
#########################################################

library(optparse)
option_list = list(
  make_option("--samples", type="character", default='all',
              help="Which samples (xxx_repX) to map (separated by commas and between single quotes)", metavar="character"),
  make_option("--steps", type="character", default='all',
              help="Which steps to perform (separated by commas and between single quotes)", metavar="character"),
  make_option("--cores", type="integer", default=5,
              help="How many cores to use?", metavar="character"),
  make_option("--doCombineReps", type="logical", default=F,
              help="Should the raw files be combined first?", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

doCombineReps <- opt$doCombineReps
if (doCombineReps == T) { system("Rscript ../_scripts/bin/master-scripts/LCAP_pooling.R", wait = T) }

SAMPLES <- if (opt$samples == 'all') { unique(gsub('.combined.r..fq.gz', '', list.files(path = '_raw_data/', pattern = 'rep..combined.r'))) } else { unlist(strsplit(opt$samples, ',')) }
which.steps <- if (opt$steps == 'all') { unlist(strsplit('trim,aln,tracks,stats', ',')) } else { unlist(strsplit(opt$steps, ',')) }
cores <- opt$cores

#########################################################
#                                                       #
# Master function                                       #
#                                                       #
#########################################################

lcap_raw2bw <- function(sample, list.samples, steps, l = 20, ncore = 8, which.genome = '../_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa', which.rrna = '../_data/ce11/WS253_ce10.rRNA_broad.bed', which.blacklist = '../_data/ce11/ce11_blacklist.bed', which.genome.chrom.sizes = '../_data/ce11/ce11.chrom.sizes') {

    suppressWarnings(suppressMessages(library(rtracklayer)))
    suppressWarnings(suppressMessages(library(GenomicAlignments)))
    suppressWarnings(suppressMessages(library(Rsamtools)))
    suppressWarnings(suppressMessages(library(BSgenome)))

    write(paste0('\n\n\n--------- Script ran @ ', date(), ' ---------'), file = list.samples[[sample]]$log, append = TRUE)
    write(paste0('--------- Now starting mapping of sample ', sample, ' ---------'), file = list.samples[[sample]]$log, append = TRUE)
    write(paste0('The sample that is going to be mapped is: \n', sample, '\nThe steps that are going to be exectuted are: \n', paste(steps, collapse = ', '), '\n'), file = list.samples[[sample]]$log, append = TRUE)

    if ('trim' %in% steps) {
        write('\nTrimming reads...', file = list.samples[[sample]]$log, append = TRUE)
            trim.reads(sample, list.samples, l = 20)
    } else { write('\nSkipping reads trimming...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('aln' %in% steps) {
        write('\nAligning reads...', file = list.samples[[sample]]$log, append = TRUE)
            lcap.pe.bwa_aln(sample, list.samples, ncore = ncore, genome = which.genome, rrna = which.rrna, blacklist = which.blacklist)
    } else { write('\nSkipping reads alignment...', file = list.samples[[sample]]$log, append = TRUE) }

    if ('tracks' %in% steps) {
        write('\nGenerating alignments tracks...', file = list.samples[[sample]]$log, append = TRUE)
            pe.bam2filled.bw(sample, list.samples, genome.chrom.sizes = which.genome.chrom.sizes)
    } else { write('\nSkipping alignments tracks...', file = list.samples[[sample]]$log, append = TRUE) }

}

#########################################################
#                                                       #
# Other functions                                       #
#                                                       #
#########################################################

# Trim paired reads to length n
trim.reads <- function(sample, list.samples, l = 20) {

    cmd <- sprintf('zcat -d -c %s | fastx_trimmer -z -l %i -Q 33 -o %s', list.samples[[sample]]$combined.r1, l, list.samples[[sample]]$trimmed.r1)
    system(cmd, wait = TRUE)

    cmd <- sprintf('zcat -d -c %s | fastx_trimmer -z -l %i -Q 33 -o %s', list.samples[[sample]]$combined.r2, l, list.samples[[sample]]$trimmed.r2)
    system(cmd, wait = TRUE)
}

# Align reads to ce11 and clean bam file (remove chrM, rRNA and blacklist regions)
lcap.pe.bwa_aln <- function(sample, list.samples = list.samples, genome = '../_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa', rrna = '../_data/ce11/WS253_ce10.rRNA_broad.bed', blacklist = '../_data/ce11/ce11_blacklist.bed', ncore = 8) {

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

# Generate FPM-normalized filled fwd and rev tracks (as well as )
pe.bam2filled.bw <- function(sample, list.samples = list.samples, genome.chrom.sizes = '../_data/ce11/ce11.chrom.sizes') {

    message('>>> Generating coverage tracks')
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

    message('>>> Generating FPM-scaled coverage tracks')
    NBPAIRS <- as.numeric(system(sprintf("samtools view -F 16 -c %s", list.samples[[sample]]$q10), intern = T))
    message('>>>>>> Scaling factor for fwd strand:', 1e6/NBPAIRS)
    cmd1 <- sprintf("samtools view -b -u %s | bedtools genomecov -ibam stdin -bg -scale %f -pc -strand - > %s", list.samples[[sample]]$q10, 1e6/NBPAIRS, paste0(sample, '.fwd.bg'))
    cmd2 <- sprintf("sort -k1,1 -k2,2n %s -o %s", paste0(sample, '.fwd.bg'), paste0(sample, '.fwd.bg'))
    cmd3 <- sprintf("bedGraphToBigWig %s %s %s", paste0(sample, '.fwd.bg'), genome.chrom.sizes, list.samples[[sample]]$fwd.FPM.bw)
    system(cmd1, wait = T)
    system(cmd2, wait = T)
    system(cmd3, wait = T)
    NBPAIRS <- as.numeric(system(sprintf("samtools view -F 32 -c %s", list.samples[[sample]]$q10), intern = T))
    message('>>>>>> Scaling factor for fwd strand:', 1e6/NBPAIRS)
    cmd1 <- sprintf("samtools view -b -u %s | bedtools genomecov -ibam stdin -bg -scale %f -pc -strand + > %s", list.samples[[sample]]$q10, 1e6/NBPAIRS, paste0(sample, '.rev.bg'))
    cmd2 <- sprintf("sort -k1,1 -k2,2n %s -o %s", paste0(sample, '.rev.bg'), paste0(sample, '.rev.bg'))
    cmd3 <- sprintf("bedGraphToBigWig %s %s %s", paste0(sample, '.rev.bg'), genome.chrom.sizes, list.samples[[sample]]$rev.FPM.bw)
    system(cmd1, wait = T)
    system(cmd2, wait = T)
    system(cmd3, wait = T)

    message('>>> Generating zscore-scaled coverage tracks')
    cov <- import.bw(list.samples[[sample]]$fwd.bw, as='RleList')
        ucov <- unlist(cov)
        mi <- mean(ucov)
        mu <- sd(ucov)
        zsc <- (cov-mi)/mu
    export.bw(zsc, list.samples[[sample]]$fwd.ZSCORE.bw)
    cov <- import.bw(list.samples[[sample]]$rev.bw, as='RleList')
        ucov <- unlist(cov)
        mi <- mean(ucov)
        mu <- sd(ucov)
        zsc <- (cov-mi)/mu
    export.bw(zsc, list.samples[[sample]]$rev.ZSCORE.bw)

}

# Get some bam stats for all bam files (from JADBtools ; must be done outside of master function)
bamStats <- function(list.samples) {

    files <- unlist(lapply(list.samples, function(x) {'[['(x,'aligned')}))

    load('../_data/rrnamodel.rda')

    out <- list()

    for (f in files) {
        
        message('\t>> Bamstats for: ', f)
        what <- c("rname", "strand", "pos", "qwidth", "mapq", "flag", "mpos", "isize")
        a <- scanBam(f, param = ScanBamParam(what = what), asMates = T)[[1]]
        flags <- bamFlagAsBitMatrix(a$flag)
        align <- !is.na(a$pos)
        grng <- GRanges(
            seqnames = a$rname[align],
            ranges = IRanges(a$pos[align], width = a$qwidth[align]),
            strand = a$strand[align],
            qual = a$mapq[align]
        )
        map_pe <- a[[1]] # All individual reads
        rm_unmapped <- grng[rowSums(flags[align,1:2]) == 2] # All individual reads mapping the genome in proper pairs (-f3)
        rm_chrM <- rm_unmapped[seqnames(rm_unmapped) != 'chrM'] # All individual reads mapping the genome in proper pairs (-f3) but not on chrM 
        rm_rRNA <- rm_chrM[!(rm_chrM %over% rtracklayer::import.bed('../_data/ce11/WS253_ce10.rRNA_broad.bed'))] # All individual reads mapping the genome in proper pairs (-f3) but not on chrM nor on rRNA clusters nor on blacklisted regions
        rm_blacklist <- rm_rRNA[!(rm_rRNA %over% rtracklayer::import.bed('../_data/ce11/ce11_blacklist.bed'))] # All individual reads mapping the genome in proper pairs (-f3) but not on chrM nor on rRNA clusters nor on blacklisted regions
        q10 <- rm_blacklist[rm_blacklist$qual >= 10L] # All individual reads mapping the genome in proper pairs (-f3) but not on chrM nor on rRNA clusters nor on blacklisted regions, with a auqlity score > 10
        
        out[f] <- sprintf(
            'all=%.2fM, -unpaired=%.2fM [%.0f%% of total], -chrM=%.2fM [%.0f%% of total], -rRNA=%.2fM [%.0f%% of total], -blacklist=%.2fM [%.0f%% of total], q10=%.2fM [%.0f%% of total]', 
            length(map_pe)/10^6, length(rm_unmapped)/10^6, (length(rm_unmapped)/length(map_pe))*100, length(rm_chrM)/10^6, (length(rm_chrM)/length(map_pe))*100, length(rm_rRNA)/10^6, (length(rm_rRNA)/length(map_pe))*100, length(rm_blacklist)/10^6, (length(rm_blacklist)/length(map_pe))*100, length(q10)/10^6, (length(q10)/length(map_pe))*100
        )

        df <- as.data.frame(matrix(NA, ncol = 2, nrow = length(out)), byrow=T)
        df[,1] <- names(out)
        df[,2] <- unlist(out)
        
    }
    
    return(df)

}

# Clean up temporary files
cleanTemp <- function(sample) {

    cmd <- sprintf("rm %s %s %s %s %s", paste0('*', sample, '*bg*'), paste0('*', sample, '*sai*'), paste0('*', sample, '*tmp*'), paste0('_raw_data/*', sample, '*tmp*'), paste0('_raw_data/*', sample, '*trimmed*'))
    system(cmd, wait = T)

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
if (!dir.exists('_bw-files/cov_fwd_tracks/')) {dir.create('_bw-files/cov_fwd_tracks/')}
if (!dir.exists('_bw-files/cov_rev_tracks/')) {dir.create('_bw-files/cov_rev_tracks/')}
if (!dir.exists('_bw-files/FPM_fwd_tracks/')) {dir.create('_bw-files/FPM_fwd_tracks/')}
if (!dir.exists('_bw-files/FPM_rev_tracks/')) {dir.create('_bw-files/FPM_rev_tracks/')}
if (!dir.exists('_bw-files/ZSCORE_fwd_tracks/')) {dir.create('_bw-files/ZSCORE_fwd_tracks/')}
if (!dir.exists('_bw-files/ZSCORE_rev_tracks/')) {dir.create('_bw-files/ZSCORE_rev_tracks/')}

list.samples <- list()
for (sample in SAMPLES) {
    dir.create(paste0('_bam-files/', sample), showWarnings = F)
    list.samples[[sample]] <- list(
        log = paste0('_log-files/', sample, '.log'),
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
        fwd.bw = paste0('_bw-files/cov_fwd_tracks/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist^q10.cov.fwd.bw'),
        rev.bw = paste0('_bw-files/cov_rev_tracks/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist^q10.cov.rev.bw'),
        fwd.FPM.bw = paste0('_bw-files/FPM_fwd_tracks/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist^q10.FPM.fwd.bw'),
        rev.FPM.bw = paste0('_bw-files/FPM_rev_tracks/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist^q10.FPM.rev.bw'),
        fwd.ZSCORE.bw = paste0('_bw-files/ZSCORE_fwd_tracks/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist^q10.ZSCORE.fwd.bw'),
        rev.ZSCORE.bw = paste0('_bw-files/ZSCORE_rev_tracks/', sample, '.combined^trimmed^map_pe^rm_unmapped^rm_chrM^rm_rRNA^rm_blacklist^q10.ZSCORE.rev.bw')
    )
}

library(parallel)
mclapply(SAMPLES, lcap_raw2bw, list.samples = list.samples, steps = which.steps, l = 20, ncore = 4, which.genome = '../_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa', which.rrna = '../_data/ce11/WS253_ce10.rRNA_broad.bed', which.blacklist = '../_data/ce11/ce11_blacklist.bed', which.genome.chrom.sizes = '../_data/ce11/ce11.chrom.sizes', mc.cores=cores)

# Then finish up the script (stats, cleanup)
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
        message('> > >   Merging ', sample, '...')
        for (STRAND in c('fwd', 'rev')) {
            files <- grep(paste0(sample, '_rep'), list.files(paste0('_bw-files/cov_', STRAND, '_tracks'), full.names = T), value = T)
            cmd1 <- sprintf('bigWigMerge %s %s %s', files[1], files[2], paste0(sample, '.tmp.bg'))
            cmd2 <- sprintf('bedGraphToBigWig %s %s %s', paste0(sample, '.tmp.bg'), '../_data/ce11/ce11.chrom.sizes', paste0('_bw-files/', sample, '_combined.', STRAND, '.cov.bw'))
            system(cmd1, wait = T)
            system(cmd2, wait = T)
            files <- grep(paste0(sample, '_rep'), list.files(paste0('_bw-files/FPM_', STRAND, '_tracks'), full.names = T), value = T)
            cmd1 <- sprintf('bigWigMerge %s %s %s', files[1], files[2], paste0(sample, '.tmp.bg'))
            cmd2 <- sprintf('bedGraphToBigWig %s %s %s', paste0(sample, '.tmp.bg'), '../_data/ce11/ce11.chrom.sizes', paste0('_bw-files/', sample, '_combined.', STRAND, '.FPM.bw'))
            system(cmd1, wait = T)
            system(cmd2, wait = T)
            files <- grep(paste0(sample, '_rep'), list.files(paste0('_bw-files/ZSCORE_', STRAND, '_tracks'), full.names = T), value = T)
            cmd1 <- sprintf('bigWigMerge %s %s %s', files[1], files[2], paste0(sample, '.tmp.bg'))
            cmd2 <- sprintf('bedGraphToBigWig %s %s %s', paste0(sample, '.tmp.bg'), '../_data/ce11/ce11.chrom.sizes', paste0('_bw-files/', sample, '_combined.', STRAND, '.ZSCORE.bw'))
            system(cmd1, wait = T)
            system(cmd2, wait = T)
        }
    }
    system('rm *.tmp.bg', wait = T)
}

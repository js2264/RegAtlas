checkInsertSize <- function(
    pairs, 
    n_reads = 1000000, 
    genome = '~/_data/ce11/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa', 
    blacklist = '~/_data/ce11/ce10_blacklist.bed',
    ncore = 10
)
{
    dir.create('insert_sizes_sandbox/', showWarnings = FALSE)
    setwd('insert_sizes_sandbox/')
    l <- lapply(1:nrow(pairs), function(K) {
        read1_file <- pairs$read1[K]
        read2_file <- pairs$read2[K]
        system(paste0('gunzip -k -c ../', read1_file, ' | head -n ', formatC(n_reads, format = 'fg'), ' > read1.fq'))
        system(paste0('gunzip -k -c ../', read2_file, ' | head -n ', formatC(n_reads, format = 'fg'), ' > read2.fq'))
        system(sprintf("bwa mem -t %i -T %i %s %s %s | samtools view -b -@ %i - | samtools sort -@ %i -m 44294967296 -O sam -T %s -n - | samtools fixmate -O sam - - | samtools sort -@ %i -m 44294967296 -O bam -T %s - > %s", ncore, 0, genome, 'read1.fq', 'read2.fq', ncore, ncore, 'read1.fq.tmp', ncore, 'read1.fq.tmp', 'out.bam'))
        system(sprintf('samtools index -@ %i %s ; samtools view -@ %i -b %s chrI chrII chrIII chrIV chrV chrX > %s ; rm %s', ncore, 'out.bam', ncore, 'out.bam', 'out^rmChrM.bam', 'out.bam.bai'))
        system(sprintf('samtools view -@ %i -b -L %s -U %s %s > /dev/null', ncore, blacklist, 'out^rmChrM^rmBlacklist.bam', 'out^rmChrM.bam'))
        system(sprintf('samtools view -@ %i -b -q 10 %s > %s', ncore, 'out^rmChrM^rmBlacklist.bam', 'out^rmChrM^rmBlacklist^q10.bam'))
        bam <- Rsamtools::BamFile("out^rmChrM^rmBlacklist^q10.bam", yieldSize = 50000000, asMates = TRUE)
        params <- Rsamtools::ScanBamParam(
            what = c("qname", "mapq", "isize", "flag"), 
            flag = Rsamtools::scanBamFlag(
                isPaired = TRUE, 
                isProperPair = TRUE, 
                isSecondaryAlignment = FALSE,
                isUnmappedQuery = FALSE,
                isNotPassingQualityControls = FALSE,
                isSupplementaryAlignment = FALSE
            )
        )
        a <- GenomicAlignments::readGAlignmentPairs(bam, param = params)
        g <- GenomicAlignments::granges(a)
        res <- list(
            "average" = mean(GenomicRanges::width(g)), 
            "stdev" = sd(GenomicRanges::width(g)) 
        )
        return(res)
    })
    df <- do.call(rbind, l)
    setwd('../')
    unlink('insert_sizes_sandbox/', recursive = TRUE)
    return(df)
}

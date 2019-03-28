#!/usr/bin/env Rscript

## This script contains commands to perform DESeq2 from bam files.
## The master function 'bam2DE()' can be used from another Rscript.
## It counts reads (using summarizeOverlaps for ATAC and featureCounts for lcap)
## It then performs DESeq2 analysis, on a chosen desgin (all combinations of tissues)
## It saves results in _DESeq2-files/ folder

suppressWarnings(suppressMessages(library(rtracklayer)))
suppressWarnings(suppressMessages(library(GenomicAlignments)))
suppressWarnings(suppressMessages(library(Rsamtools)))
suppressWarnings(suppressMessages(library(BSgenome)))
suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(parallel)))

#########################################################
#                                                       #
# Read arguments (samples and steps to perform)         #
#                                                       #
#########################################################

###

FC = 3
FDR = 0.01

##

#########################################################
#                                                       #
# Master function                                       #
#                                                       #
#########################################################

bam2DE <- function(sampleTable.path, which.save.path, which.results.path, which.model, which.design, is.stranded, is.paired, filter, value.to.compute) {

    message('\nLoading sample metadata...')
        sampleTable <- read.csv(sampleTable.path)
        message(paste0('Sample table correctly loaded\n', nrow(sampleTable), ' bam files pre-loaded'))
        if (length(filter) > 0) {
            print(filter)
            sampleTable <- sampleTable[sampleTable[filter[1]] == filter[2],]
            sampleTable['which.design'] <- sampleTable[which.design]
            message(paste0('\tApplying filters ', filter[2], ' on column ', filter[1], ': sample table filtered down to ', nrow(sampleTable), ' bam files'))
            which.contrasts = getContrasts(sampleTable, which.design)
            message(paste0('\t', length(which.contrasts),' contrasts found using design: ', which.design))
        }

    message('\nChecking files...')
        if (checkBamFiles(sampleTable) == T) {
            message('Bam files loading...')
            b.files <- BamFileList(as.character(sampleTable$bamReads))
            message(paste0(length(b.files), ' bam files correctly loaded'))
        }

    message('\nChoosing model...')
        model <- chooseModel(which.model)

    message('\nCounting reads...')
        e <- summarizeBAMs(b.files, model, which.model, sampleTable, sampleTable.path)

    message(paste0('\nGetting design informations from chosen variable: ', which.design, '...'))
        colData(e) <- DataFrame(sampleTable)
        colData(e)$chosen.design <- as.factor(unlist(as.data.frame(colData(e)[which.design])))
        colnames(e) <- gsub("\\.map.*", "", names(b.files))
        featureLength <- if (grepl('LCAP', which.results.path)) { read.table('_results-files/LCAP.featureCounts/LCAP.featureCounts.tsv', header = T, comment.char = '#', row.names = "Geneid")$Length } else if (grepl('ATAC', which.results.path)) { rep(150, times = nrow(e)) }
        mode <- ifelse(grepl('LCAP', sampleTable.path), 'one_vs_others', 'one_vs_one')

    message('\nRunning DEseq2...')
        results <- doDiffExpr(e, contrasts = which.contrasts, model = model, results.path = which.results.path, featureLength, mode, value.to.compute)

    message('\nSaving results...')
        save.image(which.save.path)
        message(paste0('\tDESeq2 analysis saved in ', which.save.path))
        if (file.exists('.tmp.DESeq2.RData')) {file.remove('.tmp.DESeq2.RData')}

    return(results)

}

#########################################################
#                                                       #
# Other functions                                       #
#                                                       #
#########################################################

# Check that all the bam files exist
checkBamFiles <- function(sampleTable) {
    if ( any(!file.exists(as.character(sampleTable$bamReads))) ) {
        stop('\nSome or all the bam files are missing. Stopping script.\n\n-----------------------')
    } else { T }
}

# Choose which model to use to count reads (load from object or bed file)
chooseModel <- function(model = NULL) {

    require(rtracklayer)

    if (class(model)[1] %in% c("GRangesList","GRanges")) {
        message(paste0('\nModel loaded: ', deparse(substitute(model))))
        which.model <- model
    } else if (is.null(model)) {
        answer <- readline(prompt="Do you want to load default gene model (_data/gnmodel_Przemek_20170531.rds) ? (y/n) :")
        if (answer == "y") {
            which.model <- readRDS('~/shared/data/gnmodel_Przemek_20170531.rds')
        } else {
            message('\nChoose which bed file to use for regions to compare: ')
            which.model <- file.choose()
        }
    } else if (is.character(model)) {
        if (model == 'gnmodel') {
            which.model <- readRDS('~/shared/data/gnmodel_Przemek_20170531.rds')
        } else if (grepl('.bed', model)) {
            message(paste0('\nReading model from bed file: ', model))
            which.model <- model
        } else if (grepl('.gtf', model)) {
            message(paste0('\nReading model from gtf file: ', model))
            which.model <- model
        }
    }
    if (class(which.model)[1] == "character") {if (file.exists(which.model)) {
        if (grepl(".bed", which.model)) {
            x <- read.table(which.model, stringsAsFactors = F)
            if (class(x[,1]) == 'character' & class(x[,2]) == 'integer' & class(x[,3]) == 'integer') {
                library(rtracklayer)
                message('\nBed file is valid and loading\n')
                write.table(file = 'tmp.model.bed', x[,c(1:3)], sep="\t", col.names=F, row.names=F, quote=F)
                which.model <- import(con='tmp.model.bed')
                if (file.exists('tmp.model.bed')) file.remove('tmp.model.bed')
            }
        } else if (grepl(".gtf", which.model)) {
            message('\nGtf file is valid and loading\n')
            which.model <- import(which.model)
        }
    } else { stop("Model file not found. Stopping script.") }}
    if (is.null(names(which.model))) {
        if (grepl(".bed", model)) { names(which.model) <- paste(as.character(seqnames(which.model)), start(ranges(which.model)), end(ranges(which.model)), sep = "_") }
        if (grepl(".gtf", model) & "gene_id" %in% names(mcols(which.model))) { names(which.model) <- mcols(which.model)$gene_id }
    }
    if (class(which.model)[1] == 'GRangesList' | class(which.model)[1] == 'GRanges') { message('Model successfully loaded') ; return(which.model) } else { stop('\nModel is invalid. Stopping script.\n\n-----------------------') }
}

# Count reads from bam files, either in single- or paired-end mode
summarizeBAMs <- function(files, model, which.model, sampleTable, sampleTable.path) {
    if (grepl('LCAP', sampleTable.path)) {
        require(GenomicFeatures)
        message('Lcap data is detected. Proceeding to counts using featureCounts.\n')
        dir.create("_results-files/LCAP.featureCounts/", showWarnings = F)
        cmd <- sprintf('featureCounts -t gene -s 2 -Q 10 -T 30 -p -a %s -o %s %s', which.model, paste('_results-files/LCAP.featureCounts/LCAP.featureCounts.tsv'), paste(sampleTable$bamReads, collapse = ' '))
        system(cmd, wait = T)
        featurecounts <- read.table('_results-files/LCAP.featureCounts/LCAP.featureCounts.tsv', header = T, comment.char = '#', row.names = "Geneid")
        colnames(featurecounts)[-c(1:5)] <- gsub(".combined.*", "", names(files))
        counts.summary <- DESeqDataSetFromMatrix(countData = featurecounts[,-c(1:5)], colData = DataFrame(sampleTable), design = ~which.design )
        genes.jj <- import('~/shared/data/ce11.genes.gtf')
        mcols(counts.summary) <- genes.jj[match(row.names(counts.summary), genes.jj$gene_id),]
        #rowRanges(counts.summary) <- exonsBy(makeTxDbFromGFF(which.model), by = 'gene')
    } else if (grepl('ATAC', sampleTable.path)) {
        message('ATAC data is detected. Proceeding to counts using summarizeOverlaps [summarizeOverlaps(model, BamFileList(files), ignore.strand = T, singleEnd = T, fragments = F)]\n')
        counts.summary <- summarizeOverlaps(model, BamFileList(files), ignore.strand = T, singleEnd = T, fragments = F)
    }
    return(counts.summary)
}

# Get list of contrasts
getContrasts <- function(sampleTable, which.design) {
    contrasts.all <- as.character(unlist(unique(sampleTable[which.design])))
    list.contrasts <- combn(contrasts.all, 2, simplify=F)
    list.contrasts <- lapply(list.contrasts, function(x) c("chosen.design", x[2], x[1]))

    return(list.contrasts)
}

# Check significant FC (l2FC > FC or < -FC and corresponding padj < FDR)
checkSignificance <- function(ROW) {
    which.l2fc <- which(results3[ROW,grepl("l2FC", colnames(results3))] > FC | results3[ROW,grepl("l2FC", colnames(results3))] < -FC)
    which.padj <- which(results3[ROW,grepl("padj", colnames(results3))] < FDR)
    is.sign <- any(which.l2fc == which.padj)
    which.sign <- colnames(results3[ROW,grepl("l2FC", colnames(results3))])[which.padj[which.l2fc == which.padj]]

    return(list(is.sign, which.sign))
}

# TPM function (from slowkow github)
counts_to_tpm <- function(counts, featureLength) {

  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))

  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength # This is taking all the REs in account (each RE is 150bp)
  }))

  # Exclude genes with efflength less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]

  # Compute TPM (process one column at a time) [this function switches to log then exp for cpmputation ease]
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(COL) {
    rate = log(counts[,COL]) - log(effLen[,COL])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

# Get DE and TPM using DESeq2, then publish summary
doDiffExpr <- function(e, contrasts, model, results.path, featureLength, mode, value.to.compute) {

    if (mode == 'one_vs_others') {

        # Perform DE analysis for each
        designs <- as.character(unique(colData(e)$which.design))
        out.l2fc <- data.frame(matrix(nrow=nrow(e), ncol=length(designs))) ; colnames(out.l2fc) <- paste0(designs, '.l2FC')
        out.padj <- data.frame(matrix(nrow=nrow(e), ncol=length(designs))) ; colnames(out.padj) <- paste0(designs, '.padj')
        
        for (DESIGN in designs) {
            message('\tPerforming DE analysis for design ', which(designs == DESIGN), ' / ', length(designs))
            colData(e)$chosen.design <- factor(ifelse(colData(e)$which.design == DESIGN, DESIGN, paste0('non.', DESIGN)))
            dds <- DESeqDataSet(e, design = ~ chosen.design)
            dds <- DESeq(dds)
            #dds <- DESeq(dds, betaPrior = F)
            ddr <- results(dds)
            out.l2fc[,paste0(DESIGN, '.l2FC')] <- -ddr$log2FoldChange
            out.padj[,paste0(DESIGN, '.padj')] <- ddr$padj
        }
        
        # Reprocess all the data with usual design values ('one_vs_one') to get appropriate gene expression estimates
        colData(e)$chosen.design <- as.factor(unlist(as.data.frame(colData(e)$which.design)))
        dds <- DESeqDataSet(e, design = ~ chosen.design)
        dds <- DESeq(dds)
        
    } else if (mode == 'one_vs_one') {
    
        # Perform DE analysis
        message('\tGenerating dds object')
        dds <- DESeqDataSet(e, design = ~ chosen.design)
        dds <- DESeq(dds)
        
        # Get all results
        message('\nFetching data for each contrast...')
        contrasts.names <- unlist(lapply(contrasts, function(x) paste0(x[2], '_vs_', x[3])))
        out.l2fc <- data.frame(matrix(nrow=nrow(dds), ncol=length(contrasts.names))) ; colnames(out.l2fc) <- paste0(contrasts.names, '.l2FC')
        out.padj <- data.frame(matrix(nrow=nrow(dds), ncol=length(contrasts.names))) ; colnames(out.padj) <- paste0(contrasts.names, '.padj')
        for (i in 1:length(contrasts)) {
            res <- results(dds, contrast = contrasts[[i]])
            metadata(res)$alpha <- FDR
            out <- as.data.frame(res)
            out.l2fc[,i] <- out[,2]
            out.padj[,i] <- out[,6]
        }

    }

    # Output TPM normalized counts
    message('\nComputing expression estimates...')
    TPM <- counts_to_tpm(counts(dds), featureLength)
    colnames(TPM) <- e$SampleName
    results2 <- matrix(unlist(lapply( unique(gsub("_rep.", "", colnames(TPM))), function(x) { rowMeans(matrix(TPM[,grep(x, colnames(TPM))], nrow = nrow(TPM))) } )), nrow = nrow(TPM))
    colnames(results2) <- unique(gsub("_rep.", "", colnames(TPM)))
    write.table(results2, gsub(".txt", paste0(".counts.TPM.txt"), results.path), quote = F, col.names = T, row.names = F, sep = '\t')
    TPM <- assay(rlog(dds, blind = FALSE))
    colnames(TPM) <- e$SampleName
    results2 <- matrix(unlist(lapply( unique(gsub("_rep.", "", colnames(TPM))), function(x) { rowMeans(matrix(TPM[,grep(x, colnames(TPM))], nrow = nrow(TPM))) } )), nrow = nrow(TPM))
    colnames(results2) <- unique(gsub("_rep.", "", colnames(TPM)))
    write.table(results2, gsub(".txt", paste0(".counts.RLOG.txt"), results.path), quote = F, col.names = T, row.names = F, sep = '\t')
    TPM <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
    colnames(TPM) <- e$SampleName
    results2 <- matrix(unlist(lapply( unique(gsub("_rep.", "", colnames(TPM))), function(x) { rowMeans(matrix(TPM[,grep(x, colnames(TPM))], nrow = nrow(TPM))) } )), nrow = nrow(TPM))
    colnames(results2) <- unique(gsub("_rep.", "", colnames(TPM)))
    write.table(results2, gsub(".txt", paste0(".counts.VST.txt"), results.path), quote = F, col.names = T, row.names = F, sep = '\t')
    if (value.to.compute == 'TPM') {
        TPM <- counts_to_tpm(counts(dds), featureLength)
    } else if (value.to.compute == 'rlog') {
        TPM <- assay(rlog(dds, blind = FALSE))
    } else if (value.to.compute == 'vst') {
        TPM <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
    }
    colnames(TPM) <- paste0(value.to.compute, '_', e$SampleName)
    results <- cbind(
        #WormBaseID = names(model),
        #gene_name = elementMetadata(model)$geneName,
        ID = row.names(dds),
        TPM,
        out.l2fc,
        out.padj
    )
    
    message('\nGenerating results...')
    results2 <- matrix(unlist(lapply( unique(gsub("_rep.", "", colnames(TPM))), function(x) { rowMeans(matrix(TPM[,grep(x, colnames(TPM))], nrow = nrow(TPM))) } )), nrow = nrow(TPM))
    colnames(results2) <- unique(gsub("_rep.", "", colnames(TPM)))
    results3 <- results
    results3[is.na(results3)] <- 0
    results3$sign.l2FC <- apply(results3[,grepl("l2FC", colnames(results3))], 1, function(x) any(x < -FC | x > FC))
    results3$sign.padj <- apply(results3[,grepl("padj", colnames(results3))], 1, function(x) any(x < FDR))

    # Save results
    message('\nWriting results...')
    write.table(results3, gsub(".txt", ".full.txt", results.path), quote = F, col.names = T, row.names = F, sep = '\t')
    write.table(results2, gsub(".txt", ".counts.txt", results.path), quote = F, col.names = T, row.names = F, sep = '\t')
    write.table(results, results.path, quote = F, col.names = T, row.names = F, sep = '\t')

    return(results)

}

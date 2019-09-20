# Gets blocks information from meme output (https://rdrr.io/github/frenkiboy/MyLib/src/R/MEME.R)
parse.MEME.output <- function(path, e.thresh = 0.01, hits.thresh = 20){
    require(TFBSTools)
    connection <- file(path, open='')
    file <- scan(connection, what='list', sep='\n', blank.lines.skip = F)
    close(connection)
    
    # Get the start / end positions of all motifs
    blocks.start <- grep('BLOCKS', file) + 3
    blocks.end   <- grep('^//', file) - 1
    ranges.start <- grep('sorted by position p-value', file) + 4
    ranges.end   <- grep('block diagrams', file) - 4
    probs.start  <- grep('position-specific probability matrix', file) + 3
    probs.end    <- grep('regular expression', file) - 4
    
    if(length(blocks.start) == 0)
    return(list(blocks=NULL, header=NULL))
    
    # Get primary info
    header.ind <- grep('^MOTIF', file)
    header <- data.frame(do.call(rbind, strsplit(file[header.ind],'\\s+'))[,c(3,6,9,15)], stringsAsFactors = F)
    colnames(header) <- c('name','width','hits','Eval')
    options(scipen = 10)
    header[,c(2:4)] <- apply(header[,c(2:4)], 2, as.numeric)
    header$name <- paste0(basename(path) %>% gsub('_meme.txt', '_' ,.), header$name)
    
    # Get blocks
    l.blocks <- lapply(seq_along(header$name), function(i) {
        block <- file[blocks.start[i]:blocks.end[i]]
        structure(
            unlist(lapply(strsplit(block, split='\\s+'), '[', 4)),
            e.val = header$Eval[i],
            nsites = header$hits[i],
            class = 'block',
            pfm = matrix(),
            pwm = matrix(),
            prior = vector()
        )
    }) %>% setNames(header$name)
    
    # Get ranges
    l.ranges <- lapply(seq_along(header$name), function(i) {
        if (all(sapply(strsplit(file[ranges.start[i]:ranges.end[i]], split='\\s+'), '[', 2) %in% c('+', '-', '*', '.'))) {
            range <- cbind(
                sapply(strsplit(file[ranges.start[i]:ranges.end[i]], split='\\s+'), '[', 1) %>% gsub('\\(.*', '', .) %>% str_split(':|-', simplify = T),
                sapply(strsplit(file[ranges.start[i]:ranges.end[i]], split='\\s+'), '[', 2)
            )
            structure(
                GRanges(seqnames = range[,1], ranges = IRanges(as.numeric(range[,2]), end = as.numeric(range[,3])), strand = range[,4]),
                class = 'GRanges'
            )
        } 
        else {
            range <- cbind(
                sapply(strsplit(file[ranges.start[i]:ranges.end[i]], split='\\s+'), '[', 1) %>% gsub('\\(.*', '', .) %>% str_split(':|-', simplify = T)
            )
            structure(
                GRanges(seqnames = range[,1], ranges = IRanges(as.numeric(range[,2]), end = as.numeric(range[,3])), strand = '*'),
                class = 'GRanges'
            )
        }
    }) %>% setNames(header$name)
    
    # Get probs
    l.probs <- lapply(seq_along(header$name), function(i) {
        prob <- file[probs.start[i]:probs.end[i]] %>% str_split(' s+')
        pwm <- BlocksToPFM(l.blocks[[i]])
        toPWM(
            PFMatrix(
                profileMatrix = as.matrix(pwm),
                ID = names(l.blocks)[i], 
                name = names(l.blocks)[i], 
                strand = '*'
            ),
            type='prob'
        )
    }) %>% setNames(header$name) %>% do.call(PWMatrixList, .)
    
    # Subset significant motifs
    l.blocks <- l.blocks[header$hits > hits.thresh & header$Eval < e.thresh]
    l.ranges <- l.ranges[header$hits > hits.thresh & header$Eval < e.thresh]
    l.probs <- l.probs[header$hits > hits.thresh & header$Eval < e.thresh]
    
    # return results
    res <- list(
        'header' = header,
        'blocks' = l.blocks, 
        'ranges' = l.ranges,
        'PWMs' = l.probs
    )
    return(res)
}

# Parse blocks information for all the meme.txt files located in 'path' (https://rdrr.io/github/frenkiboy/MyLib/src/R/MEME.R)
parseMEMEblocks <- function(path, out='PWM', type='prob'){
    
    require(TFBSTools)
    
    if(length(path) == 1){
        meme.files <- list.files(path, recursive=TRUE, pattern='.txt', full.names=TRUE)
    } else{
        meme.files <- path
    }
    if(length(meme.files) == 0)
    stop('There are no imput.files')
    
    # Read all meme files
    blocks <- lapply(meme.files, MemeBlocks)
    names(blocks) <- basename(meme.files) %>% gsub('_meme.txt', '', .)
    
    # Convert to PWMs
    pwml = list()
    for (i in 1:length(blocks)) {
        bl = blocks[[i]]
        if(is.null(bl$header))
        next()
        bl$pfm = lapply(bl$blocks, BlocksToPFM)
        sampname = names(blocks)[[i]]
        pwms = lapply(names(bl$pfm), function(x) {
            toPWM(PFMatrix(profileMatrix=as.matrix(bl$pfm[[x]]),
            ID  = paste(sampname, x, sep='.'), 
            name= paste(sampname, x, sep='.'), 
            strand='*'), type='prob')
        })
        names(pwms) <- sapply(pwms, function(PWM) PWM@ID)
        pwml[[sampname]] <- pwms
    }
    pwmu = unlist(pwml)
    names(pwmu) <- sapply(pwml, names) %>% unlist
    
    res <- PWMatrixList()
    for (k in seq_along(pwmu)) {
        res[[k]] <- pwmu[[k]]
    }
    names(res) <- names(pwmu)
    
    return(res)
    
}

# takes a meme class object and returns a PFM
BlocksToPFM <- function(char){
    d = apply(do.call(rbind, strsplit(char, split='')),2, function(x)table(factor(x, levels=c('A','C','G','T'))))
    return(d)
}

# Map a list of motifs to the genome and return GRanges (and save bed file) (https://rdrr.io/github/frenkiboy/MyLib/src/R/MEME.R)
scan_Genome <- function(
    matlist,
    genome.name,
    genome = NULL,
    chrs   = NULL,
    outpath,
    min.score = '80%',
    ncores    = 40,
    remove    = FALSE
    
){
    
    require(TFBSTools)
    
    library(doMC)
    suppressPackageStartupMessages(library(GenomicAlignments))
    library(TFBSTools)
    library(stringr)
    
    if(is.null(genome))
    stop('Genome object is not assigned')
    
    if(is.null(chrs))
    chrs = names(genome)
    
    #source(file.path(lib.path, 'FileLoader.R'), local=TRUE)
    gname = str_replace(genome.name,'^.+\\.','')
    path_out_scan_genome = file.path(outpath, gname)
    dir.create(path_out_scan_genome, showWarnings=FALSE)
    
    gl = lapply(chrs, function(x)genome[[x]])
    names(gl) = chrs
    gl = DNAStringSet(gl)
    
    registerDoMC(ncores)
    donefiles = as.character(list.files(path_out_scan_genome))
    if(remove)
    donefiles = vector()
    
    foreach(m = 1:length(matlist), .errorhandling='remove')%dopar%{
        print(m)
        mat = matlist[[m]]
        name = names(matlist)[[m]]
        outname = paste(name, 'motif_search', str_replace(min.score,'%','pct'), 'bam', sep='.')
        if((length(donefiles)==0) || (!outname %in% donefiles)){
            
            print(name)
            hits.p  = searchSeq(mat, gl, strand='+', min.score=min.score)
            hits.m  = searchSeq(mat, gl, strand='-', min.score=min.score)
            ghits.p = try(unlist(GRangesList( lapply(hits.p, function(x)as(x,'GRanges')))))
            ghits.m = try(unlist(GRangesList( lapply(hits.m, function(x)as(x,'GRanges')))))
            ghits = GRanges()
            if(!is.null(ghits.p))
            ghits = c(ghits, ghits.p)
            if(!is.null(ghits.m))
            ghits = c(ghits, ghits.m)
            
            seqlengths(ghits) = seqlengths(genome)[seqlevels(ghits)]
            genome(ghits) = genome.name
            
            rtracklayer::export(simplifyGRanges(ghits), con=gsub(".bam", ".bed", file.path(path_out_scan_genome, outname)))
            
            return(ghits)
        }
    }
}

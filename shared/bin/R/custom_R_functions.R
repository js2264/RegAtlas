# ---------------------------------------------------------------------------- #
# ---------------- R ENVIRONMENT FUNCTIONS ---------------- #
# ---------------------------------------------------------------------------- #
# Function to get the size of all objects in the environment
getEnvironmentSize <- function() {
  size = 0
  sizes = c()
  df = data.frame(Name=character(0), Size=character(0), Class=character(0), stringsAsFactors=F)
  for (x in ls(envir=.GlobalEnv)) {
    thisSize = object.size(get(x))
    size = size + thisSize
    sizes[x] = thisSize
    df[x,] = c(x, utils:::format.object_size(thisSize, 'auto'), as.character(class(x))[1])
  }
  NROW=nrow(df)
  vec=sizes[rev(order(sizes))]
  vec=vec[vec > 100000]
  df=df[names(vec),]
  row.names(df)=c()
  print(df)
  message(paste0("There are   ", NROW, "   objects in the environment."))
  message("The total workspace size is   ",appendLF = F); print(size, units='auto')
}
# Improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
# Make a .bak object of an existing object, i.e. `XX.bak <- XX`
bak <- function(X) { 
    assign(paste0(deparse(substitute(X)), ".bak"), envir = parent.frame(), X)
    warning("Use this function with care")
    message(paste0(deparse(substitute(X)), " object copied to ", deparse(substitute(X)), ".bak"))
}
# Source files in bin/ folders
.sourceDirs <- function(folders = c('~/shared/bin/R/', '~/Documents/PhD/__Bioinfo/shared/bin/R/')) {
    for (folder in folders) {
        if (any(grepl('*.R', list.files(folder)))) {
            files <- paste0(folder, '/', list.files(folder, pattern = '.R'))
            for (file in files) {
                source(file)
            }
        }
    }
}
s <- base::summary
h <- utils::head
ht <- function(d, n = 2) rbind(head(d,n),tail(d,n))
hh <- function(d) if(class(d)=="matrix"|class(d)=="data.frame") d[1:5,1:5]
read.cb <- function(...) {
  ismac <- Sys.info()[1]=="Darwin"
  if (!ismac) read.table(file="clipboard", ...)
  else read.table(pipe("pbpaste"), ...)
}
sort.per.clusters <- function(df, clusters) {
    df2 <- df[0,]
    for (cluster in unique(clusters)) {
        tmp <- df[clusters == cluster, ]
        which.max.col <- which.max(colSums(tmp))
        tmp <- tmp[order(tmp[,which.max.col], decreasing = T),]
        df2 <- rbind(df2, tmp)
    }
    colnames(df2) <- colnames(df)
    return(df2)
}
# custom print.data.frame 
print.data.frame <- function(x) {print(knitr::kable(x))}

# ---------------------------------------------------------------------------- #
# ---------------- DATA FORMATTING FUNCTIONS
# ---------------------------------------------------------------------------- #
# Transform a named list into a long-format data.frame
fromListToLongDataframe <- function(list, names = c("name", "value")) {
    require(tidyverse)
    list %<>% 
        enframe() %>%
        unnest() %>% 
        setColNames(names)
    return(list)
}
# Function to remove NAs from a vector
na.remove <- function(x) {
    which.na <- is.na(x)
    return(x[!which.na])
}
# Function to replace NAs from a vector by a given value
na.replace <- function(x, value) {
    which.na <- is.na(x)
    x[which.na] <- value
    return(x)
}
# Function to remove a specific value from a vector
value.remove <- function(x, val) {
    which.val <- x %in% val
    return(x[!which.val])
}
# Function to replace a specific value from a vector
value.replace <- function(x, val, newval) {
    which.val <- x == val
    x[which.val] <- newval
    return(x)
}
vec.split <- function(x, split) {
    vec <- unlist(strsplit(x, split))
    return(vec)
}
# Shift vector by n elements (add extra zeros)
shiftVector <- function(vec, shift, replacement = 0) {
    if (shift >= 0) {
        first <- 1 + shift
        c(
            vec[first:length(vec)],
            rep(replacement, shift) 
        )
    } 
    else {
        last <- length(vec) + shift
        c(
            rep(replacement, abs(shift)), 
            vec[1:last]
        )
    }
}
# Function to increment 1
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
# Make counts table from 2 factors
makeCounts <- function(vec1, vec2, lev.vec1 = NULL, lev.vec2 = NULL) {
    lev.vec1 <- if (length(lev.vec1) == 0) { if (is.factor(vec1)) {levels(vec1)} else { unique(vec1) } } else { lev.vec1 }
    lev.vec2 <- if (length(lev.vec2) == 0) { if (is.factor(vec2)) {levels(vec2)} else { unique(vec2) } } else { lev.vec2 }

    tab <- as.data.frame(matrix(NA, ncol = length(lev.vec1), nrow = length(lev.vec2))) ; colnames(tab) <- lev.vec1 ; row.names(tab) <- lev.vec2
    for (VAL.1 in lev.vec1) {
        tab[,VAL.1] <- table(factor(vec2[vec1 == VAL.1], levels = lev.vec2))
    }

    return(tab)
}
# Function to compute pct in a vector
pct <- function(vec) {
    t <- table(vec)
    return( round(t / sum(t) * 100, 3) )
}
pct.scale <- function(vec) {
    return(vec / sum(vec) * 100)
}
pct.cols <- function(df) {
    df2 <- apply(df, 2, function(COL) {COL/sum(COL)*100})
    return(df2)
}
pct.rows <- function(df) {
    df2 <- apply(df, 1, function(ROW) {ROW/sum(ROW)*100})
    return(df2)
}
addPct <- function(vec) {
    pct <- round(vec / sum(vec) * 100, 1)
    char <- paste0(format(vec, big.mark = ",", scientific = FALSE), ' (', pct, '%)')
    return(char)
}
## Shortcut to lapply(LIST, FUN)
heads <- function(x) {
    lapply(x, head)
}
maxs <- function(x) {
    max(unlist(lapply(x, max)))
}
mins <- function(x) {
    min(unlist(lapply(x, min)))
}
medians <- function(x) {
    sapply(x, median)
}
dims <- function(x) {
    min(unlist(lapply(x, dim)))
}
# which.median, inspired from which.min
which.median <- function(x) {
    if (length(x) %% 2 != 0) {
        which(x == median(x))
    } else if (length(x) %% 2 == 0) {
        a = sort(x)[c(length(x)/2, length(x)/2+1)]
        c(which(x == a[1]), which(x == a[2]))
    }
}
# Bin a vector in N BINS
bin.vector <- function(vec, nbins) {
    res <- as.numeric(cut(vec, breaks = seq(min(vec), max(vec), length.out = nbins + 1), include.lowest = T))
    return(res)
}
# %notin% function
`%notin%` <- Negate(`%in%`)
# change pdf() defaults
pdf10 <- function(...) {
    pdf(..., width = 10, height = 10)
}
pdftmp <- function(width = 10, height = 10) {
    pdf('tmp.pdf', width = width, heigh = height)
}
dev10 <- function(width = 10, height = 10) {
    dev.new(width = width, height = height)
}
# %in.which% function : for a vector of names, output in which element of a given list the name is
`%inwhich%` <- function(names, list) {
    t <- reshape2::melt(list)
    l <- t[match(names, t$value), ]$L1
    return(l)
}
# Function to compare 2 lists (2-way venn without plotting)
compare2lists <- function(vec1, vec2, name.vec1 = NULL, name.vec2 = NULL, do.plot = F, return.pct = T) {
    if(!do.plot) {
        all <- length(unique(c(vec1, vec2)))
        a <- list(elements_in_a = length(vec1), elements_shared_with_b = paste0(sum(vec1 %in% vec2), ' (', round(sum(vec1 %in% vec2)/length(vec1)*100, 2), ' %)'), elements_only_in_a = paste0(sum(vec1 %notin% vec2), ' (', round(sum(vec1 %notin% vec2)/length(vec1)*100, 2), ' %)'))
        b <- list(elements_in_b = length(vec2), elements_shared_with_a = paste0(sum(vec2 %in% vec1), ' (', round(sum(vec2 %in% vec1)/length(vec2)*100, 2), ' %)'), elements_only_in_b = paste0(sum(vec2 %notin% vec1), ' (', round(sum(vec2 %notin% vec1)/length(vec2)*100, 2), ' %)'))
        if (return.pct) return(list(all = all, a = a, b = b))
    } else {
        all <- length(unique(c(vec1, vec2)))
        a <- list(elements_in_a = length(vec1), elements_shared_with_b = sum(vec1 %in% vec2), elements_only_in_a = sum(vec1 %notin% vec2))
        b <- list(elements_in_b = length(vec2), elements_shared_with_a = sum(vec2 %in% vec1), elements_only_in_b = sum(vec2 %notin% vec1))
        if (return.pct) return(list(all = all, a = a, b = b))
    }
    return(list("only_a" = vec1[vec1 %notin% vec2], "only_b" = vec2[vec2 %notin% vec1], "shared" = vec1[vec1 %in% vec2]))
}
# Function to split a vector in 5
vsplit <- function(v, n) {
    l = length(v)
    r = l/n
    return(lapply(1:n, function(i) {
        s = max(1, round(r*(i-1))+1)
        e = min(l, round(r*i))
        return(v[s:e])
    }))
}
# Function to clamp a data.frame
clamp <- function(df, lims = c(-99999, 99999)) {
    df2 <- as.data.frame(df)
    df2[df2 < lims[1]] <- lims[1]
    df2[df2 > lims[2]] <- lims[2]
    return(df2)
}
# Check outliers
is_outlier <- function(x) {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
# Set row and col names
setRowNames <- function(df, names) {
    row.names(df) <- names
    return(df)
}
setColNames <- function(df, names) {
    colnames(df) <- names
    return(df)
}
# Switch from counts to cases formats (rectangular matrix of counts -> long df for ggplots2)
countsToCases <- function(x, countcol = "value") {
    x2 <- reshape2::melt(x)
    x2$row <- factor(rep(row.names(x), ncol(x)), levels = row.names(x))
    x2 <- x2[, c('variable', 'row', 'value')]
    return(x2)
}
countsToLongFormat <- countsToCases
dataframeToLongFormat <- function(df) {
    df2 <- data.frame(col.value = character(), row.value = character(), value = numeric(), stringsAsFactors = FALSE)
    i <- 0
    for (COL in colnames(df)) {
        for (ROW in row.names(df)) {
            i <- i + 1
            df2[i, 'col.value'] <- COL
            df2[i, 'row.value'] <- ROW
            df2[i, 'value'] <- df[ROW, COL]
        }
    }
    df2 <- mutate_if(df2, is.character,as.factor)
    return(df2)
}
namedListToLongFormat <- function(x) {
    lapply(names(x), function(NAME) {
        L <- x[[NAME]]
        if (is.null(ncol(L))) {
            data.frame(value = L, name = rep(NAME, length(L)))
        } 
        else {
            data.frame(L, name = rep(NAME, nrow(L)))
        }
    }) %>% do.call(rbind, .)
}

# ---------------------------------------------------------------------------- #
# ---------------- BED FORMAT-RELATED FUNCTIONS
# ---------------------------------------------------------------------------- #
# Call bedtools intersect from R
bedIntersectR <- function(functionstring = "bedtools intersect", bed1, bed2, opt.string = "") {
  #create temp files
  a.file = tempfile()
  b.file = tempfile()
  out = tempfile()
  options(scipen = 99) # not to use scientific notation when writing out

  #write bed formatted dataframes to tempfile
  write.table(bed1, file = a.file, quote = F, sep = "\t", col.names = F, row.names = F)
  write.table(bed2, file = b.file, quote = F, sep = "\t", col.names = F, row.names = F)

  # create the command string and call the command using system()
  command = paste(functionstring, opt.string, "-a", a.file, "-b", b.file, ">", out, sep=" ")
  cat(command,"\n")
  try(system(command))

  res = read.table(out, header = F)
  unlink(a.file)
  unlink(b.file)
  unlink(out)

  return(res)
}
# Write bed file from dataframe
writeBed <- function(df, file, col = 1:3) {
  
  # Check if data is a data table
  if (!is.data.frame(df[,col])){stop("Data is not formatted as a data.frame. Stopping function.")}
  
  # Write bed file
  write.table(df[,col], file = file, quote = F, sep = "\t", col.names = F, row.names = F)

}
# Merge overlapping peaks of 2 grange objetcs (or two files with columns 1:3 being chrom-start:stop)
mergeOverlappingPeaks <- function(granges1, granges2 = NULL) {
    require(GenomicRanges)
    if(is.null(granges2)) {
        if(is.character(granges1)) {
            granges1 <- read.table(granges1, stringsAsFactors = F)[,1:3]
            colnames(granges1) <- c('chrom', 'start', 'stop')
            granges1$start <- granges1$start+1
            granges1 <- makeGRangesFromDataFrame(granges1)
        }
        message('Only 1 file provided. Converting into GRanges..')
        return(granges1)
    }
    if(is.character(granges1) & is.character(granges2)) {
        granges1 <- read.table(granges1, stringsAsFactors = F)[,1:3]
        granges2 <- read.table(granges2, stringsAsFactors = F)[,1:3]
        colnames(granges1) <- c('chrom', 'start', 'stop')
        colnames(granges2) <- c('chrom', 'start', 'stop')
        granges1$start <- granges1$start+1
        granges2$start <- granges2$start+1
        granges1 <- makeGRangesFromDataFrame(granges1)
        granges2 <- makeGRangesFromDataFrame(granges2)
    }
    grange <- reduce(c(granges1[granges1 %over% granges2], granges2[granges2 %over% granges1]))
    return(grange)
}
# Get jaccard score from 2 sets of peaks
get.jaccard.from.peaks <- function(peaks.A, peaks.B, granges) {
    granges$has.A <- granges %over% peaks.A
    granges$has.B <- granges %over% peaks.B
    p <- sum(granges$has.A & granges$has.B)
    q <- sum(granges$has.A & !granges$has.B)
    r <- sum(!granges$has.A & granges$has.B)
    jaccard <- (p / (p+q+r))
    return(jaccard)
}

# ---------------------------------------------------------------------------- #
# ---------------- FUNCTIONS TO GET SOME GENE INFOS
# ---------------------------------------------------------------------------- #
WB2name <- function(x) {
    IDs <- convID$gene_name[match(x, convID$WormBaseID)]
    return(IDs)
}
name2WB <- function(x) {
    IDs <- convID$WormBaseID[match(x, convID$gene_name)]
    return(IDs)
}
getWB <- function(GENE) {
    system(paste0("open https://wormbase.org/species/c_elegans/gene/", getGeneInfos(GENE, save=F, verbose = F)$Gene.info[1]))
}
# Fecth gene infos from WormBase
fetchWBinfos <- function(GENE) {
    
    if (grepl('WBGene', GENE)) { query <- WB2name(GENE) } else { query <- GENE }
    req <- httr::GET(paste("http://www.wormbase.org/species/c_elegans/gene/", query, sep = ""), config = httr::content_type_json())
    descr <- content(req)$results[[1]][['description']]
    return(descr)
    
}
# Get expression values from cao supp. table 03
getCao03 <- function(vec) {
  if (sum(ls() == "cao03") == 0) {load(".test-plotting.RData")}
  vec=gsub(" ", "", vec)
  if (any(grepl(",", vec))) {vec=unlist(strsplit(vec, ","))}
  mat=cao03[name2WB(vec),]
  row.names(mat)=WB2name(row.names(mat))
  print(mat)
}
getCao06 <- function(vec) {
  if (sum(ls() == "cao06.2") == 0) {load(".test-plotting.RData")}
  vec=gsub(" ", "", vec)
  if (any(grepl(",", vec))) {vec=unlist(strsplit(vec, ","))}
  mat=cao06.2[cao06.2$symbol %in% vec,]
  print(mat)
}
# Function to get a lot of informations on an individual gene
getGeneInfos <- function(GENES, verbose = T, saveTXT = F, exportResult = F) {
    suppressWarnings(suppressMessages(library(GenomicRanges)))

    if (any(unlist(lapply(list(all, genes.gtf, cao03, LCAP, max.tissue.df.LCAP, order.tissues, SUMM), is.null)))) {stop('Some objects are missing. Aborting.')}

    for (GENE in GENES) {
        # Convert gene in WormBaseID and get locus name
        if (!grepl('WBGene', GENE)) { WBID <- name2WB(GENE) } else { WBID <- GENE }
        locusID <- WB2name(WBID)
        biotype <- genes.gtf[WBID]$gene_biotype
        # Get gene expression infos from cao and tissue.spe LCAP
        CAO.GENE <- cao03[WBID,]
        LCAP.GENE <- LCAP_normalized[WBID,]
        if (nrow(CAO.GENE) == 0) {stop('Gene not found in Cao data. Aborting.')}
        if (nrow(LCAP.GENE) == 0) {stop('Gene not found in LCAP data. Aborting.')}
        l2mean.CAO.GENE <- log2((CAO.GENE+1)/rowMeans(CAO.GENE+1))
        l2mean.LCAP.GENE <- log2((LCAP.GENE+1)/rowMeans(LCAP.GENE+1))
        # Get gene expression dynamics from developmental datasets
        LCAPdev.GENE <- LCAPdev[WBID,]
        l2mean.LCAPdev.GENE <- log2((LCAPdev.GENE+1)/rowMeans(LCAPdev.GENE+1))
        # Get FC over other tissues
        max.tissue.df.LCAP.GENE <- max.tissue.df.LCAP[WBID,]
        tissue <- as.vector(max.tissue.df.LCAP.GENE[,'which.tissues'])
        tau <- round(unname(genes.gtf[WBID]$TauScore), 2)
        FCs <- round(sapply(order.tissues[1:5], function(TISSUE) {LCAP.GENE[, TISSUE] / mean(unlist(LCAP.GENE[, colnames(LCAP.GENE) != TISSUE]))}), 3)
        # Get REs associated with gene (which ones, associated tissues and dev. dynamics)
        which <- grep(WBID, all$WormBaseID)
        REs <- all[which,]
        if (nrow(REs) == 0) {
            message('!!! WARNING !!! No associated REs are found !!!')
            REs.coords <- NULL ; REs.tau <- 0 ; REs.tissues <- NULL ; REs.TFs <- NULL ; REs.TFs.nb <- NULL ; ATAC.GENE <- 0 ; ATACdev.GENE <- 0 ; l2mean.ATAC.GENE <- 0 ; l2mean.ATACdev.GENE <- 0 
        } else {
            REs.coords <- REs[,c(1:3, 11)]
            row.names(REs.coords) <- paste0('RE_', seq(1:nrow(REs.coords)))
            REs.tau <- REs$TauScore
            REs.tissues <- order.tissues[REs$clustersATAC.tissues]
            REs.TFs <- SUMM[which,]
            REs.TFs.nb <- REs.TFs$TFnb
            ATAC.GENE <- ATAC[which,]
            #ATACdev.GENE <- ATACdev[row.names(ATACdev) %in% all$coords[which],]
            l2mean.ATAC.GENE <- round(log2((ATAC.GENE+1)/rowMeans(ATAC.GENE+1)),2)
            row.names(l2mean.ATAC.GENE) <- row.names(REs.coords)
            colnames(l2mean.ATAC.GENE) <- colnames(l2mean.LCAP.GENE)
            #l2mean.ATACdev.GENE <- round(log2((ATACdev.GENE+1)/rowMeans(ATACdev.GENE+1)),2)
            #row.names(l2mean.ATACdev.GENE) <- row.names(REs.coords)
            #colnames(l2mean.ATACdev.GENE) <- colnames(l2mean.LCAPdev.GENE)
        }

        # Build results object
        res <- list()
        res[['Gene.info']] <- c(WormBaseID = WBID, locus = locusID, Biotype = biotype, Tau.score = tau, Enriched.tissue = tissue)
        res[['Gene.expr.TPM']] <- round(rbind(Cao = CAO.GENE, tiss.spe.LCAP = LCAP.GENE),2)
        res[['Gene.expr.l2mean']] <- round(rbind(Cao = l2mean.CAO.GENE, tiss.spe.LCAP = l2mean.LCAP.GENE),2)
        res[['Gene.expr.dev.TPM']] <- round(rbind(LCAPdev = LCAPdev.GENE),2)
        res[['Gene.expr.dev.l2mean']] <- round(rbind(LCAPdev = l2mean.LCAPdev.GENE),2)
        res[['Associated.REs']] <- cbind(REs.coords, RE.Tau.score = round(REs.tau,2), REs.associated.tissues = REs.tissues)
        #res[['Associated.REs.dev']] <- round(ATACdev.GENE,2)
        #res[['Associated.REs.dev.l2mean']] <- l2mean.ATACdev.GENE
        res[['Associated.REs.tissue']] <- round(ATAC.GENE,2)
        res[['Associated.REs.tissue.l2mean']] <- l2mean.ATAC.GENE
        res[['Associated.REs.TFs.infos']] <- cbind(REs.coords, RE.Tau.score = round(REs.tau,2), RE.associated.tissues = REs.tissues)

        # Return results
        if (verbose == T) {
            width.bak <- options()$width
            options("width" = 200)
            cat('\n', rep('=', times = 52+nchar(locusID)), '\n', rep('>', 14),'   ', WBID, ' --- ', locusID, '  ', rep('<', 14),'\n', rep('=', times = 52+nchar(locusID)), sep = "", fill = T)

            cat('\nGene infos\n', rep('-', times = 52+nchar(locusID)), sep = "", fill = T)
                cat(paste0(capture.output(res[['Gene.info']]), collapse = "\n"), sep = "", fill = T)
            cat('\n\nGene expression\n', rep('-', times = 52+nchar(locusID)), sep = "", fill = T)
                cat('\n>> Gene tissue-specific enrichments (FC vs. other tissues)', sep = "", fill = T)
                    cat(paste0(capture.output(FCs), collapse = "\n"), sep = "", fill = T)
                cat('\n>> Gene TPM (tissue-specific)', sep = "", fill = T)
                    cat(paste0(capture.output(res[['Gene.expr.TPM']]), collapse = "\n"), sep = "", fill = T)
                cat('\n>> Gene log2 mean-centered TPM (tissue-specific)', sep = "", fill = T)
                    cat(paste0(capture.output(res[['Gene.expr.l2mean']]), collapse = "\n"), sep = "", fill = T)
                cat('\n>> Gene TPM (dev.)', sep = "", fill = T)
                    cat(paste0(capture.output(res[['Gene.expr.dev.TPM']]), collapse = "\n"), sep = "", fill = T)
                cat('\n>> Gene log2 mean-centered TPM (dev.)', sep = "", fill = T)
                    cat(paste0(capture.output(res[['Gene.expr.dev.l2mean']]), collapse = "\n"), sep = "", fill = T)
            cat('\n\nAssociated regulatory elements\n', rep('-', times = 52+nchar(locusID)), sep = "", fill = T)
                cat('\n>> REs coordinates', sep = "", fill = T)
                    cat(paste0(capture.output(res[['Associated.REs']]), collapse = "\n"), sep = "", fill = T)
                cat('\n>> REs log2 mean-centered accessibility (tissue-specific)', sep = "", fill = T)
                    cat(paste0(capture.output(res[['Associated.REs.tissue']]), collapse = "\n"), sep = "", fill = T)
                #cat('\n>> REs log2 mean-centered accessibility (dev)', sep = "", fill = T)
                #    cat(paste0(capture.output(res[['Associated.REs.dev']]), collapse = "\n"), sep = "", fill = T)

            cat('\n', rep('=', times = 52+nchar(locusID)), '\n', sep = "", fill = T)
            options("width" = width.bak)
        }

        if (saveTXT != F) {
            width.bak <- options()$width
            options("width" = 200)

            txt.file <- ifelse(is.logical(saveTXT), file(paste0(locusID, '-summary.txt'), open="wt"), as.character(saveTXT))
            cat('\n', rep('=', times = 52+nchar(locusID)), '\n', rep('>', 14),'   ', WBID, ' --- ', locusID, '   ', rep('<', 14),'\n', rep('=', times = 52+nchar(locusID)), sep = "", fill = T, file = txt.file)

            cat('\nGene infos\n', rep('-', times = 52+nchar(locusID)), sep = "", fill = T, file = txt.file, append = T)
                cat(paste0(capture.output(res[['Gene.info']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
            cat('\n\nGene expression\n', rep('-', times = 52+nchar(locusID)), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> Gene tissue-specific enrichments (FC vs. other tissues)', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(FCs), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> Gene TPM (tissue-specific)', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(res[['Gene.expr.TPM']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> Gene log2 mean-centered TPM (tissue-specific)', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(res[['Gene.expr.l2mean']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> Gene TPM (dev.)', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(res[['Gene.expr.dev.TPM']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> Gene log2 mean-centered TPM (dev.)', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(res[['Gene.expr.dev.l2mean']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
            cat('\n\nAssociated regulatory elements\n', rep('-', times = 52+nchar(locusID)), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> REs coordinates', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(res[['Associated.REs']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> REs log2 mean-centered accessibility (tissue-specific)', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(res[['Associated.REs.tissue']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
                #cat('\n>> REs log2 mean-centered accessibility (dev)', sep = "", fill = T, file = txt.file, append = T)
                #    cat(paste0(capture.output(res[['Associated.REs.dev']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)

            cat('\n', rep('=', times = 52+nchar(locusID)), '\n', sep = "", fill = T, file = txt.file, append = T)

            options("width" = width.bak)

        }

    }

    if (length(GENES) == 1 & exportResult == T) {
        return(res)
    }
}
ggi <- getGeneInfos

# ---------------------------------------------------------------------------- #
# ---------------- FUNCTIONS TO PLOT SOME INFOS LINKED TO GENES
# ---------------------------------------------------------------------------- #
# To get a heatmap of tissue expression, for a vector of genes
getHeatmapOfCaoExpr <- function(vec, name.main = "", which.column = NULL) {
    require(RColorBrewer)
    require(gplots)
    #if (sum(ls() == "cao03") == 0) {load(".test-plotting.RData")}
    cao03b <- cao03[,order.tissues[1:5]]
    if (all(grepl("WB", vec))) {WBID <- vec} else {WBID <- name2WB(vec)}
    WBID <- WBID[WBID %in% row.names(cao03)]
    CAO.mat <- log2( as.matrix(cao03b[WBID,] +1) )
    CAO.mat[is.na(CAO.mat)] <- 0
    if (!is.null(which.column)) { CAO.mat <- CAO.mat[order(CAO.mat[,which.column]),] }

    if (all(grepl("WB", vec))) {row.names(CAO.mat)=WB2name(WBID)} else {row.names(CAO.mat)=vec}
    #bins <- c(0, 1, quantile(c(na.omit(c(CAO.mat)[c(CAO.mat) >= 1])), probs=(1:100)/100)[-which.max(quantile(c(na.omit(c(CAO.mat)[c(CAO.mat) >= 1])), probs=(1:100)/100))])
    bins <- seq(0, max(CAO.mat), length.out = 100)
    cols <- colorRampPalette(c("white", "orange", "darkred"))(length(bins)-1)

    heatmap_parmfrow(CAO.mat,
          margins = c(9,4,8,6),
          col = cols,
          breaks = bins,
          main = name.main,
          key.lab = "log2(TPM)",
          cexCol = 0.4,
          cexRow = 1.2,
          ylab = "",
          scaleUP = 1.2,
          nticks = 15,
          doClust = F,
          plotSep = F
    )
}
getVioplotxOfCaoExpr <- function(vec, name.main = "") {
    require(RColorBrewer)
    require(vioplotx)
    #if (sum(ls() == "cao03") == 0) {load(".test-plotting.RData")}
    cao03b <- cao03[,order.tissues[1:5]]
    if (all(grepl("WB", vec))) {WBID <- vec} else {WBID <- name2WB(vec)}
    WBID <- WBID[WBID %in% row.names(cao03)]
    CAO.mat <- log2( as.matrix(cao03b[WBID,] +1) )
    CAO.mat[is.na(CAO.mat)] <- 0
    if (all(grepl("WB", vec))) {row.names(CAO.mat)=WB2name(WBID)} else {row.names(CAO.mat)=vec}
    vioplotx(CAO.mat[,1], CAO.mat[,2], CAO.mat[,3], CAO.mat[,4], CAO.mat[,5], col = color, main = name.main, ylab = "log2(TPM)")
}

# ---------------------------------------------------------------------------- #
# ---------------- PLOTTING FUNCTIONS
# ---------------------------------------------------------------------------- #
# Own heatmap function (doesn't work perfectly...)
heatmap_parmfrow <- function (x, margins = c(5, 3, 6.5, 5), outmargins = c(2, 2, 2, 2),
    labRow = NULL, labCol = NULL,
    cexRow = 0.2 + 1/log10(nrow(x)), cexCol = 0.2 + 1/log10(ncol(x)), cex.main = 1.2,
    col = heat.colors, breaks = NULL,
    main = NULL, key.lab = NULL, ylab = NULL, xlab = NULL, scaleUP = 1.2, scaleUP.bis = NULL,
    nticks = 10, doClust = F, plotSep = F, colSep = '#00000064', plotKey = T, printCell = F, CellCex = 1) {

    hm <- list(
        x = x,
        margins = margins,
        labRow = labRow,
        labCol = labCol,
        cexRow = cexRow,
        cexCol = cexCol,
        cex.main = cex.main,
        col = col,
        breaks = breaks,
        main = main,
        key.lab = key.lab,
        ylab = ylab,
        scaleUP = scaleUP,
        scaleUP.bis = scaleUP.bis,
        nticks = nticks,
        doClust = doClust,
        plotSep = plotSep,
        colSep = colSep
    )

    # If labels are null, set them to col/row names
    if (is.null(labRow)) { labRow <- row.names(x) } else if (labRow[1] == 'none') { labRow <- rep("", length(row.names(x))) }
    if (is.null(labCol)) { labCol <- colnames(x) } else if (labCol[1] == 'none') { labCol <- rep("", length(colnames(x))) }

    # If x is not a matrix, stop function
    if (class(x) != 'matrix') {
        stop('...Input needs to be a matrix...')
    }

    # scaleUP.bis is the factor to increase the heatmap height, to draw the TOP horizontal line of the color key above the actual heatmap
    # scaleUP is the factor to increase the heatmap height, to draw the BOTTOM horizontal line of the color key above the actual heatmap
    if (is.null(scaleUP.bis)) {scaleUP.bis <- scaleUP - 0.1}

    # Required number of rows/columns (before transposition!!)
    nr <- nrow(x)
    nc <- ncol(x)

    # If the breaks are not defined: do 16 breaks (or the number of colors -1 if the colors are given)
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        stop('...Breaks missing...')
    }

    # Clamp values of x to the breaks
    x[x < min(breaks)] <- min(breaks)
    x[x > max(breaks)] <- max(breaks)

    # If the color argument is a function (typically 'colorRampPalette'), compute the actual colors based on how many breaks there are
    if (class(col) == "function") {
        col <- col(length(breaks) - 1)
    } else if (length(col) != (length(breaks) - 1)) {
        stop('...length(colors) needs to be length(breaks)-1...')
    }

    # Define the key labels and where they should be positionned
    labels <- unique(round(seq(min(breaks), max(breaks), length.out=nticks)))
    labels.at <- seq(0, nc, length.out = length(labels))

    # If clustering is desired, do it
    if (doClust) {
        x <- x[hclust(dist(x, method = 'max'))$order,]
    }

    # Save values
    values <- as.character(hm$x)
    values[values == 0] = ""
    if (!is.logical(printCell)) {
        values <- as.character(printCell)
        values[values == 0] = ""
    }

    # Transpose the matrix (because image() is used to plot the matrix) and add some extra empty columns (which are going to be top rows when plotted) for plotting space for the color key
    x.bak <- x
    x <- t(x)
    x <- cbind(x, matrix(rep(NA, nc * (round(nr * scaleUP)-(nr))), nrow = nc))

    # Plot heatmap with empty space on top
    par(mar = margins)
    image(1L:nc, 1L:(round(nr*scaleUP)), x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, round(nr*scaleUP)),
        axes = FALSE, xlab = "", ylab = "",
        col = col, breaks = breaks)
    axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
    axis(4, 1:nr, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow, srt=180)

    # Plot frame
    polygon(0.5 + c(0,0,nc,nc), 0.5 + c(0,nr,nr,0), col=NA, border="black")

    # Plot cells separators (if wanted)
    if (plotSep) {
        for (i in seq(0, nc)) { polygon(0.5 + c(i,i,nc,nc), 0.5 + c(0,nr,nr,0), col=NA, border=colSep) }
        for (i in seq(0, nr)) { polygon(0.5 + c(0,0,nc,nc), 0.5 + c(0,i,i,0), col=NA, border=colSep) }
    }

    # Plot legend
    if (plotKey) {
        poly <- vector(mode="list", length(col))
        val <- seq(0, nc, length.out = length(breaks))
        for(i in seq(poly)){ poly[[i]] <- 0.5 + c(val[i], val[i], val[i+1], val[i+1]) } ## Define coordinates of each rectangle of the key
        for(i in seq(poly)){ polygon(poly[[i]], c(round(nr*scaleUP), round(nr*scaleUP.bis), round(nr*scaleUP.bis), round(nr*scaleUP)), col = col[i], border = NA) } ## Plot each rectangle of the color key
        polygon(0.5 + c(0,0,nc,nc), c(round(nr*scaleUP.bis), round(nr*scaleUP),round(nr*scaleUP), round(nr*scaleUP.bis)), col = NA, border = "black")
        axis(3,at = labels.at+0.5, labels = labels, las = 2, cex.axis = 0.5, pos = round(nr*scaleUP))
        #axis(3,at=seq(0, nc, length.out=nticks), labels = labels.at, las=2, cex.axis=0.5, pos=round(nr*scaleUP))
    }

    # Plot labels
    if (!is.null(ylab)) { mtext(ylab, side = 4, line = margins[4]-2, srt = 180) }
    if (!is.null(xlab)) { mtext(xlab, side = 1, line = margins[1]-2) }
    if (!is.null(key.lab)) { mtext(key.lab, side = 3, line = 2, cex=0.8) }
    mtext(main, side = 3, line = 4, cex = cex.main)

    # Plot cell values
    if (is.logical(printCell) & printCell == T) { text(x = rep(1:nrow(x.bak), each = ncol(x.bak)), y = 1:ncol(x.bak), values, cex = CellCex) }
    if (!is.logical(printCell)) {text(x = rep(1:nrow(x.bak), each = ncol(x.bak)), y = 1:ncol(x.bak), values, cex = CellCex)}

    return(hm)

}
# Re-defining boxplot function
boxplot.2 <- function(..., whisklty = 'dotted', whiskcol = 'grey30', whiskldy = 2, staplecol = 'white', outline = F, notch = T, frame = F, medlty = 0, medpch = 21) {
    boxplot(..., whisklty = whisklty, whiskcol = whiskcol, whiskldy = whiskldy, staplecol = staplecol, outline = outline, notch = notch, frame = frame, medlty = medlty, medpch = medpch)
}
# Function to plot scale of heatmaps plotted using image() function ; the scale uses the breaks and col values used for image()
image.scale <- function(breaks, col = heat.colors(12), z= NULL, zlim = NULL, axis.pos=1, add.axis=TRUE, ...) {
     if(!missing(breaks)){
      if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
     }
     if(missing(breaks) & !is.null(zlim)){
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
     }
     if(missing(breaks) & is.null(zlim)){
      zlim <- range(z, na.rm=TRUE)
      zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
      zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
     }
     poly <- vector(mode="list", length(col))
     for(i in seq(poly)){
      poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
     }
     if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
     if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
     plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)
     for(i in seq(poly)){
      if(axis.pos %in% c(1,3)){
       polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
      }
      if(axis.pos %in% c(2,4)){
       polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
      }
     }
     box()
     if(add.axis) {axis(axis.pos)}
}
# On top of a plot already existing (typically boxplot(x)), add median lines and conf. interval
plotMedSE <- function(x, functionToUse = 'mean', colorCI = rgb(0.8, 0.8, 0.8, alpha = 0.3)) {
    center <- apply(x, 2, function(y) {do.call(functionToUse, list(y))})
    sd <- apply(x, 2, sd)
    n <- nrow(x)
    error <- qnorm(0.975)*sd/sqrt(n)
    top <- center + error
    bottom <- center - error
    x.coords <- 1:ncol(x)

    # Plot mean
    lines(center)

    # Plot confidence interval
    polygon(c(rev(x.coords), x.coords), c(rev(top), bottom), col = colorCI, border = NA)
    lines(x.coords, top, lty = 'dashed', col = 'red')
    lines(x.coords, bottom, lty = 'dashed', col = 'red')

}
## Plot overlay of densities for a list of vectors
plotCumFreqs <- function(list, colors = RColorBrewer::brewer.pal(name = 'Set2', n = 8), xlim = NULL, ylim = NULL, sort = T, pdf = NULL, add0 = T, translate = F, ...) {
        
    ylim <- if(is.null(ylim)) {c(0, 1)} else {ylim}
    xlim <- if(is.null(xlim)) {c(0, max(unlist(lapply(list, max))))} else {xlim}
    if (translate) {
        bak(ylim)
        bak(xlim)
        xlim <- ylim.bak
        ylim <- xlim.bak
    }
    if(!is.null(pdf)) pdf(pdf)
    
    plot(
        NA, 
        ylim = ylim, 
        xlim = xlim,
        bty = 'l',
        ...
    )
    if (sort) {
        list <- lapply(list, sort)
    }
    abline(h = c(0, 1), lty = 2, col = 'grey90')
    if (sort) {
        lapply(1:length(list), function(K) {
            x <- if (add0) {c(0, list[[K]])} else {list[[K]]}
            y <- if (add0) {c(0, ecdf(list[[K]])(list[[K]]))} else {ecdf(list[[K]])(list[[K]])}
            if (translate) {
                bak(x)
                bak(y)
                x <- y.bak
                y <- x.bak
            }
            par(new = T)
            plot(
                x = x, 
                y = y, 
                ylim = ylim,
                xlim = xlim,
                ylab = "", 
                xlab = "", 
                main = "", 
                col = colors[K],
                cex = 0.5, 
                bty = 'n', 
                type = 'l', 
                lwd = 2, 
                axes = F
            )
        })
    } else {
        lapply(1:length(list), function(K) {
            par(new = T)
            x <- if (add0) {c(0, 1:length(list[[K]]))} else {1:length(list[[K]])}
            y <- if (add0) {c(0, cumsum(list[[K]]))} else {cumsum(list[[K]])}
            if (translate) {
                bak(x)
                bak(y)
                x <- y.bak
                y <- x.bak
            }
            points(
                x = x, 
                y = y, 
                ylim = ylim,
                xlim = xlim,
                ylab = "", 
                xlab = "", 
                main = "", 
                col = colors[[K]],
                cex = 0.5, 
                bty = 'n', 
                type = 'l', 
                lwd = 2, 
            )
        })
    }
    if(!is.null(pdf)) dev.off()

}
plotDensities <- function(LIST, colors = RColorBrewer::brewer.pal(name = 'Set2', n = 8), smooth.window = 10, plot.medians = T, plot = T, pdf = NULL, ylim = NULL, xlim = NULL, plot.separately = F, lwd.curves = 1, ...) {
    
    if (length(unique(lengths(LIST))) > 1) {
        # Get absolute minimum and maximum
        MIN <- mins(LIST)
        MAX <- maxs(LIST)
        if (is.null(xlim)) xlim = c(MIN, MAX)
        
        # Compute densities along a single vector
        LIST.2 <- lapply(LIST, function(L) {
            h <- hist(L, breaks = seq(MIN, MAX+1, 1), plot = F)$counts
            h <- h / sum(h)
        })
        if (is.null(ylim)) ylim = c(0, maxs(LIST.2))
    } else {
        LIST.2 <- LIST
        MIN <- 0
        MAX <- unique(lengths(LIST.2))
        if (is.null(xlim)) xlim = c(0, unique(lengths(LIST.2)))
        if (is.null(ylim)) ylim = c(0, maxs(LIST.2))
    }
    
    # Compute smothed densities
    LIST.2 <- lapply(LIST.2, zoo::rollmean, smooth.window)
    
    # Proceed to plotting
    if (plot == T) {
        if(!is.null(pdf)) pdf(pdf)
        
        if (plot.separately == F) {
            plot(NA, ylim = ylim, xlim = xlim, ...)
            lapply(1:length(LIST.2), function(K) {
                par(new = T)
                plot(
                    x = head(MIN:MAX, -smooth.window) + floor(smooth.window / 2),
                    y = LIST.2[[K]][1:length(head(MIN:MAX, -smooth.window) + floor(smooth.window / 2))], 
                    ylim = ylim,
                    xlim = xlim,
                    ylab = "", 
                    xlab = "", 
                    main = "", 
                    axes = F, 
                    col = colors[K], 
                    lwd = lwd.curves, 
                    type = 'l'
                )
            })
            # Plot median lines
            if (plot.medians) {
                lapply(1:length(LIST), function(K) {
                    abline(v = median(LIST[[K]]), lty = 2, col = colors[K], lwd = 0.5)
                    par(xpd = T)
                    text(median(LIST[[K]]), x = median(LIST[[K]]), y = - 0.075 * maxs(LIST.2), col = colors[K])
                    par(xpd = F)
                })
            }
        } else {
            lapply(1:length(LIST.2), function(K) {
                plot(NA, ylim = ylim, xlim = xlim, ...)
                par(new = T)
                plot(
                    x = head(MIN:MAX, -smooth.window) + floor(smooth.window / 2),
                    y = LIST.2[[K]][1:length(head(MIN:MAX, -smooth.window) + floor(smooth.window / 2))], 
                    ylim = ylim,
                    xlim = xlim,
                    ylab = "", 
                    xlab = "", 
                    main = "", 
                    axes = F, 
                    col = colors[K], 
                    lwd = lwd.curves, 
                    type = 'l'
                )
                # Plot median lines
                if (plot.medians) {
                        abline(v = median(LIST[[K]]), lty = 2, col = colors[K], lwd = 0.5)
                        par(xpd = T)
                        text(median(LIST[[K]]), x = median(LIST[[K]]), y = - 0.075 * maxs(LIST.2), col = colors[K])
                        par(xpd = F)
                }
            })
        }
        if(!is.null(pdf)) dev.off()
    }
    
    invisible(LIST.2)
    
}
# Function to get a transparent color
transpCol <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100-percent)*255/100,
                 names = name)
    return(t.col)
}
# Do a N-way venn-diagram from a list (not customizable too much)
plot.venn <- function(list, return = 'mat', ...) {
    
    require(venneuler)
    require(eulerr)
    if (is.null(names(list))) list <- setNames(list, 1:4)
    a <- plot(euler(list), ...)
    prep.mat <- function(list) {
        
        factors <- factor(unique(unlist(list)))
        M <- sapply(list, function(L) {factors %in% L})
        row.names(M) <- factors
        return(M)
        
    }
    mat <- prep.mat(list)
    grid::plotViewport()
    grid::viewport()
    a
    if (return == 'plot') {
        return(a)
    } else {
        return(mat)
    }
    
}
# Plot 2-ways Venn diagram
plot.2way.Venn <- function(vec1, vec2, names = c('a', 'b'), col = c('grey20', 'grey80'), plot.values = T, do.plot = T, ...)
{
    # Compute venn
    vec1 <- as.character(vec1)
    vec2 <- as.character(vec2)
    m <- data.frame(elements = c(vec1, vec2), sets = c(rep(names[1], length(vec1)), rep(names[2], length(vec2))))
    venn <- venneuler::venneuler(m)
    
    if (do.plot) {
        # Calculate coordinates for x-positions of each label
        venn$labels <- c("", "")
        L <- sum(venn$diameters/2)-max(venn$centers[,'x'])+min(venn$centers[,'x'])
        r1 <- venn$diameters[which.min(venn$centers[,'x'])] / 2 
        r2 <- venn$diameters[which.max(venn$centers[,'x'])] / 2 
        G <- min(venn$centers[,'x']) - r1
        p1 <- (2 * r1 - L ) / 2 + G
        p2 <- 2 * r1 + (2 * r2 - L ) / 2 + G
        pmiddle <- (2 * r1) - (L / 2)+ G
        
        # Plot circles
        plot(venn, col = col, ...)
        
        # Add labels and counts
        par(xpd =T)
        text(x = venn$centers[,'x'], y = (venn$centers[,'y'] + venn$diameters/2), c(paste0(names[1], "\n(n=", length(vec1) ,")"), paste0(names[2], "\n(n=", length(vec2) ,")")), pos = 3, cex = 0.8)
        if (plot.values) {
            text(x = c(p1, pmiddle, p2), y = c(0.5, 0.2, 0.5), c(paste0("n=", length(vec2) - sum(vec1 %in% vec2)), paste0("n=", sum(vec1 %in% vec2)), paste0(paste0("n=", length(vec1) - sum(vec1 %in% vec2)))), pos = 3, cex = 0.8)
        }
        par(xpd = F)
    }
    
    shared <- as.character(m$elements[duplicated(m$elements)])
    in_a <- as.character(m$elements[(m$elements %notin% shared) & (m$sets == names[1])])
    in_b <- as.character(m$elements[(m$elements %notin% shared) & (m$sets == names[2])])

    invisible(list(
        'shared' = shared,
        'only.in.1' = in_a, 
        'only.in.2' = in_b
    ))
    
}
# Get a reads length graph from a list of bam files
plotFragmentsLength <- function(pdf, dir = "_bam-files", files = NULL, color = brewer.pal("Spectral", n = 11), labels = NULL) 
{
    
    require(Rsamtools)
    
    # Compute size of bam files
    if (is.null(files)) { FILES <- grep("\\^rm_chrM\\^rm_blacklist\\^q10.bam", list.files(dir, full.names = T), value = T) } else { FILES <- files}
    SIZES <- list()
    for (FILE in FILES) {
        message(FILE)
        what <- c("rname", "strand", "pos", "isize")
        a <- scanBam(FILE, param = ScanBamParam(what = what))[[1]]
        SIZES[[FILE]] <- abs(a$isize) %>% na.remove()
    }
    
    # Compute histogram density plots
    HIST.LIST <- list()
    for (FILE in FILES) {
        sizes <- SIZES[[FILE]]
        sizes[sizes > 600] <- NA
        HIST.LIST[[FILE]] <- hist(sizes %>% na.remove(), plot = F, breaks = seq(0, 600, length.out = 600))$counts
        HIST.LIST[[FILE]] <- HIST.LIST[[FILE]]/max(HIST.LIST[[FILE]])
    }
    
    # Plot density curves
    pdf(pdf, width = 5, height = 5)
    k = 0
    for (FILE in FILES) {
        k <- k + 1
        plot(HIST.LIST[[FILE]],
            xlim = c(0, 600), ylim = c(0,1),
            col = color[k],
            type = 'l',
            lwd = 0.5,
            main = FILE,
            xlab = "", 
            ylab = "", 
            bty = "n"
        )
        polygon(c(0:599, 0), c(HIST.LIST[[FILE]][1:599], 0, 0), col = paste0(color[k], "64"))
    }
    dev.off()

}
## Function to generate motif occurence track using a sliding window
plotMotifBW <- function(MOTIF, WINDOW = 50, FASTA_PATH = "~/shared/sequences/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa", BW_OUTPUT = NULL, STRETCH.MOTIF = T) {
    library(GenomicRanges)
    if (is.null(BW_OUTPUT)) {
        BW_OUTPUT <- paste0(MOTIF, "-occurence_", WINDOW, "bp-sw.bw")
    }
    message("--- Plotting motif << ", MOTIF, " >> on ", basename(FASTA_PATH), " using a sliding window of ", WINDOW, "bp. ---")
    # Reading fasta file 
    fasta <- Biostrings::readBStringSet(FASTA_PATH)
    # Mapping motif occurences and center it
    if (STRETCH.MOTIF == T) {
        counts_motif <- shift(reduce(GRanges(vmatchPattern(MOTIF, fasta))), WINDOW / 2)
    } else {
        counts_motif <- shift(GRanges(vmatchPattern(MOTIF, fasta)), WINDOW / 2)
    }
    # Computing coverage of motif
    coverage_motif <- coverage(counts_motif)
    # Compute sliding window
    sliding_sum <- GRanges(runmean(coverage_motif, WINDOW))
    message("-- Writing motif occurence track to: ", BW_OUTPUT)
    rtracklayer::export.bw(sliding_sum, BW_OUTPUT)
}
plotProfile <- function(granges, bam.file) {
    what <- c("rname", "strand", "pos", "qwidth", "mapq", "flag", "mpos", "isize")
    a <- scanBam(bam.file, param = ScanBamParam(what = what), asMates = T)[[1]]
    aligned.reads <- !is.na(a$pos) # Only aligned reads
    qual.cutoff <- a$mapq >= 10L & !(is.na(a$mapq)) # q10 quality cutoff
    subset <- aligned.reads & paired.reads & qual.cutoff
    frag.granges <- GRanges(
        seqnames = a$rname[subset],
        ranges = IRanges(start = ifelse(a$pos[subset] < a$mpos[subset], a$pos[subset], a$mpos[subset]), width = abs(a$isize[subset]) )
    )
    cov <- coverage(frag.granges)
    RleList(lapply(frag.granges, rangeCoverage, cov))
}
# Plot BW signal over a set of GRanges
plot.bw.over.granges <- function(granges, bw.file = '', bw.from.granges = NULL, FLANK = c(500, 500), stranded = T, use.mean = T, BIN = 10, plotEE = T, YLIM = NULL, XLAB = NULL, YLAB = NULL, COL = NULL, ...) {
    
    require(miscset)
    
    # Deconvolute and resize GRanges 
    if (!is.null(FLANK)) {
        granges <- deconvolve.bidirectional.proms(granges)
        granges <- resize(granges, fix = 'end', width = FLANK[1])
        granges <- resize(granges, fix = 'start', width = FLANK[1] + FLANK[2])
    } else {
        granges <- deconvolve.bidirectional.proms(granges)
    }
    
    # Compute coverage over granges
    if (is.null(bw.from.granges)) {
        score <- get.granges.scores.from.bw(granges, bw.file, FUN = 'cov')
    } else {
        if (class(bw.from.granges) == 'GRanges') {
            score <- coverage(bw.from.granges)[granges]
        } else if (class(bw.from.granges) == 'SimpleRleList') {
            score <- bw.from.granges[granges]
        }
    }
    message('-- Track imported')
    score <- suppressWarnings(suppressMessages(matrix(as.vector(unlist(score)), nrow = length(granges), byrow = T)))
    score.flipped <- t(apply(score, 1, rev))
    if (stranded) {
        score <- matrix(sapply(1:nrow(score), function(K) {if((as.vector(strand(granges)) == '-')[K]) {score.flipped[K,]} else {score[K,]}}), nrow = length(granges), byrow = T)
    }
    message('-- Scores computed')

    # Compute specific metrics to plot
    medians <- apply(score, 2, median)
    means <- apply(score, 2, mean)
    smoothed_line <- if (use.mean) {zoo::rollmean(means, BIN)} else {zoo::rollmean(medians, BIN)}
    stderror <- zoo::rollmean(apply(score, 2, function(n) { sd(n, na.rm=TRUE) / sqrt( sum(!is.na(n)) ) }), BIN)
    conint  <- zoo::rollmean(apply(score, 2, function (n) { qt(0.975,sum(!is.na(n)))*sd(n,na.rm=TRUE)/sqrt(sum(!is.na(n))) }), BIN)
    topEE <- smoothed_line + stderror
    bottomEE <- smoothed_line - stderror
    message('-- Stats computed')

    # Get other parameters
    YLIM <- if (!is.null(YLIM)) {YLIM} else if (plotEE) {c(min(bottomEE), max(topEE))} else {c(min(smoothed_line), max(smoothed_line))}
    XLAB <- ifelse(!is.null(XLAB), XLAB, 'Distance from beginning of GRanges')
    YLAB <- ifelse(!is.null(YLAB), YLAB, 'Score')
    COL <- ifelse(!is.null(COL), COL, '#000000')
    COLTRANSP <- paste0(COL, '1A')
    
    # Proceed to plotting
    plot(
        smoothed_line, 
        type = 'l', 
        ylim = YLIM, 
        xlab = XLAB, ylab = YLAB, 
        col = COL,
        ...
    )
    
    # Add error-estimate transparent polygon
    if (plotEE) {
        polygon(
            x = c(seq_along(smoothed_line), rev(seq_along(smoothed_line))), 
            y = c(topEE, rev(bottomEE)), 
            ylim = YLIM, 
            col = COLTRANSP, 
            border = NA
        )
    }
    
}
plot.cov.line <- function(score, use.mean = T, BIN = 10, plotEE = T, YLIM = NULL, XLAB = NULL, YLAB = NULL, COL = NULL, verbose = F, return.extrema = F, ...) {
    
    # Compute specific metrics to plot
    medians <- apply(score, 2, median)
    means <- apply(score, 2, mean)
    smoothed_line <- if (use.mean) {zoo::rollmean(means, BIN)} else {zoo::rollmean(medians, BIN)}
    stderror <- zoo::rollmean(apply(score, 2, function(n) { sd(n, na.rm=TRUE) / sqrt( sum(!is.na(n)) ) }), BIN)
    conint  <- zoo::rollmean(apply(score, 2, function (n) { qt(0.975,sum(!is.na(n)))*sd(n,na.rm=TRUE)/sqrt(sum(!is.na(n))) }), BIN)
    topEE <- smoothed_line + conint
    bottomEE <- smoothed_line - conint
    MIN <- min(bottomEE) * 0.9
    MAX <- max(topEE) * 1.1
    if (verbose) message('.. Stats computed')
    
    if (return.extrema) {
        return(c(MIN, MAX))
    } else {

        # Get other parameters
        YLIM <- if (!is.null(YLIM)) {YLIM} else if (plotEE) {c(min(bottomEE), max(topEE))} else {c(min(smoothed_line), max(smoothed_line))}
        XLAB <- ifelse(!is.null(XLAB), XLAB, 'Distance from beginning of GRanges')
        YLAB <- ifelse(!is.null(YLAB), YLAB, 'Score')
        COL <- ifelse(!is.null(COL), COL, '#000000')
        COLTRANSP <- paste0(COL, '1A')
        
        # Proceed to plotting
        plot(
            smoothed_line, 
            type = 'l', 
            ylim = YLIM, 
            xlab = XLAB, ylab = YLAB, 
            col = COL,
            ...
        )
        
        # Add error-estimate transparent polygon
        if (plotEE) {
            polygon(
                x = c(seq_along(smoothed_line), rev(seq_along(smoothed_line))), 
                y = c(topEE, rev(bottomEE)), 
                ylim = YLIM, 
                col = COLTRANSP, 
                border = NA
            )
        }
        
        return(smoothed_line)
    }
    
} # To use with get.cov.matrix output # This one takes a rectangular matrix already computed (from get.cov.matrix)
# Compute a V-matrix of a fragments bam file over a given set of granges
compute.Vmat <- function(reads.granges, target.granges, XAXIS.CUTOFF, YAXIS.CUTOFF) {
    
    center.targets <- GRanges(
        seqnames = seqnames(target.granges), 
        IRanges(
            start = start(target.granges) + floor((end(target.granges) - start(target.granges)) / 2), 
            width = 1
        )
    )
    extended.targets <- center.targets + 2000
    center.reads <- GRanges(
        seqnames = seqnames(reads.granges), 
        IRanges(
            start = start(reads.granges) + floor((end(reads.granges) - start(reads.granges)) / 2),
             width = 1
         )
     )
    granges <- reads.granges[center.reads %over% extended.targets]
    center.reads.subset <- center.reads[center.reads %over% extended.targets]
    center.target.per.read <- center.targets[subjectHits(distanceToNearest(center.reads.subset, center.targets))]
    
    dists <- factor(as.character(start(center.target.per.read) - start(center.reads.subset)), levels = -2000:2000)
    widths <- factor(as.character(width(granges)), levels = 0:YAXIS.CUTOFF)
    
    tab <- table(dists, widths)
    tab <- tab[seq(round(nrow(tab) / 2) - XAXIS.CUTOFF, round(nrow(tab) / 2) + XAXIS.CUTOFF, 1), 1:YAXIS.CUTOFF]
    
    return(tab)
}
# Plot distance tree of a matrix 
plotDistTree <- function(matrix, rows = NULL) {
    if (!is.null(rows))
        matrix <- matrix[rows,]
    cors <- as.data.frame(matrix(0, ncol = ncol(matrix), nrow = ncol(matrix)))    
    for (ROW in 1:ncol(matrix)) {
        for (COLUMN in 1:ncol(matrix)) {
            row <- matrix[, ROW]
            col <- matrix[, COLUMN]
            cors[ROW, COLUMN] <- cor.test(row, col)$estimate
        }
    }
    row.names(cors) <- colnames(matrix)
    colnames(cors) <- colnames(matrix)
    plot(as.dendrogram(hclust(dist(cors, method = 'euc'))), horiz= T, xlab = 'Euclid. distance', main = 'Comparison with Cao')
    text(x = max(hclust(dist(cors, method = 'euc'))[[2]]) * c(1, 0.85, 0.7, 0.55, 0.4), y = 0.5, sapply(apply(rbind(paste0(order.tissues[1:5]), paste0('CAO_', order.tissues[1:5])), 2, function(VECS) VECS) %>% as.data.frame(stringsAsFactors = F) %>% as.list, function(COORDS) {cors[COORDS[[1]], COORDS[[2]]]}) %>% round(3))
}

# ---------------------------------------------------------------------------- #
# ---------------- GRANGES-RELATED FUNCTIONS
# ---------------------------------------------------------------------------- #
# Bam-file related functions
read.bam <- function(f, get.unique.reads = T, get.paired = F, filter.lengths = c(0, 0)) {
    
    # Import bam file
    message('>>> Import bam file ', f)
    require(Rsamtools)
    what <- c("rname", "strand", "pos", "qwidth", "mapq", "flag", "mpos", "isize")
    a <- scanBam(f, param = ScanBamParam(what = what), asMates = T)[[1]]
    
    # Subset valid reads
    aligned.reads <- !is.na(a$pos) # Only aligned reads
    paired.reads <- !is.na(a$mpos) & !is.na(a$pos - a$mpos > 0) & (a$pos - a$mpos) != 0 & abs(a$isize) > 0 # Only reads with a paired mate and an absolute fragment length > 0
    qual.cutoff <- a$mapq >= 10L & !(is.na(a$mapq)) # q10 quality cutoff
    if (get.paired) {
        subset <- aligned.reads & paired.reads & qual.cutoff
        frag.granges <- GRanges(
            seqnames = a$rname[subset],
            ranges = IRanges(start = ifelse(a$pos[subset] < a$mpos[subset], a$pos[subset], a$mpos[subset]), width = abs(a$isize[subset]) ),
            strand = a$strand[subset]
        )
    } else {
        subset <- aligned.reads & qual.cutoff
        frag.granges <- GRanges(
            seqnames = a$rname[subset],
            ranges = IRanges(start = a$pos[subset], width = a$qwidth[subset] ),
            strand = a$strand[subset]
        )
    }
    
    # Filter imported reads / fragments
    if (get.unique.reads) { frag.granges <- frag.granges[!duplicated(frag.granges)] }
    if (filter.lengths[1] > 0) { frag.granges <- frag.granges[width(frag.granges) >= filter.lengths[1]] }
    if (filter.lengths[2] > 0) { frag.granges <- frag.granges[width(frag.granges) <= filter.lengths[2]] }
    
    # Return filtered reads/fragments as a GRanges
    return(frag.granges)
    
}
subset.bam <- function(frags, min.length = 0, max.length = 0, subset.granges = NULL) {
    if (min.length > 0) { frags <- frags[width(frags) >= min.length] }
    if (max.length > 0) { frags <- frags[width(frags) <= max.length] }
    if (length(subset.granges) > 0) { frags <- frags[frags %over% subset.granges] }
    return(frags)
}
# Plot average of bam file over GRanges (quick way to do a seqplot)
rangeCoverage <- function(gr, cvg) {
    gr <- reduce(gr)
    ans <- unlist(cvg[gr], use.names=FALSE)
    if (runValue(strand(gr))[[1L]] == "-")
        ans <- rev(ans)
    ans
}
# Erase all metadata from a GRanges object
simplifyGRanges <- function(granges) {
    
    granges2 <- granges
    mcols(granges2) <- NULL
    return(granges2)

}
# Plot venn diagram of a list of genes against all classes from Cao data
checkAgainstCao <- function(genes, against = 'cao', pdf = NULL, label = NULL, color = NULL, new.dev = T, main = "", par.mfrow = c(2, 3), ...) {
    if (!exists("cao06.2")) {
        cao06=read.table("~/_data/Cao2016_TableS6.txt", header=T)
        cao06$max.tissue=as.character(cao06$max.tissue)
        cao06[cao06$max.tissue == "Gonad", "max.tissue"]="Germline"
        cao06[cao06$max.tissue == "Intestine", "max.tissue"]="Intest."
        cao06[cao06$max.tissue == "Hypodermis", "max.tissue"]="Hypod."
        cao06[cao06$max.tissue == "Body_wall_muscle", "max.tissue"]="Muscle"
        cao06.2=cao06[cao06$ratio >= 3 & cao06$qval <= 0.05,]
    }
    hypod.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Hypod.']
    neurons.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Neurons']
    germline.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Germline']
    muscle.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Muscle']
    intest.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Intest.']
    pharynx.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Pharynx']
    glia.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Glia']
    ubiq.genes.cao <- apply(cao03, 1, function(ROW) {all(ROW >= 5)}) %>% which() %>% names() %>% '['(. %notin% c(hypod.genes.cao, germline.genes.cao, neurons.genes.cao, muscle.genes.cao, intest.genes.cao, pharynx.genes.cao, glia.genes.cao))
    others <- names(genes.gtf)[genes.gtf$is.prot.cod & names(genes.gtf) %notin% cao06.2$gene_id & names(genes.gtf) %in% cao06$gene_id]
    cao_genes <- reshape2::melt(list(
        "Germline" = germline.genes.cao, 
        "Neurons" = neurons.genes.cao, 
        "Muscle" = muscle.genes.cao, 
        "Hypod." = hypod.genes.cao, 
        "Intest." = intest.genes.cao, 
        "Pharynx" = pharynx.genes.cao, 
        "Glia" = glia.genes.cao, 
        "Ubiq." = ubiq.genes.cao, 
        "Others?" = others)
    )
    
    if(!is.null(pdf)) {
        pdf(pdf, width = 10, height = 10)
    } else if (new.dev) {
        dev.new(width = 20, height = 14)
    }
    par(mfrow = par.mfrow)
    if (against == 'cao') {
        cao06 <- read.table("~/_data/Cao2016_TableS6.txt", header=T)
        cao06$max.tissue <- as.character(cao06$max.tissue)
        cao06[cao06$max.tissue == "Gonad", "max.tissue"]="Germline"
        cao06[cao06$max.tissue == "Intestine", "max.tissue"]="Intest."
        cao06[cao06$max.tissue == "Hypodermis", "max.tissue"]="Hypod."
        cao06[cao06$max.tissue == "Body_wall_muscle", "max.tissue"]="Muscle"
        cao06.2 <- cao06[cao06$ratio >= 3 & cao06$qval <= 0.05,]
        hypod.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Hypod.']
        neurons.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Neurons']
        germline.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Germline']
        muscle.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Muscle']
        intest.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Intest.']
        pharynx.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Pharynx']
        glia.genes.cao <- cao06.2$gene_id[cao06.2$max.tissue == 'Glia']
        ubiq.genes.cao <- apply(cao03, 1, function(ROW) {all(ROW >= 5)}) %>% which() %>% names() %>% '['(. %notin% c(hypod.genes.cao, germline.genes.cao, neurons.genes.cao, muscle.genes.cao, intest.genes.cao, pharynx.genes.cao, glia.genes.cao))
        others <- names(genes.gtf)[genes.gtf$is.prot.cod & names(genes.gtf) %notin% cao06.2$gene_id & names(genes.gtf) %in% cao06$gene_id]
        cao_genes <- reshape2::melt(list(
            "Germline" = germline.genes.cao, 
            "Neurons" = neurons.genes.cao, 
            "Muscle" = muscle.genes.cao, 
            "Hypod." = hypod.genes.cao, 
            "Intest." = intest.genes.cao, 
            "Pharynx" = pharynx.genes.cao, 
            "Glia" = glia.genes.cao, 
            "Ubiq." = ubiq.genes.cao, 
            "Others?" = others)
        )
        
        if(!is.null(label)) {LABEL = label} else {LABEL = 'input list'}
        if(!is.null(color)) {COLOR = color} else {COLOR = 'grey'}
        
        lapply(
            c("Germline", "Neurons", "Muscle", "Hypod.", "Intest.", "Ubiq."), 
            function(TISSUE) {
                plot.2way.Venn(
                    genes, 
                    cao_genes[cao_genes$L1 == TISSUE,]$value,
                    names = c(LABEL, paste0(TISSUE, '-enriched (Cao)')),
                    col = c(
                        COLOR, 
                        colortransp[c(1:5, 33)][which(c("Germline", "Neurons", "Muscle", "Hypod.", "Intest.", "Ubiq.") %in% TISSUE)]
                    )
                )
            }
        )
        mtext(side = 3, line = 1, outer = T, text = 'Comparison with Cao data')
    } else if (against == 'lcap') {
        NULL
    }
    if(!is.null(pdf)) {
        dev.off()
    }

}
# Shuffle a grange object
shuffle.GRanges <- function(granges, opt.shuffle = '-chrom -noOverlapping', excl.itself = T, excl.blacklist = '~/shared/ce11/ce10_blacklist.bed', other.excl = c(), genome = '~/shared/ce11/ce11.chrom.sizes') {
    
    export.bed(simplifyGRanges(granges), con = 'tmp.bed')
    excl.granges <- GRanges()
    if (excl.itself) {
        excl.granges <- c(excl.granges, simplifyGRanges(granges))
    } 
    if (length(excl.blacklist) > 0) {
        excl.granges <- c(excl.granges, simplifyGRanges(import.bed(excl.blacklist)))
    }
    if (length(other.excl) > 0) {
        excl.granges <- c(excl.granges, simplifyGRanges(import.bed(other.excl)))
    }
    export.bed(excl.granges, con = 'excl.bed')
    system(sprintf("shuffleBed %s -excl excl.bed -i %s -g %s | sort -k1,1 -k2.2n > tmp2.bed", opt.shuffle, 'tmp.bed', genome))
    t <- import.bed('tmp2.bed')
    system('rm tmp.bed  tmp2.bed excl.bed')
    return(t)
    
}
# Import scores from bigwig over a GRange object (using UCSC toolkit)
get.granges.scores.from.bw <- function(granges, bw.file, FUN = 'max') {
    
    names(granges) <- as.character(1:length(granges))
    bed.file <- tempfile()    
    rtracklayer::export.bed(granges, con = bed.file)
    if (FUN == 'sum') {
        system(sprintf("~/bin/UCSC_scripts/bigWigAverageOverBed %s %s out.tab", bw.file, bed.file))
        scores <- read.table("out.tab", sep = '\t', header = F)[,4]
    } else if (FUN == 'mean0'){
        system(sprintf("~/bin/UCSC_scripts/bigWigAverageOverBed %s %s out.tab", bw.file, bed.file))
        scores <- read.table("out.tab", sep = '\t', header = F)[,5]
    } else if (FUN == 'mean'){
        system(sprintf("~/bin/UCSC_scripts/bigWigAverageOverBed %s %s out.tab", bw.file, bed.file))
        scores <- read.table("out.tab", sep = '\t', header = F)[,6]
    } else if (FUN == 'max'){
        system(sprintf("~/bin/UCSC_scripts/bigWigAverageOverBed %s %s out.tab -minMax", bw.file, bed.file))
        scores <- read.table("out.tab", sep = '\t', header = F)[,8]
    } else if (FUN == 'cov') {
        cov <- rtracklayer::import.bw(bw.file, as = 'Rle')
        scores <- cov[granges]
    }
    if (file.exists('out.tab')) {system('rm out.tab')}
    return(scores)
}
# Align granges.proms to their TSS
tss.align <- function(granges, WIDTH.bef, WIDTH.aft, do.not.deconvolve = F) {
    
    if (any(strand(granges) == '*')) {
        if (do.not.deconvolve == F) {
            message('.. Splitting bidirectional promoters in 2...')
            granges <- deconvolve.bidirectional.proms(granges)
        } else {
            stop('Some granges are bidirectional. Set do.not.deconvolve argument to F. Aborting.')
        }
    }
    
    ranges(granges) <- IRanges(
        start = ifelse(as.vector(strand(granges)) == '+', (granges$TSS.fwd - WIDTH.bef), (granges$TSS.rev - WIDTH.aft + 1)),
        width = WIDTH.aft + WIDTH.bef,
        names = names(ranges(granges))
    )
    return(granges)
    
}
deconvolve.bidirectional.proms <- function(granges) {

    unid <- granges[strand(granges) == '+' | strand(granges) == '-']
    bid <- granges[strand(granges) == '*']
    bid.fwd <- bid
    strand(bid.fwd) <- '+'
    bid.rev <- bid
    strand(bid.rev) <- '-'
    granges_shifted <- sort(c(unid, bid.fwd, bid.rev), ignore.strand = T)
    
    # Get single WormBase names
    #mcols(granges)[,c('locus', 'WormBaseID')] %>% as_tibble() %>% dplyr::mutate('strand' = as.vector(strand(granges))) %>% dplyr::mutate(uniqueWormBaseID = strsplit(as.character(WormBaseID),',')) %>% unnest(uniqueWormBaseID) -> t
    #granges_shifted$singleWormBaseID <- t[(grepl(',', t$WormBaseID) & t$strand == '*') | (t$strand != '*' & !duplicated(t[,1:2])), 'uniqueWormBaseID'] %>% unlist()
    
    return(granges_shifted)
}
# Get rectangular matrix of coverage over a set of granges, from a bw file
get.cov.matrix <- function(granges, bw.file, norm = 'none', center = F, flank = NULL, stranded = T, split.bid.proms = T, verbose = T) {
    
    # Import bw file
    if (verbose) message('.. Importing bw.file...')
    if (class(bw.file)[[1]] == "character") {
        scores <- import(bw.file, as = 'Rle')
    } else if (class(bw.file) == "UCSCData" | class(bw.file)[[1]] == "SimpleRleList" | class(bw.file)[[1]] == "CompressedRleList") {
        scores <- bw.file
    }
    
    # Deconvolute and resize GRanges 
    if (any(strand(granges) == '*')) {
        if (stranded == T & split.bid.proms == T) {
            if (verbose) message('.. Bi-directional segments found and split...')
            granges <- deconvolve.bidirectional.proms(granges)
        } else {
            if (verbose) message('.. Bi-directional promoters found but not split...')
        }
    } else {
        if (verbose) message('.. No bi-directional segments found. Continuing...')
    }
        
    if (length(unique(width(granges))) == 1) { # Cases where all the GRanges have the same width
        
        if (center == T) {
            if (verbose) message('.. Centering GRanges...')
            granges <- resize(granges, width = 1, fix = 'center')
        } else {
            granges <- granges
        }
        if (length(flank) > 0) {
            if (verbose) message('.. Resizing GRanges...')
            granges <- resize(granges, fix = 'end', width = flank[1])
            granges <- resize(granges, fix = 'start', width = flank[1] + flank[2])
        }
        
        # Subset scores over GRanges
        scores.subset <- scores[granges]
        
        # Turn it into a rectangular matrix and flip reverse strand scores
        if (verbose) message('.. Converting scores into matrix...')
        scores.subset <- suppressWarnings(suppressMessages(matrix(as.vector(unlist(scores.subset)), nrow = length(granges), byrow = T)))
        scores.subset.flipped <- t(apply(scores.subset, 1, rev))
        if (stranded) {
            scores.subset <- matrix(sapply(1:nrow(scores.subset), function(K) {if((as.vector(strand(granges)) == '-')[K]) {scores.subset.flipped[K,]} else {scores.subset[K,]}}), nrow = length(granges), byrow = T)
        }
        
    } else { # Cases where all the GRanges do *NOT* have the same width
        
        widths.granges <- width(granges)
        # Get flanking ends
        if (length(flank) > 0) {
            if (verbose) message('.. Resizing GRanges...')
            granges.start.flank <- shift(resize(granges, fix = ifelse(strand(granges) == '+', 'start', 'end'), width = flank[1]), - flank[1])
            granges.end.flank <- shift(resize(granges, fix = ifelse(strand(granges) == '+', 'end', 'start'), width = flank[1]), - flank[1])
        }
        # Subset scores over GRanges
        scores.subset <- scores[granges][widths.granges > 20]
        # Bring all the scores to the same length
        if (verbose) message('.. Shrinking vectors...')
        scores.subset.shrinked <- mclapply(
            scores.subset, 
            function(score) {apply(zoo::rollapply(seq(1, length(score), length.out = 101), 2, function(RANGE) {c(floor(RANGE[1]),floor(RANGE[2]))}), 1, function(Ks) {mean(score[Ks[1]:Ks[2]])})},
            mc.cores = 50
        )
        
        # Collate the flanking matrices to the central matrix and turn it into a rectangular matrix
        if (verbose) message('.. Converting scores into matrix...')
        scores.start.flank <- suppressWarnings(suppressMessages(matrix(as.vector(unlist(scores[granges.start.flank][widths.granges > 20])), nrow = length(granges.start.flank), byrow = T)))
        scores.end.flank <- suppressWarnings(suppressMessages(matrix(as.vector(unlist(scores[granges.end.flank][widths.granges > 20])), nrow = length(granges.end.flank), byrow = T)))
        scores.subset.shrinked <- suppressWarnings(suppressMessages(matrix(as.vector(unlist(scores.subset.shrinked)), nrow = length(scores.subset), byrow = T)))
        scores.subset <- cbind(scores.start.flank, scores.subset.shrinked, scores.end.flank)
        
        # Flip reverse strand scores
        scores.subset.flipped <- t(apply(scores.subset, 1, rev))
        if (stranded) {
            scores.subset <- matrix(sapply(1:nrow(scores.subset), function(K) {if((as.vector(strand(granges)) == '-')[K]) {scores.subset.flipped[K,]} else {scores.subset[K,]}}), nrow = length(granges), byrow = T)
        }

    }
    
    # Normalize matrix
    if (norm == 'zscore') {
        scores.subset <- scale(scores.subset)
    } else if (norm == 'log2') {
        scores.subset <- log2(scores.subset)
    }
    
    # Return matrix
    return(scores.subset)
    
}
# Import list of bigwig files
import.list.bw <- function(folder = '_bw-files/', suffix = '_YA_ATAC_combined.bw', as = 'Rle') {
    bw.files <- list.files(folder, pattern = suffix, full.names = T)
    list <- mclapply(bw.files, function(bw.file) {
        message('Importing ', bw.file, '...')
        import(bw.file, as = as)
    }, mc.cores = min(10, length(bw.files)))
    names(list) <- gsub('Gonad', 'Germline', gsub('.*/', '', gsub(suffix, '', bw.files)))
    return(list)
}
scale.list.bw <- function(list.bw) {
    l <- mclapply(list.bw, function(bw.file) {
        message('Get Z-score...')
        cov <- bw.file
        ucov <- unlist(cov)
        mi <- mean(ucov)
        mu <- sd(ucov)
        zsc <- (cov-mi)/mu
        return(zsc)
    }, mc.cores = length(list.bw))
    return(l)
}
# Wrapper to plot seqplot lines
wrapper.plot.cov.lines <- function(
    list.granges, # The regions of interest
    list.bw.files, # Can be either a bw.file or an alreday imported Rle
    list.COL = NULL,
    split.bid.proms = T,
    center = F,
    flank = NULL,
    stranded = T,
    use.mean = T,
    BIN = 10,
    plotEE = T, 
    auto.scale = c(0.05, 0.95), 
    scale.max = T, 
    YLIM = NULL, 
    XLIM = NULL,
    XLAB = NULL,
    YLAB = NULL,
    verbose = T,
    by.granges = T,
    plot.central = T,
    plot.legend = F,
    ...) {
        
    if (class(list.granges) == "GRanges" & (class(list.bw.files) == 'character' | class(list.bw.files) == "SimpleRleList")) { # Case where there is only 1 bigWig file and 1 set of granges
        
        # Define starting parameters
        XLAB <- if (is.null(XLAB)) {'Segments coordinates'} else {XLAB}
        YLAB <- ifelse(is.null(YLAB), 'Score', YLAB)
        # Define colors
        if (!is.null(list.COL)) {
            list.COL <- rep(list.COL, 5)
        } else {
            list.COL <- RColorBrewer::brewer.pal('Set2', n = 8)
        }
        granges <- list.granges
        bw.file <- list.bw.files
        if (verbose) message('>>> Plotting granges')
        # Determining YLIM 
        if (verbose) message('.. Determining Y boundaries...')
        if (!is.null(auto.scale)) {
            max <- quantile(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), auto.scale[2])
            max <- quantile(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), auto.scale[1])
            YLIM <- c(min, max)
        } else if (scale.max) {
            max <- plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), verbose = F, return.extrema = T)[2]
            min <- plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F),  verbose = F, return.extrema = T)[1]
            YLIM <- c(min, max)
        } else if (!is.null(YLIM)) {
            YLIM <- YLIM
        } else {
            stop('Please provide some Y limits. Aborting.')
        }
        # Initiate plot
        if (verbose) message('.. Initiating graph...')
        if (is.null(XLIM)) XLIM <- c(1, flank[1] + flank[2])
        plot(
            NA, 
            bty = 'l',
            xlim = XLIM,
            ylim = YLIM,
            xlab = ifelse(length(XLAB) == 1, XLAB, XLAB[k.g]), 
            ylab = YLAB, 
            axes = F, 
            ...
        )
        box(bty='l')
        axis(2)
        #axis(1, at = seq(XLIM[1], XLIM[2], length.out = 3), labels = c(-unique(width(granges) / 2), '0', unique(width(granges) / 2)))
        # Add central line
        if (plot.central) abline(v = seq(XLIM[1], XLIM[2], length.out = 3)[2], lty = 3)
        # Plot track
        if (verbose) message('.. Computing values...')
        mat <- get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F)
        par(new = T)
        if (verbose) message('.. Plot track...')
        plot.cov.line(mat, use.mean, BIN, plotEE, YLIM, xlim = XLIM, axes = F, XLAB = '', YLAB = '', COL = list.COL[1])
    
    } else if (length(list.granges) > 1 & length(list.bw.files) > 1) { # Case where there are several bigWig file and 1 set of several granges 
        
        # Check that there is 1 unique set of GRanges in list.granges
        if (class(list.granges[[1]]) != "GRanges") {stop('Please provide a unique set of GRanges in list.granges.\nIf a plot showing the coverage of a *single* track over multiple GRanges is desired, please provide a *single* bw.file in list.bw.files')}
        
        # Define starting parameters
        k.g <- 0
        k.b <- 0
        XLAB <- if (is.null(XLAB)) {'Segments coordinates'} else {XLAB}
        YLAB <- ifelse(is.null(YLAB), 'Score', YLAB)
        # Define colors
        if (!is.null(list.COL)) {
            list.COL <- rep(list.COL, 5)
        } else {
            list.COL <- RColorBrewer::brewer.pal('Set2', n = 8)
        }
        
        if (by.granges) {
            
            # Loop through Granges
            for (k.g in 1:length(list.granges)) {
                # Message
                granges <- list.granges[[k.g]]
                if (verbose) message('>>> Plotting granges ', k.g, '/', length(list.granges))
                # Determining YLIM 
                if (verbose) message('.. Determining Y boundaries...')
                if (!is.null(auto.scale)) {
                    max <- max(sapply(list.bw.files, function(bw.file) {quantile(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), auto.scale[2])}))
                    min <- min(sapply(list.bw.files, function(bw.file) {quantile(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), auto.scale[1])}))
                    YLIM <- c(min, max)
                } else if (scale.max) {
                    max <- max(sapply(list.bw.files, function(bw.file) {plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[2]}))
                    min <- min(sapply(list.bw.files, function(bw.file) {plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[1]}))
                    YLIM <- c(min, max)
                } else if (!is.null(YLIM)) {
                    YLIM <- YLIM
                } else {
                    stop('Please provide some Y limits. Aborting.')
                }
                # Initiate plot
                if (is.null(XLIM)) XLIM <- c(1, (unique(width(granges)) - BIN))
                plot(
                    NA, 
                    bty = 'l',
                    xlim = XLIM,
                    ylim = YLIM,
                    xlab = ifelse(length(XLAB) == 1, XLAB, XLAB[k.g]), 
                    ylab = YLAB, 
                    axes = F, 
                    ...
                )
                box(bty='l')
                axis(2)
                axis(1, at = seq(XLIM[1], XLIM[2], length.out = 3), labels = c(-unique(width(granges) / 2), '0', unique(width(granges) / 2)))
                # Add central line
                if (plot.central) abline(v = seq(XLIM[1], XLIM[2], length.out = 3)[2], lty = 3)
                # Looping through bw.files
                for (k.b in 1:length(list.bw.files)) {
                    # Message
                    bw.file <- list.bw.files[[k.b]]
                    if (verbose) message('.. Plotting bigwig ', k.b, '/', length(list.bw.files), '...')
                    # Plot
                    mat <- get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F)
                    par(new = T)
                    plot.cov.line(mat, use.mean, BIN, plotEE, YLIM, xlim = c(1,(ncol(mat) - BIN)), axes = F, XLAB = '', YLAB = '', COL = list.COL[k.b])
                }
            }
        
        } else if (by.granges == F) {
        
            # Loop through BW
            for (k.b in 1:length(list.bw.files)) {
                
                # Message
                bw.file <- list.bw.files[[k.b]]
                if (verbose) message('>>> Plotting bw.file ', k.b, '/', length(list.bw.files))
                # Determining YLIM 
                if (verbose) message('.. Determining Y boundaries...')
                if (!is.null(auto.scale)) {
                    max <- max(sapply(list.granges, function(granges) {quantile(get.cov.matrix(granges, bw.file, center, flank, stranded, split.bid.proms, verbose = F), auto.scale[2])}))
                    min <- min(sapply(list.granges, function(granges) {quantile(get.cov.matrix(granges, bw.file, center, flank, stranded, split.bid.proms, verbose = F), auto.scale[1])}))
                    YLIM <- c(min, max)
                } else if (scale.max) {
                    max <- max(sapply(list.granges, function(granges) {plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[2]}))
                    min <- min(sapply(list.granges, function(granges) {plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[1]}))
                    YLIM <- c(min, max)
                } else if (!is.null(YLIM)) {
                    YLIM <- YLIM
                } else {
                    stop('Please provide some Y limits. Aborting.')
                }
                # Initiate plot
                if (is.null(XLIM)) XLIM <- c(1, (unique(width(list.granges[[1]])) - BIN))
                plot(
                    NA, 
                    bty = 'l',
                    xlim = XLIM,
                    ylim = YLIM,
                    xlab = ifelse(length(XLAB) == 1, XLAB, XLAB[k.b]), 
                    ylab = YLAB, 
                    axes = F, 
                    ...
                )
                box(bty='l')
                axis(2)
                axis(1, at = seq(XLIM[1], XLIM[2], length.out = 3), labels = c(-unique(width(list.granges[[1]]) / 2), '0', unique(width(list.granges[[1]]) / 2)))
                # Add central line
                if (plot.central) abline(v = seq(XLIM[1], XLIM[2], length.out = 3)[2], lty = 3)
                # Looping through granges
                for (k.g in 1:length(list.granges)) {
                    # Message
                    granges <- list.granges[[k.g]]
                    if (verbose) message('.. Plotting granges ', k.g, '/', length(list.granges), '...')
                    # Plot
                    mat <- get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F)
                    par(new = T)
                    plot.cov.line(mat, use.mean, BIN, plotEE, YLIM, xlim = c(1,(ncol(mat) - BIN)), axes = F, XLAB = '', YLAB = '', COL = list.COL[k.g])
                }
            }
            
        }
    
    } else if (length(list.granges) > 1 & length(list.bw.files) == 1 & class(list.granges[[1]]) == "GRanges") { # Case where there is only 1 bigWig file and 1 set of granges
        
        # Define starting parameters
        k.g <- 0
        XLAB <- if (is.null(XLAB)) {'Segments coordinates'} else {XLAB}
        YLAB <- ifelse(is.null(YLAB), 'Score', YLAB)
        bw.file <- list.bw.files[[1]]
        # Define colors
        if (!is.null(list.COL)) {
            list.COL <- rep(list.COL, 5)
        } else {
            list.COL <- RColorBrewer::brewer.pal('Set2', n = 8)
        }
        # Define YLIM 
        if (verbose) message('.. Determining Y boundaries...')
        if (!is.null(auto.scale)) {
            max <- max(sapply(list.granges, function(granges) {quantile(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), auto.scale[2])}))
            min <- min(sapply(list.granges, function(granges) {quantile(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), auto.scale[1])}))
            YLIM <- c(min, max)
        } else if (scale.max) {
            max <- max(sapply(list.granges, function(granges) {plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[2]}))
            min <- min(sapply(list.granges, function(granges) {plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[1]}))
            YLIM <- c(min, max)
        } else if (!is.null(YLIM)) {
            YLIM <- YLIM
        } else {
            stop('Please provide some valid Y limits. Aborting.')
        }
        # Initiate plot
        if (is.null(XLIM)) XLIM <- c(1, (unique(width(list.granges[[1]])) - BIN))
        plot(
            NA, 
            bty = 'l',
            xlim = XLIM,
            ylim = YLIM,
            xlab = ifelse(length(XLAB) > 1, 'Center of REs', XLAB), 
            ylab = YLAB, 
            axes = F, 
            ...
        )
        box(bty='l')
        axis(2)
        axis(1, at = seq(XLIM[1], XLIM[2], length.out = 3), labels = c(-unique(width(list.granges[[1]]) / 2), '0', unique(width(list.granges[[1]]) / 2)))
        # Add central line
        if (plot.central) abline(v = seq(XLIM[1], XLIM[2], length.out = 3)[2], lty = 3)
        # Loop through granges
        for (k.g in 1:length(list.granges)) {
            # Message
            granges <- list.granges[[k.g]]
            if (verbose) message('.. Plotting granges ', k.g, '/', length(list.granges), '...')
            
            # Plot
            mat <- get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F)
            par(new = T)
            plot.cov.line(mat, use.mean, BIN, plotEE, YLIM, xlim = c(1,(ncol(mat) - BIN)), axes = F, XLAB = '', YLAB = '', COL = list.COL[k.g])
            # End
        }
        # Plot legend
        if (plot.legend == T) {
            par(new = T)
            legend('topleft', legend = names(list.granges), fill = list.COL, bty = 'n')
        }
    
    } else if (length(list.granges) > 1 & length(list.bw.files) == 1 & class(list.granges[[1]]) == "CompressedGRangesList") { # Case where there is 1 bigWig file and several sets of several granges 
        
        # Define starting parameters
        list.granges.bak <- list.granges
        k.sets.g <- 0
        XLAB <- if (is.null(XLAB)) {'Segments coordinates'} else {XLAB}
        YLAB <- ifelse(is.null(YLAB), 'Score', YLAB)
        bw.file <- list.bw.files[[1]]
        # Define colors
        if (!is.null(list.COL)) {
            list.COL <- rep(list.COL, 5)
        } else {
            list.COL <- RColorBrewer::brewer.pal('Set2', n = 8)
        }
        # Loop through the different sets of GRanges
        for (k.sets.g in 1:length(list.granges)) { 
            list.granges <- list.granges.bak[[k.sets.g]]
            # Define YLIM 
            if (verbose) message('.. Determining Y boundaries...')
            if (!is.null(auto.scale)) {
                max <- max(sapply(list.granges, function(granges) {quantile(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), auto.scale[2])}))
                min <- min(sapply(list.granges, function(granges) {quantile(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), auto.scale[1])}))
                YLIM <- c(min, max)
            } else if (scale.max) {
                max <- max(sapply(list.granges, function(granges) {plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[2]}))
                min <- min(sapply(list.granges, function(granges) {plot.cov.line(get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F), return.extrema = T)[1]}))
                YLIM <- c(min, max)
            } else if (!is.null(YLIM)) {
                YLIM <- YLIM
            } else {
                stop('Please provide some valid Y limits. Aborting.')
            }
            # Initiate plot
            if (is.null(XLIM)) XLIM <- c(1, (unique(width(list.granges[[1]])) - BIN))
            plot(
                NA, 
                bty = 'l',
                xlim = XLIM,
                ylim = YLIM,
                xlab = ifelse(length(XLAB) == 1, XLAB, XLAB[k.b]), 
                ylab = YLAB, 
                axes = F, 
                ...
            )
            box(bty='l')
            axis(2)
            axis(1, at = seq(XLIM[1], XLIM[2], length.out = 3), labels = c(-unique(width(list.granges[[1]]) / 2), '0', unique(width(list.granges[[1]]) / 2)))
            # Add central line
            if (plot.central) abline(v = seq(XLIM[1], XLIM[2], length.out = 3)[2], lty = 3)
            # Loop through bw.files
            for (k.g in 1:length(list.granges)) {
                # Message
                granges <- list.granges[[k.g]]
                if (verbose) message('.. Plotting granges ', k.g, '/', length(list.granges), '...')
                
                # Plot
                mat <- get.cov.matrix(granges, bw.file, center = center, flank = flank, stranded = stranded, split.bid.proms = split.bid.proms, verbose = F)
                par(new = T)
                plot.cov.line(mat, use.mean, BIN, plotEE, YLIM, xlim = c(1,(ncol(mat) - BIN)), axes = F, XLAB = '', YLAB = '', COL = list.COL[k.g])
                # End
            }
        }
        
    }
    # End of loops
    message('All graphs successfully plotted')
}
# Add a seq column to granges 
with.seq <- function(granges, GENOME_FASTA = Rsamtools::FaFile("~/shared/sequences/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"), revcomp = F) {
    if (revcomp == T) {
        granges$seq <- Biostrings::reverseComplement(Biostrings::getSeq(GENOME_FASTA, granges))
    } 
    else {
        granges$seq <- Biostrings::getSeq(GENOME_FASTA, granges)
    }
    return(granges)
} 
extend.both.sides <- function(granges, width.bef, width.aft) {
    
    granges.bak <- granges
    granges <- resize(granges, width = width(granges.bak) + width.aft, fix = 'start')
    granges <- resize(granges, width = width.bef + width(granges.bak) + width.aft, fix = 'end')
    
    return(granges)
    
}

# ---------------------------------------------------------------------------- #
# ---------------- NUCLEOTIDE PERIODICITY FUNCTIONS
# ---------------------------------------------------------------------------- #
get.periodicity_ <- function(granges, DINUC = 'TT', doAlign = T, width.before = 0, width.after = 0, RANGE.FOR.SPECTRUM = 0:200, FREQ = 0.10, intervals = seq(2, 20, 1), doPlots = T, pdf = 'tmp.pdf') {
    
    #message('---- By default, this function automatically aligns the given promoters at their TSS (-50+250) and returns the distances between TTs and the associated FFT spectrum ----')
    if (doAlign) {
        granges <- with.seq(extend.both.sides(tss.align(granges = deconvolve.bidirectional.proms(granges), WIDTH.bef = 0, WIDTH.aft = 1), width.before, width.after))
    } else {
        granges <- with.seq(extend.both.sides(deconvolve.bidirectional.proms(granges), width.before, width.after))
    }
    dists <- get.pairwise.dists_(granges, MOTIF = DINUC, MERGE_BLOCKS = F, cores = 5)
    spectra <- get.spectra_(dists, RANGE.FOR.SPECTRUM, FREQ, return.entire.vec = T, plot = doPlots)
    mtm <- spectra$spec[unlist(lapply(c(1/(2:20)), function (FREQ) {idx <- which.min(abs(FREQ - spectra$freq)) }))]
    names(mtm) <- as.character((2:20))
    ratios <- get.ratios_(spectra, FREQ, intervals)
    if (doPlots) {
        pdf(pdf)
        plot.dists(list("Set of Proms" = dists), type = 'l', ylab = 'Counts', xlab = paste0('Distance between ', DINUC))
        plot.spectra(list("Set of Proms" = spectra))
        plot(x = 2:20, y = ratios, col = color.tissues[which(order.tissues == TISSUE)], lwd = 3, type = 'l')
        dev.off()
    }
    
    invisible(list(
        "dists" = dists, 
        "spectra" = spectra, 
        "PSD" = mtm,
        "PSD@10" = spectra$spec[which.min(abs(FREQ - spectra$freq))],
        "SNR" = ratios
    ))
    
}
get.pairwise.dists_ <- function(granges, MOTIF, cores = 4, N_MISMATCH = 0, MERGE_BLOCKS = F) {
    
    seqs <- granges$seq
    dists <- unlist(mclapply(1:length(seqs), function(k) {
        seq <- seqs[k]
        if (MERGE_BLOCKS == T) {
            c(dist(start(reduce(vmatchPattern(MOTIF, seq, max.mismatch = N_MISMATCH, fixed = F)[[1]]))))
        } else {
            c(dist(start(vmatchPattern(MOTIF, seq, max.mismatch = N_MISMATCH, fixed = F))[[1]]))
        }
    }, mc.cores = cores))
    return(dists)

}
get.spectra_ <- function(dists, RANGE.FOR.SPECTRUM = 1:100, FREQ = 0.10, plot = F, return.entire.vec = T) {
    require(multitaper)
    if(max(dists) < tail(RANGE.FOR.SPECTRUM, 1)) {
        message('Range (', head(RANGE.FOR.SPECTRUM, 1), ':', tail(RANGE.FOR.SPECTRUM, 1), ') is wider than the current vector length. Shortening the range from ', head(RANGE.FOR.SPECTRUM, 1), ' to ', max(dists), '...')
        RANGE.FOR.SPECTRUM <- head(RANGE.FOR.SPECTRUM, 1):max(dists)
    }
    if (plot) {
        spectra <- spec.mtm(hist(dists, breaks = seq(1, max(dists)+1, 1), plot = F)$counts[RANGE.FOR.SPECTRUM], k = 10, nw = 5.0, nFFT = "default", Ftest = TRUE, jackknife = FALSE, maxAdaptiveIterations = 100, plot = TRUE, na.action = na.fail)
    } else {
        spectra <- spec.mtm(hist(dists, breaks = seq(1, max(dists)+1, 1), plot = F)$counts[RANGE.FOR.SPECTRUM], k = 10, nw = 5.0, nFFT = "default", Ftest = TRUE, jackknife = FALSE, maxAdaptiveIterations = 100, plot = F, na.action = na.fail)
    }
    if (return.entire.vec) {
        return(spectra)
    } else {
        return(spectra$spec[which.min(abs(spectra$freq - FREQ))])
    }
}
get.ratios_ <- function(spectra, FREQ = FREQ, intervals = seq(2, 20, 1), log2 = F) {
    mtm <- spectra$spec[unlist(lapply(c(1/(intervals)), function (FREQ) {idx <- which.min(abs(FREQ - spectra$freq)) }))]
    if (log2) {
        ratios <- data.frame(
            freq = intervals, 
            strength = sapply(1:length(mtm), function(x) {log2(mtm[x]/mean(mtm[-x]))})
        )
    } else {
        ratios <- data.frame(
            freq = intervals, 
            strength = sapply(1:length(mtm), function(x) {mtm[x]/mean(mtm[-x])})
        )
    }
    return(ratios)
}

message(' >>>>>  Custom functions successfuly loaded. <<<<<<')



## Functions to switch from WB gene to gene_name

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

#########################################################
#########################################################
#########################################################


# PLOTTING FUNCTIONS

#1. To get a heatmap of tissue expression, for a vector of genes

getHeatmapOfCaoExpr <- function(vec, name.main = "", which.column = NULL) {
    # require(RColorBrewer)
    # require(gplots)
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
    # require(RColorBrewer)
    # require(vioplotx)
    #if (sum(ls() == "cao03") == 0) {load(".test-plotting.RData")}
    cao03b <- cao03[,order.tissues[1:5]]
    if (all(grepl("WB", vec))) {WBID <- vec} else {WBID <- name2WB(vec)}
    WBID <- WBID[WBID %in% row.names(cao03)]
    CAO.mat <- log2( as.matrix(cao03b[WBID,] +1) )
    CAO.mat[is.na(CAO.mat)] <- 0
    if (all(grepl("WB", vec))) {row.names(CAO.mat)=WB2name(WBID)} else {row.names(CAO.mat)=vec}
    vioplotx(CAO.mat[,1], CAO.mat[,2], CAO.mat[,3], CAO.mat[,4], CAO.mat[,5], col = color, main = name.main, ylab = "log2(TPM)")
}

#2. Function to plot scale of heatmaps plotted using image() function ; the scale uses the breaks and col values used for image()

image.scale <- function(z, zlim, col = heat.colors(12), breaks, axis.pos=1, add.axis=TRUE, ...) {
     if(!missing(breaks)){
      if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
     }
     if(missing(breaks) & !missing(zlim)){
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
     }
     if(missing(breaks) & missing(zlim)){
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

#3. On top of a plot already existing (typically boxplot(x)), add median lines and conf. interval
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

#4. Plot enriched GO from a vector of genes names
plotGOs <- function(x) {

    # require(gProfileR)
    reduced_c <- gprofiler(x, organism="celegans", max_p_value=0.05, correction_method="bonferroni", hier_filtering='none')
    reduced_c <- reduced_c[order(reduced_c$p.value),]
    d <- reduced_c[reduced_c$domain %in% c("MF", "BP", "CC", "keg"),]
    d <- d[!duplicated(d[, c("p.value", "overlap.size")],),]
    d <- d[1:min(25, dim(d)[1]),]
    if (nrow(d) > 0) {
        par(las = 0, mar = c(4,4,4,1))
        plot <- barplot(height=-log10(d$p.value), col=colorGO[d$domain], horiz=T, xlim=c(0,25))
        legend("right", legend=names(colorGO), fill=colorGO, col="#00000000", pch=15, bty="n")
        text(x=rep(0.2, times=length(plot)), y=plot, d$term.name, pos=4)
        mtext(text="GO-terms",side=2,line=1,outer=FALSE)
        mtext(text="-log10(adj-p.val)",side=1,line=3,outer=FALSE, cex=0.85)
    } else { plot.new() }

}

#########################################################
#########################################################
#########################################################

# UTILITARIES FUNCTION

#1. Function to get the size of all objects in the environment
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

#2. Function to increment 1
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))


#3. Improved list of objects
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




#4. Get expression values from cao supp. table 03
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

# Clean-up folder after paired-end RNA-seq mapping (lcap_mapping.R)
cleanup_folder <- function() {
    file.remove(list.files(pattern = '*.bg'))
    file.remove(list.files(pattern = '*.bai'))
    file.remove(list.files(pattern = '*.sai'))
    file.remove(list.files(pattern = '*tmp*'))
    ifelse(!dir.exists('_raw_data/'), dir.create('_raw_data/', showWarnings = FALSE), F)
    ifelse(!dir.exists('_bw-files/'), dir.create('_bw-files/', showWarnings = FALSE), F)
    ifelse(!dir.exists('_bam-files/'), dir.create('_bam-files/', showWarnings = FALSE), F)
    ifelse(!dir.exists(paste0('_bam-files/', sample, '/')), dir.create(paste0('_bam-files/', sample, '/'), showWarnings = FALSE), F)
    ifelse(!dir.exists('_log-files/'), dir.create('_log-files/', showWarnings = FALSE), F)
    ifelse(length(list.files(pattern = '.+.bamstats')), file.rename(list.files(pattern = '.+.bamstats'), paste0('_log-files/', list.files(pattern = '.+.bamstats'))), F)
    ifelse(length(list.files(pattern = '.+.fq')), file.rename(list.files(pattern = '.+.fq'), paste0('_raw_data/', list.files(pattern = '.+.fq'))), F)
    ifelse(length(list.files(pattern = '.+.bam')), file.rename(list.files(pattern = '.+.bam'), paste0('_bam-files/', sample, '/', list.files(pattern = '.+.bam'))), F)
    ifelse(length(list.files(pattern = '.+.bw')), file.rename(list.files(pattern = '.+.bw'), paste0('_bw-files/', list.files(pattern = '.+.bw'))), F)
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
        LCAP.GENE <- LCAP[WBID,]
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
        FCs <- round(max.tissue.df.LCAP.GENE[12:16],2)
        colnames(FCs) <- gsub('ratio.|.v.others', '', colnames(FCs))
        # Get REs associated with gene (which ones, associated tissues and dev. dynamics)
        REs <- all[grep(WBID, all$WormBaseID),]
        if (nrow(REs) == 0) {message('!!! WARNING !!! No associated REs are found !!!') ; REs.coords <- NULL ; REs.tau <- 0 ; REs.tissues <- NULL ; REs.TFs <- NULL ; REs.TFs.nb <- NULL ; ATAC.GENE <- 0 ; ATACdev.GENE <- 0 ; l2mean.ATAC.GENE <- 0 ; l2mean.ATACdev.GENE <- 0 } else {
            REs.coords <- REs[,c(1:3, 7)]
            row.names(REs.coords) <- paste0('RE_', seq(1:nrow(REs.coords)))
            REs.tau <- REs$TauScore
            REs.tissues <- order.tissues[REs$clustersATAC.tissues]
            REs.TFs <- SUMM[grep(WBID, all$WormBaseID),]
            REs.TFs.nb <- REs.TFs$TFnb
            ATAC.GENE <- ATAC[grep(WBID, all$WormBaseID),]
            ATACdev.GENE <- ATACdev[grep(WBID, all$WormBaseID),]
            l2mean.ATAC.GENE <- round(log2((ATAC.GENE+1)/rowMeans(ATAC.GENE+1)),2)
            row.names(l2mean.ATAC.GENE) <- row.names(REs.coords)
            colnames(l2mean.ATAC.GENE) <- colnames(l2mean.LCAP.GENE)
            l2mean.ATACdev.GENE <- round(log2((ATACdev.GENE+1)/rowMeans(ATACdev.GENE+1)),2)
            row.names(l2mean.ATACdev.GENE) <- row.names(REs.coords)
            colnames(l2mean.ATACdev.GENE) <- colnames(l2mean.LCAPdev.GENE)
        }

        # Build results object
        res <- list()
        res[['Gene.info']] <- c(WormBaseID = WBID, locus = locusID, Biotype = biotype, Tau.score = tau, Enriched.tissue = tissue)
        res[['Gene.expr.TPM']] <- round(rbind(Cao = CAO.GENE, tiss.spe.LCAP = LCAP.GENE),2)
        res[['Gene.expr.l2mean']] <- round(rbind(Cao = l2mean.CAO.GENE, tiss.spe.LCAP = l2mean.LCAP.GENE),2)
        res[['Gene.expr.dev.TPM']] <- round(rbind(LCAPdev = LCAPdev.GENE),2)
        res[['Gene.expr.dev.l2mean']] <- round(rbind(LCAPdev = l2mean.LCAPdev.GENE),2)
        res[['Associated.REs']] <- cbind(REs.coords, RE.Tau.score = round(REs.tau,2), REs.associated.tissues = REs.tissues)
        res[['Associated.REs.dev']] <- round(ATACdev.GENE,2)
        res[['Associated.REs.dev.l2mean']] <- l2mean.ATACdev.GENE
        res[['Associated.REs.tissue']] <- round(ATAC.GENE,2)
        res[['Associated.REs.tissue.l2mean']] <- l2mean.ATAC.GENE
        res[['Associated.REs.TFs.infos']] <- cbind(REs.coords, RE.Tau.score = round(REs.tau,2), RE.associated.tissues = REs.tissues)

        # Return results
        if (verbose == T) {
            width.bak <- options()$width
            options("width" = 200)
            cat('\n', rep('=', times = 31+nchar(locusID)), '\n>>>   ', WBID, ' --- ', locusID, '   <<<\n', rep('=', times = 31+nchar(locusID)), sep = "", fill = T)

            cat('\nGene infos\n', rep('-', times = 15+nchar(locusID)), sep = "", fill = T)
                cat(paste0(capture.output(res[['Gene.info']]), collapse = "\n"), sep = "", fill = T)
            cat('\n\nGene expression\n', rep('-', times = 15+nchar(locusID)), sep = "", fill = T)
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
            cat('\n\nAssociated regulatory elements\n', rep('-', times = 15+nchar(locusID)), sep = "", fill = T)
                cat('\n>> REs coordinates', sep = "", fill = T)
                    cat(paste0(capture.output(res[['Associated.REs']]), collapse = "\n"), sep = "", fill = T)
                cat('\n>> REs log2 mean-centered accessibility (tissue-specific)', sep = "", fill = T)
                    cat(paste0(capture.output(res[['Associated.REs.tissue']]), collapse = "\n"), sep = "", fill = T)
                cat('\n>> REs log2 mean-centered accessibility (dev)', sep = "", fill = T)
                    cat(paste0(capture.output(res[['Associated.REs.dev']]), collapse = "\n"), sep = "", fill = T)

            cat('\n', rep('=', times = 31+nchar(locusID)), '\n', sep = "", fill = T)
            options("width" = width.bak)
        }

        if (saveTXT != F) {
            width.bak <- options()$width
            options("width" = 200)

            txt.file <- ifelse(is.logical(saveTXT), file(paste0(locusID, '-summary.txt'), open="wt"), as.character(saveTXT))
            cat('\n', rep('=', times = 31+nchar(locusID)), '\n>>>   ', WBID, ' --- ', locusID, '   <<<\n', rep('=', times = 31+nchar(locusID)), sep = "", fill = T, file = txt.file)

            cat('\nGene infos\n', rep('-', times = 15+nchar(locusID)), sep = "", fill = T, file = txt.file, append = T)
                cat(paste0(capture.output(res[['Gene.info']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
            cat('\n\nGene expression\n', rep('-', times = 15+nchar(locusID)), sep = "", fill = T, file = txt.file, append = T)
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
            cat('\n\nAssociated regulatory elements\n', rep('-', times = 15+nchar(locusID)), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> REs coordinates', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(res[['Associated.REs']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> REs log2 mean-centered accessibility (tissue-specific)', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(res[['Associated.REs.tissue']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)
                cat('\n>> REs log2 mean-centered accessibility (dev)', sep = "", fill = T, file = txt.file, append = T)
                    cat(paste0(capture.output(res[['Associated.REs.dev']]), collapse = "\n"), sep = "", fill = T, file = txt.file, append = T)

            cat('\n', rep('=', times = 31+nchar(locusID)), '\n', sep = "", fill = T, file = txt.file, append = T)

            options("width" = width.bak)

        }

    }

    if (length(GENES) == 1 & exportResult == T) {
        return(res)
    }
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

# Function to load libraries to get genes infos
#geneInfos <- function(GENE) {
#
#    celegans_mart <- useMart(host = "metazoa.ensembl.org", biomart = "metazoa_mart", dataset = "celegans_eg_gene")
#    getBM(mart = celegans_mart, filters = "ensembl_gene_id", values = name2WB('myo-3'), attributes = listAttributes(celegans_mart)[,1][grep('description', listAttributes(celegans_mart)[,2])])
#
#    ensembl = useDataset("celegans_gene_ensembl",mart=ensembl)
#
#}


# Re-defining boxplot function
boxplot.2 <- function(..., whisklty = 'dotted', whiskcol = 'grey30', whiskldy = 2, staplecol = 'white', outline = F, notch = T, frame = F, medlty = 0, medpch = 21) {
    boxplot(..., whisklty = whisklty, whiskcol = whiskcol, whiskldy = whiskldy, staplecol = staplecol, outline = outline, notch = notch, frame = frame, medlty = medlty, medpch = medpch)
}

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

# Make a .bak object of an existing object, i.e. `XX.bak <- XX`
bak <- function(X) { 
    assign(paste0(deparse(substitute(X)), ".bak"), envir = parent.frame(), X)
    warning("Use this function with care")
    message(paste0(deparse(substitute(X)), " object copied to ", deparse(substitute(X)), ".bak"))
}
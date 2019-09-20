# Plot enriched GO from a vector of genes names

plotGOs <- function(x, hier.filtering='none', CEX = 0.4, ...) { 

    require(gProfileR)
    reduced_c <- gprofiler(x, organism="celegans", max_p_value=0.05, correction_method="bonferroni", hier_filtering=hier.filtering)
    reduced_c <- reduced_c[order(reduced_c$p.value),]
    d <- reduced_c[reduced_c$domain %in% c("MF", "BP", "CC", "keg"),]
    d <- d[!duplicated(d[, c("p.value", "overlap.size")],),]
    d <- d[1:min(25, dim(d)[1]),]
    if (nrow(d) > 0) {
        par(las = 0, mar = c(4,4,4,1))
        plot <- barplot(height=-log10(d$p.value), col=colorGO[d$domain], horiz=T, xlim=c(0,25), ...)
        legend("right", legend=names(colorGO), fill=colorGO, col="#00000000", pch=15, bty="n")
        text(x=rep(0.2, times=length(plot)), y=plot, d$term.name, pos=4, cex = CEX)
        mtext(text="GO-terms",side=2,line=1,outer=FALSE)
        mtext(text="-log10(adj-p.val)",side=1,line=3,outer=FALSE)
    } else { plot.new() }
    
    invisible(d)
    
}

plot.gos <- function(
    GENES, 
    max.p.val = 0.05, 
    correction.method = "bonferroni", 
    hier.filtering = "none", 
    onts = c("MF", "BP", "CC", "keg"), 
    max.depth = NULL, 
    max.terms = 25,
    plot = T, 
    size.by = 'p.val',
    symbol.by = 'ont', 
    color.by = 'ont' ) {
    require(gProfileR)
    
    if (is.vector(GENES) & !is.list(GENES)) {

        genes <- GENES
        gos <- gprofiler(genes, organism = "celegans", max_p_value = max.p.val, correction_method = correction.method, hier_filtering = hier.filtering)
        gos <- gos[order(gos$p.value),]
        gos.filtered <- gos[gos$domain %in% onts,]
        gos.filtered <- gos.filtered[!duplicated(gos.filtered[, c("p.value", "overlap.size")],),]
        if (!is.null(max.depth)) {gos.filtered <- gos.filtered[gos.filtered$relative.depth <= max.depth,]}
        gos.filtered <- gos.filtered[1:min(max.terms, dim(gos.filtered)[1]),]
        
        if (nrow(gos.filtered) > 0 & plot == T) {
            gos.filtered$genes.ratios <- gos.filtered$overlap.size / gos.filtered$query.size
            gos.filtered$logpval <- -log10(gos.filtered$p.value)
            gos.filtered$counts <- gos.filtered$overlap.size
            gos.filtered$symbol <- as.numeric(c('15', '16', '17', '18')[factor(gos.filtered$domain, levels = c('MF', 'BP', 'CC', 'keg'))])
            gos.filtered$color <- colorRampPalette(c('blue', 'red'))(10)[cut(gos.filtered$counts, breaks = seq(floor(min(gos.filtered$counts)), ceiling(max(gos.filtered$counts)), length.out = 11), include.lowest = T)]
            gos.filtered$size <- as.numeric(cut(gos.filtered$genes.ratios, breaks = seq(floor(min(gos.filtered$genes.ratios/0.1)), ceiling(max(gos.filtered$genes.ratios/0.1)), 1) / 10, include.lowest = T)) / 10 * 3 + 1
            gos.filtered$plotted.val <- gos.filtered$logpval
            gos.filtered.reordered <- gos.filtered[order(gos.filtered$plotted.val, decreasing = F),]
            
            layout(matrix(1:2, ncol = 2), widths = c(6, 2))
            par(mar = c(5, 1, 1, 1))
            dotchart(
                gos.filtered.reordered$plotted.val, 
                labels = gos.filtered.reordered$term.name,
                pt.cex = gos.filtered.reordered$size, 
                pch = gos.filtered.reordered$symbol, 
                color = gos.filtered.reordered$color, 
                lcolor = NA,
                xlab = '-log10(adj. p-val.)', 
                ylab = ''
            )
            plot.new()
            par(oma = c(0,0,0,0), mar = c(0,0,0,0))
            # Symbols (GO domain)
            legend('topleft', title = "GO db", legend = c('MF', 'BP', 'CC', 'keg'), pch = c(15, 16, 17, 18), bty = 'n')
            # Size (% of GO term)
            legend(x = 0, y = 0.8, title = "% of GO term", legend = as.character(seq(floor(min(gos.filtered$genes.ratios/0.1)), ceiling(max(gos.filtered$genes.ratios/0.1)), 1) / 10), pch = 19, pt.cex = (seq(floor(min(gos.filtered$genes.ratios/0.1)), ceiling(max(gos.filtered$genes.ratios/0.1)), 1) / 10 * 3 + 1), bty = 'n')
            # Color (nb of genes)
            legend(x = 0, y = 0.5, title = "# of genes", legend = as.character(round(seq(floor(min(gos.filtered$counts)), ceiling(max(gos.filtered$counts)), length.out = 10))), col = colorRampPalette(c('blue', 'red'))(10), pch = 15, bty = 'n')
        }
        
        return(gos.filtered)
    
    }
    
    if (is.list(GENES)) {
        list.gos <- list()
        
        for (K in 1:length(GENES)) {
            genes <- GENES[[K]]
            gos <- gprofiler(genes, organism="celegans", max_p_value = max.p.val, correction_method = correction.method, hier_filtering= hier.filtering)
            gos <- gos[order(gos$p.value),]
            gos.filtered <- gos[gos$domain %in% onts,]
            gos.filtered <- gos.filtered[!duplicated(gos.filtered[, c("p.value", "overlap.size")],),]
            if (!is.null(max.depth)) {gos.filtered <- gos.filtered[gos.filtered$relative.depth <= max.depth,]}
            gos.filtered <- gos.filtered[1:min(max.terms, dim(gos.filtered)[1]),]
            list.gos[[K]] <- gos.filtered[, c("term.id", "domain", "term.size", "query.size", "overlap.size", "p.value", "recall", "precision", "term.name")]
        }
        
        return(list.gos)
    }
    
}

homemadeGOplot <- function(GENES_LIST, DICT = 'anatomy', pdf = NULL) {
    DICT_FILE <- if (DICT == 'anatomy') {
        "~/Downloads/anatomy_dict.csv"
    }

    DICT <- read.csv(DICT_FILE, header = T, stringsAsFactors = F, row.names = 1)

    #GENES_LIST <- readLines('germline.genes.list')
    DICT_SUB <- DICT[row.names(DICT) %in% GENES_LIST,]
    DICT_SUB <- DICT_SUB[,colSums(DICT_SUB) > 10]
    GO_LIST_SUB <- colnames(DICT_SUB)

    P_VALS <- c()
    ODD_VALS <- c()
    PCT_VALS <- c()
    for (GO in GO_LIST_SUB) {
        TEST <- fisher.test(rbind( c( sum(DICT_SUB[,GO]), sum(DICT[,GO]) - sum(DICT_SUB[,GO]) ), c( (sum(DICT_SUB) - sum(DICT_SUB[,GO])), (sum(DICT) - sum(DICT_SUB) - sum(DICT[,GO]) + sum(DICT_SUB[,GO])) )))
        P_VALS[GO] <- TEST$p.val
        ODD_VALS[GO] <- TEST$estimate
        PCT_VALS[GO] <- round(sum(DICT_SUB[,GO]) / (sum(DICT[,GO]) - sum(DICT_SUB[,GO])) * 100, 2)
    }
    
    # Filter significant GO terms
    FILTER <- which(p.adjust(P_VALS, method = "bonf") < 0.001 & ODD_VALS > 1 & PCT_VALS > 5)
    # Get the order of p-values
    P_VALS_ORDER <- rev(order(p.adjust(P_VALS, method = "bonf")[FILTER], decreasing = F))
    # Order the rest of the values and adjust p-values
    GOs_SORTED <- gsub("\\.", " ", gsub(".WBbt.*", "", names(P_VALS[FILTER][P_VALS_ORDER])))
    P_VALS_SORTED <- p.adjust(P_VALS, method = "bonf")[FILTER][P_VALS_ORDER]
    ODD_VALS_SORTED <- ODD_VALS[FILTER][P_VALS_ORDER]
    PCT_VALS_SORTED <- PCT_VALS[FILTER][P_VALS_ORDER]

    # Proceed to plotting results
    if (length(P_VALS_SORTED) > 0) {
        if(!is.null(pdf)) {
            pdf(pdf)
        }
        par(las = 0, mar = c(4,10,4,1))
        plot <- barplot(height = -log10(P_VALS_SORTED), horiz = T, names = "")
        grid(col = "black", ny = 1, lwd = 0.2, lty = 1)
        plot <- barplot(height = -log10(P_VALS_SORTED), horiz = T, names = "", add = T)
        par(xpd = T)
        text(x=rep(-0.2, times = length(P_VALS_SORTED)), y = plot, GOs_SORTED, pos = 2, las = 0)
        par(xpd = F)
        mtext(text = "-log10(adj-p.val)", side = 1, line = 2, outer = FALSE, cex = 0.85)
        if(!is.null(pdf)) {
            dev.off()
        }
    } else {  }
}

wrapper.clusterProfiler <- function(genes) {
    require(AnnotationHub)
    require(clusterProfiler)
    require(org.Ce.eg.db)
    hub <- AnnotationHub()
    Celegans <- hub[["AH52249"]]
    ego <- enrichGO(gene = bitr(genes, fromType = "WORMBASE", toType = c("ENTREZID"), OrgDb = org.Ce.eg.db)$ENTREZID,
        OrgDb            = org.Ce.eg.db,
        ont              = "CC",
        pAdjustMethod    = "BH",
        pvalueCutoff     = 0.01,
        qvalueCutoff     = 0.05
    )
    message("Run clusterProfiler::dotplot(x, showCategory = NULL) + ggplot2::theme_bw(base_size = 10) to see results")
    return(ego)
}

wrapper.compareCluster <- function(list) {
    require(AnnotationHub)
    require(clusterProfiler)
    require(org.Ce.eg.db)
    hub <- AnnotationHub()
    Celegans <- hub[["AH52249"]]
    genes <- lapply(list, function(genes) {
        bitr(genes, fromType = "WORMBASE", toType = c("ENTREZID"), OrgDb = org.Ce.eg.db)$ENTREZID
    })
    names(genes) <- names(list)
    df <- data.frame(genes = unlist(genes), group = lapply(names(genes), function(group) {rep(group, length(genes[[group]]))}) %>% unlist())
    res <- compareCluster(
        genes ~ group, 
        data = df, 
        fun = "enrichGO", 
        OrgDb            = org.Ce.eg.db,
        ont              = "CC",
        pAdjustMethod    = "BH",
        pvalueCutoff     = 0.01,
        qvalueCutoff     = 0.05
    )
    message("Run clusterProfiler::dotplot(x, showCategory = NULL) + ggplot2::theme_bw(base_size = 10) to see results")
    return(res)
}

wraper.gprofiler <- function(
    query = genes, 
    compare_genes = TRUE,
    organism = "celegans", 
    max_p_value = 0.05, 
    correction_method = "bonferroni", 
    sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG", "TF"), 
    plotly = TRUE
    ) {
    
    if (class(query) == "list" & compare_genes == TRUE) 
        multi_query <- TRUE
    results <- gprofiler2::gost(query = query, 
        organism = "celegans", 
        ordered_query = FALSE, 
        multi_query = compare_genes,
        significant = TRUE, 
        exclude_iea = TRUE, 
        measure_underrepresentation = FALSE, 
        evcodes = FALSE, 
        user_threshold = max_p_value, 
        correction_method = correction_method, 
        domain_scope = "annotated", 
        custom_bg = NULL, 
        numeric_ns = "", 
        sources = sources
    )
    plot <- gprofiler2::gostplot(
        results, 
        interactive = plotly
    )
    results <- list(
        "table" = results$result, 
        "plot" = plot
    )
    return(results)
}

fetchGOsDB <- function(species = c('hs', 'mm', 'dm', 'ce'), update = FALSE) {
    species.file <- switch(species,
        'hs' = "goa_human.gaf.gz",
        'mm' = "mgi.gaf.gz",
        'dm' = "fb.gaf.gz",
        'ce' = "wb.gaf.gz"
    )
    annot.file <- paste0("~/shared/annots/", species, ".GO.annots.gaf")
    if (update | !file.exists(annot.file)) {
        system(paste0("wget http://current.geneontology.org/annotations/", species.file))
        system(paste0("mv ", species.file, " ~/shared/annots/", species.file))
        system(paste0("gunzip ~/shared/annots/", species.file))
        system(paste0("mv ~/shared/annots/", gsub('.gz', '', species.file), " ", annot.file))
    }
    return(annot.file)
}

getGenesWithGO <- function(
    gos,
    species = c('hs', 'mm', 'dm', 'ce'), 
    ev.codes = c('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP', 'IBA', 'IBD', 'IKR', 'IRD', 'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA', 'TAS', 'NAS', 'IC', 'ND', 'IEA'),
    exclude.electronic.annots = TRUE, 
    update = FALSE
    ) {
    if (exclude.electronic.annots & any(grepl('IEA', ev.codes))) 
        ev.codes <- ev.codes[ev.codes != 'IEA']
    GOs.file <- fetchGOsDB(species, update)
    lines <- suppressWarnings(readLines(GOs.file)) %>% 
        grep(paste(paste0("\t", ev.codes, "\t"), collapse = '|'), ., value = T)
    genes <- sapply(gos, function(GO) {
        lines %>% 
            grep(paste0("\t", GO, "\t"), ., value = T) %>% 
            strsplit('\t') %>% 
            sapply(., '[[', 2)
    }) %>% 
        unlist()
    res <- genes[table(genes) >= length(gos)] %>% unname()
    return(res)
}

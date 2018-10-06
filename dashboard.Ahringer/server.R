#=============================================================================#
#                                                                             #
#     USAGE: THIS IS NOT A FUNCTIONAL SCRIPT.                                 #
#     It is a part of a Shiny app                                             #
#                                                                             #
#     AUTHOR: Jacques SERIZAY                                                 #
#     CREATED: 2018/07/13                                                     #
#     REVISION: 2018/10/05                                                    #
#                                                                             #
#=============================================================================#

## Define server function ----------------------------------------------------------------------------------------------------

shinyServer <- function(input, output, session) {

  # Read gene name and get gene infos
  {
    gene <- reactive ({ input$searchGene })
    output$gene <- renderText({ gene() })
    infos.gene <- reactive ({ getGeneInfos(gene(), saveTXT = F, verbose = F, exportResult = T) })

    output$geneInfos <- renderUI({
        WBID <- paste("<b>WormBase ID:</b>   ", infos.gene()$Gene.info[1])
        LOCUS <- paste("<b>Locus:</b>   ", infos.gene()$Gene.info[2])
        COORDS <- paste("<b>Coordinates:</b>   ", as.character(base::range(genes.gtf[infos.gene()$Gene.info[1]])))
        BIOTYPE <- paste("<b>Gene biotype:</b>   ", infos.gene()$Gene.info[3])
        ENRICHED <- paste("<b>Enriched in tissue:</b>   ", infos.gene()$Gene.info[5])
        HTML(paste("", WBID, LOCUS, COORDS, BIOTYPE, ENRICHED, sep = '<br/>'))
    })
  }

  # Generate button to download gene-specific text file
  {
      output$downloadINFOS <- downloadHandler(
          reactive(paste0(infos.gene()$Gene.info[2], "_summary.txt")),
          content = function(file) {
              temp <- tempfile()
              getGeneInfos(gene(), saveTXT = temp, verbose = F, exportResult = F)
              writeLines(readLines(temp), file)
          }
      )
  }
  
  # Generate button to access to Genome Browser
  {
      observeEvent(input$switchToGenome, { updateTabItems(session, "tabs", "browser") })
  }

  # Generate browser
  {
    RE.coords <- reactive ({ c(
      as.character(infos.gene()$Associated.REs[1,1]),
      as.character(min(infos.gene()$Associated.REs[,2])-1000),
      as.character(max(infos.gene()$Associated.REs[,3])+1000)
    ) })
    gene.coords <- reactive ({ c(
      paste0('chr', as.character(seqnames(genes.gtf[infos.gene()$Gene.info[1]]))),
      as.character(start(genes.gtf[infos.gene()$Gene.info[1]])-3000),
      as.character(end(genes.gtf[infos.gene()$Gene.info[1]])+3000)
    ) })
    coords <- reactive ({ c(RE.coords()[1], min(RE.coords()[2], gene.coords()[2]), max(RE.coords()[3], gene.coords()[3])) })
    #-->
    url <- renderText({ getURL(coords()[1], coords()[2], coords()[3], "1.12.5") })
    output$url <- reactive ({ getURL(coords()[1], coords()[2], coords()[3], "1.12.5") })
    output$jbrowser <- renderJbrowse({ iframeJbrowse.2(url()) })
  }

  # Generate the REs subset table
  {
    #-->
    output$REs.table <- renderDataTable({
        if (all(infos.gene()$Associated.REs == 0)) { datatable(matrix(0)) } else {
            datatable(
                setNames(cbind(infos.gene()$Associated.REs[c(1:4, 6)], round(infos.gene()$Associated.REs.dev, 1), round(infos.gene()$Associated.REs.tissue, 1)),
                c("Chr", "Start", "Stop", "Regulatory Class", "Enriched in tissue(s)", "Coverage Emb. (mixed)", "Coverage L1 (mixed)", "Coverage L2 (mixed)", "Coverage L3 (mixed)", "Coverage L4 (mixed)", "Coverage YA (mixed)", "Coverage Hypod. (YA)", "Coverage Neurons (YA)", "Coverage Gonad (YA)", "Coverage Muscle (YA)", "Coverage Intest. (YA)")),
                autoHideNavigation = T,
                rownames = F,
                style = 'bootstrap',
                extensions = c('Buttons'),
                options = list(
                    pageLength = 10,
                    dom = 'Brtip',
                    buttons = list('copy', 'print', list(
                    extend = 'collection',
                    buttons = c('csv', 'excel', 'pdf'),
                    text = 'Download')),
                    scrollX = 200,
                    autoWidth = T,
                    searching = F
                )
            )
        }
    })
  }

  # Generate the expression profiles plots
  {
      #1: Plot expr. during development
      output$Expr.plots_dev <- renderPlot({
          par(mar=c(6,5,2,5))
          plot(
              as.numeric(infos.gene()$Gene.expr.dev.TPM),
              type = 'l', lty = 1, lwd = 3, xlab = "", ylab = "Expr. (TPM) over time", bty = 'n',
              main = "Stage-specific gene expression (mixed)", xaxt = 'n'
          )
          axis(1, labels = c('Emb.', 'L1', 'L2', 'L3', 'L4', 'YA'), at = 1:6)
      })
      #2: Plot expr. / tissue
      output$Expr.plots_tis <- renderPlot({
          par(mar=c(6,5,2,5))
          barplot(
              as.numeric(infos.gene()$Gene.expr.TPM[2,]),
              col = color.tissues, ylab = "Expr. / tissue (TPM)",
              names = order.tissues[1:5],
              main = "Tissue-specific gene expression (YA)", 
              las = 3
          )
      })

  }

  # Generate tSNE plots (proms and genes)
  {

    # Plot genes
    output$tSNE.plots_genes <- renderPlot({
        par(mar=c(2,4,2,2))
        plot(tSNE.df.LCAP,
            col = CLASSES.genes,
            pch = 20,
            cex = 1.5,
            xlab = '',
            ylab = '',
            axes = T,
            xaxt='n',
            yaxt='n',
            bty = 'o',
            main = 't-SNE projection of genes',
            cex.main = 1.5
        )
        points(
            tSNE.df.LCAP[genes.gtf[genes.valid]$gene_id %in% infos.gene()$Gene.info[1],'V1'],
            tSNE.df.LCAP[genes.gtf[genes.valid]$gene_id %in% infos.gene()$Gene.info[1],'V2'],
            col = 'black',
            pch=20,
            cex = 4
        )
        mtext(line = 1, side = 1, 'tSNE dim. 1')
        mtext(line = 1, side = 2, 'tSNE dim. 2')
    })
    
    # Plot promoters
    output$tSNE.plots_proms <- renderPlot({
        par(mar=c(2,4,2,2))
        plot(tSNE.df.proms,
            col = CLASSES.proms,
            pch = 20,
            cex = 1.5,
            xlab = '',
            ylab = '',
            axes = T,
            xaxt='n',
            yaxt='n',
            bty = 'o',
            main = 't-SNE projection of promoters',
            cex.main = 1.5
        )
        if (!all(infos.gene()$Associated.REs == 0)) {
            points(
                tSNE.df.proms[all.proms.valid$chr %in% infos.gene()$Associated.REs[,1] & all.proms.valid$start %in% infos.gene()$Associated.REs[,2] & all.proms.valid$stop %in% infos.gene()$Associated.REs[,3],'V1'],
                tSNE.df.proms[all.proms.valid$chr %in% infos.gene()$Associated.REs[,1] & all.proms.valid$start %in% infos.gene()$Associated.REs[,2] & all.proms.valid$stop %in% infos.gene()$Associated.REs[,3],'V2'],
                col = 'black',
                pch=20,
                cex = 4
            )
        }
        mtext(line = 1, side = 1, 'tSNE dim. 1')
        mtext(line = 1, side = 2, 'tSNE dim. 2')
    })
    
    # Plot legend
    output$tSNE.plots_legend <- renderPlot({
        plot.new()
        par(mar=c(0,0,0,0), oma = c(0,0,0,0), xpd = NA)
        legend("left",
              legend = c(order.tissues[1:5], '2 tissues', '3 tissues', '4 tissues', order.tissues[32:34]),
              fill = c(color.tissues[1:5], color.tissues[6], color.tissues[16], color.tissues[26], color.tissues[32:34]),
              ncol = 1, cex = 1.5, bty = 'n'
        )
    })

  }

  # Get multiple gene names
  {
    multipleGenes <- reactive ({ unlist(strsplit(x = input$searchMulitpleGenePaste, split = '[\r\n]' )) %>% ifelse(!grepl('WBGene', .), name2WB(.), .) %>% .[!is.na(.)]})
    output$multipleGenes <- renderText ({ multipleGenes() })
    output$multipleGenesLength <- reactive ({ paste(length(multipleGenes()), "valid gene(s) found.") })
  }

  # Get GOs plot (enriched in multiple genes)
  {

    output$GO.plot <- renderPlot({

        lay <- layout(matrix(seq(1:2), nrow = 1), widths = c(3,1))
        # Plot GOs
        reduced_c <- gprofiler(multipleGenes(), organism="celegans", max_p_value=0.05, correction_method="bonferroni", hier_filtering='moderate')
        reduced_c <- reduced_c[order(reduced_c$p.value),]
        d <- reduced_c[reduced_c$domain %in% c("MF", "BP", "CC", "keg"),]
        d <- d[!duplicated(d[, c("p.value", "overlap.size")],),]
        d <- d[1:min(20, dim(d)[1]),]
        if (nrow(d) > 0) {
          par(las = 0, mar = c(4,20,4,1))
          plot <- barplot(height=-log10(d$p.value), col=colorGO[d$domain], horiz=T, xlim=c(0,25))
          par(xpd=T)
          #text(x=rep(-1, times = length(plot)), y = plot, paste0(substr(d$term.name, start = 1, stop = 40), ifelse(nchar(d$term.name) > 40, '...', ''), ' (n=', d$overlap.size, ')'), pos = 2)
          text(x=rep(-1, times = length(plot)), y = plot, paste0(substr(d$term.name, start = 1, stop = 40), ifelse(nchar(d$term.name) > 40, '...', '')), pos = 2)
          par(xpd=F)
          legend("topright", legend=names(colorGO), fill=colorGO, col="#00000000", pch=15, bty="n")
          mtext(text="-log10(adj-p.val)",side=1,line=3,outer=FALSE, cex=0.85)
        } else { plot.new() }
    })

  }

  # Get HMs.plot (multiple genes heatmaps of LCAP/ATAC fold-changes)
  {
    output$HMs.plot <- renderPlot({
        mat.LCAP <- log2(as.matrix(max.tissue.df.LCAP[multipleGenes(), 12:16]))
        mat.ATAC <- as.matrix(all.deconv[all.deconv$uniqueWormBaseID %in% multipleGenes(), 15:19])
        mat.ATAC2 <- mat.ATAC
        for (NCOL in 1:ncol(mat.ATAC)) { mat.ATAC2[,NCOL] <- log2((mat.ATAC[,NCOL]+1)/(rowMeans(mat.ATAC[,-NCOL]+1))) }

        lay <- layout(matrix(seq(1:2), nrow = 1))
        # Plot LCAP
        hm <- heatmap_parmfrow(
            mat.LCAP,
            margins = c(6,1,4,1),
            col = colorRampPalette(c("royalblue", "white", "darkred"))(99),
            breaks = seq(-5, 5, length.out = 100),
            main = '',
            key.lab = "Tissue-specific gene expression enrichment (log2)",
            cexCol = 1,
            cexRow = 0.4,
            ylab = "",
            scaleUP = 1.05,
            scaleUP.bis = 1.04,
            nticks = 15,
            doClust = T,
            labRow = 'none',
            labCol = order.tissues[1:5]
        )
        # Plot ATAC
        hm <- heatmap_parmfrow(
            mat.ATAC2,
            margins = c(6,1,4,1),
            col = colorRampPalette(c("royalblue", "white", "darkred"))(99),
            breaks = seq(-5, 5, length.out = 100),
            main = '',
            key.lab = "Tissue-specific RE accessibility enrichment (log2)",
            cexCol = 1,
            cexRow = 0.4,
            ylab = "",
            scaleUP = 1.05,
            scaleUP.bis = 1.04,
            nticks = 15,
            doClust = T,
            labRow = 'none',
            labCol = order.tissues[1:5]
        )

    })

  }

  # Generate tables to download
  {
      output$downloadATAC.txt <- downloadHandler( "tissue-specific.ATAC-seq.dataset.txt", content = function(file) {write.table(atac.dt, file, quote = F, row = F, col = T, sep = '\t')} )
      output$downloadATAC.gff <- downloadHandler( "tissue-specific.ATAC-seq.dataset.gff", content = function(file) {
          infos=paste0("ID=", make.unique(paste(all$gene_name, all$regulatory_class, sep = '--'), sep = '_'), ";Associated-gene-Name=", all$gene_name, ";Associated-gene-WB=", all$WormBaseID, ";CV=", round(all$cv, 3), ";domain=", all$domain)
          values <- paste0(';color=', color.tissues[max.tissue.df$which.tissues], '; =: : : : : : : : : : : : : : : : : : : : : : : : : ', ';Ranked-tissues=', apply(max.tissue.df[,5:9], 1, function(x) paste0(x, collapse = ' * ')), ';Tissue-ATACeq-TPM=', apply(max.tissue.df[,10:14], 1, function(x) paste0(round(x, 3), collapse = ' / ')), ';Consecutive-ratios=', apply(max.tissue.df[,15:18], 1, function(x) paste0(round(x, 3), collapse = ' / ')), ";Enriched-tissue.s.=", gsub('\\.', '', max.tissue.df$which.tissues))
          GFF=cbind(as.character(all$chr), rep('Ahringer-JJ-JA', times = nrow(all)), as.character(all$regulatory_class), all$start, all$stop, rep(".", times=nrow(all)), ifelse(is.null(all$strand), ".", all$strand), rep(".", times=nrow(all)), paste0(infos, values))
          write.table(rbind(c("##gffTags=on\n##displayName=Name\n##gff-version 3", rep(" ", 8)), GFF), file, quote = F, row = F, col = F, sep = '\t')
      } )
      output$downloadLCAP.txt <- downloadHandler( "tissue-specific.RNA-seq.dataset.txt", content = function(file) {write.table(lcap.dt, file, quote = F, row = F, col = T, sep = '\t')}  )
      output$downloadLCAP.gff <- downloadHandler( "tissue-specific.RNA-seq.dataset.gff", content = function(file) {
          
            names(genes.gtf) <- genes.gtf$gene_name
            cv <- genes.gtf$cv

            infos=paste0(
                "ID=", genes.gtf$gene_id,
                ";Name=", genes.gtf$gene_name,
                ";CV=", cv,
                ';color=', color.tissues[c(1:5,33:35)][max.tissue.df.LCAP$which.tissues],
                '; =: : : : : : : : : : : : : : : : : : : : : : : : : ',
                ';Ranked-tissues=', apply(max.tissue.df.LCAP[,2:6], 1, function(x) paste0(x, collapse = ' / ')),
                ';Tissue-RNAseq-TPM=', apply(max.tissue.df.LCAP[,7:12], 1, function(x) paste0(round(x, 3), collapse = ' / ')),
                ';Hypod-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Hypod..v.others'],
                ';Neurons-FCvs-mean=', max.tissue.df.LCAP[,'ratio.Neurons.v.others'],
                ';Gonad-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Gonad.v.others'],
                ';Muscle-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Muscle.v.others'],
                ';Intest-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Intest..v.others'],
                ";Enriched-in=", max.tissue.df.LCAP$which.tissues
            )

            GFF=cbind(
                as.character(as.character(seqnames(genes.gtf))),
                rep('WormBase', times = length(genes.gtf)),
                as.character(genes.gtf$gene_biotype),
                start(genes.gtf),
                end(genes.gtf),
                rep(".", length(genes.gtf)),
                as.character(strand(genes.gtf)),
                rep(".", length(genes.gtf)),
                infos
            )
            
            write.table(rbind(c("##gffTags=on\n##displayName=Name\n##gff-version 3", rep(" ", 8)), GFF), file, quote = F, row = F, col = F, sep = '\t')
            
      }  )
  }
  {
      output$atac.table <- renderDataTable({
          datatable(atac.dt,
              autoHideNavigation = T,
              class = 'cell-border compact',
              rownames = F,
              style = 'bootstrap',
              filter = list(position = 'top', clear = FALSE),
              extensions = c('Buttons'),
              options = list(
                  pageLength = 200,
                  dom = 'Brtip',
                  buttons = list('copy', 'print', list(extend = 'collection', buttons = c('csv', 'excel', 'pdf'), text = 'Download visible data'), I('colvis')),
                  scrollX = 200,
                  scrollY = 200,
                  searching = T,
                  scroller = TRUE,
                  deferRender = TRUE
              )
          )
      })

      output$lcap.table <- renderDataTable({
          datatable(lcap.dt,
              autoHideNavigation = T,
              class = 'cell-border compact',
              rownames = F,
              style = 'bootstrap',
              filter = list(position = 'top', clear = FALSE),
              extensions = c('Buttons'),
              options = list(
                  pageLength = 200,
                  dom = 'Brtip',
                  buttons = list('copy', 'print', list(extend = 'collection', buttons = c('csv', 'excel', 'pdf'), text = 'Download visible data'), I('colvis')),
                  scrollX = 200,
                  scrollY = 200,
                  searching = T,
                  scroller = TRUE,
                  deferRender = TRUE
              )
          )
      })

  }

} #EOF

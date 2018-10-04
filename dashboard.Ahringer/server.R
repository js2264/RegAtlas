#=============================================================================#
#                                                                             #
#     USAGE: THIS IS NOT A FUNCTIONAL SCRIPT.                                 #
#     It is a part of a Shiny app                                             #
#                                                                             #
#     AUTHOR: Jacques SERIZAY                                                 #
#     CREATED: 2018/07/13                                                     #
#     REVISION: ..../../..                                                    #
#                                                                             #
#=============================================================================#

## Define server function ----------------------------------------------------------------------------------------------------

shinyServer <- function(input, output) {

  # Read gene name and get gene infos
  {
    gene <- reactive ({ input$searchGene })
    output$gene <- renderText({ gene() })
    infos.gene <- reactive ({ getGeneInfos(gene(), saveTXT = F, verbose = F, exportResult = T) })

    output$geneInfos <- renderUI({
        WBID <- paste("<b>WormBase ID:</b>   ", infos.gene()$Gene.info[1])
        LOCUS <- paste("<b>Locus ID:</b>   ", infos.gene()$Gene.info[2])
        BIOTYPE <- paste("<b>Gene biotype:</b>   ", infos.gene()$Gene.info[3])
        ENRICHED <- paste("<b>Enriched in tissue:</b>   ", infos.gene()$Gene.info[5])
        HTML(paste("", WBID, LOCUS, BIOTYPE, ENRICHED, sep = '<br/>'))
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
      output$Expr.plots <- renderPlot({
          lay <- layout(matrix(seq(1:2), nrow = 1))
          #1: Plot expr. during development
          par(mar=c(4,5,2,2))
          plot(
              as.numeric(infos.gene()$Gene.expr.dev.TPM),
              type = 'l', lty = 1, lwd = 3, xlab = "", ylab = "Expr. (TPM) over time", bty = 'n',
              main = "Stage-specific gene expression (mixed)", xaxt = 'n'
          )
          axis(1, labels = c('Emb.', 'L1', 'L2', 'L3', 'L4', 'YA'), at = 1:6)
          #2: Plot expr. / tissue
          par(mar=c(4,5,2,2))
          barplot(
              as.numeric(infos.gene()$Gene.expr.TPM[2,]),
              col = color.tissues, ylab = "Expr. / tissue (TPM)",
              names = order.tissues[1:5],
              main = "Tissue-specific gene expression (YA)"
          )
      })

  }

  # Generate tSNE plots (proms and genes)
  {
    #-->
    output$tSNE.plots <- renderPlot({
      lay <- layout(matrix(seq(1:3), nrow = 1), widths = c(0.5, 0.4, 0.5))
      par(oma = c(1, 5, 1, 5))
      # Plot genes
      par(mar=c(2,2,2,2))
      plot(tSNE.df.LCAP,
           col = CLASSES.genes,
           pch = 20,
           cex = 2,
           xlab = '',
           ylab = '',
           axes = T,
           xaxt='n',
           yaxt='n',
           bty = 'o',
           main = 'Genes tSNE (RNA-seq)',
           cex.main = 2
      )
      points(
        tSNE.df.LCAP[genes.gtf[genes.valid]$gene_id %in% infos.gene()$Gene.info[1],'V1'],
        tSNE.df.LCAP[genes.gtf[genes.valid]$gene_id %in% infos.gene()$Gene.info[1],'V2'],
        col = 'black',
        pch=20,
        cex = 4
      )
      #
      # Plot legend
      plot.new()
      par(mar=c(0,0,0,0))
      legend("center",
            legend = c(order.tissues[1:5], '2 tissues', '3 tissues', '4 tissues', order.tissues[32:34]),
            fill = c(color.tissues[1:5], color.tissues[6], color.tissues[16], color.tissues[26], color.tissues[32:34]),
            ncol = 1, cex = 2, bty = 'n'
      )
      #
      # Plot proms
      par(mar=c(2,2,2,2))
      plot(tSNE.df.proms,
           col = CLASSES.proms,
           pch = 20,
           cex = 2,
           xlab = '',
           ylab = '',
           axes = T,
           xaxt='n',
           yaxt='n',
           bty = 'o',
           main = 'Promoters tSNE (ATAC-seq)',
           cex.main = 2
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
      #
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
      output$downloadATAC <- downloadHandler( "tissue-specific.ATAC-seq.dataset.txt", content = function(file) {write.table(atac.dt, file, quote = F, row = F, col = T, sep = '\t')} )
      output$downloadLCAP <- downloadHandler( "tissue-specific.RNA-seq.dataset.txt", content = function(file) {write.table(lcap.dt, file, quote = F, row = F, col = T, sep = '\t')}  )
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

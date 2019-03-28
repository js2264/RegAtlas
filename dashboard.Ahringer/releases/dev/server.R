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
    
    output$geneDescr <- renderUI({ INFOS <- HTML(paste0("<h4>Description:</h4><br/>", fetchWBinfos(infos.gene()$Gene.info[2]))) })
    
  }

  # Generate buttons to download gene-specific text file, all tracks in zip, and genes list full report
  {
      output$downloadINFOS <- downloadHandler(
          reactive(paste0(infos.gene()$Gene.info[2], "_summary.txt")),
          content = function(file) {
              temp <- tempfile()
              getGeneInfos(gene(), saveTXT = temp, verbose = F, exportResult = F)
              writeLines(readLines(temp), file)
          }
      )
      
      output$downloadBWLCAP <- downloadHandler(
        "tissue-specific-RNA-seq_bigwig-tracks.zip",
        content <- function(file) {
            file.copy("tissue-specific-RNA-seq_bigwig-tracks.zip", file)
        }
      )
      
      output$downloadBWATAC <- downloadHandler(
        "tissue-specific-ATAC-seq_bigwig-tracks.zip",
        content <- function(file) {
            file.copy("tissue-specific-ATAC-seq_bigwig-tracks.zip", file)
        }
      )
      
      output$downloadGenesListGFF <- downloadHandler(
        "gene-list_full-report.gff",
        content <- function(file) {
            
            # Generate HEADER
            HEADER <- rbind(
                c("##gffTags=on\n##displayName=Name\n##gff-version 3", rep(" ", 8)),
                c(paste0("##FILE GENERATE ON: ", date(), ".\n##NUMBER OF GENES FOUND IN THE QUERY: ", length(multipleGenes())), rep(" ", 8)),
                c("##\n##This file contains both genes and associated regulatory elements annotations. Last column contains general information as well as details of tissue-specific expression / accessibility of the genes / associated REs.", rep(" ", 8)),
                c("##This file can easily be loaded in IGV. Coordinates are ce11.", rep(" ", 8))
            )
            
            # Generate genes GFF part
            infos <- paste0(
                "ID=", genes.gtf$gene_id,
                ";Name=", genes.gtf$gene_name,
                ';color=', color.tissues[c(1:5,33:35)][max.tissue.df.LCAP$which.tissues],
                '; =: : : : : : : : : : : : : : : : : : : : : : : : : ',
                ';Ranked-tissues=', apply(max.tissue.df.LCAP[,2:6], 1, function(x) paste0(x, collapse = ' / ')),
                ';Ranked-tissue-RNAseq-TPM=', apply(max.tissue.df.LCAP[,7:12], 1, function(x) paste0(round(x, 3), collapse = ' / ')),
                ';Hypod-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Hypod..v.others'],
                ';Neurons-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Neurons.v.others'],
                ';Germline-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Germline.v.others'],
                ';Muscle-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Muscle.v.others'],
                ';Intest-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Intest..v.others'],
                ";Enriched-in=", max.tissue.df.LCAP$which.tissues
            )

            GFF_LCAP <- cbind(
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
            
            # Generate REs GFF part
            infos <- paste0(
                "ID=", make.unique(paste(all.deconv$gene_name, all.deconv$regulatory_class, sep = '--'), sep = '_'), 
                ";Associated-gene-Name=", all.deconv$gene_name, 
                ";Associated-gene-WB=", all.deconv$WormBaseID, 
                ";domain=", all.deconv$domain
            )
            
            values <- paste0(
                ';color=', color.tissues[all.deconv$which.tissues], 
                '; =: : : : : : : : : : : : : : : : : : : : : : : : : ', 
                ';Ranked-tissues=', apply(all.deconv[, grep('max.tissue$', colnames(all.deconv))], 1, function(x) paste0(x, collapse = ' / ')), 
                ';Ranked-tissue-ATACseq-TPM=', apply(all.deconv[, grep('max.tissue.cov$', colnames(all.deconv))], 1, function(x) paste0(round(x, 3), collapse = ' / ')), 
                ';Consecutive-ratios=', apply(all.deconv[, grep('ratio', colnames(all.deconv))], 1, function(x) paste0(round(x, 3), collapse = ' / ')), 
                ";Enriched-tissue.s.=", gsub('\\.', '', all.deconv$which.tissues)
            )
            
            GFF_ATAC <- cbind(
                as.character(all.deconv$chr), 
                rep('Ahringer-tissue-spe-REs', times = nrow(all.deconv)), 
                as.character(all.deconv$regulatory_class), 
                all.deconv$start, 
                all.deconv$stop, 
                rep(".", times=nrow(all.deconv)), 
                ifelse(is.null(all.deconv$strand), ".", all.deconv$strand), 
                rep(".", times=nrow(all.deconv)), 
                paste0(infos, values)
            )
            
            # Export GFF file
            write.table(rbind(HEADER, GFF_LCAP[names(genes.gtf) %in% multipleGenes(),], GFF_ATAC[all.deconv$uniqueWormBaseID %in% multipleGenes(),]), file, quote = F, row = F, col = F, sep = '\t')

        }
      )
     
  }
  
  # Generate buttons to access to Genome Browser and WormBase
  {
      
      observeEvent(input$switchToGenome, { updateTabItems(session, "tabs", "browser") })
      
      output$Link <- renderUI({
          tags$a( "WormBase Link", href = paste0("https://www.wormbase.org/species/c_elegans/gene/", infos.gene()$Gene.info[1]), target = "_blank" )
      })
      
  }
  
  # Generate browser
  {
    RE.coords <- reactive ({ c(
      as.character(infos.gene()$Associated.REs[1,1]),
      as.character(min(infos.gene()$Associated.REs[,2])-3000),
      as.character(max(infos.gene()$Associated.REs[,3])+3000)
    ) })
    gene.coords <- reactive ({ c(
      paste0('chr', as.character(seqnames(genes.gtf[infos.gene()$Gene.info[1]]))),
      as.character(start(genes.gtf[infos.gene()$Gene.info[1]])-3000),
      as.character(end(genes.gtf[infos.gene()$Gene.info[1]])+3000)
    ) })
    coords <- reactive ({ c(RE.coords()[1], min(RE.coords()[2], gene.coords()[2]), max(RE.coords()[3], gene.coords()[3])) })
    #-->
    url <- reactive ({ getURL(as.character(coords()[1]), as.numeric(coords()[2]), as.numeric(coords()[3]), "1.12.5") })
    output$jbrowser <- renderUI(
        tags$div(
            id="jbrowser", 
            style="width: 100%; height: 100%; visibility: inherit;",
            class="trewjb html-widget html-widget-output shiny-bound-output",
            div(
                style = "width: 100%; height: calc(100vh - 150px);", 
                tags$iframe(
                    style = "border: 1px solid black",
                    width = "100%",
                    height = "100%",
                    src = url()
                )
            )
        )
    )
  }
  
  # Generate the quickView bsModal 
  {
    
    output$quickResults <- renderUI({
        
        gene <- reactive ({ input$quickGene })
        infos.gene <- reactive ({ getGeneInfos(gene(), saveTXT = F, verbose = F, exportResult = T) })
        WBID <- paste("<b>WormBase ID:</b>   ", infos.gene()$Gene.info[1])
        LOCUS <- paste("<b>Locus:</b>   ", infos.gene()$Gene.info[2])
        COORDS <- paste("<b>Coordinates:</b>   ", as.character(base::range(genes.gtf[infos.gene()$Gene.info[1]])))
        BIOTYPE <- paste("<b>Gene biotype:</b>   ", infos.gene()$Gene.info[3])
        ENRICHED <- paste("<b>Enriched in tissue:</b>   ", infos.gene()$Gene.info[5])
        
        h3("Quick gene view")
        HTML(paste(
            h3("Quick gene view"), 
            WBID,
            LOCUS, 
            COORDS,
            BIOTYPE, 
            ENRICHED, 
            br(),
            h3("Gene description from WormBase"), 
            p(HTML(fetchWBinfos(infos.gene()$Gene.info[2]))),
            sep = '<br/>'
        ))
        
    })
    
    observeEvent(input$quickSearch, {
        toggleModal(session, "quickGENE", toggle = "open")
    })
    
  }
  
  # Generate the REs subset table
  {
    #-->
    
    output$REs.table <- renderDataTable({
        if (all(infos.gene()$Associated.REs == 0)) { datatable(matrix(0)) } else {
            datatable(
                setNames(cbind(infos.gene()$Associated.REs[c(1:4, 6)], round(infos.gene()$Associated.REs.dev, 1), round(infos.gene()$Associated.REs.tissue, 1)),
                c("Chr", "Start", "Stop", "Regulatory Class", "Enriched in tissue(s)", "Coverage Emb. (mixed)", "Coverage L1 (mixed)", "Coverage L2 (mixed)", "Coverage L3 (mixed)", "Coverage L4 (mixed)", "Coverage YA (mixed)", "Coverage Hypod. (YA)", "Coverage Neurons (YA)", "Coverage Germline (YA)", "Coverage Muscle (YA)", "Coverage Intest. (YA)")),
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
#  {
#
#    # Plot genes
#    output$tSNE.plots_genes <- renderPlot({
#        par(mar=c(2,4,2,2), oma = c(0,0,0,0))
#        plot(tSNE.df.LCAP,
#            col = CLASSES.genes,
#            pch = 20,
#            cex = 1.5,
#            xlab = '',
#            ylab = '',
#            axes = T,
#            xaxt='n',
#            yaxt='n',
#            bty = 'o',
#            main = 't-SNE projection of genes',
#            cex.main = 1.5
#        )
#        points(
#            tSNE.df.LCAP[genes.gtf[genes.valid]$gene_id %in% infos.gene()$Gene.info[1],'V1'],
#            tSNE.df.LCAP[genes.gtf[genes.valid]$gene_id %in% infos.gene()$Gene.info[1],'V2'],
#            col = 'black',
#            pch=20,
#            cex = 4
#        )
#        mtext(line = 1, side = 1, 'tSNE dim. 1')
#        mtext(line = 1, side = 2, 'tSNE dim. 2')
#    })
#    
#    # Plot promoters
#    output$tSNE.plots_proms <- renderPlot({
#        par(mar=c(2,4,2,2), oma = c(0,0,0,0))
#        plot(tSNE.df.proms,
#            col = CLASSES.proms,
#            pch = 20,
#            cex = 1.5,
#            xlab = '',
#            ylab = '',
#            axes = T,
#            xaxt='n',
#            yaxt='n',
#            bty = 'o',
#            main = 't-SNE projection of promoters',
#            cex.main = 1.5
#        )
#        if (!all(infos.gene()$Associated.REs == 0)) {
#            points(
#                tSNE.df.proms[all.proms.valid$chr %in% infos.gene()$Associated.REs[,1] & all.proms.valid$start %in% infos.gene()$Associated.REs[,2] & all.proms.valid$stop %in% infos.gene()$Associated.REs[,3],'V1'],
#                tSNE.df.proms[all.proms.valid$chr %in% infos.gene()$Associated.REs[,1] & all.proms.valid$start %in% infos.gene()$Associated.REs[,2] & all.proms.valid$stop %in% infos.gene()$Associated.REs[,3],'V2'],
#                col = 'black',
#                pch=20,
#                cex = 4
#            )
#        }
#        mtext(line = 1, side = 1, 'tSNE dim. 1')
#        mtext(line = 1, side = 2, 'tSNE dim. 2')
#    })
#    
#    # Plot legend
#    output$tSNE.plots_legend <- renderPlot({
#        plot.new()
#        par(mar=c(0,0,0,0), oma = c(0,0,0,0), xpd = NA)
#        legend("left",
#              legend = c(order.tissues[1:5], '2 tissues', '3 tissues', '4 tissues', order.tissues[32:34]),
#              fill = c(color.tissues[1:5], color.tissues[6], color.tissues[16], color.tissues[26], color.tissues[32:34]),
#              ncol = 1, cex = 2, bty = 'n'
#        )
#    })
#
#  }

  # Get the list of multiple genes
  {
    
    multipleGenes <- eventReactive( input$getList, { 
        
        # Get gene names from text box (allows for *)
        input.genes <- unlist(strsplit(x = gsub(" ", "", input$searchMulitpleGenePaste), split = ',|;|[\r\n]' ))
        if(any(grepl('\\*', input.genes))) { input.genes <- c(genes.gtf$gene_name[grep(paste(input.genes[grepl('\\*', input.genes)], collapse = '|'), genes.gtf$gene_name)], input.genes[!grepl('\\*', input.genes)]) }
        text.list <- unique(input.genes %>% ifelse(!grepl('WBGene', .), name2WB(.), .) %>% .[!is.na(.)])
        # Get gene names from classes of genes
        checkbox.list <- unlist(list.genes[input$checkGroupGeneClasses])
        # Get gene names from the bsModal of genes associated with tissue-specific promoters
        bsmodal.list <- row.names(mat.proms.gene)[which(apply(mat.proms.gene, 1, function(x) {all(names(x)[x == T] %in% input$checkGroupGeneClasses_2) & all(input$checkGroupGeneClasses_2 %in% names(x)[x == T])}))]
        
        genes <- unique(c(text.list, checkbox.list, bsmodal.list))
        return(genes[genes %in% names(genes.gtf)])
        
    } )
    
    observe( if (input$resetGenes > 0) {
          
          updateTextAreaInput(
              session, 
              "searchMulitpleGenePaste", 
              value = "", 
              placeholder = 'Paste here (one gene per line)'
          )
          updateCheckboxGroupInput(
              session, 
              "checkGroupGeneClasses", 
              choices = list(
                  "Hypodermis-enriched genes" = "hypod.genes", 
                  "Neurons-enriched genes" = "neurons.genes", 
                  "Germline-enriched genes" = "germline.genes", 
                  "Muscle-enriched genes" = "muscle.genes", 
                  "Intestine-enriched genes" = "intest.genes" 
              ),
              selected = NULL
          )
          updateCheckboxGroupInput(
              session, 
              "checkGroupGeneClasses_2", 
              choiceNames = c(
                  "Hypodermis-specific promoter(s)",
                  "Neurons-specific promoter(s)",
                  "Germline-specific promoter(s)",
                  "Muscle-specific promoter(s)",
                  "Intestine-specific promoter(s)",
                  "Hypodermis & Neurons-specific promoter(s)",
                  "Hypodermis & Germline-specific promoter(s)",
                  "Hypodermis & Muscle-specific promoter(s)",
                  "Hypodermis & Intestine-specific promoter(s)",
                  "Neurons & Germline-specific promoter(s)",
                  "Neurons & Muscle-specific promoter(s)",
                  "Neurons & Intestine-specific promoter(s)",
                  "Germline & Muscle-specific promoter(s)",
                  "Germline & Intestine-specific promoter(s)",
                  "Muscle & Intestine-specific promoter(s)",
                  "Hypodermis & Neurons & Germline-specific promoter(s)",
                  "Hypodermis & Neurons & Muscle-specific promoter(s)",
                  "Hypodermis & Neurons & Intestine-specific promoter(s)",
                  "Hypodermis & Germline & Muscle-specific promoter(s)",
                  "Hypodermis & Germline & Intestine-specific promoter(s)",
                  "Hypodermis & Muscle & Intestine-specific promoter(s)",
                  "Neurons & Germline & Muscle-specific promoter(s)",
                  "Neurons & Germline & Intestine-specific promoter(s)",
                  "Neurons & Muscle & Intestine-specific promoter(s)",
                  "Germline & Muscle & Intestine-specific promoter(s)",
                  "Neurons & Germline & Muscle & Intestine-specific promoter(s)",
                  "Hypodermis & Germline & Muscle & Intestine-specific promoter(s)",
                  "Hypodermis & Neurons & Germline & Intestine-specific promoter(s)",
                  "Hypodermis & Neurons & Germline & Muscle-specific promoter(s)",
                  "Soma-specific promoter(s)",
                  "Ubiquitous promoter(s)"
              ),
              choiceValues = order.tissues[c(1:27, 29, 30, 32:33)],
              selected = NULL
          )
          
    } )
    
    output$multipleGenesPromsGroupsLength <- reactive ({ paste(length(input$checkGroupGeneClasses_2), "group(s) selected.") })
    output$multipleGenesLength <- reactive ({ 
        l <- length(multipleGenes())
        if (l < 2) {
            'No genes found.\nEnter or select at least 2 genes and click on "Perform analysis" button.'
        } else {
        paste(l, "valid genes found.") 
        }
    })
    
    output$genesList <- renderUI({ HTML(paste(c("<h3>Genes query</h3>", sort(WB2name(multipleGenes()))), collapse = '<br/>')) })

  }

  # Get GOs plot (enriched in multiple genes) and button to download results
  {
    
    pathways <- reactive ({ input$checkGroupGOs })
    filteringSetting <- reactive ({ input$hierarchyFiltering })
    gProfileR_results <- eventReactive( input$runGO, { gprofiler(multipleGenes(), organism="celegans", max_p_value=0.05, correction_method="bonferroni", hier_filtering=filteringSetting()) } )
    
    output$GO.plot <- renderPlot({

        reduced_c <- gProfileR_results()
        reduced_c <- reduced_c[order(reduced_c$p.value),]
        d <- reduced_c[reduced_c$domain %in% pathways(),]
        d <- d[!duplicated(d[, c("p.value", "overlap.size")],),]
        d <- d[1:min(20, dim(d)[1]),]
        if (nrow(d) > 0) {
          par(las = 0, mar = c(4,20,4,1))
          plot <- barplot(height=-log10(d$p.value), col=colorGO[d$domain], horiz=T, main = "Enrichment of associated GOs", xlim = c(0, max(25, d$p.value)))
          par(xpd=T)
          #text(x=rep(-1, times = length(plot)), y = plot, paste0(substr(d$term.name, start = 1, stop = 40), ifelse(nchar(d$term.name) > 40, '...', ''), ' (n=', d$overlap.size, ')'), pos = 2)
          text(x=rep(-1, times = length(plot)), y = plot, paste0(substr(d$term.name, start = 1, stop = 40), ifelse(nchar(d$term.name) > 40, '...', '')), pos = 2)
          par(xpd=F)
          legend("topright", legend=names(colorGO), fill=colorGO, col="#00000000", pch=15, bty="n")
          mtext(text="-log10(adj-p.val)",side=1,line=3,outer=FALSE, cex=0.85)
        } else { plot.new() }
    })
    
    output$downloadGO <- downloadHandler(
        "gene-list_GO-summary.txt",
        content = function(file) {
            reduced_c <- gprofiler(multipleGenes(), organism="celegans", max_p_value=0.05, correction_method="bonferroni", hier_filtering='none')
            reduced_c <- reduced_c[order(reduced_c$p.value),]
            write.table(reduced_c, file, row.names = F, quote = F, sep = "\t")
        }
    )


  }

  # Get HMs.plot (multiple genes heatmaps of LCAP/ATAC fold-changes)
  {
    colorScale_LCAPdev <- reactive ({ colorRampPalette(brewer.pal(brewer.pal.info[input$colorScale_LCAPdev,1], input$colorScale_LCAPdev))(100) })
    colorScale2_LCAPdev <- reactive ({ if( input$colorScale_doRev_LCAPdev ) { rev(colorScale_LCAPdev()) } else { colorScale_LCAPdev() } })
    titleLab_LCAPdev <- reactive ({ if( input$LCAPdev_TPMZscore ) { "log2TPM" } else { "z-score" } })
    clusterFUN_LCAPdev <- reactive ({ if( input$clusterFUN_LCAPdev ) { "hclust" } else { "kmeans" } })
    titleTAB_LCAPdev <- reactive ({ paste0("heatmap_dev-RNA-seq_", titleLab_LCAPdev(),"_", ifelse(clusterFUN_LCAPdev() == 'hclust', 'hclust', paste0(input$NCLUST_LCAPdev, "-clusters")),"_genesets.txt") })
    
    colorScale_LCAP <- reactive ({ colorRampPalette(brewer.pal(brewer.pal.info[input$colorScale_LCAP,1], input$colorScale_LCAP))(100) })
    colorScale2_LCAP <- reactive ({ if( input$colorScale_doRev_LCAP ) { rev(colorScale_LCAP()) } else { colorScale_LCAP() } })
    titleLab_LCAP <- reactive ({ if( input$LCAP_TPMZscore ) { "log2TPM" } else { "z-score" } })
    clusterFUN_LCAP <- reactive ({ if( input$clusterFUN_LCAP ) { "hclust" } else { "kmeans" } })
    titleTAB_LCAP <- reactive ({ paste0("heatmap_tissues-RNA-seq_", titleLab_LCAP(),"_", ifelse(clusterFUN_LCAP() == 'hclust', 'hclust', paste0(input$NCLUST_LCAP, "-clusters")),"_genesets.txt") })

    colorScale_ATAC <- reactive ({ colorRampPalette(brewer.pal(brewer.pal.info[input$colorScale_ATAC,1], input$colorScale_ATAC))(100) })
    colorScale2_ATAC <- reactive ({ if( input$colorScale_doRev_ATAC ) { rev(colorScale_ATAC()) } else { colorScale_ATAC() } })
    titleLab_ATAC <- reactive ({ if( input$ATAC_TPMZscore ) { "log2TPM" } else { "z-score" } })
    clusterFUN_ATAC <- reactive ({ if( input$clusterFUN_ATAC ) { "hclust" } else { "kmeans" } })
    titleTAB_ATAC <- reactive ({ paste0("heatmap_tissues-ATAC-seq_", titleLab_ATAC(),"_", ifelse(clusterFUN_ATAC() == 'hclust', 'hclust', paste0(input$NCLUST_ATAC, "-clusters")),"_genesets.txt") })

    # Plot LCAPdev
    output$HMs.plot_LCAPdev <- renderPlot({
        par(mar = c(6,1,4,1))
        if (titleLab_LCAPdev() == "log2TPM") {
            mat <- log2(LCAPdev[multipleGenes(),]+1)
            breaks <- NA
        } else {
            mat <- t(apply(LCAPdev[multipleGenes(),], 1, scale)) %>% na.replace(., 0)
            colnames(mat) <- colnames(LCAPdev)
            breaks <- seq(-2, 2, length.out = 101)
        }
        row.names(mat) <- WB2name(row.names(mat))
        if (clusterFUN_LCAPdev() == 'hclust') {
            doHCLUST <- T
            ANNOTS_DF <- NULL
            ANNOTS_COL <- NULL
        } else {
            doHCLUST <- F
            NCLUST_LCAPdev <- input$NCLUST_LCAPdev
            set.seed(32) ; ANNOTS <- sort(kmeans(mat, NCLUST_LCAPdev)$cluster, decreasing = F)
            mat <- mat[names(ANNOTS),]
            ANNOTS_DF <- data.frame(Cluster = ANNOTS)
            ANNOTS_COL <- list(Cluster = rep(brewer.pal('Set1', n = 11), 2)[1:NCLUST_LCAPdev])
        }
        pheatmap(
            mat,
            color = colorScale2_LCAPdev(),
            scale = "none", 
            breaks = breaks, 
            cluster_rows = doHCLUST, 
            cluster_cols = F, 
            treeheight_row = 20,
            show_rownames = ifelse(nrow(mat) > 20, F, T), 
            main = paste0("Mixed-tissue gene expression (", titleLab_LCAPdev(), ")\nat each developmental stage"),
            annotation_row = ANNOTS_DF,
            annotation_colors = ANNOTS_COL
        )
    })

    # Plot LCAP
    output$HMs.plot_LCAP <- renderPlot({
        par(mar = c(6,1,4,1))
        if(titleLab_LCAP() == "log2TPM") {
            mat <- log2(LCAP[multipleGenes(),]+1)
            breaks <- NA
        } else {
            mat <- t(apply(LCAP[multipleGenes(),], 1, scale)) %>% na.replace(., 0)
            colnames(mat) <- order.tissues[1:5]
            breaks <- seq(-2, 2, length.out = 101)
        }
        row.names(mat) <- WB2name(row.names(mat))
        if (clusterFUN_LCAP() == 'hclust') {
            doHCLUST <- T
            ANNOTS_DF <- NULL
            ANNOTS_COL <- NULL
        } else {
            doHCLUST <- F
            NCLUST_LCAP <- input$NCLUST_LCAP
            set.seed(32) ; ANNOTS <- sort(kmeans(mat, NCLUST_LCAP)$cluster, decreasing = F)
            mat <- mat[names(ANNOTS),]
            ANNOTS_DF <- data.frame(Cluster = ANNOTS)
            ANNOTS_COL <- list(Cluster = rep(brewer.pal('Set1', n = 11), 2)[1:NCLUST_LCAP])
        }
        pheatmap(
            mat,
            color = colorScale2_LCAP(),
            scale = "none", 
            breaks = breaks, 
            cluster_rows = doHCLUST, 
            cluster_cols = F, 
            treeheight_row = 20,
            show_rownames = ifelse(nrow(mat) > 20, F, T), 
            main = paste0("Gene expression (", titleLab_LCAP(), ")\nin each tissue"),
            annotation_row = ANNOTS_DF,
            annotation_colors = ANNOTS_COL
        )
    })
    
    # Plot ATAC
    output$HMs.plot_ATAC <- renderPlot({
        par(mar = c(6,1,4,1))
        if(titleLab_ATAC() == "log2TPM") {
            mat <- log2(all.deconv[all.deconv$uniqueWormBaseID %in% multipleGenes(), grep(paste(order.tissues, collapse = '|'), colnames(all.deconv))] + 1)
            breaks <- NA
        } else {
            mat <- t(apply(all.deconv[all.deconv$uniqueWormBaseID %in% multipleGenes(), grep(paste(order.tissues, collapse = '|'), colnames(all.deconv))], 1, scale)) %>% na.replace(., 0)
            colnames(mat) <- order.tissues[1:5]
            breaks <- seq(-2, 2, length.out = 101)
        }
        if (clusterFUN_ATAC() == 'hclust') {
            doHCLUST <- T
            ANNOTS_DF <- NULL
            ANNOTS_COL <- NULL
        } else {
            doHCLUST <- F
            NCLUST_ATAC <- input$NCLUST_ATAC
            set.seed(32) ; ANNOTS <- sort(kmeans(mat, NCLUST_ATAC)$cluster, decreasing = F)
            mat <- mat[names(ANNOTS),]
            ANNOTS_DF <- data.frame(Cluster = ANNOTS)
            ANNOTS_COL <- list(Cluster = rep(brewer.pal('Set1', n = 11), 2)[1:NCLUST_ATAC])
        }
        pheatmap(
            mat,
            color = colorScale2_ATAC(),
            scale = "none", 
            breaks = breaks, 
            cluster_rows = doHCLUST, 
            cluster_cols = F, 
            treeheight_row = 20, 
            show_rownames = F,
            main = paste0("Associated REs accessibility (", titleLab_ATAC(), ")\nin each tissue"),
            annotation_row = ANNOTS_DF,
            annotation_colors = ANNOTS_COL
        )
    })
    
    ## Generate HM tables to download
    output$downloadHM_LCAPdev <- downloadHandler( 
        titleTAB_LCAPdev, 
        content = function(file) {
            if(titleLab_LCAPdev() == "log2TPM") {
                mat <- log2(LCAPdev[multipleGenes(),]+1)
            } else {
                mat <- t(apply(LCAPdev[multipleGenes(),], 1, scale)) %>% na.replace(., 0)
                colnames(mat) <- colnames(LCAPdev)
            }
            row.names(mat) <- WB2name(row.names(mat))
            if (clusterFUN_LCAPdev() == 'kmeans') {
                set.seed(32) ; VEC <- sort(kmeans(mat, input$NCLUST_LCAPdev)$cluster, decreasing = F)
                mat <- cbind(mat[names(VEC),], clusterNB = VEC)
            }
            write.table(round(mat, 3), file, quote = F, row = T, col = T, sep = '\t')
        }
    )
    output$downloadHM_LCAP <- downloadHandler( 
        titleTAB_LCAP, 
        content = function(file) {
            if(titleLab_LCAP() == "log2TPM") {
                mat <- log2(LCAP[multipleGenes(),]+1)
            } else {
                mat <- t(apply(LCAP[multipleGenes(),], 1, scale)) %>% na.replace(., 0)
                colnames(mat) <- colnames(LCAP)
            }
            row.names(mat) <- WB2name(row.names(mat))
            if (clusterFUN_LCAP() == 'kmeans') {
                set.seed(32) ; VEC <- sort(kmeans(mat, input$NCLUST_LCAP)$cluster, decreasing = F)
                mat <- cbind(mat[names(VEC),], clusterNB = VEC)
            }
            write.table(round(mat, 3), file, quote = F, row = T, col = T, sep = '\t')
        }
    )
    output$downloadHM_ATAC <- downloadHandler( 
        titleTAB_ATAC, 
        content = function(file) {
            VEC0 <- all.deconv$uniqueWormBaseID %in% multipleGenes()
            if(titleLab_ATAC() == "log2TPM") {
                mat <- log2(all.deconv[VEC0, grep(paste(order.tissues, collapse = '|'), colnames(all.deconv))] + 1)
            } else {
                mat <- t(apply(all.deconv[VEC0, grep(paste(order.tissues, collapse = '|'), colnames(all.deconv))], 1, scale)) %>% na.replace(., 0)
                colnames(mat) <- order.tissues[1:5]
            }
            row.names(mat) <- make.unique(all.deconv[VEC0,]$locus)
            if (clusterFUN_ATAC() == 'kmeans') {
                set.seed(32) ; VEC <- sort(kmeans(mat, input$NCLUST_ATAC)$cluster, decreasing = F)
                mat <- cbind(mat[names(VEC),], clusterNB = VEC)
            }
            write.table(round(mat, 3), file, quote = F, row = T, col = T, sep = '\t')
        }
    )

  }
  
  # Get venn Diagrams
  {
      
      output$Venn.Germline <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[[1]], multipleGenes(), names = c("Germline-enrich.\ngenes", "Genes query"), col = c(color.tissues[1], 'grey50')) })
      output$Venn.Neurons <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[[2]], multipleGenes(), names = c("Neurons-enrich.\ngenes", "Genes query"), col = c(color.tissues[2], 'grey50')) })
      output$Venn.Muscle <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[[3]], multipleGenes(), names = c("Muscle-enrich.\ngenes", "Genes query"), col = c(color.tissues[3], 'grey50')) })
      output$Venn.Hypod <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[[4]], multipleGenes(), names = c("Hypod.-enrich.\ngenes", "Genes query"), col = c(color.tissues[4], 'grey50')) })
      output$Venn.Intest <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[[5]], multipleGenes(), names = c("Intest.-enrich.\ngenes", "Genes query"), col = c(color.tissues[5], 'grey50')) })
      
  }
  
  # Generate tables to display and button to download them
  {
      output$downloadATAC.txt <- downloadHandler( "tissue-specific.ATAC-seq.dataset.txt", content = function(file) {write.table(atac.dt, file, quote = F, row = F, col = T, sep = '\t')} )
      output$downloadATAC.gff <- downloadHandler( "tissue-specific.ATAC-seq.dataset.gff", content = function(file) {
          infos=paste0("ID=", make.unique(paste(all$gene_name, all$regulatory_class, sep = '--'), sep = '_'), ";Associated-gene-Name=", all$gene_name, ";Associated-gene-WB=", all$WormBaseID, ";domain=", all$domain)
          values <- paste0(';color=', color.tissues[max.tissue.df$which.tissues], '; =: : : : : : : : : : : : : : : : : : : : : : : : : ', ';Ranked-tissues=', apply(max.tissue.df[,5:9], 1, function(x) paste0(x, collapse = ' / ')), ';Ranked-tissue-ATACseq-TPM=', apply(max.tissue.df[,10:14], 1, function(x) paste0(round(x, 3), collapse = ' / ')), ';Consecutive-ratios=', apply(max.tissue.df[,15:18], 1, function(x) paste0(round(x, 3), collapse = ' / ')), ";Enriched-tissue.s.=", gsub('\\.', '', max.tissue.df$which.tissues))
          GFF=cbind(as.character(all$chr), rep('Ahringer-tissue-spe-REs', times = nrow(all)), as.character(all$regulatory_class), all$start, all$stop, rep(".", times=nrow(all)), ifelse(is.null(all$strand), ".", all$strand), rep(".", times=nrow(all)), paste0(infos, values))
          write.table(rbind(c("##gffTags=on\n##displayName=Name\n##gff-version 3", rep(" ", 8)), GFF), file, quote = F, row = F, col = F, sep = '\t')
      } )
      output$downloadLCAP.txt <- downloadHandler( "tissue-specific.RNA-seq.dataset.txt", content = function(file) {write.table(lcap.dt, file, quote = F, row = F, col = T, sep = '\t')}  )
      output$downloadLCAP.gff <- downloadHandler( "tissue-specific.RNA-seq.dataset.gff", content = function(file) {
          
            names(genes.gtf) <- genes.gtf$gene_name

            infos=paste0(
                "ID=", genes.gtf$gene_id,
                ";Name=", genes.gtf$gene_name,
                ';color=', color.tissues[c(1:5,33:35)][max.tissue.df.LCAP$which.tissues],
                '; =: : : : : : : : : : : : : : : : : : : : : : : : : ',
                ';Ranked-tissues=', apply(max.tissue.df.LCAP[,2:6], 1, function(x) paste0(x, collapse = ' / ')),
                ';Ranked-tissue-RNAseq-TPM=', apply(max.tissue.df.LCAP[,7:12], 1, function(x) paste0(round(x, 3), collapse = ' / ')),
                ';Hypod-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Hypod..v.others'],
                ';Neurons-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Neurons.v.others'],
                ';Germline-FC-vs-mean=', max.tissue.df.LCAP[,'ratio.Germline.v.others'],
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

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

## Define server function ------------------------------------------------------

shinyServer <- function(input, output, session) {
    
    # Read gene name and get gene infos
    {
        gene <- reactive({ input$searchGene })
        output$gene <- renderText({ gene() })
        infos.gene <- reactive ({ getGeneInfos(gene(), saveTXT = F, verbose = F, exportResult = T) })
        # Get the gene info / description
        output$geneInfos <- renderUI({
            WBID <- paste("<b>WormBase ID:</b>   ", infos.gene()$Gene.info[['WormBaseID']])
            LOCUS <- paste("<b>Locus:</b>   ", infos.gene()$Gene.info[['locus']])
            COORDS <- paste("<b>Coordinates:</b>   ", as.character(base::range(genes.gtf[infos.gene()$Gene.info[['WormBaseID']]])) %>% gsub('.{2}$', '', .))
            STRAND <- paste("<b>Orientation:</b>   ", as.character(base::range(genes.gtf[infos.gene()$Gene.info[['WormBaseID']]])) %>% substrRight(1))
            BIOTYPE <- paste("<b>Gene biotype:</b>   ", infos.gene()$Gene.info[['Biotype']])
            ENRICHED <- paste("<b>Gene class:</b>   ", infos.gene()$Gene.info[['Enriched.tissue']])
            HTML(paste(WBID, LOCUS, COORDS, STRAND, BIOTYPE, ENRICHED, sep = '<br/>'))
        })
        output$geneDescr <- renderUI({ INFOS <- HTML(fetchWBinfos(infos.gene()$Gene.info[2])) })
        # Get the gene expression graphs
        output$Expr.plots_dev <- renderPlot({
            par(mar=c(6,5,4,5))
            vec <- unlist(infos.gene()[['Gene.expr.dev.TPM']])
            if (is.null(vec)) vec <- rep(0, 5)
            barplot(
                vec,
                col = 'grey', 
                names = c('Emb.', 'L1', 'L2', 'L3', 'L4', 'YA'),
                xlab = "", ylab = "log2 TPM", 
                main = "Developmental gene expression", 
                las = 3
            )
        })
        output$Expr.plots_tis <- renderPlot({
            par(mar=c(6,5,4,5))
            vec <- unlist(infos.gene()[['Gene.expr.TPM']]['tiss.spe.LCAP',])
            if (is.null(vec)) vec <- rep(0, 5)
            barplot(
                vec,
                col = color.tissues, ylab = "log2 TPM",
                names = order.tissues[1:5],
                main = "Background-normalized\ntissue-specific gene expression (YA)", 
                las = 3
            )
        })
        addTooltip(session, "Expr.plots_tis", "Important", content = paste0(
            "<p>This plot shows gene expression in TPMs (Transcripts Per Millions) from nuclear RNA-seq, 
            with background RNA contamination corrected by DSA. 
            This normalization improves the quantification of gene expression for tissue-specific genes.
            Hence, these values are recommended when looking at a tissue-specific gene.</p>"
        ), trigger = 'hover')
        output$Expr.plots_tis_uncorrected <- renderPlot({
            par(mar=c(6,5,4,5))
            vec <- unlist(LCAP[infos.gene()$Gene.info[['WormBaseID']],])
            if (is.null(vec)) vec <- rep(0, 5)
            barplot(
                vec,
                col = color.tissues, ylab = "log2 TPM",
                names = order.tissues[1:5],
                main = "Tissue-specific gene expression (YA)", 
                las = 3
            )
        })
        addTooltip(session, "Expr.plots_tis_uncorrected", "Important", content = paste0(
            "<p>This plot shows gene expression in TPMs (Transcripts Per Millions) from nuclear RNA-seq. 
            These values are recommended when looking at a broadly expressed gene (i.e. Ubiquitous).</p>"
        ), trigger = 'hover')
        # Trigger the display of the information
        output$displayPanelsTab1 <- renderText({ ifelse(infos.gene()$valid, 1, 0) })
        outputOptions(output, "displayPanelsTab1", suspendWhenHidden = FALSE)
    }
    
    # Generate buttons to download gene-specific text file, all tracks in zip, 
    # and genes list full report
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
            "tissue-specific-RNAseq_bigwig-tracks.tar.gz",
            content = function(file) {
                download.file("http://ahringerlab.com/public/RNAseq_tracks.tar.gz", file)
            }
        )
        
        output$downloadBWATAC <- downloadHandler(
            "tissue-specific-ATAC-seq_bigwig-tracks.tar.gz",
            content = function(file) {
                download.file("http://ahringerlab.com/public/ATACseq_tracks.tar.gz", file)
            }
        )
        
        output$downloadGenesListGFF <- downloadHandler(
            "gene-list_full-report.gff",
            content = function(file) {
                
                # Generate HEADER
                HEADER <- rbind(
                    c("##gffTags=on\n##displayName=Name\n##gff-version 3", rep(" ", 8)),
                    c(paste0("##FILE GENERATE ON: ", date(), ".\n##NUMBER OF GENES FOUND IN THE QUERY: ", length(multipleGenes())), rep(" ", 8)),
                    c("##\n##This file contains both genes and associated regulatory elements annotations. Last column contains general information as well as details of tissue-specific expression / accessibility of the genes / associated REs.", rep(" ", 8)),
                    c("##This file can easily be loaded in IGV. Coordinates are ce11.", rep(" ", 8))
                )
                
                # Generate REs GFF part
                infos <- paste0(
                    "ID=", make.unique(paste(all$gene_name, all$regulatory_class, sep = '--'), sep = '_'), 
                    ";Associated-gene-Name=", all$gene_name, 
                    ";Associated-gene-WB=", all$WormBaseID, 
                    ";domain=", all$domain
                )
                values <- paste0(
                    ';color=', color.tissues[max.tissue.df$which.tissues], 
                    '; =: : : : : : : : : : : : : : : : : : : : : : : : : ', 
                    ';Ranked-tissues=', apply(max.tissue.df[,5:9], 1, function(x) paste0(x, collapse = ' / ')), 
                    ';Ranked-tissue-ATACseq-TPM=', apply(max.tissue.df[,10:14], 1, function(x) paste0(round(x, 3), collapse = ' / ')), 
                    ';Consecutive-ratios=', apply(max.tissue.df[,15:18], 1, function(x) paste0(round(x, 3), collapse = ' / ')), 
                    ";Tissue-specificity=", gsub('\\.', '', max.tissue.df$which.tissues)
                )
                GFF_ATAC <- cbind(
                    as.character(all$chr), 
                    rep('Ahringer-tissue-spe-REs', 
                    times = nrow(all)), 
                    as.character(all$regulatory_class), 
                    all$start, 
                    all$stop, 
                    rep(".", times=nrow(all)), 
                    ifelse(is.null(all$strand), ".", all$strand), 
                    rep(".", times=nrow(all)), 
                    paste0(infos, values)
                )
                
                # Generate genes GFF part
                infos <- paste0(
                    "ID=", genes.gtf$gene_id,
                    ";domain=", genes.gtf$domain
                )
                values <- paste0(
                    ';color=', c(color.tissues[1], "#081d89", "#081d89", color.tissues[2:36])[max.tissue.df.LCAP$which.tissues], 
                    '; =: : : : : : : : : : : : : : : : : : : : : : : : : ',
                    ';Ranked-tissues=', apply(max.tissue.df.LCAP[,5:9], 1, function(x) paste0(x, collapse = ' * ')), 
                    ';Ranked-tissue-RNAseq-TPM=', apply(max.tissue.df.LCAP[,grep('cov', colnames(max.tissue.df.LCAP))], 1, function(x) paste0(round(x, 3), collapse = ' * ')), 
                    ';Tissue-zscore=', apply(round(t(apply(max.tissue.df.LCAP[,grepl('max.tissue.cov', colnames(max.tissue.df.LCAP))], 1, scale)), 2), 1, function(x) { paste(x, collapse = ' * ') } ),
                    ';Consecutive-ratios=', apply(max.tissue.df.LCAP[,15:18], 1, function(x) paste0(round(x, 3), collapse = ' / ')), 
                    ";Tissue-specificity=", gsub('\\.', '', max.tissue.df.LCAP$which.tissues)
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
                    paste0(infos, values)
                )
                
                # Export GFF file
                write.table(rbind(HEADER, GFF_LCAP[names(genes.gtf) %in% multipleGenes(),], GFF_ATAC[all$WormBaseID %in% multipleGenes(),]), file, quote = F, row = F, col = F, sep = '\t')
                
            }
        )
        
    }
    
    # Generate buttons to access to Genome Browser and WormBase
    {
        observeEvent(input$switchToGenome, { 
            updateTabItems(session, "tabs", "browser")
        })
        output$Link <- renderUI({
            url <-paste0("https://www.wormbase.org/species/c_elegans/gene/", infos.gene()$Gene.info[['WormBaseID']])
            HTML(paste0('<button class="action-button bttn bttn-fill bttn-sm bttn-primary bttn-no-outline" id="Link" type="button"><i class="fa fa-book"></i> <a href="', url, '" target="_blank">Go to Wormbase entry</a> </button>'))
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
                    style = "width: 100%; height: calc(100vh - 100px);", 
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
    
    # Generate browser pop-up when clicking on the sidebar
    {
        observe({
            if(input$tabs == 'browser' & counter$cnt < 2) {
                showModal(
                    modalDialog(
                        p('Welcome to our embedded genome browser. This browser is built using ', a(href = 'https://jbrowse.org/', target = "_blank", "JBrowse"), '.'),
                        p('The genome version currently in use is WBcel235/ce11.'),
                        p('You can access the browser outside of our Shiny app by clicking ', a(href = 'http://ahringerlab.com/JBrowse-1.12.5/index.html?data=data%2Fjson%2Fce11&amp;menu=1&amp;nav=1&amp;tracklist=1&amp;overview=1', target = "_blank", "here"), '.'),
                        hr(),
                        h4("Quick tips:"),
                        p("- Local files can be uploaded temporarily and anonymously to this browser using the button in the 'Track' tab."),
                        p("- Each track can be individually re-scaled. Click on the dropdown button appearing when hovering over the track name for more options."),
                        p("- Sessions can be shared with other using the 'Share' button located in the top right corner."),
                        p("- Tissue annotations are available for REs and genes. Click on any element to display more information."),
                        br(),
                        h4('Important:'),
                        h4('More tracks are available in the "Select tracks" tab located on the left side of the genome-browser!'),
                        easyClose = TRUE, 
                        fade = FALSE
                    )
                )
            }
        })
        counter <- reactiveValues(cnt = 0)
        observeEvent(input$tabs, { 
            if(input$tabs == 'browser') {
                counter$cnt <- counter$cnt + 1
            } 
        })
        
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
            ENRICHED <- paste("<b>Class:</b>   ", infos.gene()$Gene.info[4])
            
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
        output$REs.table <- renderDataTable({
            if (all(infos.gene()$Associated.REs == 0)) { datatable(matrix(0)) } else {
                datatable(
                    setNames(
                        cbind(
                            infos.gene()$Associated.REs[c(1:4, 6)],
                            round(infos.gene()$Associated.REs.tissue, 1)
                        ),
                        c("Chr", "Start", "Stop", "Regulatory Class", "Class", "Cov. Germline (YA)", "Cov. Neurons (YA)", "Cov. Muscle (YA)", "Cov. Hypod. (YA)", "Cov. Intest. (YA)")
                    ),
                    autoHideNavigation = T,
                    rownames = F,
                    style = 'bootstrap',
                    extensions = c('Buttons'),
                    options = list(
                        pageLength = 10,
                        dom = 'Brtip',
                        buttons = list('copy', 'print', list(
                            extend = 'collection',
                            buttons = c('csv'),
                            text = 'Download')),
                            scrollX = 200,
                            autoWidth = T,
                            searching = F
                        )
                    )
                }
            })
        }
        
    # Get the list of multiple genes
    {
        # Get the full list of genes
        multipleGenes <- eventReactive( 
            input$getList, 
            { 
                # Get gene names from text box (allows for *)
                input.genes <- unlist(strsplit(x = gsub(" ", "", input$searchMulitpleGenePaste), split = ',|;|[\r\n]' ))
                if(any(grepl('\\*', input.genes))) { input.genes <- c(genes.gtf$gene_name[grep(paste(input.genes[grepl('\\*', input.genes)], collapse = '|'), genes.gtf$gene_name)], input.genes[!grepl('\\*', input.genes)]) }
                text_list <- unique(input.genes %>% ifelse(!grepl('WBGene', .), name2WB(.), .) %>% .[!is.na(.)])
                # Get gene names from the bsModal (list of examples of gene lists)
                bsmodal_list <- unlist(
                    list(
                        "Germline-specific genes" = germline.genes, 
                        "Neurons-specific genes" = neurons.genes, 
                        "Muscle-specific genes" = muscle.genes, 
                        "Hypodermis-specific genes" = hypod.genes, 
                        "Intestine-specific genes" = intest.genes, 
                        "Soma-specific genes" = names(genes.gtf)[genes.gtf$which.tissues == 'Soma'], 
                        "Ubiquitous genes" = names(genes.gtf)[genes.gtf$which.tissues == 'Ubiq.'], 
                        "Genes with germline-specific promoter(s)" = row.names(mat.proms.gene)[row.names(mat.proms.gene) %in% names(genes.gtf) & mat.proms.gene$Germline],
                        "Genes with neurons-specific promoter(s)" = row.names(mat.proms.gene)[row.names(mat.proms.gene) %in% names(genes.gtf) & mat.proms.gene$Neurons],
                        "Genes with muscle-specific promoter(s)" = row.names(mat.proms.gene)[row.names(mat.proms.gene) %in% names(genes.gtf) & mat.proms.gene$Muscle],
                        "Genes with hypodermis-specific promoter(s)" = row.names(mat.proms.gene)[row.names(mat.proms.gene) %in% names(genes.gtf) & mat.proms.gene$Hypod.],
                        "Genes with intestine-specific promoter(s)" = row.names(mat.proms.gene)[row.names(mat.proms.gene) %in% names(genes.gtf) & mat.proms.gene$Intest.],
                        "Genes with Neurons & Muscle-specific promoter(s)" = row.names(mat.proms.gene)[row.names(mat.proms.gene) %in% names(genes.gtf) & mat.proms.gene$Neurons_Muscle],
                        "Genes with Hypodermis & Intestine-specific promoter(s)" = row.names(mat.proms.gene)[row.names(mat.proms.gene) %in% names(genes.gtf) & mat.proms.gene$Hypod._Intest.],
                        "Genes with soma-specific promoter(s)" = row.names(mat.proms.gene)[row.names(mat.proms.gene) %in% names(genes.gtf) & mat.proms.gene$Soma],
                        "Genes with ubiquitous promoter(s)" = row.names(mat.proms.gene)[row.names(mat.proms.gene) %in% names(genes.gtf) & mat.proms.gene$Ubiq.]
                    )[as.numeric(input$checkGroupGeneClasses)]
                )
                # Fuse lists
                genes <- unique(c(text_list, bsmodal_list))
                return(genes[genes %in% names(genes.gtf)])
            } 
        )
        # Reset gene lists when reset button clicked
        observeEvent(input$resetGenes, {
            updateTextAreaInput(
                session, 
                "searchMulitpleGenePaste", 
                value = "", 
                placeholder = 'Paste here (one gene per line)'
            )
            updateCheckboxGroupInput(
                session, 
                "checkGroupGeneClasses", 
                choiceNames = c(
                    "Germline-specific genes",
                    "Neurons-specific genes",
                    "Muscle-specific genes",
                    "Hypodermis-specific genes",
                    "Intestine-specific genes",
                    "Soma-specific genes",
                    "Ubiquitous genes",
                    "Genes with germline-specific promoter(s)",
                    "Genes with meurons-specific promoter(s)",
                    "Genes with muscle-specific promoter(s)",
                    "Genes with hypodermis-specific promoter(s)",
                    "Genes with Intestine-specific promoter(s)",
                    "Genes with Neurons & Muscle-specific promoter(s)",
                    "Genes with Hypodermis & Intestine-specific promoter(s)",
                    "Genes with soma-specific promoter(s)",
                    "Genes with ubiquitous promoter(s)"
                ),
                choiceValues = 1:16,
                selected = NULL
            )
        } )
        # Text infos
        output$multipleGenesPromsGroupsLength <- reactive ({ paste(length(input$checkGroupGeneClasses), "group(s) selected.") })
        output$multipleGenesLength <- reactive({ 
            l <- length(multipleGenes())
            if (l == 0) {
                'No input genes.\nEnter or select at least 2 genes and click on "Perform analysis" button.'
            } else if (l == 1) {
                'No genes found.\nEnter or select at least 2 genes and click on "Perform analysis" button.'
            } else {
                paste(l, "valid genes found.") 
            }
        })
        output$genesList <- renderUI({ HTML(paste(c("<h3>Genes query</h3>", sort(WB2name(multipleGenes()))), collapse = '<br/>')) })
        # Trigger the display of the information
        output$displayPanelsTab2 <- renderText({ ifelse(length(multipleGenes()) > 0, 1, 0) })
        outputOptions(output, "displayPanelsTab2", suspendWhenHidden = FALSE)
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
                par(las = 0, mar = c(4,30,4,1))
                plot <- barplot(height=-log10(d$p.value), col=colorGO[d$domain], horiz=T, main = "Enrichment of associated GOs\n(gProfileR)", xlim = c(0, max(-log10(d$p.value))))
                par(xpd=T)
                #text(x=rep(-1, times = length(plot)), y = plot, paste0(substr(d$term.name, start = 1, stop = 40), ifelse(nchar(d$term.name) > 40, '...', ''), ' (n=', d$overlap.size, ')'), pos = 2)
                text(x=rep(-1/25*max(-log10(d$p.value)), times = length(plot)), y = plot, paste0(substr(d$term.name, start = 1, stop = 60), ifelse(nchar(d$term.name) > 60, '...', '')), pos = 2)
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
        
        output$Venn.Germline <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[['germline.genes']], multipleGenes(), names = c("Germline-specific\ngenes", "Query"), col = c(color.tissues[1], rgb(0.85, 0.85, 0.85, 0.5))) })
        output$Venn.Neurons <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[['neurons.genes']], multipleGenes(), names = c("Neurons-specific\ngenes", "Query"), col = c(color.tissues[2], rgb(0.85, 0.85, 0.85, 0.5))) })
        output$Venn.Muscle <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[['muscle.genes']], multipleGenes(), names = c("Muscle-specific\ngenes", "Query"), col = c(color.tissues[3], rgb(0.85, 0.85, 0.85, 0.5))) })
        output$Venn.Hypod <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[['hypod.genes']], multipleGenes(), names = c("Hypod.-specific\ngenes", "Query"), col = c(color.tissues[4], rgb(0.85, 0.85, 0.85, 0.5))) })
        output$Venn.Intest <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[['intest.genes']], multipleGenes(), names = c("Intest.-specific\ngenes", "Query"), col = c(color.tissues[5], rgb(0.85, 0.85, 0.85, 0.5))) })
        output$Venn.Ubiq <- renderPlot({ par(mar = c(0,0,0,0), oma = c(0,0,0,0)) ; a <- plot.2way.Venn(list.genes[['ubiq.genes']], multipleGenes(), names = c("Ubiq.-specific\ngenes", "Query"), col = c(color.tissues[33], rgb(0.85, 0.85, 0.85, 0.5))) })
        
    }
    
    # Generate tables to display and button to download them
    {
        output$downloadATAC.txt <- downloadHandler( "tissue-specific.ATAC-seq.dataset.txt", content = function(file) {write.table(atac.dt, file, quote = F, row = F, col = T, sep = '\t')} )
        output$downloadATAC.gff <- downloadHandler( "tissue-specific.ATAC-seq.dataset.gff", content = function(file) {
            infos <- paste0(
                "ID=", make.unique(paste(all$gene_name, all$regulatory_class, sep = '--'), sep = '_'), 
                ";Associated-gene-Name=", all$gene_name, 
                ";Associated-gene-WB=", all$WormBaseID, 
                ";domain=", all$domain
            )
            values <- paste0(
                ';color=', color.tissues[max.tissue.df$which.tissues], 
                '; =: : : : : : : : : : : : : : : : : : : : : : : : : ', 
                ';Ranked-tissues=', apply(max.tissue.df[,5:9], 1, function(x) paste0(x, collapse = ' / ')), 
                ';Ranked-tissue-ATACseq-TPM=', apply(max.tissue.df[,10:14], 1, function(x) paste0(round(x, 3), collapse = ' / ')), 
                ';Consecutive-ratios=', apply(max.tissue.df[,15:18], 1, function(x) paste0(round(x, 3), collapse = ' / ')), 
                ";Tissue-specificity=", gsub('\\.', '', max.tissue.df$which.tissues)
            )
            GFF <- cbind(
                as.character(all$chr), 
                rep('Ahringer-tissue-spe-REs', 
                times = nrow(all)), 
                as.character(all$regulatory_class), 
                all$start, 
                all$stop, 
                rep(".", times=nrow(all)), 
                ifelse(is.null(all$strand), ".", all$strand), 
                rep(".", times=nrow(all)), 
                paste0(infos, values)
            )
            write.table(rbind(c("##gffTags=on\n##displayName=Name\n##gff-version 3", rep(" ", 8)), GFF), file, quote = F, row = F, col = F, sep = '\t')
        } )
        output$downloadLCAP.txt <- downloadHandler( "tissue-specific.RNA-seq.dataset.txt", content = function(file) {write.table(lcap.dt, file, quote = F, row = F, col = T, sep = '\t')}  )
        output$downloadLCAP.gff <- downloadHandler( "tissue-specific.RNA-seq.dataset.gff", content = function(file) {
            infos <- paste0(
                "ID=", genes.gtf$gene_id,
                ";domain=", genes.gtf$domain
            )
            values <- paste0(
                ';color=', c(color.tissues[1], "#081d89", "#081d89", color.tissues[2:36])[max.tissue.df.LCAP$which.tissues], 
                '; =: : : : : : : : : : : : : : : : : : : : : : : : : ',
                ';Ranked-tissues=', apply(max.tissue.df.LCAP[,5:9], 1, function(x) paste0(x, collapse = ' * ')), 
                ';Ranked-tissue-RNAseq-TPM=', apply(max.tissue.df.LCAP[,grep('cov', colnames(max.tissue.df.LCAP))], 1, function(x) paste0(round(x, 3), collapse = ' * ')), 
                ';Tissue-zscore=', apply(round(t(apply(max.tissue.df.LCAP[,grepl('max.tissue.cov', colnames(max.tissue.df.LCAP))], 1, scale)), 2), 1, function(x) { paste(x, collapse = ' * ') } ),
                ';Consecutive-ratios=', apply(max.tissue.df.LCAP[,15:18], 1, function(x) paste0(round(x, 3), collapse = ' / ')), 
                ";Tissue-specificity=", gsub('\\.', '', max.tissue.df.LCAP$which.tissues)
            )
            GFF <- cbind(
                as.character(as.character(seqnames(genes.gtf))),
                rep('WormBase', times = length(genes.gtf)),
                as.character(genes.gtf$gene_biotype),
                start(genes.gtf),
                end(genes.gtf),
                rep(".", length(genes.gtf)),
                as.character(strand(genes.gtf)),
                rep(".", length(genes.gtf)),
                paste0(infos, values)
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
    
    # Prevent timeout (https://support.dominodatalab.com/hc/en-us/articles/360015932932-Increasing-the-timeout-for-Shiny-Server)
    output$keepAlive <- renderText({
        req(input$count) 
        paste("Shiny app is live.")
    })
    
} #EOF

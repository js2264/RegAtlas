#=============================================================================#
#                                                                             #
#     USAGE: THIS IS NOT A FUNCTIONAL SCRIPT.                                 #
#     It is a part of a Shiny app                                             #
#                                                                             #
#     AUTHOR: Jacques SERIZAY                                                 #
#     CREATED: 2018/07/13                                                     #
#     REVISION: 2018/10/05                                                    #
#                                                                             #
#     Must be in ./dev/ folder                                                #
#                                                                             #
#=============================================================================#

## Load data and libraries -----------------------------------------------------

load('../dashboard.Ahringer/releases/dev/data/minimal-data.RData')
source('../dashboard.Ahringer/releases/dev/R/custom_R_functions.R')
require(shiny)
require(shinydashboard)
require(shinyWidgets)
require(shinycssloaders)
require(shinyBS)
require(magrittr)
require(GenomicRanges)
require(devtools)
require(DT)
require(urltools)
require(htmltools)
require(httr)
require(gProfileR)
require(pheatmap)
require(d3heatmap)
require(RColorBrewer)
require(apputils)
require(venneuler)

# Define variables
colorGO <- c("#52c23164","#20109564","#ef690664", "white")
names(colorGO) <- c("MF", "BP", "CC", "kegg")
all.deconv <- cbind(
    all, 
    ATAC, 
    max.tissue.df[, grepl('max.tissue|ratio', colnames(max.tissue.df))]
) %>% 
    dplyr::mutate(uniqueWormBaseID = strsplit(as.character(WormBaseID),',')) %>% 
    tidyr::unnest(uniqueWormBaseID)
atac.dt <- cbind(
    all[, grepl('chr|start|stop|gene_name|regulatory_class|domain|which.tissues', colnames(all))], 
    round(ATAC, 3), 
    max.tissue.df[, grepl('max.tissue$', colnames(max.tissue.df))]
)
row.names(atac.dt) <- all$coords
colnames(atac.dt)[colnames(atac.dt) == 'gene_name'] <- "geneID"
colnames(atac.dt)[grep(paste(order.tissues, collapse = '|'), colnames(atac.dt))] <- paste0(
    order.tissues[1:length(grep(paste(order.tissues, collapse = '|'), colnames(atac.dt)))], 
    '_TPM'
)
atac.dt$regulatory_class <- factor(atac.dt$regulatory_class)
atac.dt$domain <- factor(atac.dt$domain)
lcap.dt <- cbind(
    as.data.frame(genes.gtf)[,c(1:3,5,10,11,13,16,18)], 
    round(LCAP, 3),
    max.tissue.df.LCAP[, grepl('max.tissue$', colnames(max.tissue.df.LCAP))]
)
colnames(lcap.dt)[1:3] <- c('chr', 'start', 'stop')
colnames(lcap.dt)[colnames(lcap.dt) == 'gene_id'] <- "WormBaseID"
colnames(lcap.dt)[colnames(lcap.dt) == 'gene_name'] <- "geneID"
lcap.dt$gene_biotype <- factor(lcap.dt$gene_biotype)
lcap.dt$domain <- factor(lcap.dt$domain)
hypod.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Hypod.']
neurons.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Neurons']
germline.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Germline']
muscle.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Muscle']
intest.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Intest.']
list.genes <- list(hypod.genes, neurons.genes, germline.genes, muscle.genes, intest.genes)
names(list.genes) <- c("germline.genes", "neurons.genes", "muscle.genes", "hypod.genes", 'intest.genes')
colors.decimals <- c(0.10395294, 0.3906374, 0.1192665, 0.14261010, 0.14226132, 0.13421772)
link <- "http://ahringerlab.com/JBrowse-1.12.5/index.html?data=data%2Fjson%2Fce11&loc=chrIII&tracks=genes%2Cregulatory_elements%2Chypod.atac%2Cneurons.atac%2Cgonad.atac%2Cmuscle.atac%2Cintest.atac%2Chypod.lcap.fwd%2Cneurons.lcap.fwd%2Cgonad.lcap.fwd%2Cmuscle.lcap.fwd%2Cintest.lcap.fwd%2Chypod.lcap.rev%2Cneurons.lcap.rev%2Cgonad.lcap.rev%2Cmuscle.lcap.rev%2Cintest.lcap.rev%2Ctranscripts&highlight=&menu=1&nav=1&tracklist=1&overview=1"
NCLUST_LCAPdev <- 5
NCLUST_LCAP <- 5
NCLUST_ATAC <- 5
version <- gsub("../dashboard.Ahringer_", "", list.dirs('..', recursive = F)[length(list.dirs('..', recursive = F)) - 1])

# Function to generate an URL to visit jserizay.site/JBrowse
getURL <- function (chr, start, end, release = "1.12.5", 
                    tracks = c("genes", "regulatory_elements", "hypod.atac", "neurons.atac", "gonad.atac", "muscle.atac", "intest.atac", "hypod.lcap.fwd", "neurons.lcap.fwd", "gonad.lcap.fwd", "muscle.lcap.fwd", "intest.lcap.fwd", "hypod.lcap.rev", "neurons.lcap.rev", "gonad.lcap.rev", "muscle.lcap.rev", "intest.lcap.rev", "transcripts"),
                    show_menu = TRUE, show_navigation = TRUE, show_tracklist = TRUE, show_overview = TRUE)
{
  baseurl <- paste0("http://ahringerlab.com/JBrowse-", release, "/index.html")
  range <- if (missing(start) || missing(end)) {""} else { paste0("%3A", paste0( round(start - 0.25 * (end - start + 1)), "..", round(end + 0.25 * (end - start + 1)) )) }
  tracks <- paste(unique(tracks), collapse = "%2C")
  menu <- if (show_menu) {"&menu=1"} else {"&menu=0"}
  navigation <- if (show_navigation) {"&nav=1"} else {"&nav=0"}
  tracklist <- if (show_tracklist) {"&tracklist=1"} else {"&tracklist=0"}
  overview <- if (show_overview) {"&overview=1"} else {"&overview=0"}
  url <- param_set(baseurl, key = "data", value = "data%2Fjson%2Fce11")
  url <- param_set(url, key = "loc", value = paste0(chr, range))
  url <- param_set(url, key = "tracks", value = tracks)
  return(paste0(url, "&highlight=", navigation, tracklist, overview))
}

# Function to input several genes separated by newlines
textareaInput <- function(inputId, label, value = "", placeholder = "", rows = 2)
{
  tagList(
    div(strong(label), style="margin-top: 5px;"),
    tags$style(type="text/css", "textarea {width:100%; margin-top: 5px;}"),
    tags$textarea(id = inputId, placeholder = placeholder, rows = rows, value))
}

# Function to bypass idle time of Shiny Server
inactivity <- "function idleTimer() {
  var t = setTimeout(logout, 500000000);
  window.onmousemove = resetTimer; // catches mouse movements
  window.onmousedown = resetTimer; // catches mouse movements
  window.onclick = resetTimer;     // catches mouse clicks
  window.onscroll = resetTimer;    // catches scrolling
  window.onkeypress = resetTimer;  //catches keyboard actions

  function logout() {
    window.close();  //close the window
  }

  function resetTimer() {
    clearTimeout(t);
    t = setTimeout(logout, 500000000);  // time is in milliseconds (1000 is 1 second)
  }
}
idleTimer();"


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
    output$REs.table <- renderDataTable({
        if (all(infos.gene()$Associated.REs == 0)) { datatable(matrix(0)) } else {
            datatable(
                setNames(
                    cbind(
                        infos.gene()$Associated.REs[c(1:4, 6)],
                        round(infos.gene()$Associated.REs.tissue, 1)
                    ),
                    c("Chr", "Start", "Stop", "Regulatory Class", "Enriched in tissue(s)", "Cov. Germline (YA)", "Cov. Neurons (YA)", "Cov. Muscle (YA)", "Cov. Hypod. (YA)", "Cov. Intest. (YA)")
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

## Define each dashbord tab individually ---------------------------------------------------------------------------

# Single-gene entry
TAB1 <- tabItem(
    tabName = 'genelookup',
    fluidPage(
        ## Row 1: Gene entry
        fluidRow(
            column(width = 8, {
                searchInput(
                    inputId = "searchGene",
                    label = "Enter a single gene name (WormBase ID or locus ID):",
                    value = 'hlh-1',
                    btnSearch = icon("search", lib = "glyphicon"),
                    btnReset = icon("remove", lib = "glyphicon"),
                    width = "95%"
                )
            } )
        ),
        fluidRow(
          column(width = 4, { 
              fluidRow(
                  htmlOutput("geneInfos", height = '200px') %>% withSpinner(type = 6, color = "#421500", size = 0.5),
                  br(),
                  actionBttn("WBdescr", label = "Get gene description", icon = icon("search", lib = "font-awesome"), size = "sm", style = "minimal"),
                  bsModal("WBDESCR", "WormBase gene description", "WBdescr", size = "large", { htmlOutput("geneDescr") %>% withSpinner(type = 6, color = "#421500", size = 0.5) } )
              )
          } ),
          column(width = 2, { 
            fluidRow(
              br(),
              downloadBttn("downloadINFOS", label = "Download full gene report (.txt file)", size = "sm", style = "fill", color = "primary"),
              br(),
              br(),
              actionBttn("switchToGenome", label = "Go to Genome Browser", icon = icon("area-chart", lib = "font-awesome"), size = "sm", style = "fill", color = "primary"),
              br(),
              br(),
              actionBttn("WBlink", label = "View gene in WormBase", icon = icon("book", lib = "font-awesome"), size = "sm", style = "fill", color = "primary"),
              bsModal("openWB", "WormBase gene entry", "WBlink", size = "large", htmlOutput("Link"))
            )
          } )
        ),
        br(),
        hr(),
        br(),
        
        ## Row 2: OUTPUT LCAPdev and LCAP-tissues graphs ,as well as gene infos.
        h4("Temporal and spatial gene expression profiles"),
        br(),
        fluidRow(
            column(width = 5, { plotOutput("Expr.plots_dev", height = '300px')  }),
            column(width = 5, { plotOutput("Expr.plots_tis", height = '300px')  })
        ),
        br(),
        hr(),
        br(),
        
        ## Row 3: OUTPUT table of associated REs table
        h3('Table of associated regulatory elements (REs)'),
        fluidRow( dataTableOutput("REs.table") ),
        br(),
        hr(),
        br(),

        ## Row 4: OUTPUT tSNE plots (ATAC / LCAP)
        #h4('t-SNE plots of promoters (left) and genes (right)'),
        #br(),
        #h5("Black dots represent the gene and its associated promoters."),
        #br(),
        #fluidRow(
        #    column(width = 5, { plotOutput("tSNE.plots_genes")  }),
        #    column(width = 2, { plotOutput("tSNE.plots_legend")  }),
        #    column(width = 5, { plotOutput("tSNE.plots_proms")  })
        #),
        #br(),
        #hr(),
        br()
    )
)

# Multiple-genes entry
TAB2 <- tabItem(
    tabName = 'geneslookup',
    fluidPage(

        ## Row 1: Multiple genes entry
        fluidRow(
            column(width = 2, {
                textAreaInput(
                    inputId = "searchMulitpleGenePaste",
                    label = h4("Paste gene names:"),
                    value = "",
                    placeholder = 'Paste here (one gene per line)\n\n* can be used to find patterns (e.g. nhr-*)',
                    rows = 8
                )
            }), 
            column(width = 2, h4("And / Or")),
            column(width = 2, { 
                tags$div(class = "multicol2", checkboxGroupInput("checkGroupGeneClasses", label = h4("Select a class of genes to query:"), selected = "neurons.genes", choices = list(
                    "Hypodermis-enriched genes" = "hypod.genes", 
                    "Neurons-enriched genes" = "neurons.genes", 
                    "Germline-enriched genes" = "germline.genes", 
                    "Muscle-enriched genes" = "muscle.genes", 
                    "Intestine-enriched genes" = "intest.genes"
                )))
            }),
            column(width = 2, h4("And / Or")),
            column(width = 4, {
                fluidRow(
                    actionBttn("MoreChoicesGenesList", label = "Choose genes based on their promoters", icon = icon("filter", lib = "font-awesome"), size = "sm", style = "fill", color = "primary"),
                    textOutput("multipleGenesPromsGroupsLength")
                )
            } ),
            bsModal(
                id = "moreChoices", 
                title = h3("Select groups of genes with only:"), 
                trigger = "MoreChoicesGenesList", 
                size = "large", 
                fluidPage(
                    checkboxGroupInput(
                        inputId = "checkGroupGeneClasses_2", 
                        label = "", 
                        selected = NULL, 
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
                            "Gonad & Muscle & Intestine-specific promoter(s)",
                            "Neurons & Gonad & Muscle & Intestine-specific promoter(s)",
                            "Hypodermis & Gonad & Muscle & Intestine-specific promoter(s)",
                            "Hypodermis & Neurons & Gonad & Intestine-specific promoter(s)",
                            "Hypodermis & Neurons & Gonad & Muscle-specific promoter(s)",
                            "Soma-specific promoter(s)",
                            "Ubiquitous promoter(s)"
                        ),
                        choiceValues = order.tissues[c(1:27, 29, 30, 32:33)]
                    ),
                    br()
                )
            )
        ),
        fluidRow(
            fluidRow(
                column(width = 4, ""),
                column(width = 4, { actionBttn("getList", label = "Perform analysis", icon = icon("search", lib = "font-awesome"), style = 'bordered', color = "primary", block = T) } ),
                column(width = 4, "")
            ),
            fluidRow(
                column(width = 4, ""),
                column(width = 2, { textOutput("multipleGenesLength") %>% withSpinner(type = 6, color = "#421500", size = 0.5) }),
                column(width = 2, { actionBttn("resetGenes", label = "Reset genes query", size = "xs", style = "fill") }),
                column(width = 4, "")
            ),
            br()
        ),
        br(),
        br(),
        br(), 
        downloadBttn("downloadGenesListGFF", label = "Download detailed report of input genes and associated REs (GFF file, IGV friendly)", size = "sm", style = "fill", color = "primary", block = F),
        br(),
        hr(),
        br(),

        ## Row 2a: OUTPUT Venn diagrms
        fluidRow(
            h4("Intersection of genes query with tissue-enriched genes"),
            column(width = 2, { plotOutput("Venn.Germline")  }),
            column(width = 2, { plotOutput("Venn.Neurons")  }),
            column(width = 2, { plotOutput("Venn.Muscle")  }),
            column(width = 2, { plotOutput("Venn.Hypod")  }),
            column(width = 2, { plotOutput("Venn.Intest")  })
        ),
        br(),
        hr(),
        br(),

        ## Row 2b: OUTPUT HMs.plot
        fluidRow(
            column(width = 2, { selectizeInput("colorScale_LCAPdev", "Choose a color scale (dev. RNA-seq): ", choices = rownames(RColorBrewer::brewer.pal.info), selected = "Spectral", multiple = FALSE) } ),
            column(width = 2, { dropdownButton(
                circle = TRUE, status = "warning", icon = icon("gear"),
                tooltip = tooltipOptions(title = "Options for developmental RNA-seq heatmap"),
                tags$h4("Options for developmental RNA-seq heatmap"),
                switchInput(
                    "clusterFUN_LCAPdev", 
                    label = "Type of clustering", 
                    value = F, 
                    onLabel = "Hierarchical", 
                    offLabel = "K-means", 
                    onStatus = 'success', 
                    offStatus = 'primary', 
                    labelWidth = 100, 
                    size = 'small', 
                    inline = T
                ),
                numericInput("NCLUST_LCAPdev", label = "Number of k-means clusters", value = 4),
                switchInput(
                    "colorScale_doRev_LCAPdev", 
                    label = "Reverse color scale?", 
                    value = T, 
                    onLabel = "Yes", 
                    offLabel = "No", 
                    onStatus = 'success', 
                    offStatus = 'primary', 
                    labelWidth = 150, 
                    width = "400px", 
                    size = 'small', 
                    inline = T
                ),
                switchInput(
                    "LCAPdev_TPMZscore", 
                    label = "Which values to plot (dev. RNA-seq)?", 
                    value = F, 
                    onLabel = "TPM", 
                    offLabel = "Z-score", 
                    onStatus = 'primary', 
                    offStatus = 'warning',
                    labelWidth = 200, 
                    width = "400px", 
                    size = 'small', 
                    inline = T
                )
            ) } ),
            column(width = 2, { selectizeInput("colorScale_LCAP", "Choose a color scale (RNA-seq): ", choices = rownames(RColorBrewer::brewer.pal.info), selected = "Spectral", multiple = FALSE) } ),
            column(width = 2, { dropdownButton(
                circle = TRUE, status = "warning", icon = icon("gear"),
                tooltip = tooltipOptions(title = "Options for tissue RNA-seq heatmap"),
                tags$h4("Options for tissue RNA-seq heatmap"),
                switchInput(
                    "clusterFUN_LCAP", 
                    label = "Type of clustering", 
                    value = T, 
                    onLabel = "Hierarchical", 
                    offLabel = "K-means", 
                    onStatus = 'success', 
                    offStatus = 'primary', 
                    labelWidth = 100, 
                    size = 'small', 
                    inline = T
                ),
                numericInput("NCLUST_LCAP", label = "Number of k-means clusters", value = 5),
                switchInput(
                    "colorScale_doRev_LCAP", 
                    label = "Reverse color scale?", 
                    value = T, 
                    onLabel = "Yes", 
                    offLabel = "No", 
                    onStatus = 'success', 
                    offStatus = 'primary', 
                    labelWidth = 150, 
                    width = "400px", 
                    size = 'small', 
                    inline = T
                ),
                switchInput(
                    "LCAP_TPMZscore", 
                    label = "Which values to plot (RNA-seq)?", 
                    value = F, 
                    onLabel = "TPM", 
                    offLabel = "Z-score", 
                    onStatus = 'primary', 
                    offStatus = 'warning',
                    labelWidth = 200, 
                    width = "400px", 
                    size = 'small', 
                    inline = T
                )
            ) } ),
            column(width = 2, { selectizeInput("colorScale_ATAC", "Choose a color scale (ATAC-seq): ", choices = rownames(RColorBrewer::brewer.pal.info), selected = "YlOrBr", multiple = FALSE) } ),
            column(width = 2, { dropdownButton(
                circle = TRUE, status = "warning", icon = icon("gear"), right = T,
                tooltip = tooltipOptions(title = "Options for tissue ATAC-seq heatmap"),
                tags$h4("Options for tissue ATAC-seq heatmap"),
                switchInput(
                    "clusterFUN_ATAC", 
                    label = "Type of clustering", 
                    value = T, 
                    onLabel = "Hierarchical", 
                    offLabel = "K-means", 
                    onStatus = 'success', 
                    offStatus = 'primary', 
                    labelWidth = 100, 
                    size = 'small', 
                    inline = T
                ),
                numericInput("NCLUST_ATAC", label = "Number of k-means clusters", value = 5),
                switchInput(
                    "colorScale_doRev_ATAC", 
                    label = "Reverse color scale?", 
                    value = F, 
                    onLabel = "Yes", 
                    offLabel = "No", 
                    onStatus = 'success', 
                    offStatus = 'primary', 
                    labelWidth = 150, 
                    width = "400px", 
                    size = 'small', 
                    inline = T
                ),
                switchInput(
                    "ATAC_TPMZscore", 
                    label = "Which values to plot (ATAC-seq)?", 
                    value = T, 
                    onLabel = "TPM", 
                    offLabel = "Z-score", 
                    onStatus = 'primary', 
                    offStatus = 'warning',
                    labelWidth = 200, 
                    width = "400px", 
                    size = 'small', 
                    inline = T
                )
            ) } )
        ),
        fluidRow(
            column(width = 4, { plotOutput("HMs.plot_LCAPdev")  }),
            column(width = 4, { plotOutput("HMs.plot_LCAP")  }),
            column(width = 4, { plotOutput("HMs.plot_ATAC")  })
        ),
        br(),
        fluidRow(
            column(width = 4, { downloadBttn("downloadHM_LCAPdev", label = "Download heatmap", size = "sm", style = "fill", color = "primary", block = F) }),
            column(width = 4, { downloadBttn("downloadHM_LCAP", label = "Download heatmap", size = "sm", style = "fill", color = "primary", block = F) }),
            column(width = 4, { downloadBttn("downloadHM_ATAC", label = "Download heatmap", size = "sm", style = "fill", color = "primary", block = F) })
        ),
        br(),
        hr(),
        br(),

        ## Row 3: OUTPUT GO.plot
        fluidRow(
            column(width = 2, {
                fluidRow(
                    actionBttn("runGO", label = "Perform GO analysis", icon = icon("search", lib = "font-awesome"), style = 'bordered', color = "primary"), 
                    br(),
                    br(),
                    dropdownButton(
                        tags$h4("Select GO databses"),
                        checkboxGroupInput("checkGroupGOs", label = "", 
                            choices = list("MF (Molecular Function)" = "MF", "BP (Biological Process)" = "BP", "CC (Cellular Component)" = "CC", "kegg (KEGG pathway)" = "keg"),
                            selected = c("MF", "BP", "CC", "keg")
                        ),
                        br(),
                        tags$h4("Apply filtering?"),
                        selectInput(
                            inputId = 'hierarchyFiltering', 
                            label = '', 
                            choices = c("none", "moderate", "strong"), 
                            selected = "moderate"
                        ),
                        circle = TRUE, status = "warning", icon = icon("gear"),
                        tooltip = tooltipOptions(title = "Set filtering")
                    )
                )
            } ),
            column(width = 7, { plotOutput("GO.plot") %>% withSpinner(type = 6, color = "#421500", size = 0.5) } ),
            column(width = 2, { fluidRow(
                br(),
                downloadBttn("downloadGO", label = "Download full GO report (.txt file)", size = "sm", style = "fill", color = "primary", block = T)
            ) } )
        ),
        br(),
        hr(),
        br(),

        ## Row 4: DISPLAY GENES LIST
        column( width = 2, { htmlOutput("genesList") } ),
        br()
    )
)

# Genome browser
TAB3 <- tabItem(
    tabName = 'browser',
    fluidPage( 
        fluidRow(
            column(width = 6, downloadBttn("downloadBWLCAP", label = "Download all the tissue-specific RNA-seq tracks (bigwig format)", size = "sm", style = "fill", color = "primary", block = F) ),
            column(width = 6, downloadBttn("downloadBWATAC", label = "Download all the tissue-specific ATAC-seq tracks (bigwig format)", size = "sm", style = "fill", color = "primary", block = F) )
        ),
        br(),
        fluidRow( { 
            #JbrowseOutput("jbrowser", height = "100%") ### THIS ONE IS GOOD!!! (IT DOESN'T MATTER)
            htmlOutput("jbrowser")
        } )
    )
)

# Download data
TAB4 <- tabItem(
    tabName = 'download',
    fluidPage(
        # Download buttons
        h2("Download our tissue-specific datasets"),
        fluidRow(
            column(width = 5, { fluidRow(
                downloadBttn("downloadATAC.txt", label = "Download the entire tissue-specific ATAC-seq dataset (txt format, Excel friendly)", size = "sm", style = "fill", color = "primary", block = T),
                br(),
                br(),
                downloadBttn("downloadATAC.gff", label = "Download the entire tissue-specific ATAC-seq dataset (GFF format, IGV friendly)", size = "sm", style = "fill", color = "primary", block = T)
            ) } ),
            column(width = 1, p(" ")),
            column(width = 5, { fluidRow(
                downloadBttn("downloadLCAP.txt", label = "Download the entire tissue-specific RNA-seq dataset (txt format, Excel friendly)", size = "sm", style = "fill", color = "primary", block = T),
                br(),
                br(),
                downloadBttn("downloadLCAP.gff", label = "Download the entire tissue-specific RNA-seq dataset (GFF format, IGV friendly)", size = "sm", style = "fill", color = "primary", block = T)
            ) } )
        ),
        hr(),

        # Browse ATAC table
        h4("Navigate ATAC-seq data"),
        fluidRow( column(width = 12, { dataTableOutput("atac.table") } ) ),
        hr(),

        # Browse LCAP table
        h4("Navigate RNA-seq data"),
        fluidRow( column(width = 12, { dataTableOutput("lcap.table") } ) )
    )
)

# Contact us
TAB5 <- tabItem(
    tabName = 'contact',
    fluidPage(
        fluidRow(
            column(width = 6, {
                HTML(
                    '
                    <div class="card">
                      <img src="http://ahringerlab.com/assets/img/ahringer-group-2017.jpg" alt="Lab" style="height: 284.16px">
                      <h1>Ahringer Lab</h1>
                      <br/>
                      <p class="cardtitle"> </p>
                      <p>Gurdon Institute, UK</p>
                      <p>Cambridge University, UK</p>
                      </br>
                      <br/>
                      <p><buttoncard><a style="text-decoration: none;font-size: 22px; color: white;" href="http://www.ahringer.group.gurdon.cam.ac.uk/" target="_blank">Visit us</a></button></p>
                    </div>
                    '
                )
            }),
            column(width = 6, {
                HTML(
                    '
                    <div class="card">
                      <img src="http://ahringerlab.com/assets/img/JS.jpg" alt="JS" style="height: 284.16px">
                      <h1>Jacques Serizay</h1>
                      <a itemprop="sameAs" href="https://github.com/js2264" target="_blank">
                          <span class="fa-stack fa-lg">
                              <i class="fa fa-circle fa-stack-2x" style="color:black"></i>
                              <i class="fa fa-github fa-stack-1x fa-inverse"></i>
                          </span>
                      </a>
                      <a itemprop="sameAs" href="https://www.linkedin.com/in/jacques-serizay-55103460/" target="_blank">
                          <span class="fa-stack fa-lg">
                              <i class="fa fa-circle fa-stack-2x" style="color:black"></i>
                              <i class="fa fa-linkedin fa-stack-1x fa-inverse"></i>
                          </span>
                      </a>
                      <a itemprop="sameAs" href="https://scholar.google.co.uk/citations?user=e5QTBIAAAAAJ" target="_blank">
                          <span class="fa-stack fa-lg">
                              <i class="fa fa-circle fa-stack-2x" style="color:black"></i>
                              <i class="fa fa-google fa-stack-1x fa-inverse"></i>
                          </span>
                      </a>
                      <br/>
                      <br/>
                      <p class="cardtitle">PhD candidate</p>
                      <p class="cardtitle">Developer/Maintenance of JABrowse</p>
                      <br/>
                      <p>Gurdon Institute, UK</p>
                      <p>Cambridge University, UK</p>
                      <br/>
                      <br/>
                      <p><buttoncard><a style="text-decoration: none;font-size: 22px; color: white;" href="mailto:js2264@cam.ac.uk" target="_blank">Contact</a></button></p>
                    </div>
                    '
                )
            })
        )
    )
)

## Finalise UI ----------------------------------------------------------------------------------------------------

SIDEBAR <- sidebarMenu(
	id = "tabs", 
    menuItem("Look-up gene", tabName = "genelookup", icon = icon("ellipsis-h", lib = "font-awesome")),
    menuItem("Look-up multiple genes", tabName = "geneslookup", icon = icon("ellipsis-h", lib = "font-awesome")),
    menuItem("Genome browser", tabName = "browser", icon = icon("area-chart", lib = "font-awesome")),
    menuItem("Explore/Download datasets", tabName = "download", icon = icon("download", lib = "font-awesome")),
	menuItem("Contact us", tabName = "contact", icon = icon("envelope-open", lib = "font-awesome")),
    sidebarSearchForm(textId = "quickGene", buttonId = "quickSearch", label = "Quick gene search..."), 
    tags$footer(
        img(src = "http://ahringerlab.com/assets/img/sidebar-img_150x150.png", alt = "", style = "
            color: #b8c7ce;
            padding: 0px 40px 30px 40px;
            z-index: 1000;
        "),
        br(),
        HTML(paste("Ahringer lab -", icon("copyright", lib = "font-awesome"), "2019")), 
        align = "left",
        style = "
            position: fixed;
            bottom: 0;
            text-align: center;
            color: #b8c7ce;
            padding: 0px 0px 10px 0px;
            left: 0px;
            z-index: 1000;
        "
    )
)

BODY <- tabItems(TAB1, TAB2)

shinyUI <- dashboardPage(
    dashboardHeader(
        title = paste0("Ahringer lab C. elegans tissue-specific database ", version), titleWidth = 600,
        tags$li(
            a(
                href = 'http://ahringerlab.com',
                img(src = 'http://ahringerlab.com/assets/img/favicon.ico', title = "Go back to main page", height = "30px"),
                style = "padding-top:10px; padding-bottom:10px;"
            ),
            class = "dropdown"
        )
    ),
    dashboardSidebar(SIDEBAR),
    dashboardBody(
        tags$script(inactivity),
        tags$head(includeCSS("../dashboard.Ahringer/releases/dev/assets/custom.css")),
        tags$head(tags$link(rel = "shortcut icon", href = "http://ahringerlab.com/assets/img/favicon.ico")),
        tags$head(tags$style(HTML('.left-side, .main-sidebar {position: fixed;} .navbar {position: fixed; right: 0px; left: 0px;} .main-header .logo {position: fixed} .content {padding: 65px 15px 15px 15px;}'))),        
        tags$head(tags$style(HTML('.logo .sidebar-toogle .navbar-navbar-static-top {position: fixed;}'))),
        tags$head(tags$style(HTML('.skin-blue .main-header .logo {background-color: #333;} .skin-blue .main-header .logo:hover {background-color: #333;}'))),
        tags$head(tags$style(HTML('.skin-blue .main-header .navbar {background-color: #333;} .skin-blue .main-header .navbar {background-color: #333;}'))),
        tags$head(tags$style(HTML('.skin-blue .main-header .navbar .sidebar-toggle {background-color: #333;} .skin-blue .main-header .navbar .sidebar-toggle:hover {background-color: #444;}'))),
        tags$head(tags$style(HTML('.skin-blue .left-side, .skin-blue .main-sidebar, .skin-blue .wrapper {background-color: #333;}'))),
        tags$head(tags$style(HTML('.main-header .sidebar-toggle {size:100px}'))),
        tags$head(tags$style(HTML('.main-header .sidebar-toggle:before {size:100px; content: "\\f0d9"}'))),
        tags$head(tags$style(HTML('.main-header .sidebar-toggle:after {size:100px; content: "\\f0da"}'))),
        tags$head(tags$style(HTML("hr {border-top: 1px dashed #b7b7b7;}"))),
        tags$head(tags$style(HTML(".btn {border-radius: 30px;} .btn:hover {transform: scale(1.05);}"))),
        tags$head(tags$style(HTML(".bttn-fill.bttn-primary {background: #ddd;color: #333;}.bttn-fill.bttn-sm {padding: 4px 10px;font-size: 16px;font-family: inherit;}.bttn {border-radius: 0px;border-color: #555;}"))),
        tags$head(tags$style(HTML(".bttn-bordered.bttn-primary {background: #fff;color: #600;border-radius: 5px; border-color: #600} .bttn-bordered.bttn-primary:hover, .bttn-bordered.bttn-primary:focus {background: #600;color: #fff;}"))),
        tags$head(tags$style(HTML(".bttn-minimal.bttn-default {background: #fff;color: #333;border-color: #fff;}"))),
        tags$head(tags$style(HTML("a {color: #333} a:hover, a:focus, a:active, a:visited {color: #1d89ff;}"))),
        tags$head(tags$style(HTML(".multicol2 {-webkit-column-count: 1; /* Chrome, Safari, Opera */ -moz-column-count: 1; /* Firefox */ column-count: 1;}"))),
        tags$head(tags$style(HTML(".multicol5 {-webkit-column-count: 5; /* Chrome, Safari, Opera */ -moz-column-count: 5; /* Firefox */ column-count: 5;}"))),
        tags$head(tags$style(HTML(".content-wrapper, .right-side { background-color: #FFF;}"))),
        tags$head(tags$style(HTML(".card { width: 80%;}"))),
        tags$head(tags$style(HTML(".img-round { width: 80%;}"))),
        BODY, 
        bsModal(
            id = "quickGENE", 
            title = "Quick gene view",
            trigger = '', 
            size = "large", 
            htmlOutput("quickResults")
        )
    )
)

## 
## 
shinyApp(ui = shinyUI, server = shinyServer)
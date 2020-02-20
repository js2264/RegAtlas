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

# Home tab
HOME <- tabItem(
    tabName = 'home',
    fluidPage(
        HTML("<h1 style = 'text-align: center;'><b>Welcome to the Ahringer lab <em>C. elegans</em> regulatory atlas</b></h1>"), 
        HTML("<h3 style = 'text-align: center;'>Explore and mine <em>C. elegans</em> developmental and tissue-specific regulatory elements and genes<h3>"), 
        br(), 
        br(), 
        fluidRow(
            column(width = 1),
            column(width = 2, style = "cursor: pointer", {
                uiOutput("img_link_genelookup")
            }),
            column(width = 2, style = "cursor: pointer", {
                uiOutput("img_link_geneslookup")
            }),
            column(id = 'browser_click', width = 2, style = "cursor: pointer", {
                uiOutput("img_link_browser")
            }),
            column(width = 2, style = "cursor: pointer", {
                uiOutput("img_link_download")
            }),
            column(width = 2, style = "cursor: pointer", {
                uiOutput("img_link_infos")
            }),
            column(width = 1)
        )
    )
)

# Single-gene entry
TAB1 <- tabItem(
    add_busy_bar(color = "red", height = '10px'),
    tabName = 'genelookup',
    fluidPage(
        
        ## Row 0: Gene entry
        fluidRow(
            h3('Single gene information'), 
            column(width = 3, {
                searchInput(
                    inputId = "searchGene",
                    label = "Look up a gene:",
                    # value = 'unc-120',
                    placeholder = 'e.g. WBGene00001251 or hlh-1', 
                    btnSearch = icon("search", lib = "glyphicon"),
                    btnReset = icon("remove", lib = "glyphicon"),
                    width = "95%"
                )
            } )
        )
    ),
    conditionalPanel(
        condition = "output.displayPanelsTab1 == '0'", 
        fluidPage(
            ## WRONG ENTRY ERROR MESSAGE
            fluidRow(
                column(width = 12, {
                    h4("Gene not found. Please enter a valid gene locus (e.g. hlh-1) or WormBaseID (e.g. WBGene00001948)", style = "color: FireBrick ; font-weight: bold")
                })
            )
        )
    ),
    conditionalPanel(
        condition = "output.displayPanelsTab1 == '1'", 
        fluidPage(
            
            ## Row 1: Gene infos
            fluidRow(
                column(width = 3, { 
                    fluidRow(
                        h4('Quick gene information:'), 
                        br(), 
                        htmlOutput("geneInfos", height = '200px') %>% withSpinner(proxy.height = '200px', type = 6, color = "#421500", size = 0.5)
                    )
                } ),
                column(width = 5, { 
                    fluidRow(
                        h4('Wormbase description:'), 
                        br(), 
                        htmlOutput("geneDescr") %>% withSpinner(proxy.height = '200px', type = 6, color = "#421500", size = 0.5)
                    )
                } ), 
                column(width = 1),
                column(width = 3, { 
                    fluidRow(
                        h4('Additional resources:'), 
                        br(),
                        downloadBttn("downloadINFOS", label = "Download full gene report (.txt file)", size = "sm", style = "simple", color = "primary"),
                        br(),
                        br(),
                        actionBttn("switchToGenome", label = "Go to Genome Browser", icon = icon("area-chart", lib = "font-awesome"), size = "sm", style = "simple", color = "primary"),
                        br(),
                        br(),
                        uiOutput("Link")
                    )
                } )
            ),
            br(),
            hr(),
            br(),
            
            ## Row 2: OUTPUT LCAPdev and LCAP-tissues graphs ,as well as gene infos.
            fluidRow(column(width = 12,
                div(style="display: inline-block;",tags$h4('Temporal and spatial gene expression profiles    ')),
                div(style="display: inline-block; width: 50px"),
                div(style="display: inline-block; width: 300px", popify(
                    actionBttn(
                        inputId = NA, 
                        color = "royal", 
                        icon = icon("question"), 
                        style = 'simple'
                    ), 
                    options=list(container="body"),
                    placement = 'bottom', 
                    "<p><b>These plots show gene expression in TPMs (Transcripts Per Millions) from nuclear RNA-seq.</b><br><br>Temporal gene expression values (left) are obtained from Janes et al., 2018.<br><br>Tissue-specific gene expression values (from YA, middle and right) are obtained from Serizay et al. (submitted).<br><br>The right plot shows gene expression values with background RNA contamination corrected by DSA. This normalization improves the quantification of gene expression for tissue-specific genes, and these values are recommended when looking at a tissue-specific gene. Otherwise, the uncorrected values (middle plot) are recommended for broadly expressed genes. </br></br>See the Information tab for more details.</p>"
                ))
            )),
            br(),
            br(),
            fluidRow(
                tipify(column(
                  width = 4,
                  plotOutput("Expr.plots_dev", height = '300px')
                ), placement = 'top', "<p>Gene expression across development, in TPMs (Transcripts Per Million) [Janes et al., 2018]</p>"),
                tipify(column(
                    width = 4,
                    plotOutput("Expr.plots_tis_uncorrected", height = '300px')
                ), placement = 'top', "<p>Gene expression in YA tissues, in TPMs (Transcripts Per Million) [Serizay et al. (submitted)]</p>"),
                tipify(column(
                    width = 4,
                    plotOutput("Expr.plots_tis", height = '300px')
                ), placement = 'top', "<p>Gene expression in YA tissues, in TPMs (Transcripts Per Million) with background RNA contamination corrected by DSA [Serizay et al. (submitted)]</p>")
            ),
            br(),
            hr(),
            br(),
            
            ## Row 3: OUTPUT table of associated REs table
            h4('Associated accessible elements and tissue-specific accessibility values'),
            br(),
            fluidRow( dataTableOutput("REs.table") ),
            br(),
            hr(),
            br(),
            br()
        )
    )
)

# Multiple-genes entry
TAB2 <- tabItem(
    tabName = 'geneslookup',
    fluidPage(
        ## Row 1: Multiple genes entry
        fluidRow(
            h3('Input a list of genes to determine: '), 
            h4('(1) their overlap with the gene expression classes defined in Serizay et al. (submitted);', style = 'padding-left: 50px;'),
            h4('(2) their expression in different young adult (YA) tissues;', style = 'padding-left: 50px;'),
            h4('(3) the tissue-specific accessibility of the associated accessible chromatin loci;', style = 'padding-left: 50px;'),
            h4('(4) the GO terms enriched in the gene list.', style = 'padding-left: 50px; padding-bottom: 30px;'),
            column(width = 3, fluidRow(
                textAreaInput(
                    inputId = "searchMulitpleGenePaste",
                    label = "Paste gene names:",
                    value = "",
                    placeholder = 'Paste here (one gene per line)\n\n* can be used to find patterns (e.g. nhr-*)',
                    rows = 10
                ), 
                actionBttn("MoreChoicesGenesList", label = "Examples of gene lists", icon = icon("filter", lib = "font-awesome"), size = "sm", style = "simple", color = "primary"),
                bsModal(
                    id = "moreChoices", 
                    title = h3("Select examples of gene lists:"), 
                    trigger = "MoreChoicesGenesList", 
                    size = "large", 
                    fluidPage(
                        checkboxGroupInput(
                            inputId = "checkGroupGeneClasses", 
                            label = "",
                            width = '70%', 
                            choiceNames = c(
                                "Germline-specific genes",
                                "Neurons-specific genes",
                                "Muscle-specific genes",
                                "Hypodermis-specific genes",
                                "Intestine-specific genes",
                                "Soma-specific genes",
                                "Ubiquitous genes",
                                "Genes with germline-specific promoter(s)",
                                "Genes with neurons-specific promoter(s)",
                                "Genes with muscle-specific promoter(s)",
                                "Genes with hypodermis-specific promoter(s)",
                                "Genes with intestine-specific promoter(s)",
                                "Genes with Neurons & Muscle-specific promoter(s)",
                                "Genes with Hypodermis & Intestine-specific promoter(s)",
                                "Genes with soma-specific promoter(s)",
                                "Genes with ubiquitous promoter(s)"
                            ),
                            choiceValues = 1:16
                        )
                    )
                ), 
                textOutput("multipleGenesPromsGroupsLength")
            )), 
            column(width = 3, style = "padding-left: 40px;", fluidRow(
                h4(' '), 
                br(), 
                br(), 
                actionBttn("getList", label = "Perform analysis", icon = icon("search", lib = "font-awesome"), style = 'bordered', color = "primary", block = FALSE), 
                p(style = "color: grey;", 'Long lists will take a while to be processed'),
                br(), 
                textOutput("multipleGenesLength"),
                actionBttn("resetGenes", label = "Reset query genes", size = "xs", style = "simple", color = 'primary')
            ))
        )
    ),
    conditionalPanel(
        condition = "output.displayPanelsTab2 == '1'", 
        fluidPage(
            br(),
            hr(),
            br(), 
            downloadBttn("downloadGenesListGFF", label = "Download annotated GFF file of input genes and associated accessible sites (GFF file, IGV friendly)", size = "sm", style = "simple", color = "primary", block = F),
            downloadBttn("downloadGenesListTXT", label = "Download detailed text report of input genes and associated accessible sites (txt file, Excel friendly)", size = "sm", style = "simple", color = "primary", block = F),
            br(),
            br(),
            hr(),
            br(),
            
            ## Row 2a: OUTPUT intersection diagrams
            fluidRow(
                h4("Intersection of query genes with adult gene expression classes from Serizay et al. (submitted)"),
                column(width = 12, { plotOutput("intersection_bars", height = '400px') %>% withSpinner(type = 6, color = "#421500", size = 0.5) })
            ),
            br(),
            hr(),
            br(),
            
            ## Row 2b: OUTPUT HMs.plot
            fluidRow(
                column(width = 4, { plotOutput("HMs.plot_LCAPdev", height = '500px') }),
                column(width = 4, { plotOutput("HMs.plot_LCAP", height = '500px') }),
                column(width = 4, { plotOutput("HMs.plot_ATAC", height = '500px') })
            ),
            br(),
            fluidRow(
                column(width = 1, { dropdownButton(
                    circle = TRUE, status = "warning", icon = icon("gear"),
                    tooltip = tooltipOptions(title = "Options for developmental RNA-seq heatmap"),
                    tags$h4("Options for developmental RNA-seq heatmap"),
                    selectizeInput(
                        "colorScale_LCAPdev", 
                        "Choose a color scale (dev. RNA-seq): ", 
                        choices = rownames(RColorBrewer::brewer.pal.info), 
                        selected = "Spectral", 
                        multiple = FALSE
                    ), 
                    switchInput(
                        "clusterFUN_LCAPdev", 
                        label = "Type of clustering", 
                        value = TRUE, 
                        onLabel = "Hierarchical", 
                        offLabel = "K-means", 
                        onStatus = 'success', 
                        offStatus = 'primary', 
                        labelWidth = 100, 
                        size = 'small', 
                        inline = TRUE
                    ),
                    numericInput("NCLUST_LCAPdev", label = "Number of k-means clusters", value = 4),
                    switchInput(
                        "colorScale_doRev_LCAPdev", 
                        label = "Reverse color scale?", 
                        value = TRUE, 
                        onLabel = "Yes", 
                        offLabel = "No", 
                        onStatus = 'success', 
                        offStatus = 'primary', 
                        labelWidth = 150, 
                        width = "400px", 
                        size = 'small', 
                        inline = TRUE
                    ),
                    switchInput(
                        "LCAPdev_TPMZscore", 
                        label = "Which values to plot (dev. RNA-seq)?", 
                        value = TRUE, 
                        onLabel = "TPM", 
                        offLabel = "Z-score", 
                        onStatus = 'primary', 
                        offStatus = 'warning',
                        labelWidth = 200, 
                        width = "400px", 
                        size = 'small', 
                        inline = TRUE
                    )
                ) } ),
                column(width = 3, { 
                    HTML('<button class="action-button bttn bttn-simple bttn-sm bttn-primary bttn-no-outline narrow" id="downloadHM_LCAPdev_bttn" type="button">
                        <i class="fa fa-download"></i>
                        <a id="downloadHM_LCAPdev" class="shiny-download-link" href="" target="_blank" download>Download heatmap</a>
                    </button>')
                }),
                column(width = 1, { dropdownButton(
                    circle = TRUE, status = "warning", icon = icon("gear"),
                    tooltip = tooltipOptions(title = "Options for tissue RNA-seq heatmap"),
                    tags$h4("Options for tissue RNA-seq heatmap"),
                    selectizeInput(
                        "colorScale_LCAP", 
                        "Choose a color scale (RNA-seq): ", 
                        choices = rownames(RColorBrewer::brewer.pal.info), 
                        selected = "Spectral", 
                        multiple = FALSE
                    ), 
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
                column(width = 3, { 
                    HTML('<button class="action-button bttn bttn-simple bttn-sm bttn-primary bttn-no-outline narrow" id="downloadHM_LCAPdev_bttn" type="button">
                        <i class="fa fa-download"></i>
                        <a id="downloadHM_LCAP" class="shiny-download-link" href="" target="_blank" download>Download heatmap</a>
                    </button>')
                }),
                column(width = 1, { dropdownButton(
                    circle = TRUE, status = "warning", icon = icon("gear"), right = T,
                    tooltip = tooltipOptions(title = "Options for tissue ATAC-seq heatmap"),
                    tags$h4("Options for tissue ATAC-seq heatmap"),
                    selectizeInput(
                        "colorScale_ATAC", 
                        "Choose a color scale (ATAC-seq): ", 
                        choices = rownames(RColorBrewer::brewer.pal.info), 
                        selected = "YlOrBr", 
                        multiple = FALSE
                    ), 
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
                ) } ),
                column(width = 3, { 
                    HTML('<button class="action-button bttn bttn-simple bttn-sm bttn-primary bttn-no-outline narrow" id="downloadHM_LCAPdev_bttn" type="button">
                        <i class="fa fa-download"></i>
                        <a id="downloadHM_ATAC" class="shiny-download-link" href="" target="_blank" download>Download heatmap</a>
                    </button>')
                })
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
                            checkboxGroupInput(
                                "checkGroupGOs", label = "", 
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
                        ),
                        br(), 
                        br(), 
                        downloadBttn("downloadGO", label = "Download full GO report (.txt file)", size = "sm", style = "simple", color = "primary", block = T)
                    )
                } ),
                column(width = 10, { plotOutput("GO.plot") } )
            ),
            br(),
            hr(),
            br(),
            
            ## Row 4: DISPLAY GENES LIST
            column( width = 2, { htmlOutput("genesList") } ),
            br()
        )
    )
)

# Genome browser
TAB3 <- tabItem(
    tabName = 'browser',
    fluidPage( 
        fluidRow(
            htmlOutput("jbrowser")
        )
    )
)

# Download data
TAB4 <- tabItem(
    tabName = 'download',
    fluidPage(
        fluidRow(
            column(3, h4("Download full tissue-specific datasets", style = 'display: inline-block')), 
            column(7, dropdownButton(
                    downloadBttn("downloadATAC.txt", label = "Download tissue-specific ATAC-seq dataset (txt format, Excel friendly)", size = "sm", style = "simple", color = "primary", block = TRUE), 
                    downloadBttn("downloadATAC.gff", label = "Download tissue-specific ATAC-seq dataset (GFF format, IGV friendly)", size = "sm", style = "simple", color = "primary", block = TRUE), 
                    downloadBttn("downloadBWATAC", label = "Download all zipped tissue-specific ATAC-seq tracks (bigWig format)", size = "sm", style = "simple", color = "primary", block = TRUE), 
                    downloadBttn("downloadLCAP.txt", label = "Download tissue-specific RNA-seq dataset (txt format, Excel friendly)", size = "sm", style = "simple", color = "primary", block = TRUE), 
                    downloadBttn("downloadLCAP.gff", label = "Download tissue-specific RNA-seq dataset (GFF format, IGV friendly)", size = "sm", style = "simple", color = "primary", block = TRUE),
                    downloadBttn("downloadBWLCAP", label = "Download all zipped tissue-specific RNA-seq tracks (bigWig format)", size = "sm", style = "simple", color = "primary", block = TRUE),
                    status = "warning", 
                    icon = icon("download")
                )
            )
        ),
        hr(),
        
        # Browse ATAC table
        h4("Search ATAC-seq data"),
        fluidRow( column(width = 12, { dataTableOutput("atac.table") } ) ),
        hr(),
        
        # Browse LCAP table
        h4("Search RNA-seq data"),
        fluidRow( column(width = 12, { dataTableOutput("lcap.table") } ) )
    )
)

# Information about the data
TAB5 <- tabItem(
    tabName = 'infos',
    fluidPage(
        fluidRow(
            column(width = 10, {
                HTML(
                    '
                    <p>
                    <h2>General Information</h2>
                        <h3>Data availability</h3>
                        The tissue-specific accessibility and expression data displayed in this application are from Serizay et al, submitted; raw data are available from GEO accession number GSE141213. Processed data are downloadable from the "Explore/Download datasets" of this application.
                        Developmental time course accessibility and expression data are from Janes et al, 2018 (DOI: 10.7554/eLife.37344e).
                        Tissue-specific expresssion at the L2 stage are from single cell sequencing data in Cao et al, 2017 (DOI: 10.1126/science.aam8940).
                        <h3>Genome version</h3>
                            Gene and regulatory element coordinates are in the WBcel235/ce11 version of the genome. Genome sequence and gene annotations were obtained from <a href = ftp://ftp.ensembl.org/pub/release-92/>Ensembl release 92</a>.</br>
                    </br>
                    </br>
                    <h2>Additional information</h2>
                        <h3>Background RNA contamination</h3>
                            We found that RNA from some highly expressed tissue-specific genes was unexpectedly detected in all sorted tissue samples. This appears to be due to contamination by cytoplasmic mature RNA released during nuclear isolation. This RNA, released from all of the tissues of the animal, may have adhered to the outside nuclei, and thus be carried over into all samples. We found that the carry-over is primarily detectable for highly expressed genes (e.g, muscle myosin <i>unc-54</i>).
                        <h3>Background removal</h3>
                            To estimate and remove background contamination in our RNA-seq samples, we used the DSA 1.0 package [<a href = https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-89>Zhong et al., BMC Bioinformatics 2013</a>] in R with default parameters and the linear regression method. The corrected gene expression values are displayed in the "Look-up gene" tab (third barplot).</br>
                            The DSA package estimates tissue-specific gene expression from mixed samples (e.g. complex tissues or samples with contamination) using profiles of "pure" genes known to be expressed in a single tissue. DSA normalisation improves profiles for tissue-specific genes, but introduce tissue biases for broadly expressed genes.</br>
                            </br>
                            The following genes were used as "pure" tissue-specific genes for the DSA-based background correction:</br>
                                - Germline: maph-1.2, prom-1, zim-3, htp-2, xnd-1; </br>
                                - Neuron: ncx-4, tsp-7, oig-8, Y106G6G.6, Y106G6G.2, C35E7.2, F32B4.5, ttr-39; </br>
                                - Muscle: B0379.1, srp-3, ttr-36, Y97E10AR.2, C13C4.7, tsp-11, lev-8; </br>
                                - Hypodermis: fip-5, osm-8, ptr-16, F36H9.5, nipi-4, R07E3.7; </br>
                                - Intestine: Y106G6H.1, F46A8.11, ugt-6, T02B5.3.</br>
                            </br>
                    <h2>App development</h2>
                        This app was built in R 3.5.1 "Feather Spray", with Shiny 1.1.0. The source code of this app is available on <a href = https://github.com/js2264/RegAtlas/>Github</a>.
                    </p>
                    </br>
                    </br>
                    '
                )
            })
        ), 
        fluidRow(column(width = 5, h2('Contact information'))),
        fluidRow(
            column(width = 3, {
                HTML(
                    '
                    <h3>Jacques Serizay</h3>
                    <br/>
                    <p class="cardtitle">Developer of RegAtlas</p>
                    <br/>
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
                    <a itemprop="sameAs" href="mailto:js2264@cam.ac.uk" target="_blank">
                    <span class="fa-stack fa-lg">
                    <i class="fa fa-circle fa-stack-2x" style="color:black"></i>
                    <i class="fa fa-envelope fa-stack-1x fa-inverse"></i>
                    </span>
                    </a>
                    '
                )
            }),
            column(width = 5, {
                HTML(
                    '
                    <h3>Ahringer lab</h3>
                    <br/>
                    <p class="cardtitle">Gurdon Institute, University of Cambridge, UK</p>
                    <br/>
                    <a itemprop="sameAs" href="http://www.ahringer.group.gurdon.cam.ac.uk/" target="_blank">
                    <span class="fa-stack fa-lg">
                    <i class="fa fa-circle fa-stack-2x" style="color:black"></i>
                    <i class="fa fa-home fa-stack-1x fa-inverse"></i>
                    </span>
                    </a>
                    '
                )
            })
        ),
        br(), 
        br(), 
        tags$head(
            HTML(
              "
              <script>
              var socket_timeout_interval
              var n = 0
              $(document).on('shiny:connected', function(event) {
              socket_timeout_interval = setInterval(function(){
              Shiny.onInputChange('count', n++)
              }, 15000)
              });
              $(document).on('shiny:disconnected', function(event) {
              clearInterval(socket_timeout_interval)
              });
              </script>
              "
            )
        ), 
        textOutput("keepAlive")
    )
)

# router <- make_router(
#   route("home", HOME),
#   route("genelookup", TAB1),
#   route("geneslookup", TAB2),
#   route("browser", TAB3),
#   route("download", TAB4),
#   route("infos", TAB5),
#   route("contact", TAB6)
# )

## Finalise UI ----------------------------------------------------------------------------------------------------

SIDEBAR <- sidebarMenu(
    id = "tabs", 
    menuItem("Home", tabName = "home", icon = icon("home", lib = "font-awesome")),
    menuItem("Look-up gene", tabName = "genelookup", icon = icon("search", lib = "font-awesome")),
    menuItem("Gene set analyses", tabName = "geneslookup", icon = icon("ellipsis-h", lib = "font-awesome")),
    menuItem("Genome browser", tabName = "browser", icon = icon("area-chart", lib = "font-awesome")),
    menuItem("Explore & download data", tabName = "download", icon = icon("download", lib = "font-awesome")),
    menuItem("Information & contact", tabName = "infos", icon = icon("info", lib = "font-awesome")),
    # sidebarSearchForm(textId = "quickGene", buttonId = "quickSearch", label = "Quick gene search..."), 
    tags$footer(
        img(src = "http://ahringerlab.com/assets/img/sidebar-img_150x150.png", alt = "", style = "
        color: #b8c7ce;
        padding: 0px 40px 30px 40px;
        z-index: 1000;
        "),
        br(),
        HTML(paste("Ahringer lab -", icon("copyright", lib = "font-awesome"), "2020")), 
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

BODY <- tabItems(HOME, TAB1, TAB2, TAB3, TAB4, TAB5)

shinyUI <- dashboardPage(
    dashboardHeader(
        title = HTML(paste0("Ahringer lab C. elegans regulatory atlas ", version)), titleWidth = 600,
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
        # ### Startup modal 
        # bsModal(
        #     id = 'startupModal', 
        #     title = NULL, 
        #     trigger = '', 
        #     size = 'large', 
        #     HTML('<img src="http://ahringerlab.com/assets/img/startup_picture.png" title="" style="max-width: 100%; max-height: 100%; display: block; ">')
        # ),
        # Quick search modal
        # bsModal(
        #     id = "quickGENE", 
        #     title = "Quick gene view",
        #     trigger = '', 
        #     size = "large", 
        #     htmlOutput("quickResults")
        # ), 
        tags$script(inactivity),
        tags$head(includeCSS("assets/custom.css")),
        tags$head(tags$link(rel = "shortcut icon", href = "http://ahringerlab.com/assets/img/favicon.ico")),
        tags$head(tags$style(HTML('.left-side, .main-sidebar {position: fixed;} .navbar {position: fixed; right: 0px; left: 0px;} .main-header .logo {position: fixed} .content {padding: 65px 15px 15px 15px;}'))),
        tags$head(tags$style(HTML('.logo .sidebar-toogle .navbar-navbar-static-top {position: fixed;}'))),
        tags$head(tags$style(HTML('.skin-blue .main-header .logo {background-color: #333;} .skin-blue .main-header .logo:hover {background-color: #333;}'))),
        tags$head(tags$style(HTML('.skin-blue .main-header .navbar {background-color: #333;} .skin-blue .main-header .navbar {background-color: #333;}'))),
        tags$head(tags$style(HTML('.skin-blue .main-header .navbar .sidebar-toggle {background-color: #333;} .skin-blue .main-header .navbar .sidebar-toggle:hover {background-color: #444;}'))),
        tags$head(tags$style(HTML('.skin-blue .left-side, .skin-blue .main-sidebar, .skin-blue .wrapper {background-color: #333;}'))),
        tags$head(tags$style(HTML('.main-header .sidebar-toggle {size:100px}'))),
        tags$head(tags$style(HTML('.main-header .sidebar-toggle:before {size:100px; content: ""}'))),
        tags$head(tags$style(HTML('.main-header .sidebar-toggle:after {size:100px; content: ""}'))),
        tags$head(tags$style(HTML("hr {border-top: 1px dashed #b7b7b7;}"))),
        tags$head(tags$style(HTML(".btn {border-radius: 30px;} .btn:hover {transform: scale(1.05);}"))),
        tags$head(tags$style(HTML(".bttn-fill.bttn-primary {background: #ddd;color: #333;}.bttn-fill.bttn-sm {padding: 4px 10px;font-size: 16px;font-family: inherit;}.bttn {border-radius: 0px;border-color: #555;}"))),
        tags$head(tags$style(HTML(".bttn-simple.bttn-primary {background: #ddd;color: #333;}.bttn-bordered.bttn-primary:hover, .bttn-bordered.bttn-primary:focus {background: #eee;color: #333;}.bttn-simple.bttn-sm {padding: 4px 10px;font-size: 16px;font-family: inherit;}.bttn {border-radius: 0px;border-color: #555;}"))),
        tags$head(tags$style(HTML(".bttn-bordered.bttn-primary {background: #fff;color: #600;border-radius: 5px; border-color: #600} .bttn-bordered.bttn-primary:hover, .bttn-bordered.bttn-primary:focus {background: #600;color: #fff;}"))),
        tags$head(tags$style(HTML(".bttn-minimal.bttn-default {background: #fff;color: #333;border-color: #fff;}"))),
        tags$head(tags$style(HTML("a {color: #333} a:hover, a:focus, a:active, a:visited {color: #1d89ff;}"))),
        tags$head(tags$style(HTML(".multicol2 {-webkit-column-count: 1; /* Chrome, Safari, Opera */ -moz-column-count: 1; /* Firefox */ column-count: 1;}"))),
        tags$head(tags$style(HTML(".multicol5 {-webkit-column-count: 5; /* Chrome, Safari, Opera */ -moz-column-count: 5; /* Firefox */ column-count: 5;}"))),
        tags$head(tags$style(HTML(".content-wrapper, .right-side { background-color: #FFF;}"))),
        tags$head(tags$style(HTML(".card { width: 80%;}"))),
        tags$head(tags$style(HTML(".img-round { width: 80%;}"))),
        tags$head(tags$style(HTML(".narrow { width: 120px; }"))),
        tags$head(tags$style(HTML(".popover{max-width: 600px;}"))),
        tags$head(tags$style("#startupModal .modal-dialog{ width: 90vw; height: 90vh}")),
        tags$head(tags$style("#startupModal .modal-content{ background-color: #fff; }")),
        tags$head(tags$style("#startupModal .modal-header{ height: 0px; visibility: hidden; padding: 0px; }")),
        tags$script(HTML("
            var openTab = function(tabName){
                $('a', $('.sidebar')).each(function() {
                    if(this.getAttribute('data-value') == tabName) {
                        this.click()
                    };
                });
            }
        ")), 
        BODY
    )
)

## 
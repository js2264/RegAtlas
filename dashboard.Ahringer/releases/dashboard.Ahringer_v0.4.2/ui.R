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
        ## Row 0: Gene entry
        fluidRow(
            column(width = 3, {
                searchInput(
                    inputId = "searchGene",
                    label = "Look up a gene:",
                    value = 'hlh-1',
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
                    h4("Gene not found. Please enter a valid gene locus (e.g. hlh-1) or WormBaseID (e.g. WBGene00001948)", style = "padding-left: 40px; color: FireBrick ; font-weight: bold")
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
            h4("Temporal and spatial gene expression profiles"),
            br(),
            br(),
            fluidRow(
                column(width = 1), 
                column(width = 5, { plotOutput("Expr.plots_dev", height = '300px') }) %>% withSpinner(),
                column(width = 5, { plotOutput("Expr.plots_tis", height = '300px') }) %>% withSpinner()
            ),
            br(),
            hr(),
            br(),
            
            ## Row 3: OUTPUT table of associated REs table
            h3('Table of associated regulatory elements (REs)'),
            br(),
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
)

# Multiple-genes entry
TAB2 <- tabItem(
    tabName = 'geneslookup',
    fluidPage(
        
        ## Row 1: Multiple genes entry
        fluidRow(
            column(width = 3, fluidRow(
                textAreaInput(
                    inputId = "searchMulitpleGenePaste",
                    label = h4("Paste gene names:"),
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
                p(style = "color: grey;", 'List longer than 500 genes may take a while to be processed'),
                br(), 
                textOutput("multipleGenesLength"),
                actionBttn("resetGenes", label = "Reset genes query", size = "xs", style = "simple", color = 'primary')
            ))
        )
    ),
    conditionalPanel(
        condition = "output.displayPanelsTab2 == '1'", 
        fluidPage(
            br(),
            hr(),
            br(), 
            downloadBttn("downloadGenesListGFF", label = "Download detailed report of input genes and associated REs (GFF file, IGV friendly)", size = "sm", style = "simple", color = "primary", block = F),
            br(),
            br(),
            hr(),
            br(),
            
            ## Row 2a: OUTPUT Venn diagrams
            fluidRow(
                h4("Intersection of genes query with tissue-enriched genes"),
                column(width = 2, { plotOutput("Venn.Germline")%>% withSpinner(type = 6, color = "#421500", size = 0.5)  }),
                column(width = 2, { plotOutput("Venn.Neurons")%>% withSpinner(type = 6, color = "#421500", size = 0.5)  }),
                column(width = 2, { plotOutput("Venn.Muscle")%>% withSpinner(type = 6, color = "#421500", size = 0.5)  }),
                column(width = 2, { plotOutput("Venn.Hypod")%>% withSpinner(type = 6, color = "#421500", size = 0.5)  }),
                column(width = 2, { plotOutput("Venn.Intest")%>% withSpinner(type = 6, color = "#421500", size = 0.5)  })
            ),
            br(),
            hr(),
            br(),
            
            ## Row 2b: OUTPUT HMs.plot
            fluidRow(
                column(width = 4, { plotOutput("HMs.plot_LCAPdev")%>% withSpinner(type = 6, color = "#421500", size = 0.5)  }),
                column(width = 4, { plotOutput("HMs.plot_LCAP")%>% withSpinner(type = 6, color = "#421500", size = 0.5)  }),
                column(width = 4, { plotOutput("HMs.plot_ATAC")%>% withSpinner(type = 6, color = "#421500", size = 0.5)  })
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
                column(width = 10, { plotOutput("GO.plot") %>% withSpinner(type = 6, color = "#421500", size = 0.5) } )
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
        fluidRow(
            column(3, h4("Download our tissue-specific datasets", style = 'display: inline-block')), 
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
        h4("Navigate ATAC-seq data"),
        fluidRow( column(width = 12, { dataTableOutput("atac.table") } ) ),
        hr(),
        
        # Browse LCAP table
        h4("Navigate RNA-seq data"),
        fluidRow( column(width = 12, { dataTableOutput("lcap.table") } ) )
    )
)

# Informations about the data
TAB5 <- tabItem(
    tabName = 'infos',
    fluidPage(
        fluidRow(
            column(width = 12, {
                HTML(
                    '
                    <h4>Informations</h4>
                    <p>The data made available here has been published by the Ahringer lab (see list 
                    <a href = http://www.ahringer.group.gurdon.cam.ac.uk/publications.html>here</a>
                    ).</br>
                    All the gene annotations and regulatory elements annotations are performed using 
                    ce11 genome.</br>
                    For detailed Material and Methods, please refer to the corresponding publication(s)
                    [<a href = http://www.ahringer.group.gurdon.cam.ac.uk/publications.html>here</a>].</br>
                    To get access to all the scripts used in our studies, please refer to our GitHub 
                    repositories 
                    <a href = https://github.com/js2264>here</a>, 
                    <a href = https://github.com/jurgjn>here</a> and 
                    <a href = https://github.com/Przemol/>here</a>. 
                    Feel free to contact us directly for any request!
                    </p>
                    '
                )
            })
        ), 
        br(), 
        br(), 
        br(), 
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

# Contact us
TAB6 <- tabItem(
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
                    <img src="http://ahringerlab.com/assets/img/JS.jpg" alt="JS" style="height: 147px; margin-top: 20px">
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
    menuItem("Informations", tabName = "infos", icon = icon("info", lib = "font-awesome")),
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

BODY <- tabItems(TAB1, TAB2, TAB3, TAB4, TAB5, TAB6)

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
        tags$head(includeCSS("assets/custom.css")),
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
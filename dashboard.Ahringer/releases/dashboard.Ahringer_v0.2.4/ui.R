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
        h4('t-SNE plots of promoters (left) and genes (right)'),
        br(),
        h5("Black dots represent the gene and its associated promoters."),
        br(),
        fluidRow(
            column(width = 5, { plotOutput("tSNE.plots_genes")  }),
            column(width = 2, { plotOutput("tSNE.plots_legend")  }),
            column(width = 5, { plotOutput("tSNE.plots_proms")  })
        ),
        br(),
        hr(),
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
                column(width = 3, ""),
                column(width = 6, { actionBttn("getList", label = "Perform analysis", icon = icon("search", lib = "font-awesome"), style = 'bordered', color = "primary", block = T) }),
                column(width = 3, "")
            ),
            fluidRow(
                column(width = 3, ""),
                column(width = 3, { textOutput("multipleGenesLength") }),
                column(width = 3, { actionBttn("resetGenes", label = "Reset genes query", size = "xs", style = "fill", ) }),
                column(width = 3, "")
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
            column(width = 2, { plotOutput("Venn.Hypod")  }),
            column(width = 2, { plotOutput("Venn.Neurons")  }),
            column(width = 2, { plotOutput("Venn.Germline")  }),
            column(width = 2, { plotOutput("Venn.Muscle")  }),
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
                numericInput("NCLUST_LCAPdev", label = "Number of k-means clusters", value = 5),
                switchInput(
                    "colorScale_doRev_LCAPdev", 
                    label = "Reverse color scale?", 
                    value = T, 
                    onLabel = "Yes", 
                    offLabel = "No", 
                    onStatus = 'success', 
                    offStatus = 'primary', 
                    labelWidth = 100, 
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
                    labelWidth = 100, 
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
                    labelWidth = 100, 
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
        column( width = 2, { htmlOutput("genesList") } ),
        br()
    )
)

# Genome browser
TAB3 <- tabItem(
    tabName = 'browser',
    fluidPage( 
        fluidRow(
            column(width = 5, downloadBttn("downloadBWLCAP", label = "Download all the tissue-specific RNA-seq tracks (bigwig format)", size = "sm", style = "fill", color = "primary", block = T) ),
            column(width = 5, downloadBttn("downloadBWATAC", label = "Download all the tissue-specific ATAC-seq tracks (bigwig format)", size = "sm", style = "fill", color = "primary", block = T) )
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
    	h2("Contact us"),
        fluidRow(
            column(width = 4, {
                HTML(
                    paste(
                        h3("Lab:\n\n"),'<br/>',
                        h5("Ahringer lab\n\n"),'<br/>',
                        h5("Gurdon Institute\n\n"),'<br/>',
                        h5("University of Cambridge, UK")
                    )
                )
            }),
            column(width = 4, {
                HTML(
                    paste(
                        h3("Developer/Maintenance:\n\n"),'<br/>',
                        h5("Jacques Serizay\n\n"),'<br/>',
                        h5("jserizay.site\n\n"),'<br/>',
                        h5("js2264 -at- cam.ac.uk")
                    )
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
    menuItem("Download datasets", tabName = "download", icon = icon("download", lib = "font-awesome")),
	menuItem("Contact us", tabName = "contact", icon = icon("envelope-open", lib = "font-awesome")),
    sidebarSearchForm(textId = "quickGene", buttonId = "quickSearch", label = "Quick gene search..."), 
    tags$footer(
        img(src = "http://tispelegans.site/www_JABrowse/sidebar-img_150x150.png", alt = "", style = "
            color: #b8c7ce;
            padding: 0px 40px 30px 40px;
            z-index: 1000;
        "),
        br(),
        HTML(paste("Ahringer lab -", icon("copyright", lib = "font-awesome"), "2018")), 
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

BODY <- tabItems(TAB1, TAB2, TAB3, TAB4, TAB5)

shinyUI <- dashboardPage(
    dashboardHeader(title = "Ahringer lab C. elegans tissue-specific database", titleWidth = 450),
    dashboardSidebar(SIDEBAR),
    dashboardBody(
        tags$script(inactivity),
        tags$head(tags$link(rel = "shortcut icon", href = "http://tispelegans.site/www_JABrowse/favicon.ico")),
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
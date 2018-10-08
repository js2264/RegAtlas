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
        })
    ),
    fluidRow( 
      column(width = 4, { htmlOutput("geneInfos", height = '200px') }),
      column(width = 2, { 
        fluidRow(
          br(), 
          downloadButton("downloadINFOS", label = "Download full gene report (.txt file)"),
          br(),
          br(),
          actionButton("switchToGenome", label = "Go to Genome Browser", icon = icon("area-chart", lib = "font-awesome"))
        )
      })
    ),
    br(),
    hr(),
    br(),

    ## Row 2: OUTPUT LCAPdev and LCAP-tissues graphs ,as well as gene infos.
    h4("Temporal and spatial gene expression profiles"),
    br(),
    fluidRow(
        column(width = 5, { plotOutput("Expr.plots_dev", height = '300px') }),
        column(width = 5, { plotOutput("Expr.plots_tis", height = '300px') })
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
        column(width = 5, { plotOutput("tSNE.plots_genes") }),
        column(width = 2, { plotOutput("tSNE.plots_legend") }),
        column(width = 5, { plotOutput("tSNE.plots_proms") })
    ),
    br(),
    hr(),
    br()
)

# Multiple-genes entry
TAB2 <- tabItem(
	tabName = 'geneslookup',

    ## Row 1: Multiple genes entry
    fluidRow(
        column(width = 2, {
            textAreaInput(
                inputId = "searchMulitpleGenePaste",
                label = "Paste gene names:",
                value = paste(germline.genes, collapse = '\n'),
                placeholder = 'Paste here (one gene per line)',
                rows = 8
            )
        })
    ),
    textOutput("multipleGenesLength"),
    h6("The toy example contains genes enriched in germline compared to other tissues"),
    hr(),
    br(),

    ## Row 2: OUTPUT GO.plot
    h4("Enriched Gene Ontology terms (gProfileR)"),
    fluidRow( column(width = 10, { plotOutput("GO.plot") } ) ),
    br(),
    hr(),
    br(),

    ## Row 3: OUTPUT HMs.plot
    h4("Tissue-specific enrichment of query genes and associated REs"),
    fluidRow( column(width = 6, { plotOutput("HMs.plot") } ) ),
    br(),
    hr(),
    br()

)

# Genome browser
TAB3 <- tabItem(
    fluidPage( fluidRow( { 
        JbrowseOutput("jbrowser") 
    } ) ), 
    tabName = 'browser'
)

# Download data
TAB4 <- tabItem(
    tabName = 'download',

    # Download buttons
    h2("Download our tissue-specific datasets"),
    fluidRow(
        column(width = 5, { fluidRow(
            downloadButton("downloadATAC.txt", label = "Download the entire tissue-specific ATAC-seq dataset (txt format, Excel friendly)"),
            br(),
            br(),
            downloadButton("downloadATAC.gff", label = "Download the entire tissue-specific ATAC-seq dataset (GFF format, IGV friendly)")            
        ) } ),
        column(width = 5, { fluidRow(
            downloadButton("downloadLCAP.txt", label = "Download the entire tissue-specific RNA-seq dataset (txt format, Excel friendly)"), 
            br(),
            br(),
            downloadButton("downloadLCAP.gff", label = "Download the entire tissue-specific RNA-seq dataset (GFF format, IGV friendly)")
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

# Contact us
TAB5 <- tabItem(
	tabName = 'contact',
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

## Finalise UI ----------------------------------------------------------------------------------------------------

SIDEBAR <- sidebarMenu(
	id = "tabs", 
    menuItem("Look-up gene", tabName = "genelookup", icon = icon("ellipsis-h", lib = "font-awesome")),
    menuItem("Look-up multiple genes", tabName = "geneslookup", icon = icon("ellipsis-h", lib = "font-awesome")),
    menuItem("Genome browser", tabName = "browser", icon = icon("area-chart", lib = "font-awesome")),
    menuItem("Download datasets", tabName = "download", icon = icon("download", lib = "font-awesome")),
	menuItem("Contact us", tabName = "contact", icon = icon("envelope-open", lib = "font-awesome"))
)

BODY <- tabItems(TAB1, TAB2, TAB3, TAB4, TAB5)

shinyUI <- dashboardPage(
    dashboardHeader(title = "Ahringer lab C. elegans tissue-specific database", titleWidth = 450),
    dashboardSidebar(SIDEBAR),
    dashboardBody(
        tags$head(tags$style(HTML('.left-side, .main-sidebar {background-image: url("http://tispelegans.site/img/sidebar-img_230x700.jpg"); background-size: 230px;}'))),
        # tags$head(tags$style(HTML('.content {background-image: url("http://tispelegans.site/img/body-img_2000x3000.jpg"); background-size: 98%;}'))),
        tags$head(tags$style(HTML('.skin-blue .main-header .logo {background-color: #333;} .skin-blue .main-header .logo:hover {background-color: #333;}'))),
        tags$head(tags$style(HTML('.skin-blue .main-header .navbar {background-color: #333;} .skin-blue .main-header .navbar {background-color: #333;}'))),
        tags$head(tags$style(HTML('.skin-blue .main-header .navbar .sidebar-toggle {background-color: #333;} .skin-blue .main-header .navbar .sidebar-toggle:hover {background-color: #444;}'))),
        tags$head(tags$style(HTML('.skin-blue .left-side, .skin-blue .main-sidebar, .skin-blue .wrapper {background-color: #333;}'))),
        tags$head(tags$style(HTML('.main-header .sidebar-toggle {size:100px}'))),
        tags$head(tags$style(HTML('.main-header .sidebar-toggle:before {size:100px; content: "\\f0d9"}'))),
        tags$head(tags$style(HTML('.main-header .sidebar-toggle:after {size:100px; content: "\\f0da"}'))),
        tags$head(tags$style(HTML("hr {border-top: 1px dashed #b7b7b7;}"))),
        tags$head(tags$style(HTML(".btn {border-radius: 30px;} .btn:hover {transform: scale(1.05); background-color: #bebebe}"))),
        BODY
    )
)

## To run the app locally ( DO NOT INCLUDE WHEN DEPLOYING ) ------------------------------------------------------

#shinyApp(ui = shinyUI, server = shinyServer)

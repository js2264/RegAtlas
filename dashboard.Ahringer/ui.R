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

## Define each dashbord tab individually ---------------------------------------------------------------------------

# Single-gene entry
TAB1 <- tabItem(
    tabName = 'genelookup',

    ## Row 1: Gene entry
    fluidRow(
        column(width = 6, {
            searchInput(
                inputId = "searchGene",
                label = "Enter a single gene name (WormBase ID or locus ID):",
                value = 'hlh-1',
                btnSearch = icon("search", lib = "glyphicon"),
                btnReset = icon("remove", lib = "glyphicon"),
                width = "95%"
            )
        }),
        column(width = 4, { htmlOutput("geneInfos", height = '200px') })
    ),
    fluidRow( column(width = 3, { downloadButton("downloadINFOS", label = "Download gene results (.txt file)") } ) ),
    hr(),
    br(),

    ## Row 2: OUTPUT LCAPdev and LCAP-tissues graphs ,as well as gene infos.
    h4("Temporal and spatial gene expression profiles"),
    fluidRow( column(width = 8, { plotOutput("Expr.plots", height = '300px') }) ),
    br(),
    hr(),
    br(),

    ## Row 3: OUTPUT tSNE plots (ATAC / LCAP)
    h4('t-SNE plots of promoters (left) and genes (right)'),
    h5("Black dots represent the gene and its associated promoters."),
    fluidRow( column(width = 10, { plotOutput("tSNE.plots") } ) ),
    br(),
    hr(),
    br(),

    ## Row 4: OUTPUT table of associated REs table
    h3('Table of associated regulatory elements (REs)'),
    fluidRow( dataTableOutput("REs.table") ),
    br(),
    hr(),
    br()
)

# Genome browser
TAB2 <- tabItem(
    fluidPage( fluidRow( { 
        JbrowseOutput("jbrowser") 
    } ) ), 
    tabName = 'browser'
)

# Multiple-genes entry
TAB3 <- tabItem(
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

# Download data
TAB4 <- tabItem(
    tabName = 'download',

    # Download buttons
    h2("Download our tissue-specific datasets"),
    fluidRow(
        column(width = 4, { downloadButton("downloadATAC", label = "Download the entire tissue-specific ATAC-seq dataset") } ),
        column(width = 4, { downloadButton("downloadLCAP", label = "Download the entire tissue-specific RNA-seq dataset") } )
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
	menuItem("Look-up gene", tabName = "genelookup", icon = icon("ellipsis-h", lib = "font-awesome")),
    menuItem("Genome browser", tabName = "browser", icon = icon("area-chart", lib = "font-awesome")),
    menuItem("Look-up multiple genes", tabName = "geneslookup", icon = icon("ellipsis-h", lib = "font-awesome")),
    menuItem("Download datasets", tabName = "download", icon = icon("download", lib = "font-awesome")),
	menuItem("Contact us", tabName = "contact", icon = icon("envelope-open", lib = "font-awesome"))
)

BODY <- tabItems(TAB1, TAB2, TAB3, TAB4, TAB5)

shinyUI <- dashboardPage(
    dashboardHeader(title = "Ahringer lab"),
    dashboardSidebar(SIDEBAR),
    dashboardBody(
        tags$head(tags$style(HTML('
          .main-header .sidebar-toggle:before {
            content: "\\f0d9";}'))),
        tags$head(tags$style(HTML('
          .main-header .sidebar-toggle:after {
            content: "\\f0da";}'))),
        tags$head( tags$style(HTML("hr {border-top: 1px solid #000000;}")) ),
        BODY
    )
)

## This runs the app locally ( DO NOT INCLUDE WHEN DEPLOYING ) ------------------------------------------------------

#shinyApp(ui = shinyUI, server = shinyServer)

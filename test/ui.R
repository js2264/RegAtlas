
# Multiple-genes entry
shinyUI <- fluidPage(
    
        ## Row 1: Multiple genes entry
        fluidRow(
            column(width = 2, {
                textAreaInput(
                    inputId = "searchMulitpleGenePaste",
                    label = "Paste gene names:",
                    value = "",
                    placeholder = 'Paste here (one gene per line)',
                    rows = 8
                )
            }), 
            column(width = 1, h4("... Or ...")),
            column(width = 2, { 
                tags$div(class = "multicol2", checkboxGroupInput("checkGroupGeneClasses", label = h4("Select the classes of genes to query"), selected = "muscle.genes", choices = list(
                    "Hypodermis-enriched genes" = "hypod.genes", 
                    "Neurons-enriched genes" = "neurons.genes", 
                    "Germline-enriched genes" = "germline.genes", 
                    "Muscle-enriched genes" = "muscle.genes", 
                    "Intestine-enriched genes" = "intest.genes"
                )))
            }),
            column(width = 1, h4("... Or ...")),
            column(width = 4, {
                actionBttn("MoreChoicesGenesList", label = "Choose genes based on their promoters", icon = icon("filter", lib = "font-awesome"), size = "sm", style = "fill", color = "primary")
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
                        selected = "Neurons_Intest.", 
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
                    br(),
                    br(),
                    br(),
                    helpText("To be able to tick any checkbox, all the options in the main page must be unticked."),
                    br()
                )
            )
        ),
        
        fluidRow(
            column(width = 1, { textOutput("multipleGenesLength") %>% withSpinner(type = 6, color = "#421500") }),
            column(width = 2, { actionBttn("resetGenes", label = "Reset genes query", size = "xs", style = "fill") })
        ),
        br(),
        h6("Plots in this tab may take a while to load..."),
        br(), 
        downloadBttn("downloadGenesListGFF", label = "Download detailed report of input genes and associated REs (GFF file, IGV friendly)", size = "sm", style = "fill", color = "primary", block = T),
        br(),
        hr(),
        br()

)

## 
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

  # Get the list of multiple genes
  {
    multipleGenes <- reactive ({ 
        
        
        text.list <- unique(unlist(strsplit(x = input$searchMulitpleGenePaste, split = '[\r\n]' )) %>% ifelse(!grepl('WBGene', .), name2WB(.), .) %>% .[!is.na(.)])
        checkbox.list <- unlist(list.genes[input$checkGroupGeneClasses])
        bsmodal.list <- c()
        
        genes <- unique(c(text.list, checkbox.list, bsmodal.list))
        return(genes[genes %in% names(genes.gtf)])
        
    })
    
    observe( if(input$resetGenes > 0){
        
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
                "Gonad & Muscle & Intestine-specific promoter(s)",
                "Neurons & Gonad & Muscle & Intestine-specific promoter(s)",
                "Hypodermis & Gonad & Muscle & Intestine-specific promoter(s)",
                "Hypodermis & Neurons & Gonad & Intestine-specific promoter(s)",
                "Hypodermis & Neurons & Gonad & Muscle-specific promoter(s)",
                "Soma-specific promoter(s)",
                "Ubiquitous promoter(s)"
            ),
            choiceValues = order.tissues[c(1:27, 29, 30, 32:33)],
            selected = NULL
        )
        
    } )

    output$multipleGenesLength <- reactive ({ paste(length(multipleGenes()), "valid gene(s) found.") })
    
  }

} #EOF

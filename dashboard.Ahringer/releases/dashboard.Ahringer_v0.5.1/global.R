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

## Load data and libraries ----------------------------------------------------------------------------------------------------

load('data/minimal-data.RData')
source('R/custom_R_functions.R')
source('R/helpers.R')
suppressMessages(require(shiny))
suppressMessages(require(shinybusy))
suppressMessages(require(shinydashboard))
suppressMessages(require(shinyWidgets))
suppressMessages(require(shinycssloaders))
suppressMessages(require(shinyBS))
suppressMessages(require(DT))
suppressMessages(require(urltools))
suppressMessages(require(htmltools))
suppressMessages(require(httr))
suppressMessages(require(apputils))

# Define variables
colorGO <- c("#52c23164","#20109564","#ef690664", "white")
names(colorGO) <- c("MF", "BP", "CC", "kegg")
colors.decimals <- c(0.10395294, 0.3906374, 0.1192665, 0.14261010, 0.14226132, 0.13421772)
link <- "http://ahringerlab.com/JBrowse-1.12.5/index.html?data=data%2Fjson%2Fce11&loc=chrIII&tracks=genes%2Cregulatory_elements%2Cgonad.atac%2Cneurons.atac%2Cmuscle.atac%2Chypod.atac%2Cintest.atac%2Cgonad.lcap.fwd%2Cneurons.lcap.fwd%2Cmuscle.lcap.fwd%2Chypod.lcap.fwd%2Cintest.lcap.fwd%2Cgonad.lcap.rev%2Cneurons.lcap.rev%2Cmuscle.lcap.rev%2Chypod.lcap.rev%2Cintest.lcap.rev%2Ctranscripts&highlight=&menu=1&nav=1&tracklist=1&overview=1"
NCLUST_LCAPdev <- 5
NCLUST_LCAP <- 5
NCLUST_ATAC <- 5
version <- gsub("../dashboard.Ahringer_", "", list.dirs('..', recursive = F)[length(list.dirs('..', recursive = F)) - 1])
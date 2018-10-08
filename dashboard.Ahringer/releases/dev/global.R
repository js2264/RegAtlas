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

require(shiny)
require(shinydashboard)
require(shinyWidgets)
require(shinyBS)
require(shinycssloaders)
require(devtools)
require(DT)
require(urltools)
require(htmltools)
require(httr)
require(raster)
require(magrittr)
require(gProfileR)
require(dplyr)
require(tidyr)
require(pheatmap)
require(RColorBrewer)
if (!require(trewjb)) { devtools::install_github("Marlin-Na/trewjb") ; require(trewjb)}
if (!require(GenomicRanges)) { source("https://bioconductor.org/biocLite.R") ; biocLite("GenomicRanges") ; require(GenomicRanges)}
options("repos" = BiocInstaller::biocinstallRepos())
load('data/tSNE.minimal.RData')
source('bin/useful_R_functions.R')
source('bin/heatmap_parmfrrow.R')

# Define variables
colorGO=c("#52c23164","#20109564","#ef690664", "white")
names(colorGO)=c("MF", "BP", "CC", "kegg")
all.deconv <- cbind(all, ATAC, max.tissue.df[,5:19]) %>% mutate(uniqueWormBaseID = strsplit(as.character(WormBaseID),',')) %>% unnest(uniqueWormBaseID)
atac.dt <- cbind(all[,c(1:5,7,10)], which.tissues = max.tissue.df$which.tissues, round(ATAC, 3), max.tissue.df[,5:9])
row.names(atac.dt) <- all$coords
colnames(atac.dt)[5] <- "geneID"
colnames(atac.dt)[9:13] <- paste0(colnames(atac.dt)[9:13], '_TPM')
atac.dt$regulatory_class <- factor(atac.dt$regulatory_class)
atac.dt$domain <- factor(atac.dt$domain)
lcap.dt <- cbind(as.data.frame(genes.gtf)[,c(1:3,5,10,11,13,16,18,19)], which.tissues = genes.gtf$which.tissues, round(LCAP, 3), max.tissue.df.LCAP[,2:6])
colnames(lcap.dt)[1:3] <- colnames(atac.dt)[1:3]
colnames(lcap.dt)[5] <- "WormBaseID"
colnames(lcap.dt)[6] <- "geneID"
lcap.dt$gene_biotype <- factor(lcap.dt$gene_biotype)
lcap.dt$domain <- factor(lcap.dt$domain)
cv <- genes.gtf$cv
hypod.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Hypod.']
neurons.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Neurons']
germline.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Gonad']
muscle.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Muscle']
intest.genes <- row.names(max.tissue.df.LCAP)[max.tissue.df.LCAP$which.tissues == 'Intest.']
list.genes <- list(hypod.genes, neurons.genes, germline.genes, muscle.genes, intest.genes)
names(list.genes) <- c("hypod.genes", "neurons.genes", "germline.genes", "muscle.genes", 'intest.genes')

# Function to generate an URL to visit jserizay.site/JBrowse
getURL <- function (chr, start, end, release = "1.12.5", 
                    tracks = c("genes", "regulatory_elements", "hypod.atac", "neurons.atac", "gonad.atac", "muscle.atac", "intest.atac", "hypod.lcap.fwd", "neurons.lcap.fwd", "gonad.lcap.fwd", "muscle.lcap.fwd", "intest.lcap.fwd", "hypod.lcap.rev", "neurons.lcap.rev", "gonad.lcap.rev", "muscle.lcap.rev", "intest.lcap.rev", "transcripts"),
                    show_menu = TRUE, show_navigation = TRUE, show_tracklist = TRUE, show_overview = TRUE)
{
  baseurl <- paste0("http://tispelegans.site/JBrowse-", release, "/index.html")
  range <- if (missing(start) || missing(end)) {""} else {paste0("%3A", parseRange(start = as.numeric(start), end = as.numeric(end), resizeFactor = 1.5))}
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

# Re-WRITE iframeJBrowse to get functional height setting
iframeJbrowse.2 <- function (link, width = NULL, height = "900px", elementId = NULL)
{
    div_style = paste0("width: 100%; height: ", height, ";")
    iframe_html <- htmltools::tags$iframe(src = link, style = "border: 1px solid black", width = "100%", height = "100%", "Sorry, your browser does not support iframe.")
    div_html <- htmltools::tags$div(iframe_html, style = div_style)
    htmlwidgets::createWidget(name = "trewjb", list(html = as.character(div_html)), package = "trewjb", elementId = elementId)
}

# Function to input several genes separated by newlines
textareaInput <- function(inputId, label, value = "", placeholder = "", rows = 2)
{
  tagList(
    div(strong(label), style="margin-top: 5px;"),
    tags$style(type="text/css", "textarea {width:100%; margin-top: 5px;}"),
    tags$textarea(id = inputId, placeholder = placeholder, rows = rows, value))
}

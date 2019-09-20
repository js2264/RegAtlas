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
    all[, c('chr','start','stop','gene_name','regulatory_class','domain','which.tissues')], 
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

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

# Function to extract last character of a string
substrRight <- function(x, n) {
    substr(x, (nchar(x)-n+1), nchar(x))
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

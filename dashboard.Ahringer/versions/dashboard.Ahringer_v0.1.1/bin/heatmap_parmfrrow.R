heatmap_parmfrow <- function (x, margins = c(5, 3, 6.5, 5), outmargins = c(2, 2, 2, 2),
    labRow = NULL, labCol = NULL,
    cexRow = 0.2 + 1/log10(nrow(x)), cexCol = 0.2 + 1/log10(ncol(x)), cex.main = 1.2,
    col = heat.colors, breaks = NULL,
    main = NULL, key.lab = NULL, ylab = NULL, xlab = NULL, scaleUP = 1.2, scaleUP.bis = NULL,
    nticks = 10, doClust = F, plotSep = F, colSep = '#00000064', plotKey = T, printCell = F, CellCex = 1)

{

hm <- list(
    x = x,
    margins = margins,
    labRow = labRow,
    labCol = labCol,
    cexRow = cexRow,
    cexCol = cexCol,
    cex.main = cex.main,
    col = col,
    breaks = breaks,
    main = main,
    key.lab = key.lab,
    ylab = ylab,
    scaleUP = scaleUP,
    scaleUP.bis = scaleUP.bis,
    nticks = nticks,
    doClust = doClust,
    plotSep = plotSep,
    colSep = colSep
)

# If labels are null, set them to col/row names
if (is.null(labRow)) { labRow <- row.names(x) } else if (labRow[1] == 'none') { labRow <- rep("", length(row.names(x))) }
if (is.null(labCol)) { labCol <- colnames(x) } else if (labCol[1] == 'none') { labCol <- rep("", length(colnames(x))) }

# If x is not a matrix, stop function
if (class(x) != 'matrix') {
    stop('...Input needs to be a matrix...')
}

# scaleUP.bis is the factor to increase the heatmap height, to draw the TOP horizontal line of the color key above the actual heatmap
# scaleUP is the factor to increase the heatmap height, to draw the BOTTOM horizontal line of the color key above the actual heatmap
if (is.null(scaleUP.bis)) {scaleUP.bis <- scaleUP - 0.1}

# Required number of rows/columns (before transposition!!)
nr <- nrow(x)
nc <- ncol(x)

# If the breaks are not defined: do 16 breaks (or the number of colors -1 if the colors are given)
if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    stop('...Breaks missing...')
}

# Clamp values of x to the breaks
x[x < min(breaks)] <- min(breaks)
x[x > max(breaks)] <- max(breaks)

# If the color argument is a function (typically 'colorRampPalette'), compute the actual colors based on how many breaks there are
if (class(col) == "function") {
    col <- col(length(breaks) - 1)
} else if (length(col) != (length(breaks) - 1)) {
    stop('...length(colors) needs to be length(breaks)-1...')
}

# Define the key labels and where they should be positionned
labels <- unique(round(seq(min(breaks), max(breaks), length.out=nticks)))
labels.at <- seq(0, nc, length.out = length(labels))

# If clustering is desired, do it
if (doClust) {
    x <- x[hclust(dist(x, method = 'max'))$order,]
}

# Save values
values <- as.character(hm$x)
values[values == 0] = ""

# Transpose the matrix (because image() is used to plot the matrix) and add some extra empty columns (which are going to be top rows when plotted) for plotting space for the color key
x <- t(x)
x <- cbind(x, matrix(rep(NA, nc * (round(nr * scaleUP)-(nr))), nrow = nc))

# Plot heatmap with empty space on top
par(mar = margins)
image(1L:nc, 1L:(round(nr*scaleUP)), x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, round(nr*scaleUP)),
    axes = FALSE, xlab = "", ylab = "",
    col = col, breaks = breaks)
axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
axis(4, 1:nr, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow, srt=180)

# Plot frame
polygon(0.5 + c(0,0,nc,nc), 0.5 + c(0,nr,nr,0), col=NA, border="black")

# Plot cells separators (if wanted)
if (plotSep) {
    for (i in seq(0, nc)) { polygon(0.5 + c(i,i,nc,nc), 0.5 + c(0,nr,nr,0), col=NA, border=colSep) }
    for (i in seq(0, nr)) { polygon(0.5 + c(0,0,nc,nc), 0.5 + c(0,i,i,0), col=NA, border=colSep) }
}

# Plot legend
if (plotKey) {
    poly <- vector(mode="list", length(col))
    val <- seq(0, nc, length.out = length(breaks))
    for(i in seq(poly)){ poly[[i]] <- 0.5 + c(val[i], val[i], val[i+1], val[i+1]) } ## Define coordinates of each rectangle of the key
    for(i in seq(poly)){ polygon(poly[[i]], c(round(nr*scaleUP), round(nr*scaleUP.bis), round(nr*scaleUP.bis), round(nr*scaleUP)), col = col[i], border = NA) } ## Plot each rectangle of the color key
    polygon(0.5 + c(0,0,nc,nc), c(round(nr*scaleUP.bis), round(nr*scaleUP),round(nr*scaleUP), round(nr*scaleUP.bis)), col = NA, border = "black")
    axis(3,at = labels.at+0.5, labels = labels, las = 2, cex.axis = 0.5, pos = round(nr*scaleUP))
    #axis(3,at=seq(0, nc, length.out=nticks), labels = labels.at, las=2, cex.axis=0.5, pos=round(nr*scaleUP))
}

# Plot labels
if (!is.null(ylab)) { mtext(ylab, side = 4, line = margins[4]-2, srt = 180) }
if (!is.null(xlab)) { mtext(xlab, side = 1, line = margins[1]-2) }
if (!is.null(key.lab)) { mtext(key.lab, side = 3, line = 2, cex=0.8) }
mtext(main, side = 3, line = 4, cex = cex.main)

# Plot cell values
if (is.logical(printCell) & printCell == T) { text(x = rep(1:nrow(x), each = ncol(x)), y = 1:ncol(x), values, cex = CellCex) }
if (!is.logical(printCell)) {text(x = rep(1:nrow(x), each = ncol(x)), y = 1:ncol(x), printCell, cex = CellCex)}

return(hm)

}

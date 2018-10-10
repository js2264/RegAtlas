# Changelog

## To do list
* Add loading screen 
* Make better contact cards
* Add bookmark options

------------ 

## [Unreleased]


------------ 


## [0.2.4] - 2018-10-10
##### Added
- Patterns for multiple genes entries
- Foot copyright and logo in the sidebar
- Settings for heatmaps now in a dropdown button
- Button in multi-genes to run all analysis AFTER selection of genes of interest
- WB description in quickView
- Button to switch between k-means and hclustering

##### Changed
- Improved string split for multiple genes entries
- Resized elements for better displays in lower resolutions
- Fixed navbar and sidebar to the screen
- Fixed browser height, now fit to any screen size

##### Deleted
- Libraries that were not needed (mostly trewjb)

## [0.2.3] - 2018-10-08
##### Added
- Added Euler diagrams

##### Changed
- Switched to Neurons-enriched genes for default in multiple genes list view

## [0.2.2] - 2018-10-08
##### Added
- Quick View button
- Added filtering option in GO analysis
- Added copyright
- In the multiple genes view, add temporal expression heatmap

## [0.2.1] - 2018-10-08
##### Added
- Gene description fetched from WB, when clicking on a button

##### Changed
- Modified directories architecture
- Modified procedure to upload the app
- Style of actionBttns (to launch GO analysis and get gene description)

## [0.2.0] - 2018-10-08
##### Added
- Heavily enhanced the heatmaps plotted in multi-genes tab (possibility to choose which colors, which units to use, reverse colors)
- Button to download GO results
- Link to WormBase
- Button to download the tracks bigwig files
- In the multiple genes view, there is now a button to download a GFF file of the query results (rbind-ing the genes and the REs in 1 file)
- Spinner wheel in each tab is displayed while loading the tab
- You now have to click on a button to run GO analysis. This shortens the loading time of the tab.
- Reset button to erase all the queried genes
- Selection boxes to select tissue-enriched genes or genes associated with tissue-specific promoter(s) 

##### Changed
- Switched to pheatmap to plot heatmaps
- Updated button style, using Bttns instead of regular buttons

## [0.1.1] - 2018-10-06
##### Added
- APP AES. CHANGES: 
	- Added Coordinates to geneInfos
	- Added a button linking to genome Browser
	- Add a background image in the dashboard sidebar
- Possibility to download a gff file of the annotations / genes

##### Changed
- APP AES. CHANGES: 
	- Changed theme to black
	- Changed order of tabs (now 1 gene ; multiple genes ; genome browser ; data download ; contact)
	- Moved geneInfos below the gene entry space (to see it better when querying smth)
	- Turn labels in plots to vertical
	- Moved t-SNE plots after REs table
	- Each plot is plotted individually now, so that it displays correctly
	- Changed style of clickable buttons
- Updated getGeneInfos function to work even when no REs are associated

## [0.1.0] - 2018-10-04
##### Added
- A working app (v0.1.0)
- This CHANGELOG file
- README contains explanations about how to set up a Digital Ocean Droplet server and to configure it to server Shiny Apps

##### Changed
- APP AES. CHANGES: 
	- HEADER: Icon to collapse the sidebar is now left+side carets
	- TAB1: Moved tSNE legend on top of the graphs
	- TAB2: Genome browser: height is adjusted to 850px
	- TAB2: Genome browser: title has been removed
	- TAB2: Zoom out on genes when querying a gene in the browser
- BROWSER AES. CHANGES: 
	- Re-add the “add track” button
	- Updated few plugins
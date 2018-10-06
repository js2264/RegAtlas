# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## To do list
* Button to download the tracks bigwig files


## [Unreleased]
The unreleased version of the app is contained in the parrent directory (`dashboard.Ahringer/`).
It will be released once synchronizing with the server (tispelegans.site = 167.99.196.115).



## [0.1.1] - 2018-10-05
##### Added
- APP AES. CHANGES: 
	- Added Coordinates to geneInfos
	- Added a button linking to genome Browser
	- Add a background image in the dashboard sidebar
- Added the possibility to download a bed file of the annotations / genes

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
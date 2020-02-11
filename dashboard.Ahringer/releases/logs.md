# Changelog

## To do list
- Browser: 
	- [ ] Display color key
	- [ ] Change the genes model tracks to see the gene names at the top
	- [ ] Rename tracks without 'specific' and adding 'YA' to RNA-seq
	- [ ] Possibility to save Genome browser Settings
	- [ ] Add heteroK tracks
	- [ ] bookmark options
- Update the annotation tracks in the browser

------------ 

## [Unreleased]

------------ 

## [0.5.2] - 2020-02-11
- added clickable links to each vignette in home page
- added support for SSL in the genome browser tab

## [0.5.1] - 2020-02-10
- Added home tab
- Speed up the loading
- Corrected the "regulatory_class" column in the atac table (exploring tab)
- Updated lab photo and icons

## [0.5.0] - 2020-02-09
- Moved many loading functions to server file. Notably, GenomicRanges is loaded automatically, and only when needed.
- put a reactiveEvent for single gene searches. 
- Landing page info
- loading buffer indicator

## [0.4.7] - 2020-02-08
##### Changed
- Corrected columns of REs table in first tab
- Removed "print" option from the REs table in first tab
- Updated single gene txt report
- Fixed regex match in multiple gene queries
- Changed height of heatmaps in multiple gene queries

## [0.4.6] - 2020-01-20
##### Changed
- Lots of small tweaks. Notably the ggi function now outputs a tidier results txt file.

## [0.4.5] - 2020-01-20
##### UI changes
- Changed gene report through ggi function (simplified)
- Added the un-normalized gene expression values

## [0.4.4] - 2019-11-28
##### Added 
##### UI changes
- Barplot for developmental LCAP
- Added "Orientation" to the quick gene infos
- Changed 2nd tab label
- Changed 1st tab icon
- Tidied up the text output (ggi function)
- Added pop-up for genome browser 
- Added blurb for methods
- Added Venn for ubiquitous genes
- Changed some labels in single-gene tab

## [0.4.2] - 2019-10-20
##### Added 
- Error message when looking for unknown gene [TAB 1]

## [0.4.1] - 2019-09-21
##### Added 
- Prevent timeout of the app. 

## [0.4.0] - 2019-09-20
##### Added 
- Tab 1 : Flexible match of genes (not case sensitive, can remove weird characters)
- Reference gProfileR
- Added Informations tab

##### UI changes
- Removed the download buttons in explore/download data, put them in a dropdown
- Also added the bigwig tracks download buttons to the same dropdown
- Reduced size of my photo
- Tab 1: remodelled the display
- Tab 2: Remove 2nd / 3rd column, replace by a dropdown list as "example". 

##### Fixed
- Fix all download buttons
- Fixed the error messages (tab 1) when genes don't exist

##### Deleted
- Fixed TAB 1 REs table issue: removed the developmental TPM values for REs

## [0.3.1] - 2019-03-29
##### Added
- A landing page has been designed!!
- Changed order of columns in REs table, and move tissue data before dev data
- Renamed column names in the REs table
- Right side headbar: added icon to go back to ahringerlab.com

##### Changed
- Because the ATACdev df only has the ATAC scores for 42245 peaks and not the new annotations, we can't recover dev. ATAC scores. One fix is to add the name of the REs (locus field) as the ATACdev row.names, so that the genes can be feched.  
- Links to img changed to assets/img/
- Updated copyright to 2019

## [0.3.0] - 2019-03-28
##### Changed
- New domain name: ahringerlab.com. Acquired on March 28th 2019 for 5years (provider: Namecheap). 
- Background color of the body switched to real white (#FFF)
- Changed photo for my contact card
- Increased size of photos in contact cards
- ADDED LICENSE GNU v3

## [0.2.4] - 2018-10-10
##### Added
- Patterns for multiple genes entries
- Foot copyright and logo in the sidebar
- Settings for heatmaps now in a dropdown button
- Button in multi-genes to run all analysis AFTER selection of genes of interest
- WB description in quickView
- Button to switch between k-means and hclustering
- Spin while plotting multiple genes heatmaps
- Added version # in the title
- Button to download heatmaps results
- Contact cards

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
This folder contains one folder per Shiny app.  
The attached script (`prepare.env.sh`) is used to:
1. Import some required data (i.e. `_tissuesTSNE/.tSNE.RData` and some R functions)
2. Create a minimal self-contained data file (`./${APP_PATH}/bin/data/tSNE.minimal.RData`).
3. Eventually, it tries to upload the app automatically, by sourcing global.R then deploying with
```r
rsconnect::deployApp()
```

*Remarque:*  
Within each application folder, 3 files are important:  
—> ui.R defines the UI of the app  
—> server.R defines the backend app computation part  
—> global.R is ran each time the app wakes up, to load required dependencies (i.e. libraries, data files, etc).  

!!! WARNING !!!  
`${APP_PATH}` is hard-coded: it is the name of the app being developed.  
!!! WARNING !!!  

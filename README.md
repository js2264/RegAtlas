# Steps to set up a Dashboard Shiny app:

__PURPOSE:__ This markdown contains instructions to create and upload a Shiny app to present results obtained during my PhD.  
__Dependencies:__ `app_wrapper.sh` is a wrapper of 3 functions (located in `bin/`) which aim to upload an existing Shiny App folder. It will first attempt to fetch foreign data from GI server. The rest of the project is self-contained.  

To deploy any app, run:  
```
bash bin/app_wrapper.sh -f=APP.NAME
```
This will execute the following commands (in the appropriate folder):  
1. ```bash bin/import.data.sh ${APP_PATH}```  
This script will automatically fetch the data from the Gurdom Institute server.

2. ```Rscript bin/get.mininaldataset.R```  
This script will reduce the size of the original data to a more manageable minimal self-contained object.

3. ```Rscript bin/deployApp.R```  
This script will deploy the application to Shiny server

This app is backed-up on Github. To push changes to Github, run:
```
git add .
git commit -m "update"
git push -u origin master
```

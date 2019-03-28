---
title: "Shiny App upload procedure"
date: "2019-03-28"
author: Jacques SERIZAY
output:
    html_document:
        theme: sandstone
        highlight: tango
        preserve_yaml: true
        df_print: paged
        dev: pdf
---
This code has been written by Jacques SERIZAY, for the purpose of analysing datasets generated in the Ahringer lab, at the Gurdon Institute, University of Cambridge.  

# Steps to set up a Dashboard Shiny app:
 
__PURPOSE:__ This markdown contains instructions to create and upload a Shiny app to present results obtained during my PhD on the project: *Tissue-specific characteristics of the chromatin in C. elegans*.  
__Dependencies:__ `app_wrapper.sh` is a wrapper of 3 functions (located in `bin/`) which aim to upload an existing Shiny App folder. The first step is to update existing data from Gurdon Institute server, which requires a password. The rest of the project is self-contained.  

To deploy any app, run:  
```
bash bin/app_release_wrapper.sh --folder "dashboard.Ahringer" --release-version "X.X.X"
```

This will execute the following commands (in the appropriate folder):  
1. ```bash bin/import.data.sh```  
This script will automatically fetch the required data from the Gurdom Institute server.

2. ```Rscript bin/get.mininaldataset.R```  
This script will reduce the size of the original RData file to a more manageable minimal self-contained object.

3. ```bash bin/makeNewRelease.sh - ${APP_RELEASE}```  
This script will create a new release directory with a specified release #, and then sync the app folder with 
This app is backed-up on Github. To push changes to Github, run:

```
git add .
git commit -m "update"
git push -u origin master
```

_Edited on 2019/03/28 to improve/streamline the process of uploading an app._

-------------------

#### Infos on Digital Ocean hosting procedure:
I subscribed to a DigitalOcean droplet (IP address: **167.99.196.115**). I followed the following steps to set up R on the server (mostly from https://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/):
1. Register to DO and create a droplet. Save droplet name and IP/passowrd (sent by email).  
**CURRENT IP ADDRESS: 167.99.196.115**
2. SSH to your droplet (```ssh root@167.99.196.115```), then enter your password (the one sent by email) then change it to a more convenient password. 
3. Create a new user and give yourself sudo rights: 
```
adduser jserizay ; gpasswd -a jserizay sudo
su - jserizay
```
4. Install Nginx
```
sudo apt-get update
sudo apt-get -y install nginx
```
5. Install R  
```
sudo sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list'
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
sudo apt-get update
sudo apt-get -y install r-base
```
6. Add more RAM memory
```
sudo /bin/dd if=/dev/zero of=/var/swap.1 bs=1M count=1024
sudo /sbin/mkswap /var/swap.1
sudo /sbin/swapon /var/swap.1
sudo sh -c 'echo "/var/swap.1 swap swap defaults 0 0 " >> /etc/fstab'
```
7. Install devtools
```
sudo apt-get -y install libcurl4-gnutls-dev libxml2-dev libssl-dev
sudo su - -c "R -e \"install.packages('devtools', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"devtools::install_github('daattali/shinyjs')\""
```
8. Install Rstudio and Shiny
```
sudo apt-get -y install gdebi-core
wget https://download2.rstudio.org/rstudio-server-1.1.442-amd64.deb
sudo gdebi rstudio-server-1.1.442-amd64.deb
sudo su - -c "R -e \"install.packages('shiny', repos='http://cran.rstudio.com/')\""
wget https://download3.rstudio.org/ubuntu-12.04/x86_64/shiny-server-1.5.6.875-amd64.deb
sudo gdebi shiny-server-1.5.6.875-amd64.deb
```
9. Upload Shiny app to /srv/shiny-server/ (ideally using the app_release_wrapper.sh script provided in ./bin/)
10. Install all the required packages directly from the server after launching R using ```sudo -i R```!  
<br>
**IMPORTANT: There is no more need to source global.R and to execute rsconnect::deployApp(). The app works by itself when ui.R and server.R are located in /srv/shiny-server/, as long as all the dependencies (required packages!!!) are installed.**
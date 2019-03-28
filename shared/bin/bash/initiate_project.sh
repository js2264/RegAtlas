#!/bin/bash

# Read variables
PROJECT=${1}
FOLDER="`date +%Y%m%d`_${PROJECT}"
README="README_${PROJECT}.Rmd"

# Create project folder and cp Rmd template
cd ~/
mkdir ${FOLDER}
cp ~/shared/bin/R/.Rmarkdown_template.Rmd "${FOLDER}/${README}"

# Modify Rmd fields
sed -i "s,TITLE_PROJECT,${PROJECT}," "${FOLDER}/${README}"
sed -i "s,DATE_PROJECT,`date +%Y-%m-%d`," "${FOLDER}/${README}"
sed -i "s,FOLDER_PROJECT,${FOLDER}," "${FOLDER}/${README}"

# Go to project folder
cd ${FOLDER}

###
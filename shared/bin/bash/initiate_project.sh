#!/bin/bash

# Read variables
PROJECT="${1}"
FOLDER="`date +%Y%m%d`_${PROJECT}"
README="README_${PROJECT}.Rmd"
SUFFIX="${2}"

# Create project folder and cp Rmd template
cd ~/
mkdir ${FOLDER}
cp ~/shared/bin/R/.Rmarkdown_template.Rmd "${FOLDER}/${README}"

# Modify Rmd fields
sed -i "s,%TITLE,${PROJECT},g" "${FOLDER}/${README}"
sed -i "s,%DATE,`date +%Y-%m-%d`,g" "${FOLDER}/${README}"
sed -i "s,%FOLDER,${FOLDER},g" "${FOLDER}/${README}"
sed -i "s,%SUFFIX,${SUFFIX},g" "${FOLDER}/${README}"
# Go to project folder
cd ${FOLDER}

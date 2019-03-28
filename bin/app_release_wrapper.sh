#!/usr/bin/env bash

## Define functions
function checkfolder() {
    if [ ! -d ${1} ]; then
        echo "ERROR: app folder not found."
        exit
    fi
}
function checkVersionNumber() {
    if [ -z "$1" ]; then
        echo "ERROR: No version number provided (-v X.X.X)"
        exit 0
    fi
}
function checkcurrentfolder() {
    if [ "$PWD" != "/Users/jacquesserizay/Documents/PhD/__Bioinfo/_shinyapps" ]; then
        echo "Moving to /Users/jacquesserizay/Documents/PhD/__Bioinfo/_shinyapps/"
        cd "/Users/jacquesserizay/Documents/PhD/__Bioinfo/_shinyapps/"
    fi
}
function usage() {
    echo -e "\nThis script is a wrapper to develop a Shiny app located in the directory ./APP_PATH/"
    echo -e ""
    echo -e "USAGE:\tapp_wrapper.sh -f APP_PATH"
    echo -e ""
}

## Read the argument
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -f | --folder)
            APP_FOLDER="$2"
            checkfolder ${APP_FOLDER}
            shift # past argument
            shift # past value
            ;;
        -id | --release-version)
            APP_RELEASE="$2"
            checkVersionNumber ${APP_RELEASE}
            shift # past argument
            shift # past value
            ;;
        -h | --help)
            usage
            exit
            ;;
        *)
            echo "ERROR: unknown parameter \"$key\""
            usage
            exit 1
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

#0. Go to shinyapps folder
checkcurrentfolder

#1. Import data from GI server
{
    echo -e "\n\n-- Fetching data --\n\n"
    rsync \
        --recursive \
        --copy-links \
        --times \
        --stats \
        --itemize-changes \
        --human-readable \
        'js2264@cb-head3.gurdon.private.cam.ac.uk:~/shared/' \
        './shared/'
}

#2. Create minimal dataset
{
    echo -e "\n\n-- Minimizing data set --\n\n"
    Rscript ./bin/get.minimaldataset.R
}

#3. Deploy App
{
    echo -e "\n\n-- Deploying app --\n\n"
    ## Copy dev folder to new release folder and make a symlink for "current release"
    DEV_FOLDER="${APP_FOLDER}/releases/dev"
    RELEASE_FOLDER="${APP_FOLDER}/releases/${APP_FOLDER}_v${APP_RELEASE}"
    CURRENT_FOLDER="${APP_FOLDER}/releases/current"
    rm -rf ${RELEASE_FOLDER}
    mkdir ${RELEASE_FOLDER}
    cp -rf ${DEV_FOLDER}/* ${RELEASE_FOLDER}/
    rm ${CURRENT_FOLDER}
    ln -s `basename ${RELEASE_FOLDER}`/ ${CURRENT_FOLDER}
    # Sync with Digital Ocean server
    rsync \
        --recursive \
        --links \
        --stats \
        --itemize-changes \
        --human-readable \
        ${APP_FOLDER}/ \
        "root@167.99.196.115:/srv/shiny-server/${APP_FOLDER}/"
}


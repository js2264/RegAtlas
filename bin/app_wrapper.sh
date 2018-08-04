#!/usr/bin/env bash

#=============================================================================#
#                                                                             #
#     USAGE:                                                                  #
#     bash app_wrapper.sh -f APP.NAME                                         #
#                                                                             #
#     AUTHOR: Jacques SERIZAY                                                 #
#     CREATED: 2018/07/12                                                     #
#     REVISION: 2018-08-04                                                    #
#                                                                             #
#=============================================================================#
 
## Define functions
function checkfolder() {
    if [ ! -d ${1} ]; then
        echo "ERROR: app folder not found."
        exit
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
            APP_PATH="$2"
            checkfolder ${APP_PATH}
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


## Go to the appropriate app folder
cd ${APP_PATH}/

#1. Import data from GI server
../bin/import.data.sh

#2. Create minimal dataset and save it in ${APP_PATH}/bin
Rscript ../bin/get.minimaldataset.R

#3. Deploy App
Rscript ../bin/deployApp.R

cd ..

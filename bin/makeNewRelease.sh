#!/usr/bin/env bash

#=============================================================================#
#                                                                             #
#     USAGE:                                                                  #
#     bash makeNewRelease.sh -v release#                                      #
#                                                                             #
#     AUTHOR: Jacques SERIZAY                                                 #
#     CREATED: 2018-08-08                                                     #
#     REVISION: ....-..-..                                                    #
#                                                                             #
#     Must be in ./dev/ to run this script                                    #
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

    echo -e "\nThis script is used to create a new release of dashboard.Ahringer"
    echo -e ""
    echo -e "USAGE:\tmakeNewRelease.sh -v release#"
    echo -e ""
}

## Read the argument
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -v | --version)
            VERSION="$2"
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


## Copy new release in a new numbered folder and in releases/current/
FOLDER=$(echo "../dashboard.Ahringer_v"${VERSION})
rm -rf ${FOLDER}
mkdir ${FOLDER}
cp -rf * ${FOLDER}

FOLDER="../current/"
rm -rf ${FOLDER}
mkdir ${FOLDER}
cp -rf * ${FOLDER}


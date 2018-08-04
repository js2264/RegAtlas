#!/usr/bin/env bash

echo -e "\n\n-- Fetching data --\n\n"

scp js2264@cb-head2.gurdon.private.cam.ac.uk:"_tissuesTSNE/.tSNE.RData" ../
scp js2264@cb-head2.gurdon.private.cam.ac.uk:"_scripts/bin/useful_R_functions.R" bin/
scp js2264@cb-head2.gurdon.private.cam.ac.uk:"_scripts/R-functions/heatmap_parmfrrow.R" bin/

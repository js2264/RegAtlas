#!/usr/bin/env bash

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
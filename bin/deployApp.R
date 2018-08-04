#!/usr/bin/env Rscript

message('\n\n-- Deploying app --\n\n')

source('global.R')
rsconnect::deployApp()

q('no')

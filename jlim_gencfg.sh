#! /bin/bash

Rscript -e 'library(getopt);' -e 'jlimR:::jlim.gencfg()' ARGSTART $* 

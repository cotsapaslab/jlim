#! /bin/bash

Rscript -e 'library(getopt); jlimR:::jlim.gencfg()' ARGSTART $* 

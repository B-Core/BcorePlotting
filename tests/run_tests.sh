#!/bin/bash

DIR=./

# Test makeMDSplot and related functions
Rscript $DIR/test1.R
rm $DIR/Rplots.pdf

# Test makeHeatmap and related functions
Rscript $DIR/test2.R



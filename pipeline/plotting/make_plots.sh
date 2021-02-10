#!/bin/bash

for item in `find ../combine/output/ | grep "inventory.txt.gz"`; do Rscript make_plots.R ${item} output mouse; done

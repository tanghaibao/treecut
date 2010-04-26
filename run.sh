#!/bin/bash

#treecut.py data/tree.nwk data/continuous.csv --cutoff .01
treecut.py data/tree.nwk data/discrete.csv --discrete
#treecut.py data/tree.nwk data/continuous.csv --printall

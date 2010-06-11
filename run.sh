#!/bin/bash

#/usr/bin/python treecut.py data/tree.nwk data/continuous.csv tree.svg --cutoff .01
/usr/bin/python treecut.py data/simple.nwk data/simple.csv simple.png
treecut.py data/tree.nwk data/discrete.csv --discrete discrete.png
#treecut.py data/simple.nwk data/simple.csv simple.svg
#eisen_to_newick.py data/data.gtr data/data.cdt data/data.nwk

